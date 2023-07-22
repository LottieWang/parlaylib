#include <atomic>
#include <utility>
#include <limits>

#include <parlay/primitives.h>
#include <parlay/sequence.h>


#include "helper/ligra_light.h"

// **************************************************************
// Parallel Breadth First Search (Using the Ligra interface)
// The graph is a sequence of sequences of vertex ids, representing
// the outedges for each vertex.
// Requires the transpose graph (i.e the back edges).
// Returns a sequence of sequences, with the ith element corresponding to
// all vertices at distance i (i.e. the i-th frontier during the search).
// This version uses the ligra interface.  See: helper/ligra_light.h
// Importantly it supports the forward-backwards optimization.  See:
//  Julian Shun, Guy E. Blelloch:
//  Ligra: a lightweight graph processing framework for shared memory.
//  PPOPP 2013:
// **************************************************************
using namespace parlay;

template <typename vertex, typename distance, typename bits>
class PrunedLandmarkLabeling{
  public:
    // constexpr distance INF_D = std::numeric_limits<distance>::max();
    distance INF_D = std::numeric_limits<distance>::max();
    sequence<sequence<vertex>> index_v;
    sequence<sequence<distance>> index_d;
    sequence<sequence<bits>> bitwise_indexv;
    sequence<sequence<sequence<distance>>> bitwise_indexd;

    PrunedLandmarkLabeling(vertex n){
      index_v = sequence<sequence<vertex>>(n);
      index_d = sequence<sequence<distance>>(n);
      bitwise_indexv = sequence<sequence<bits>>(4);
      bitwise_indexd = sequence<sequence<sequence<distance>>>(4);
    }
    void insert(vertex i, vertex j, distance dist){
      index_v[i].push_back(j);
      index_d[i].push_back(dist);
    }
    std::pair<sequence<sequence<vertex>>, sequence<sequence<distance>> > pack() {
      return std::make_pair(index_v, index_d);
  }
};

template <typename vertex, typename distance, typename graph>
auto Single_PrunedBFS(graph& G, vertex r,
		  sequence<vertex>& orders,
      sequence<vertex>& inv_orders,
		  PrunedLandmarkLabeling<vertex, distance, uint64_t>& L){
  vertex n = G.size();
  distance dist = 0;
  size_t n_explore = 0;
  // initial temp distance to INF
  auto delta = tabulate<std::atomic<distance>>(n, [&](vertex i){return (i==r)?  0: L.INF_D;});

  // initialize the distance from r to others
  auto delta_r = tabulate(n, [&](vertex i){return L.INF_D;});
  parallel_for(0, (L.index_v[r]).size(), [&](vertex i){
    delta_r[orders[L.index_v[r][i]]]=L.index_d[r][i];});
  
  // function to apply on each edge s->d
  auto edge_f = [&] (vertex s, vertex d) -> bool {
    if (delta[d].compare_exchange_strong(L.INF_D, dist)){
      L.insert(d, inv_orders[r], dist);
      return true;
    }
    return false;
  };

  
  auto cond_f = [&] (vertex d) {
    // if d's order is smaller than r, it must already settled.
    if (inv_orders[d] < inv_orders[r] || delta[d] <= dist) return false;
    // check whether delta[d] can already got by previous labels
    for (size_t i = 0; i<(L.index_v[d]).size(); i++){
      vertex v = orders[L.index_v[d][i]];
      if (delta_r[v] != L.INF_D && delta_r[v]+L.index_d[d][i] <= dist) return false;
    }
    return true;
  };

  // because the graph is symmetric, GT = G
  auto frontier_map = ligra::edge_map(G, G, edge_f, cond_f);
  auto frontier = ligra::vertex_subset(r);
  // printf("BFS from %u, initialize done\n", r);
  while (frontier.size() > 0) {
    n_explore+= frontier.size();
    dist++;
    frontier = frontier_map(frontier);
  }
  printf("  total vertices explored %lu, rounds %u\n", n_explore, dist);
}

template <typename vertex, typename distance, typename graph, typename bits>
auto BitwiseMulti_BFS(graph&GT, vertex offset,
      sequence<vertex>& orders,
      sequence<vertex>& inv_orders){
  vertex n = GT.size();
  vertex n_bits = sizeof(bits)*8;
  distance INF = std::numeric_limits<distance>::max();
  bits ALL = std::numeric_limits<bits>::max();
  auto visited = tabulate<bits>(n, [](vertex i){return 0;});
  auto distances = tabulate(n, [&](vertex i){return tabulate<distance>(n_bits, [INF](vertex i){return INF;});});
  auto changed = tabulate<bool>(n, [&](vertex i){return false;});
  for (int i = 0; i<n_bits;i++){
    visited[orders[offset+i]] = (bits)1<<i;
    distances[orders[offset+i]][i] = 0;
    changed[orders[offset+i]] = true;
  }
  distance dist = 0;
  // not sure whether reduce(sequence<bool>) is true/false or size_t
  while (dist == 0 || parlay::count(changed, true)!=0){
    dist ++;
    printf("  round %u, frontier %u\n", dist, parlay::count(changed, true));
    changed = tabulate(n, [&](vertex u){
      bool result=false;
      if (visited[u]!= ALL){
        parallel_for(0, GT[u].size(), [&](vertex j){
          vertex v = GT[u][j];
          if (visited[u]!= ALL && changed[v]){
            bits diff = visited[v]&~visited[u];
            if (diff != 0){
              result |= true;
              visited[u] |= visited[v];
              for (int i = 0; i<n_bits; i++){
                if ((diff & (bits)1<<i)!=0){
                  distances[u][i]=dist;
                }
              }
            }
          }
        }, 1000);
      }
      return result;
    });
  }
  return std::make_pair(visited, distances);
}

// template <typename vertex, typename graph>
// auto Multi_BFS(graph& G, sequence<vertex>& sources, ){

// }
template <class Graph, typename vertex, typename distance>
auto create_PrunedLandmarkLabeling(Graph& G) {
  vertex n = G.size();

  // auto ids = tabulate(n, [&](vertex i){return (vertex)i;});
  // auto increase_orders = internal::integer_sort(make_slice(ids), [&](vertex i){return G[i].size();});
  // printf("Check sorted\n");
  auto degrees = parlay::sort(tabulate(n, [&](vertex i){return std::make_pair(G[i].size(), i);}));
  // auto orders = parlay::map(parlay::sort(make_slice(tabulate(n, [&](vertex i){return std::make_pair(G[i].size(), i);})), std::greater<std::pair<size_t, vertex>>()),[&](auto kv){return std::get<1>(kv);});
  auto orders = tabulate(n, [&](vertex i){return std::get<1>(degrees[i]);});
  printf("get orders\n");
  auto inv_orders = sequence<vertex>(n);
  parallel_for(0, n, [&] (vertex i){ inv_orders[orders[i]] = i;});
  printf("get inverse orders\n");

  // empty le-lists
  printf("Initial Labelings\n");
  PrunedLandmarkLabeling<vertex, distance, uint64_t> L(n);
  
  
  for(vertex i = 0; i <std::min((unsigned)256,n); i+=64){
    printf("BFS round %u\n", i);
    // Single_PrunedBFS(G, r, orders, inv_orders,  L);
    auto result = BitwiseMulti_BFS<vertex, distance,Graph, uint64_t>(
      G, i,orders,inv_orders);
    L.bitwise_indexv[int(i/64)] = std::get<0>(result);
    L.bitwise_indexd[int(i/64)] = std::get<1>(result);
  }
  // for (vertex i = min(4,n); i<n; i++){

  // }

  // return L.pack();
}
