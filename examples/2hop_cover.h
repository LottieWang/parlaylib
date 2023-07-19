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
#define K_CAP 32

using namespace parlay;

template <typename vertex, typename distance>
class PrunedLandmarkLabeling{
  public:
    // constexpr distance INF_D = std::numeric_limits<distance>::max();
    distance INF_D = std::numeric_limits<distance>::max();
    sequence<sequence<vertex>> index_v;
    sequence<sequence<distance>> index_d;

    PrunedLandmarkLabeling(vertex n){
      index_v = sequence<sequence<vertex>>(n);
      index_d = sequence<sequence<distance>>(n);
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
		  PrunedLandmarkLabeling<vertex, distance>& L){
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

template <class Graph, typename vertex, typename distance>
auto create_PrunedLandmarkLabeling(Graph& G) {
  vertex n = G.size();

  // auto ids = tabulate(n, [&](vertex i){return (vertex)i;});
  // auto increase_orders = internal::integer_sort(make_slice(ids), [&](vertex i){return G[i].size();});
  // printf("Check sorted\n");
  auto degrees = tabulate(n, [&](vertex i){return std::make_pair(G[i].size(), i);});
  auto orders = parlay::map(parlay::sort(make_slice(degrees), std::greater<std::pair<size_t, vertex>>()),
          [&](auto kv){return std::get<1>(kv);});
  printf("get orders\n");
  auto inv_orders = sequence<vertex>(n);
  parallel_for(0, n, [&] (vertex i){ inv_orders[orders[i]] = i;});
  printf("get inverse orders\n");

  // empty le-lists
  printf("Initial Labelings\n");
  PrunedLandmarkLabeling<vertex, distance> L(n);
  
  
  // BFS one by one; TODO: change latter rounds to prefix-doubling
  for(vertex i = 0; i <n; i++){
    vertex r = orders[i];
    printf("BFS from %u, degree: %lu\n", r,G[r].size());
    Single_PrunedBFS(G, r, orders, inv_orders,  L);
  }
  

  return L.pack();
}
