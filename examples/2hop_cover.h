#include <atomic>
#include <utility>
#include <limits>
#include <unordered_set>
#include <bitset>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include "BFS_ligra.h"


// #include "helper/ligra_light.h"
uint32_t n_bits= 64;
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

template <typename vertex, typename distance, int kNumBitParallelRoots>
class PrunedLandmarkLabeling{
  public:
    // constexpr distance INF_D = std::numeric_limits<distance>::max();
    distance INF_D = std::numeric_limits<distance>::max();
    sequence<sequence<vertex>> index_v;
    sequence<sequence<distance>> index_d;
    // sequence<sequence<uint64_t>> bitwise_indexv;
    sequence<sequence<distance>> bitwise_indexd;

    PrunedLandmarkLabeling(vertex n){
      index_v = sequence<sequence<vertex>>(n);
      index_d = sequence<sequence<distance>>(n);
      bitwise_indexd = sequence<sequence<distance>>(n);
    }
    void insert(vertex i, vertex j, distance dist){
      index_v[i].push_back(j);
      index_d[i].push_back(dist);
    }
    std::pair<sequence<sequence<vertex>>, sequence<sequence<distance>> > pack() {
      return std::make_pair(index_v, index_d);
  }
};

template <typename vertex, typename distance, typename graph, int kNumBitParallelRoots>
auto Single_PrunedBFS(graph& G, vertex r,
		  sequence<vertex>& orders,
      sequence<vertex>& inv_orders,
		  PrunedLandmarkLabeling<vertex, distance, kNumBitParallelRoots>& L){
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


template <typename vertex, typename distance, typename graph>
  auto Single_BFS(graph& G, vertex r){
    auto result = BFS(r, G, G);
    for (vertex i = 0; i<result.size(); i++){
      printf("round %u, frontier %u\n", i+1, result[i].size());
    }
    distance INF8 = std::numeric_limits<distance>::max();
    auto distances = parlay::tabulate(G.size(), [&](vertex i){return INF8;});
    parlay::parallel_for(0, result.size(), [&](distance dist){
      parallel_for(0, result[dist].size(), [&](vertex j){
        distances[result[dist][j]]=dist;
      });
    });
    return distances;
  }

template <typename vertex, typename distance, typename graph>
  auto BFS_Pruned(graph& G, vertex r, sequence<sequence<distance>>& bitwise_dist, sequence<sequence<vertex>>& label_v, 
  sequence<sequence<distance>>& label_d,
  sequence<vertex>& orders,
  sequence<vertex>& inv_orders){
  sequence<vertex> que;
  std::unordered_set<vertex> visited;
  sequence<std::pair<vertex, distance>> distances;
  vertex que_t0=0, que_t1=0, que_h=0;
  que.push_back(r); que_h++;
  visited.insert(r);
  que_t1 = que_h;
  for (distance dist = 0; que_t0<que_h; ++dist){
    for (vertex que_i=que_t0; que_i < que_t1; ++que_i){
      vertex v = que[que_i];
      if (inv_orders[v] < inv_orders[r]){continue;}
      // Prune 1
      bool prune = false;
      for (vertex i = 0; i<bitwise_dist[0].size(); i++){
        if (prune) {break;}
        distance old_d = bitwise_dist[r][i]+bitwise_dist[v][i];
        if (old_d <= dist){
          prune=true;
          break;
        }
      }
      if (prune) {continue;}
      // Prune 2
      vertex i = 0, j= 0;
      while (i<label_v[v].size() && j < label_v[r].size()){
        if (label_v[v][i]==label_v[r][j]){
          distance old_d = label_d[v][i]+label_d[r][j];
          if (old_d<=dist){
            prune = true;
            break;
          }
          i++; j++;
        }else {
          i += label_v[v][i] < label_v[r][j] ? 1 : 0;
          j += label_v[v][i] > label_v[r][j] ? 1 : 0;
        }
      }
      if (prune) {continue;}
      // traverse
      distances.push_back(std::make_pair(v, dist));
      for (vertex i = 0; i<G[v].size(); i++){
        vertex u = G[v][i];
        if (std::get<1>(visited.insert(u))){
          que.push_back(u); que_h++;
        }
      }
    }
      que_t0 = que_t1;
      que_t1 = que_h;
  }
  return distances;
}

template <typename vertex, typename distance, typename graph>
auto BitwiseMulti_BFS(graph&GT, vertex offset,
      sequence<vertex>& orders,
      sequence<vertex>& inv_orders,
      sequence<sequence<distance>>& Ldistances){
  vertex n = GT.size();
  // vertex n_bits = sizeof(uint64_t)*8;
  distance INF = std::numeric_limits<distance>::max();
  uint64_t ALL = std::numeric_limits<uint64_t>::max();
  auto visited = tabulate<uint64_t>(n, [](vertex i){return 0;});
  auto distances = tabulate(n, [&](vertex i){return tabulate<distance>(n_bits, [INF](vertex i){return INF;});});
  auto changed = tabulate<bool>(n, [&](vertex i){return false;});
  for (int i = 0; i<n_bits; i++){
    visited[orders[offset+i]] = (uint64_t)1<<i;
    distances[orders[offset+i]][i] = 0;
    changed[orders[offset+i]] = true;
  }
  distance dist = 0;
  // not sure whether reduce(sequence<bool>) is true/false or size_t
  while (parlay::count(changed, true) !=0){
    dist ++;
    printf("  round %u, frontier %u\n", dist, parlay::count(changed, true));
    visited = tabulate(n, [&](vertex u){
      bool local_change=false;
      uint64_t visited_next = visited[u];
      if (visited_next!= ALL){
        parallel_for(0, GT[u].size(), [&](vertex j){
          vertex v = GT[u][j];
          if (visited_next!= ALL){
            visited_next |= visited[v];
          }
        }, 1000);
        uint64_t diff = visited_next &~ visited[u];
        if (diff != 0){
          local_change=true;
          for (int i = 0; i<n_bits; i++){
            if ((diff & (uint64_t)1<<i) != 0){
              distances[u][i]=dist;
            }
          }
        }
      }
      changed[u]=local_change;
      return visited_next;
    });
  }
  // May limit the performance, rewrite latter
  parlay::parallel_for(0, n, [&](size_t i){Ldistances[i].append(distances[i]);});
}

template <typename vertex, typename graph, typename distance, int kNumBitParallelRoots>
auto Multi_BFS(graph& G, vertex start, vertex end,
  sequence<vertex>& orders,
  sequence<vertex>& inv_orders,
  PrunedLandmarkLabeling<vertex, distance, kNumBitParallelRoots>& L){
  auto label_pairs = parlay::map(make_slice(orders).cut(start, end), [&](vertex r){
    sequence<vertex> que;
    std::unordered_set<vertex> visited;
    sequence<std::pair<vertex, std::pair<vertex, distance>>> distances;
    // vertex n_bits = sizeof(uint64_t)*8;
    vertex que_t0=0, que_t1=0, que_h=0;
    que.push_back(r); que_h++;
    visited.insert(r);
    que_t1 = que_h;    
    for (distance dist = 0; que_t0<que_h; ++dist){
      // printf("  source %u round %u explore %u\n", r, dist, que_t1-que_t0);
      for (vertex que_i=que_t0; que_i < que_t1; ++que_i){
        vertex v = que[que_i];
        // if (inv_orders[v] < start) continue;
        // Prune 0; not sure whether it is correct;
        if (inv_orders[v] < inv_orders[r]){
          // printf("Prune 0\n");
          continue;}
        // Prune 1
        bool prune = false;
        for (vertex i = 0; i<L.bitwise_indexd[0].size() && (!prune); i++){
          distance old_d = L.bitwise_indexd[r][i]+L.bitwise_indexd[v][i];
          if (old_d <= dist){
            prune=true;
            break;
          }
        }
        if (prune) {continue;}
        // printf("  source %u vertex %u begin prune cheking 2\n", r,v);
        vertex i = 0, j= 0;
        while (i<L.index_v[v].size() && j < L.index_v[r].size()){
          if (L.index_v[v][i]==L.index_v[r][j]){
            distance old_d = L.index_d[v][i]+L.index_d[r][j];
            if (old_d<=dist){
              prune = true;
              break;
            }
            i++; j++;
          }else {
            if (L.index_v[v][i] < L.index_v[r][j]){
              i++;
            }else{
              j++;
            }
          }
        }
        if (prune) {
          // printf("Prune 2\n");
          continue;
        }
        // traverse
        distances.push_back(std::make_pair(v, std::make_pair(inv_orders[r], dist)));
        for (vertex i = 0; i<G[v].size(); i++){
          vertex u = G[v][i];
          if (std::get<1>(visited.insert(u))){
            // que[que_h++]=u;
            que.push_back(u); que_h++;
          }
        }
      }
      que_t0 = que_t1;
      que_t1 = que_h;
    }
    return distances;
  });
  auto pairs = flatten(label_pairs);
  printf("# labels added: %lu\n", pairs.size());
  auto indexs = group_by_index(pairs, G.size());
  parallel_for(0, G.size(), [&](vertex i){
    L.index_v[i].append(map(indexs[i], [&](auto p){return std::get<0>(p);}));
    L.index_d[i].append(map(indexs[i], [&](auto p){return std::get<1>(p);}));
  });
}
template <class Graph, typename vertex, typename distance, int kNumBitParallelRoots>
auto create_PrunedLandmarkLabeling(Graph& G) {
  vertex n = G.size();

  auto degrees = parlay::sort(tabulate(n, [&](vertex i){return std::make_pair(G[i].size(), i);}));
  auto orders = tabulate(n, [&](vertex i){return std::get<1>(degrees[n-1-i]);});
  printf("get orders\n");
  auto inv_orders = sequence<vertex>(n);
  parallel_for(0, n, [&] (vertex i){ inv_orders[orders[i]] = i;});
  printf("get inverse orders\n");

  // empty le-lists
  printf("Initial Labelings\n");
  PrunedLandmarkLabeling<vertex, distance, kNumBitParallelRoots> L(n);
  
  // vertex n_bits = 64;
  for(vertex i = 0; i <kNumBitParallelRoots; i+=1){
    printf("BFS round %u\n", i);
    BitwiseMulti_BFS<vertex, distance,Graph>(
      G, i*n_bits, orders, inv_orders, L.bitwise_indexd);
    printf("bitwiseLabeling_size: %lu x %lu\n", (L.bitwise_indexd).size(), L.bitwise_indexd[0].size());    
  }
  vertex step = 10;
  double beta = 1.2;
  vertex start = kNumBitParallelRoots*n_bits;
  int round = 0;
  while (start < n){
    vertex end = std::min(start + step, n);
    printf("Multi_BFS round %u: batch_size %u\n", round, end-start);
    Multi_BFS(G, start, end, orders,inv_orders,L);
    start = end;
    step = std::ceil(step*beta);
    round++;
  }
  // for (int i = kNumBitParallelRoots*n_bits; i<n; i++){
  //   auto distances = BFS_Pruned(G, orders[i], L.bitwise_indexd, 
  //       L.index_v,L.index_d,orders, inv_orders);
  //   printf("round %d add %u labels\n", i, distances.size());
  //   parallel_for(0, distances.size(), [&](vertex j){
  //     auto p = distances[j];
  //     L.index_v[std::get<0>(p)].push_back(i);
  //     L.index_d[std::get<0>(p)].push_back(std::get<1>(p));
  //   });
  // }

  // return L.pack();
}
