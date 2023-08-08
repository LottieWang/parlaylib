#include <atomic>
#include <utility>
#include <limits>
#include <unordered_set>
#include <bitset>
#include <xmmintrin.h>
#include <vector>


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
    static const distance INF8 = 100;
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
  auto delta = tabulate<std::atomic<distance>>(n, [&](vertex i){return (i==r)?  0: L.INF8;});

  // initialize the distance from r to others
  auto delta_r = tabulate(n, [&](vertex i){return L.INF8;});
  parallel_for(0, (L.index_v[r]).size(), [&](vertex i){
    delta_r[orders[L.index_v[r][i]]]=L.index_d[r][i];});
  
  // function to apply on each edge s->d
  auto edge_f = [&] (vertex s, vertex d) -> bool {
    if (delta[d].compare_exchange_strong(L.INF8, dist)){
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
      if (delta_r[v] != L.INF8 && delta_r[v]+L.index_d[d][i] <= dist) return false;
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
    static const uint8_t INF8 = 100;
    auto distances = parlay::tabulate(G.size(), [&](vertex i){return INF8;});
    parlay::parallel_for(0, result.size(), [&](distance dist){
      parallel_for(0, result[dist].size(), [&](vertex j){
        distances[result[dist][j]]=dist;
      });
    });
    return distances;
  }

template <typename vertex, typename distance, typename graph, int kNumBitParallelRoots>
  auto Pruned_labeling(graph& G, sequence<sequence<distance>>& bitwise_dist, 
                      sequence<bool>& usd, PrunedLandmarkLabeling<vertex, distance, kNumBitParallelRoots>& L){
  vertex n = G.size();
  static const distance INF8 = 100;  
  std::vector<std::pair<std::vector<vertex>, std::vector<distance> > >
        tmp_idx(n, std::make_pair(std::vector<vertex>(1, n),
                             std::vector<distance>(1, INF8)));
  std::vector<bool> vis(n);  
  std::vector<vertex> que(n);
  std::vector<distance> dst_r(n+1, INF8); 
  long int total_size = 0; 
  vertex que_t0=0, que_t1=0, que_h=0;
  for (vertex r = 0; r<n; r++){
    if (usd[r]) continue;
    sequence<distance> &idx_r = bitwise_dist[r];
    std::pair<std::vector<vertex>, std::vector<distance>> & tmp_idx_r = tmp_idx[r];
    for (size_t i = 0; i<tmp_idx_r.first.size(); i++){
      dst_r[tmp_idx_r.first[i]]=tmp_idx_r.second[i];
    }
    int que_t0=0, que_t1=0, que_h=0;
    que[que_h++]=r;
    vis[r]=true;
    que_t1=que_h;
    long int label_size = 0;
    for (distance d = 0; que_t0<que_h; d++){
      for (vertex que_i = que_t0; que_i <que_t1; que_i++){
        vertex v = que[que_i];
        if (usd[v]) continue;
        std::pair<std::vector<vertex>, std::vector<distance>> & tmp_idx_v = tmp_idx[v];
        sequence<distance> & idx_v = bitwise_dist[v];

        // Prefetch
        // _mm_prefetch(&idx_v[0], _MM_HINT_T0);
        // _mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
        // _mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

        for (vertex i = 0;i<n_bits*kNumBitParallelRoots; i++){
          distance td = idx_r[i]+idx_v[i];
          if (td<=d) goto pruned;
        }

        for (size_t i = 0; i<tmp_idx_v.first.size(); ++i){
          vertex w = tmp_idx_v.first[i];
          distance td= tmp_idx_v.second[i]+dst_r[w];
          if (td<=d) goto pruned;
        }

        // Traverse
        label_size++;
        tmp_idx_v.first.back() = r;
        tmp_idx_v.second.back()=d;
        tmp_idx_v.first.push_back(n);
        tmp_idx_v.second.push_back(INF8);
        for (size_t i = 0; i<G[v].size(); i++){
          vertex w = G[v][i];
          if (!vis[w]){
            que[que_h++] = w;
            vis[w]=true;
          }
        }
      pruned:
        {}
      }
      que_t0 = que_t1;
      que_t1 = que_h;
    }
    for (vertex i = 0; i<que_h; i++) vis[que[i]]=false;
    for (size_t i = 0; i<tmp_idx_r.first.size(); ++i){
      dst_r[tmp_idx_r.first[i]]=INF8;
    }
    usd[r]=true;
    // #ifdef DEBUG
    printf("vertex %d add labels %d\n", r,         label_size);
    // #endif
    total_size += que_h;
  }
  for (vertex v = 0; v < n; ++v) {
    L.index_v[v]=tabulate(tmp_idx[v].first.size(), [&](vertex i){return tmp_idx[v].first[i];});
    L.index_d[v]=tabulate(tmp_idx[v].second.size(), [&](vertex i){return tmp_idx[v].second[i];});
    // tmp_idx[v].first.clear();
    // tmp_idx[v].second.clear();
  }
  long total = reduce(map(tmp_idx, [&](auto p){return p.first.size()-1;}));
  printf("average normal label size %f\n", (double)total/(double)n);
  printf("total explored nodes: %ld\n", total_size);
}

template <typename vertex, typename distance, typename graph>
auto BitwiseMulti_BFS(graph&GT, vertex offset,
      sequence<sequence<distance>>& Ldistances,
      sequence<bool>& usd){
  vertex n = GT.size();
  // vertex n_bits = sizeof(uint64_t)*8;
  static const distance INF8 = 100;
  uint64_t ALL = std::numeric_limits<uint64_t>::max();
  auto visited = tabulate<uint64_t>(n, [](vertex i){return 0;});
  auto distances = tabulate(n, [&](vertex i){return tabulate<distance>(n_bits, [INF8](vertex i){return INF8;});});
  auto changed = tabulate<bool>(n, [&](vertex i){return false;});
  for (int i = 0; i<n_bits; i++){
    visited[offset+i] = (uint64_t)1<<i;
    distances[offset+i][i] = 0;
    changed[offset+i] = true;
    usd[offset+i]=true;
  }
  distance dist = 0;
  // not sure whether reduce(sequence<bool>) is true/false or size_t
  while (parlay::count(changed, true) !=0){
    dist ++;
    // printf("  round %u, frontier %u\n", dist, parlay::count(changed, true));
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
  // printf("# labels added: %lu\n", pairs.size());
  auto indexs = group_by_index(pairs, G.size());
  parallel_for(0, G.size(), [&](vertex i){
    L.index_v[i].append(map(indexs[i], [&](auto p){return std::get<0>(p);}));
    L.index_d[i].append(map(indexs[i], [&](auto p){return std::get<1>(p);}));
  });
}
template <class Graph, typename vertex, typename distance, int kNumBitParallelRoots>
auto create_PrunedLandmarkLabeling(Graph& G) {
  vertex n = G.size();
  sequence<bool> usd(n, false);
  parlay::internal::timer t("Indexing");
  // auto degrees = parlay::sort(tabulate(n, [&](vertex i){return std::make_pair(G[i].size(), i);}));
  // auto orders = tabulate(n, [&](vertex i){return std::get<1>(degrees[n-1-i]);});
  // printf("get orders\n");
  // auto inv_orders = sequence<vertex>(n);
  // parallel_for(0, n, [&] (vertex i){ inv_orders[orders[i]] = i;});
  // printf("get inverse orders\n");
  std::vector<std::pair<float, vertex> > deg(n);
  for (vertex v = 0; v < n; ++v) {
    // We add a random value here to diffuse nearby vertices
    deg[v] = std::make_pair(G[v].size() + float(rand()) / RAND_MAX, v);
  }
  std::sort(deg.rbegin(), deg.rend());
  sequence<vertex> orders(n);
  sequence<vertex> inv_orders(n);
  for (vertex i = 0; i < n; ++i) {
    orders[i] = deg[i].second;
    inv_orders[deg[i].second]=i;
  }

  auto newG =  parlay::tabulate(n, [&] (vertex i) {
      return parlay::map(G[orders[i]], [&](vertex j){return inv_orders[j];});});  
  

  // empty le-lists
  printf("Initial Labelings\n");
  PrunedLandmarkLabeling<vertex, distance, kNumBitParallelRoots> L(n);
  
  // vertex n_bits = 64;
  t.next("loading");
  for(vertex i = 0; i <kNumBitParallelRoots; i+=1){
    // printf("BFS round %u\n", i);
    BitwiseMulti_BFS<vertex, distance,Graph>(
      newG, i*n_bits, L.bitwise_indexd, usd);
    // printf("bitwiseLabeling_size: %lu x %lu\n", (L.bitwise_indexd).size(), L.bitwise_indexd[0].size());    
  }
  t.next("bitwise index");
  // vertex step = 10;
  // double beta = 1.2;
  // vertex start = kNumBitParallelRoots*n_bits;
  // int round = 0;
  // while (start < n){
  //   vertex end = std::min(start + step, n);
  //   // printf("Multi_BFS round %u: batch_size %u\n", round, end-start);
  //   Multi_BFS(G, start, end, orders,inv_orders,L);
  //   start = end;
  //   step = std::ceil(step*beta);
  //   round++;
  // }
  // t.next("pruned index");
  Pruned_labeling<vertex, distance, Graph, kNumBitParallelRoots>(newG, L.bitwise_indexd, usd, L);
  t.next("pruned labeling");
  return L.pack();
}
