#include <atomic>
#include <utility>
#include <limits>
#include <unordered_set>
#include <bitset>
#include <xmmintrin.h>
#include <vector>


#include <parlay/primitives.h>
#include <parlay/sequence.h>
// #include "helper/ligra_light.h"
#include "BFS_ligra.h"

using vertex = uint32_t;
using distance = uint8_t;
using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
using graph = nested_seq;

struct node_info {                 // On return
  std::atomic<ulong> visited;      // S_r^0 (from Akiba, Iwata, Yoshina paper)
  std::atomic<ulong> visited_pre; // S_r^{-1}
  std::atomic<distance> d;             // P
  node_info() : visited(0), visited_pre(0), d((distance)100) {}
};

template<int kBitParallelRounds = 64>
class PrunedLandmarkLabeling {
  private:
    static const distance INF8;  // For unreachable pairs
    struct index_t {
        parlay::sequence<distance> bpspt_d;
        parlay::sequence<uint64_t> bpspt_s0;
        parlay::sequence<uint64_t> bpspt_s1;
        parlay::sequence<vertex> spt_v;
        parlay::sequence<distance> spt_d;

        index_t(size_t size) : bpspt_d(size), 
            bpspt_s0(size), bpspt_s1(size) {}
    }; 
    //   __attribute__((aligned(64)));  // Aligned for cache lines
    void BitPar_BFS(const graph &G);
    void Pruned_labeling(const graph &G);
    parlay::sequence<bool> usd;
    parlay::sequence<vertex> orders;
    parlay::sequence<index_t> index_;
 public:
    // Constructs an index from a graph. 
    // Return true if construct successfully
    bool ConstructIndex(const graph& G);

    // Returns distance vetween vertices |v| and |w| if they are connected.
    // Otherwise, returns |INT_MAX|.
    //   inline distance QueryDistance(vertex v, vertex w);
    PrunedLandmarkLabeling(){}
};

template<int kBitParallelRounds>
const distance PrunedLandmarkLabeling<kBitParallelRounds>::INF8 = 100;

template<int kBitParallelRounds>
void PrunedLandmarkLabeling<kBitParallelRounds>::
BitPar_BFS(const graph &G){
    vertex n = G.size();
    vertex r = 0;
    for (int i_bpspt = 0; i_bpspt < kBitParallelRounds; ++i_bpspt) {
        while (r < n && usd[r]) ++r;
        if (r == n) {
            for (int v = 0; v < n; ++v) index_[v].bpspt_d[i_bpspt] = INF8;
            continue;
        }
        usd[r] = true;
        parlay::sequence<node_info> nodes(n);
        nodes[r].d=(uint8_t)0;
        nodes[r].visited=1ul;
        printf("source: %d ", orders[r]);
        parlay::sequence<vertex> vertices(64);
        int ns = 0;
        vertices[ns]=r;
        auto ngb = parlay::sort(G[r]);
        for (size_t i = 0; i < G[r].size(); ++i) {
            vertex v = ngb[i];
            if (!usd[v]) {
                ns++;
                usd[v] = true;
                printf("%d ", orders[v]);
                nodes[v].visited =1ul << ns;
                vertices[ns]=v;
                if (ns == 63) break;
            }
        }
        printf("\n");
        distance round=0;
        auto edge_f = [&] (vertex u, vertex v) -> bool {
            uint64_t visit_u = nodes[u].visited.load();
            uint64_t visit_v = nodes[v].visited_pre.load();
            if ((visit_u | visit_v) != visit_v){
              nodes[v].visited_pre.fetch_or(visit_u);
              distance old_d = nodes[v].d.load();
              return (old_d!=round) && 
                nodes[v].d.compare_exchange_strong(old_d, round);
            }else{return false;}
        };
        auto cond_f = [&] (vertex v) {return !(nodes[v].visited.load()&1);};
        auto frontier_map = ligra::edge_map(G, G, edge_f, cond_f);
        // auto frontier = ligra::vertex_subset(vertices);
        auto frontier = ligra::vertex_subset<vertex>();
        frontier.add_vertices(vertices);
        while (frontier.size() > 0) {
            // printf("round %d frontier %d\n", round, frontier.size());
            round++;
            frontier = frontier_map(frontier);
            frontier.apply([&](vertex v){
              if ((nodes[v].visited_pre.load() & 1)==1){
                uint64_t tmp = nodes[v].visited.load();
                nodes[v].visited= nodes[v].visited_pre.load();
                nodes[v].visited_pre=tmp;
              }else{
                nodes[v].visited = nodes[v].visited_pre.load();
              }
            });
        }
        parlay::internal::timer t("inner ");
        t.start();
        parlay::parallel_for(0, n, [&](vertex v) {
            index_[orders[v]].bpspt_d[i_bpspt] = nodes[v].d.load();
            index_[orders[v]].bpspt_s1[i_bpspt] = nodes[v].visited_pre.load();
            index_[orders[v]].bpspt_s0[i_bpspt] = nodes[v].visited.load() &~nodes[v].visited_pre.load();
        });
        t.next("write back");
    }
}
/**/
template<int kBitParallelRounds>
void PrunedLandmarkLabeling<kBitParallelRounds>::
Pruned_labeling(const graph &G){
  vertex n = G.size();
  // parlay::sequence<std::pair<std::vector<vertex>, std::vector<distance> > >
  //       tmp_idx(n, std::make_pair(std::vector<vertex>(1, n),
  //                            std::vector<distance>(1, INF8)));
  parlay::sequence<std::pair<std::vector<vertex>, std::vector<distance> > >
        tmp_idx(n);
  parlay::sequence<bool> vis(n);  
  parlay::sequence<vertex> que(n);
  parlay::sequence<distance> dst_r(n+1, INF8); 
  long int total_size = 0; 
  vertex que_t0=0, que_t1=0, que_h=0;
  for (vertex r = 0; r<n; r++){
    if (usd[r]) continue;
    index_t &idx_r = index_[orders[r]];
    std::pair<std::vector<vertex>, std::vector<distance>> & tmp_idx_r = tmp_idx[r];
    for (size_t i = 0; i<tmp_idx[r].first.size(); i++){
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
        // std::pair<std::vector<vertex>, std::vector<distance>> & tmp_idx_v = tmp_idx[v];
        index_t & idx_v = index_[orders[v]];

        // Prefetch
        // _mm_prefetch(&idx_v[0], _MM_HINT_T0);
        // _mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
        // _mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

        for (int i = 0; i < kBitParallelRounds; ++i) {
            int td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
            if (td - 2 <= d) {
                td +=
                    (idx_r.bpspt_s1[i] & idx_v.bpspt_s1[i]) ? -2 :
                    ((idx_r.bpspt_s0[i] & idx_v.bpspt_s1[i]) |
                    (idx_r.bpspt_s1[i] & idx_v.bpspt_s0[i]))
                    ? -1 : 0;
                if (td <= d) goto pruned;
            }
        }

        for (size_t i = 0; i<tmp_idx[v].first.size(); ++i){
          vertex w = tmp_idx[v].first[i];
          distance td= tmp_idx[v].second[i]+dst_r[w];
          if (td<=d){ goto pruned;}
        }

        // Traverse
        label_size++;
        // tmp_idx[v].first.back() = r;
        // tmp_idx[v].second.back()=d;
        // tmp_idx[v].first.push_back(n);
        // tmp_idx[v].second.push_back(INF8);

        tmp_idx[v].first.push_back(r);
        tmp_idx[v].second.push_back(d);
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
    for (size_t i = 0; i<tmp_idx[r].first.size(); ++i){
        dst_r[tmp_idx[r].first[i]]=INF8;
    }
    usd[r]=true;
    // #ifdef DEBUG
    // printf("vertex %d add labels %d\n", r, label_size);
    // #endif
    total_size += que_h;
  }
  for (vertex v = 0; v < n; ++v) {
    size_t k = tmp_idx[v].first.size();
    index_[orders[v]].spt_v.reserve(k);
    index_[orders[v]].spt_d.reserve(k);
    for (int i = 0; i < k; ++i) index_[orders[v]].spt_v[i] = tmp_idx[v].first[i];
    for (int i = 0; i < k; ++i) index_[orders[v]].spt_d[i] = tmp_idx[v].second[i];
    // tmp_idx[v].first.clear();
    // tmp_idx[v].second.clear();
  }
  long total = reduce(map(tmp_idx, [&](auto p){return p.first.size();}));
  printf("average normal label size %f\n", (double)total/(double)n);
  printf("total explored nodes: %ld\n", total_size);
}


template<int kBitParallelRounds>
bool PrunedLandmarkLabeling<kBitParallelRounds>::
ConstructIndex(const graph& G) {
    parlay::internal::timer t("Indexing");
    vertex n =G.size();
    index_=parlay::sequence<index_t>(n,index_t(kBitParallelRounds));
    auto sizes = parlay::tabulate(n, [&](vertex j){return std::make_pair(G[j].size(), j);});
    sizes = parlay::sort(sizes, [&] (auto a, auto b) {return a.first > b.first;});
    // std::vector<std::pair<float, vertex> > deg(n);
    // for (vertex v = 0; v < n; ++v) {
    //   // We add a random value here to diffuse nearby vertices
    //   deg[v] = std::make_pair(G[v].size() + float(rand()) / RAND_MAX, v);
    // }
    // std::sort(deg.rbegin(), deg.rend());

    orders= parlay::sequence<vertex>(n);
    parlay::sequence<vertex> inv_orders(n);
    for (vertex i = 0; i < n; ++i) {
      orders[i] = sizes[i].second;
      inv_orders[sizes[i].second]=i;
    }

    auto newG =  parlay::tabulate(n, [&] (vertex i) {
      return parlay::map(G[orders[i]], [&](vertex j){return inv_orders[j];});});  

    t.next("loading");
    for (int repeat = 0; repeat<5; repeat++){
      usd = parlay::tabulate(n, [](vertex i){return false;});
      BitPar_BFS(newG);
      t.next("BitParallel BFS");
    }
    for (int i = 0; i<30; i++){
      std::cout << std::hex << index_[i].bpspt_s0[0] << ", " << index_[i].bpspt_s1[0] << ", " << ((unsigned int) index_[i].bpspt_d[0]) << std::endl;
    }
    // Pruned_labeling(newG);
    // t.next("Pruned Labeling");
    return true;
}