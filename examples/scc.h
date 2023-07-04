#include <atomic>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "helper/ligra_light.h"
#include "find_if.h"


constexpr auto max_label = std::numeric_limits<vertex>::max();

template <typename vertex, typename graph>
auto single_reach(vertex start, const graph& G, const graph& GT) {
  auto visited = parlay::tabulate<std::atomic<bool>>(G.size(), [&] (long i) {
    return (i==start) ? true : false; });

  auto edge_f = [&] (vertex u, vertex v) -> bool {
    bool expected = false;
    return visited[v].compare_exchange_strong(expected, true);};
  auto cond_f = [&] (vertex v) { return !visited[v];};
  auto frontier_map = ligra::edge_map(G, GT, edge_f, cond_f);

  auto frontier = ligra::vertex_subset(start);
  while (frontier.size() > 0) {
    frontier = frontier_map(frontier.to_seq());
  }
  return visited;
}


template <class Graph>
auto find_scc(Graph& G, Graph& GT){
    // Initial label & done and Trimming.
    auto label = parlay::tabulate<vertex>(G.size(), [&](size_t i){
        return G[i].size()==0||GT.size()==0 ? i: max_label;});
    auto done = parlay::tabulate<std::atomic<bool>>(G.size(), [&](size_t i){
        return G[i].size()==0||GT.size()==0 ? true: false;});
    
    // random permutation
    auto order = parlay::random_permutation((int)G.size());
    int start = ::find_if(order, [&] (vertex i) { return !done[order[i]];});
    
    // first round: single_reach
    auto forw_reach = single_reach(start, G, GT);
    auto back_reach = single_reach(start, GT, G);
    parlay::parallel_for(0, G.size(), [&](size_t i){
        done[i] = forw_reach[i] && back_reach[i]; 
        label[i] = order[start];
    });

    // multi_reach

    return label;
}