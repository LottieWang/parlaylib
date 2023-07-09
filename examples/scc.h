#include <atomic>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "helper/ligra_light.h"
#include "find_if.h"
#include "hash_table.h"



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

// TODO: design a nested-loop version of multi-search, look at (key,val) insert for all edges.

template<typename vertex, typename Fa, typename graph, 
  typename Get = ligra::identity<typename graph::value_type::value_type>>
auto edge_map_sparse(parlay::sequence<vertex> const &vertices,
  const graph& G, Fa f, Get get={}) {
    using edge = typename graph::value_type::value_type;
    auto nested_pairs = parlay::map(vertices, [&] (vertex v) {
      return delayed::map(G[v], [=] (edge e) {
        return std::pair(v, e);});});
    auto pairs = delayed::flatten(nested_pairs);
    return delayed::to_sequence(delayed::map_maybe(pairs, [&] (auto p) {
      auto [u,e] = p;
      vertex v = get(e);
      return ((f(u,v))
              ? std::make_optional(v)
              : std::nullopt);}));
  }
// using edge_map for load balance
template <typename vertex, typename label_type, typename graph>
auto multi_search_edge_map(parlay::sequence<vertex> sources, const graph& G, const graph& GT, parlay::sequence<label_type>& label, size_t m, size_t scc_offset){
  using Table = hash_table<vertex, vertex>;
  Table T(m, std::make_pair(G.size(),G.size()));
  // insert all frontier pairs to hash_table
  parlay::parallel_for(0, sources.size(), [&](size_t i){
    T.insert(sources[i], (vertex)(scc_offset+i));});

  auto edge_f = [&](vertex u, vertex v) -> bool {
    if (label[u]!=label[v]) return false;
    iter_k <vertex, vertex, Table> entry(u, T);
    entry.init();
    bool label_changed=false;
    while (entry.has_next()){
      vertex source=entry.get();
      label_changed |= T.insert(v,source);
      if (T.overfull) break;
    }
    return label_changed && T.overfull;
  };
  // // cond_f is used for backward, but not used here, can be anything
  // auto cond_f = [&] (vertex v){ return label[v]<scc_offset+sources.size();};
  // auto frontier_map = ligra::edge_map(G, GT, edge_f, cond_f);
  // auto frontier = ligra::vertex_subset(sources);
  parlay::sequence<vertex> frontier = sources;
  while (frontier.size()>0){
    frontier = edge_map_sparse(frontier, G, edge_f);
  }
  return T;
}

template <typename vertex, typename graph>
auto find_scc(graph& G, graph& GT){
    // Initial label & done and Trimming.
    using label_type = uint64_t;
    constexpr label_type TOP_BIT = static_cast<label_type>(1) << (8*sizeof(label_type)-1);
    constexpr label_type VAL_MASK = std::numeric_limits<vertex>::max();
    using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
    size_t scc_offset=1;
    auto label = parlay::tabulate<std::atomic<label_type>>(G.size(), [](size_t i){return 0;});
    
    // random permutation
    auto order = parlay::random_permutation((int)G.size());

    // Trimming
    auto trim = parlay::filter(order, [&](size_t i){
      return G[order[i]].size()==0 || GT[order[i]].size()==0;});
    parlay::parallel_for(0, trim.size(), [&](size_t i){label[trim[i]]=(scc_offset+i) | TOP_BIT;});
    scc_offset+=trim.size();
    
    // first round: single_reach
    vertex source = ::find_if(order, [&] (size_t i) { return !(label[order[i]] & TOP_BIT);});
    auto forw_reach = single_reach(source, G, GT);
    auto back_reach = single_reach(source, GT, G);
    parlay::parallel_for(0, G.size(), [&](size_t i){
        if (forw_reach[i] && back_reach[i])
          label[i] = TOP_BIT | scc_offset;
    });
    scc_offset+=1;

    // remove identidied SCCs
    auto vertices = parlay::filter(order, [&](size_t i){return !(label[order[i]] & TOP_BIT);});
    
    // multi_reach
    float beta = 1.5;
    int step = 2;
    int start = 0;
    int end;
    int n = vertices.size();
    std::cout << "n_remain is " << n << std::endl;
    while (n > 0){
      end = std::min(n, start+step);
      auto sources = parlay::filter(vertices.cut(start, end), [&](size_t i){return !(label[vertices[i]]|TOP_BIT);});
      size_t n_remain = parlay::internal::count_if_index(n-start, [&](size_t i){return !(TOP_BIT & label[vertices[start+i]]);});
      size_t m = n_remain;
      auto fwd_table = multi_search_edge_map(sources, G, GT, label, m, scc_offset);
      while (fwd_table.overfull==true){
        m=2*m;
        fwd_table = multi_search_edge_map(sources, G, GT, label , m, scc_offset);
      }
      auto bwd_table = multi_search_edge_map(sources, GT, G, label, m, scc_offset);
      while (bwd_table.overfull==true){
        m=2*m;
        bwd_table = multi_search_edge_map(sources, GT, G, label, m, scc_offset);
      }
      auto small_table = (fwd_table.size()<= bwd_table.size())? fwd_table: bwd_table;
      auto large_table = (fwd_table.size()>bwd_table.size())? fwd_table: bwd_table;

      auto map_f1 = [&](std::pair<vertex, vertex> a){
        vertex u = parlay::internal::get_key(a);
        label_type l = parlay::internal::get_val(a);
        if (large_table.contain(a)){
          parlay::write_max(&label[u], (TOP_BIT | l), std::less<label_type>());
        }else{
          parlay::write_max(&label[u], l, std::less<label_type>());
        }
      };
      small_table.map(map_f1);

      auto map_f2 = [&](std::pair<vertex, vertex> a){
        vertex u = parlay::internal::get_key(a);
        label_type l = parlay::internal::get_val(a);
        parlay::write_max(&label[u], l, std::less<label_type>());
      };
      large_table.map(map_f2);

      n-=(end-start);
      start += step;
      scc_offset+=step;
      step = step * beta;
    }

    auto pairs = parlay::delayed::tabulate(G.size(), [&] (vertex i) {
        return std::pair(label[i] & VAL_MASK, i);});
    auto lists = parlay::group_by_index(pairs, scc_offset);
    auto SCCs = parlay::pack(lists, parlay::delayed_map(lists, [=] (auto b) {return b.size()>0;}));
    // auto SCCs= parlay::filter(lists, [&](size_t i){return lists[i].size()>0;});
    // std::cout << "number of SCC: " << SCCs.size() << std::endl;
    // return SCCs;
    std::cout << "number of lists: " << SCCs.size() << std::endl;
    return SCCs;
    // return parlay::map(label, [=](label_type l){return (vertex)(l & VAL_MASK);});
}