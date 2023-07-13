#include <atomic>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "helper/ligra_light.h"
#include "find_if.h"
#include "hash_table.h"



template <typename vertex,  typename graph>
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


template<typename vertex, typename Fa, typename Cond, typename graph, 
  typename Get = ligra::identity<typename graph::value_type::value_type>>
auto edge_map_sparse(parlay::sequence<vertex> const &vertices,
  const graph& G, Fa f, Cond cond, Get get={}) {
    using edge = typename graph::value_type::value_type;
    auto nested_pairs = parlay::map(vertices, [&] (vertex v) {
      return delayed::map(G[v], [=] (edge e) {
        return std::pair(v, e);});});
    auto pairs = delayed::flatten(nested_pairs);
    return delayed::to_sequence(delayed::map_maybe(pairs, [&] (auto p) {
      auto [u,e] = p;
      vertex v = get(e);
      return (cond(v) && (f(u,v))
              ? std::make_optional(v)
              : std::nullopt);}));
  }
// using edge_map for load balance
template <typename vertex, typename label_type, typename graph>
auto multi_search_edge_map(parlay::sequence<vertex> sources, const graph& G,  parlay::sequence<std::atomic<label_type>>& label, size_t m, vertex scc_offset){
  constexpr label_type TOP_BIT = static_cast<label_type>(1) << (8*sizeof(label_type)-1);
  auto bits = parlay::tabulate<std::atomic<bool>>(G.size(), [](size_t i){return false;});
  using Table = hash_table<vertex, vertex>;
  Table T(m, std::make_pair(G.size(),G.size()));
  // insert all frontier pairs to hash_table
  parlay::parallel_for(0, sources.size(), [&](size_t i){
    T.insert(sources[i], (vertex)(scc_offset+i));});
  
  auto edge_f = [&](vertex u, vertex v) -> bool {
    if (label[u]!=label[v]) return false;
    iter_k <vertex, vertex, Table> entry(u, T);
    entry.init();
    vertex source=entry.get();
    bool label_changed=T.insert(v,source);
    while (entry.has_next()){
      source=entry.get();
      label_changed |= T.insert(v,source);
      if (T.overfull) break;
    }
    if (label_changed && !(T.overfull)){
      bool expected = false;
      return bits[v].compare_exchange_strong(expected, true);
    }else{
      return false;
    }
  };
  auto cond_f = [&] (vertex v) { return !(label[v] & TOP_BIT) && !bits[v];};

  parlay::sequence<vertex> frontier = sources;
  int round = 0;
  while (frontier.size()>0){
    // std::cout << "    round " << round << " frontier.size " << frontier.size() << std::endl;
    parlay::parallel_for(0, frontier.size(), [&](size_t i){bits[frontier[i]]=false;});
    frontier = edge_map_sparse(frontier, G, edge_f, cond_f);
    round ++;
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
    vertex scc_offset=1;
    auto label = parlay::tabulate<std::atomic<label_type>>(G.size(), [](size_t i){return 0;});
    
    // random permutation
    // auto order = parlay::random_permutation((vertex)G.size());

    // Trimming
    auto id = parlay::sequence<vertex>::from_function(G.size(), [&] (vertex i) -> vertex { return i; });
    auto trim = parlay::filter(id, [&](vertex i){
      return G[i].size()==0 || GT[i].size()==0;});
    parlay::parallel_for(0, trim.size(), [&](size_t i){label[trim[i]]=(scc_offset+i) | TOP_BIT;});
    scc_offset+=trim.size();
    std::cout << "singleton SCCs: " << trim.size() << std::endl;
    auto non_zeros = parlay::filter(id, [&](vertex i){return G[i].size()!=0 && GT[i].size() !=0;});
    // std::cout << "  n_remain: " << parlay::internal::count_if_index(G.size(), [&](size_t i){return !(label[i]&TOP_BIT);}) << std::endl;
    std::cout << "   n_remain: " << non_zeros.size() << std::endl;
    
    // First round: single_reach
    // size_t source_id = ::find_if(order, [&] (vertex i) { return (label[i] & TOP_BIT)==0;});
    auto perms = parlay::random_shuffle(non_zeros);
    vertex source = perms[0];
    auto forw_reach = single_reach(source, G, GT);
    auto back_reach = single_reach(source, GT, G);
    parlay::parallel_for(0, perms.size(), [&](size_t i){
        vertex u = perms[i];
        if (forw_reach[u] && back_reach[u]){
          label[u] = scc_offset|TOP_BIT;
        }else if (forw_reach[u] || back_reach[u]){
          label[u] = scc_offset;
        }
    });
    std::cout << "first scc size: " <<parlay::internal::count_if_index(perms.size(), [&](size_t i){return forw_reach[perms[i]] && back_reach[perms[i]];}) << std::endl;
    scc_offset+=1;

    std::cout << "  n_remain: " << parlay::internal::count_if_index(G.size(), [&](size_t i){return !(label[i]&TOP_BIT);}) << std::endl;

    // remove identidied SCCs
    auto vertices = parlay::filter(perms, [&](vertex x){return !(label[x] & TOP_BIT);});
    
    // multi_reach
    float beta = 1.5;
    int step = 2;
    int start = 0;
    int end;
    int n = vertices.size();
    std::cout << "vertices.size " << n << std::endl;
    int round = 1;
    while (start< n){
      end = std::min(n, start+step);
      std::cout << "round " << round++ << ":  ";
      auto sources = parlay::filter(make_slice(vertices).cut(start, end), [&](vertex x){return !(label[x] & TOP_BIT);});
      std::cout << "source size: " << sources.size();
      size_t n_remain = parlay::internal::count_if_index(n-start, [&](size_t i){return !(TOP_BIT & label[vertices[start+i]]);});
      std::cout << "  n_remain: " << n_remain << std::endl;
      size_t m = n_remain;
      auto fwd_table = multi_search_edge_map(sources, G, label, m, scc_offset);
      while (fwd_table.overfull==true){
        m=2*m;
        fwd_table = multi_search_edge_map(sources, G, label, m, scc_offset);
        std::cout << "resize forward table" << std::endl;
      }      
      std::cout << "  forward table size " << fwd_table.size() << std::endl;
      auto bwd_table = multi_search_edge_map(sources, GT, label, m, scc_offset);
      while (bwd_table.overfull==true){
        m=2*m;
        bwd_table = multi_search_edge_map(sources, GT, label, m, scc_offset);
        std::cout << "resize backward table" << std::endl;
      }
      std::cout << "  backward table size " << bwd_table.size() << std::endl;
      auto small_table = (fwd_table.size()<= bwd_table.size())? fwd_table: bwd_table;
      auto large_table = (fwd_table.size()>bwd_table.size())? fwd_table: bwd_table;

      // size_t found_cnt = 0;
      auto map_f1 = [&](std::pair<vertex, vertex> a){
        vertex u = parlay::internal::get_key(a);
        vertex l = parlay::internal::get_val(a);
        if (large_table.contain(a)){
          parlay::write_max(&label[u], (TOP_BIT | (label_type)l), std::less<label_type>());
          // found_cnt += 1;
        }else{
          parlay::write_max(&label[u], (label_type)l, std::less<label_type>());
        }
      };
      small_table.map(map_f1);
      // printf("\n");
      // printf("found %ld \n", found_cnt);
      auto map_f2 = [&](std::pair<vertex, vertex> a){
        vertex u = parlay::internal::get_key(a);
        vertex l = parlay::internal::get_val(a);
        parlay::write_max(&label[u], (label_type)l, std::less<label_type>());
      };
      // n_remain = parlay::internal::count_if_index(n-start, [&](size_t i){return !(TOP_BIT & label[vertices[start+i]]);});
      // std::cout << "  n_remain: " << n_remain << std::endl;
      large_table.map(map_f2);
      scc_offset+=(end-start);
      start =end;
      step = ceil(step * beta);
    }

    auto pairs = parlay::delayed::tabulate(G.size(), [&] (vertex i) {
        return std::pair(label[i] & VAL_MASK, i);});
    auto lists = parlay::group_by_index(pairs, scc_offset);
    auto SCCs = parlay::pack(lists, parlay::delayed_map(lists, [=] (auto b) {return b.size()>0;}));
    std::cout << "number of lists: " << SCCs.size() << std::endl;
    return SCCs;
}