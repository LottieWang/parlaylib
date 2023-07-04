#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

using vertex = int;
#include "scc.h"
#include "helper/graph_utils.h"

// **************************************************************
// Driver
// **************************************************************
using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
using graph = nested_seq;
using utils = graph_utils<vertex>;

int main(int argc, char* argv[]) {
  auto usage = "Usage: scc <n> || scc <filename>";
  if (argc != 2) std::cout << usage << std::endl;
  else {
    long n = 0;
    graph G, GT;
    try { n = std::stol(argv[1]); }
    catch (...) {}
    if (n == 0) {
      G = utils::read_graph_from_file(argv[1]);
      utils::print_graph_stats(G);
      GT = utils::transpose(G);
      n = G.size();
    } else {
      G = utils::rmat_graph(n, 20*n);
      GT = utils::transpose(G);
    }
    utils::print_graph_stats(G);
    parlay::sequence<vertex> label;
    parlay::internal::timer t("Time");
    for (int i=0; i < 5; i++) {
      label = find_scc(G, GT);
      t.next("SCC ");
    }

    std::cout << "num vertices visited: " << filter(label, [&](size_t i){return label[i]!=max_label;}).size() << std::endl;
  }
}
