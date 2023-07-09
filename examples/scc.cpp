#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

using vertex = unsigned int;
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
      GT = utils::transpose(G);
      n = G.size();
    } else {
      G = utils::rmat_graph(n, 20*n);
      GT = utils::transpose(G);
    }
    utils::print_graph_stats(G);
    // nested_seq SCCs;
    parlay::internal::timer t("Time");
    for (int i=0; i < 1; i++) {
      find_scc<vertex, graph>(G, GT);
      t.next("SCC ");
    }

  }
}
