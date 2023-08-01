#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "2hop_cover.h"
#include "helper/graph_utils.h"

// **************************************************************
// Driver
// **************************************************************
using vertex = uint32_t;
using distance = uint8_t;
using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
using graph = nested_seq;
using utils = graph_utils<vertex>;
using result = std::pair<sequence<sequence<vertex>>, sequence<sequence<distance>>>;

int main(int argc, char* argv[]) {
  auto usage = "Usage: 2hop_cover <n> || BFS_ligra <filename>";
  if (argc != 2) std::cout << usage << std::endl;
  else {
    long n = 0;
    graph G, GT;
    try { n = std::stol(argv[1]); }
    catch (...) {}
    if (n == 0) {
      // G = utils::read_symmetric_graph_from_file(argv[1]);
      // GT = G;
      G = utils::read_graph_from_bin(argv[1]);
      n = G.size();
    } else {
      G = utils::rmat_graph(n, 20*n);
      GT = utils::transpose(G);
    }
    utils::print_graph_stats(G);
    // result result;
    parlay::internal::timer t("Time");
    for (int i=0; i < 1; i++) {
      create_PrunedLandmarkLabeling<graph, vertex, distance, 0>(G);
      t.next("2hop_cover");
    }
    // long total = reduce(map(std::get<0>(result), parlay::size_of()));
    // std::cout << "Average label size: " << ((double) total/(double)n) << std::endl;
  }
}
