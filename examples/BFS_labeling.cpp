#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "BFS_labeling.h"
#include "helper/graph_utils.h"

using utils = graph_utils<vertex>;

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: construct_index GRAPH" << std::endl;
        exit(EXIT_FAILURE);
    }
    graph G = utils::read_graph_from_bin(argv[1]);
    // graph G = utils::read_symmetric_graph_from_file(argv[1]);
    vertex n = G.size();
    utils::print_graph_stats(G);
    PrunedLandmarkLabeling<0> pll;
    for (int i=0; i < 1; i++) {
      pll.ConstructIndex(G);
    }
}
