#include <iostream>
#include <string>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "BFS_ligra.h"
#include "helper/graph_utils.h"

// **************************************************************
// Driver
// **************************************************************
using vertex = int;
using nested_seq = parlay::sequence<parlay::sequence<vertex>>;
using graph = nested_seq;
using utils = graph_utils<vertex>;

int main(int argc, char* argv[]) {
  auto usage = "Usage: BFS_ligra <n> || BFS_ligra <filename>";
  if (argc != 2) std::cout << usage << std::endl;
  else {
    long n = 0;
    graph G, GT;
    try { n = std::stol(argv[1]); }
    catch (...) {}
    if (n == 0) {
      // G = utils::read_symmetric_graph_from_file(argv[1]);
      G = utils::read_graph_from_bin(argv[1]);
      GT = G;
      n = G.size();
    } else {
      G = utils::rmat_graph(n, 20*n);
      GT = utils::transpose(G);
    }
    utils::print_graph_stats(G);
    nested_seq result;
    parlay::internal::timer t("Time");
    parlay::sequence<vertex> sources = {14, 363, 1677, 1867, 530, 185, 230, 40, 132, 323, 529, 2776, 23, 114, 2137, 4971, 8499, 2510, 220, 30, 124, 1016, 140, 36, 1960, 231, 1961, 24, 2743, 150, 73, 2136, 717, 173, 1051, 563, 121, 0, 1530, 186, 292, 763, 157, 183, 131, 540, 489, 746, 174, 8491, 335, 4125, 27, 2902, 372, 1937, 286, 4953, 1139, 542, 4960, 70, 294, 362};
    // for (int i=0; i < 5; i++) {
    //   result = BFS(1, G, GT);
    //   t.next("BFS_ligra");
    // }
    auto distance = parlay::tabulate(n, 
        [&](vertex i){return parlay::tabulate(10, [&](vertex i){return std::numeric_limits<vertex>::max();});});
    for (size_t s = 0; s<10; s++){
      result = BFS(sources[s], G, GT);
      printf("result size %lu\n", result.size());
      for (size_t i = 0; i<result.size(); i++){
        for (size_t j =0; j<std::min(result[i].size(), (size_t)10); j++){
          printf("%d ", result[i][j]);
          distance[result[i][j]][s]=i;
        }
        printf("\n");
      }
    }
    // for (int i = 0; i<10; i++){
    //   printf("distances: ");
    //   for (int j = 0; j<64; i++){
    //     printf("%d ",distance[i][j]);
    //   }
    //   printf("\n");
    // }

    long visited = reduce(map(result, parlay::size_of()));
    std::cout << "num vertices visited: " << visited << std::endl;
  }
}
