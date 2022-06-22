#ifndef GAPBS_TC_H
#define GAPBS_TC_H

// Encourage use of gcc's parallel algorithms (for sort for relabeling)
#ifdef _OPENMP
#ifndef _GLIBCXX_PARALLEL
  #define _GLIBCXX_PARALLEL
#endif
#endif

#include <algorithm>
#include <iostream>
#include <cinttypes>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"

namespace gapbs {
namespace tc {

// heuristic to see if sufficently dense power-law graph
template <typename GraphT_>
bool WorthRelabelling(const GraphT_ &g) {
  int64_t average_degree = g.num_edges() / g.num_nodes();
  if (average_degree < 10)
    return false;
  SourcePicker<GraphT_> sp(g);
  int64_t num_samples = std::min(int64_t(1000), g.num_nodes());
  int64_t sample_total = 0;
  pvector<int64_t> samples(num_samples);
  for (int64_t trial=0; trial < num_samples; trial++) {
    samples[trial] = g.out_degree(sp.PickNext());
    sample_total += samples[trial];
  }
  std::sort(samples.begin(), samples.end());
  double sample_average = static_cast<double>(sample_total) / num_samples;
  double sample_median = samples[num_samples/2];
  return sample_average / 1.3 > sample_median;
}

template <typename GraphT_=Graph>
void PrintTriangleStats(const GraphT_ &g, size_t total_triangles) {
    std::cout << total_triangles << " triangles" << std::endl;
}

// Compares with simple serial implementation that uses std::set_intersection
template <typename GraphT_=Graph>
bool TCVerifier(const GraphT_ &g, size_t test_total) {
  #ifdef _OPENMP
  #ifdef _GLIBCXX_PARALLEL
    #undef _GLIBCXX_PARALLEL
  #endif
  #endif

  size_t total = 0;
  std::vector<NodeID> intersection;
  std::vector<NodeID> neigh_u, neigh_v;
  intersection.resize(g.num_nodes());
  for (NodeID u : g.vertices() ) {
    neigh_u.resize(g.out_degree(u));
    NodeID i = 0;
    for ( NodeID v : g.out_neigh(u) ) {
      neigh_u[i++] = v;
    }

    for (NodeID v : g.out_neigh(u)) {
        neigh_v.resize(g.out_degree(v));
        NodeID j = 0;
        for ( NodeID w : g.out_neigh(v) ) {
          neigh_v[j++] = w;
        }

        auto new_end = std::set_intersection(neigh_u.begin(),
                                             neigh_u.end(),
                                             neigh_v.begin(),
                                             neigh_v.end(),
                                             intersection.begin());
      total += std::distance(intersection.begin(), new_end);
    }
  }
  #ifdef _OPENMP
  #ifndef _GLIBCXX_PARALLEL
    #define _GLIBCXX_PARALLEL
  #endif
  #endif

  total = total / 6;  // each triangle was counted 6 times
  if (total != test_total)
      std::cout << total << " != " << test_total << std::endl;
  return total == test_total;
}
}
}

#endif