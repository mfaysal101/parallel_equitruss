#ifndef AFFOREST_H
#define AFFOREST_H

#include "platform_atomics.h"
#include "pvector.h"
#include <omp.h>


#include <iostream>
#include <random>
#include <unordered_map>
#include <iterator>

namespace gapbs {

// Place nodes u and v in same component of lower component ID
template <typename NodeT_>
void Link(NodeT_ u, NodeT_ v, pvector<NodeT_>& comp) {
  NodeT_ p1 = comp[u];
  NodeT_ p2 = comp[v];
  while (p1 != p2) {
    NodeT_ high = p1 > p2 ? p1 : p2;
    NodeT_ low = p1 + (p2 - high);
    NodeT_ p_high = comp[high];
    // Was already 'low' or succeeded in writing 'low'
    if ((p_high == low) ||
        (p_high == high && compare_and_swap(comp[high], high, low)))
      break;
    p1 = comp[comp[high]];
    p2 = comp[low];
  }
}


// Reduce depth of tree for each component to 1 by crawling up parents
template <typename GraphT_, typename NodeT_, typename RangeT_>
void Compress(const GraphT_ &g, RangeT_ &&r, gapbs::pvector<NodeT_>& comp) {
  #pragma omp parallel for
  for (auto n_iter=r.begin(); n_iter < r.end(); ++n_iter) {
    auto& n = *n_iter;
    while (comp[n] != comp[comp[n]]) {
      comp[n] = comp[comp[n]];
    }
  }
}

template <typename NodeT_, typename RangeT_>
NodeT_ SampleFrequentElement(const gapbs::pvector<NodeT_>& comp,
                             RangeT_&& r,
                             int64_t num_samples = 1024) {
  std::unordered_map<NodeT_, int> sample_counts(32);
  using kvp_type = typename std::unordered_map<NodeT_, int>::value_type;
  // Sample elements from 'comp'
  std::mt19937 gen;
  NodeT_ dist = r.end() - r.begin();
  std::uniform_int_distribution<NodeT_> distribution(0, dist-1);
  for (NodeT_ i = 0; i < num_samples; i++) {
    NodeT_ howfar = distribution(gen);
    NodeT_ n = *( r.begin() + howfar);
    sample_counts[comp[n]]++;
  }
  // Find most frequent element in samples (estimate of most frequent overall)
  auto most_frequent = std::max_element(
    sample_counts.begin(), sample_counts.end(),
    [](const kvp_type& a, const kvp_type& b) { return a.second < b.second; });
  float frac_of_graph = static_cast<float>(most_frequent->second) / num_samples;
  std::cout
    << "Skipping largest intermediate component (ID: " << most_frequent->first
    << ", approx. " << static_cast<int>(frac_of_graph * 100)
    << "% of the graph)" << std::endl;
  return most_frequent->first;
}

template <typename GraphT_, typename NodeT_, typename RangeT_>
void Afforest(GraphT_ &g, RangeT_ &&range, gapbs::pvector<NodeT_>& comp,
              int32_t neighbor_rounds = 2) {
  // Initialize each node to a single-node self-pointing tree
  #pragma omp parallel for
  for (auto n_iter=range.begin(); n_iter < range.end(); ++n_iter) {
    auto& n = *n_iter;
    comp[n] = n;
  }

  // Process a sparse sampled subgraph first for approximating components.
  // Sample by processing a fixed number of neighbors for each node (see paper)
  for (int r = 0; r < neighbor_rounds; ++r) {
  #pragma omp parallel for
    for (auto u_iter=range.begin(); u_iter < range.end(); ++u_iter) {
    auto& u = *u_iter;
      for (NodeT_ v : g.out_neigh(u, r)) {
        // Link at most one time if neighbor available at offset r
        Link(u, v, comp);
        break;
      }
    }
    Compress(g, range, comp);
  }

  // Sample 'comp' to find the most frequent element -- due to prior
  // compression, this value represents the largest intermediate component
  NodeT_ c = SampleFrequentElement(comp, range);

  // Final 'link' phase over remaining edges (excluding largest component)
  if (!g.directed()) {
    #pragma omp parallel for
    for (auto u_iter=range.begin(); u_iter < range.end(); ++u_iter) {
      auto& u = *u_iter;
      // Skip processing nodes in the largest component
      if (comp[u] == c)
        continue;
      // Skip over part of neighborhood (determined by neighbor_rounds)
      for (NodeT_ v : g.out_neigh(u, neighbor_rounds)) {
        Link(u, v, comp);
      }
    }
  } else {
    #pragma omp parallel for schedule(dynamic, 16384)
    for (auto u_iter=range.begin(); u_iter < range.end(); ++u_iter) {
      auto& u = *u_iter;
      if (comp[u] == c)
        continue;
      for (NodeT_ v : g.out_neigh(u, neighbor_rounds)) {
        Link(u, v, comp);
      }
      // To support directed graphs, process reverse graph completely
      for (NodeT_ v : g.in_neigh(u)) {
        Link(u, v, comp);
      }
    }
  }
  // Finally, 'compress' for final convergence
  Compress(g, range, comp);
}
}
#endif