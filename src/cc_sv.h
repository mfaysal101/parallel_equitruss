#ifndef CC_SV_H
#define CC_SV_H

#include "pvector.h"
#include <omp.h>

#include <iostream>

namespace gapbs {

// The hooking condition (comp_u < comp_v) may not coincide with the edge's
// direction, so we use a min-max swap such that lower component IDs propagate
// independent of the edge's direction.
template <class GraphT_, class NodeT_>
pvector<NodeT_> ShiloachVishkin(const GraphT_ &g) {
  pvector<NodeT_> comp(g.num_nodes());
  #pragma omp parallel for
  for (NodeT_ n=0; n < g.num_nodes(); n++)
    comp[n] = n;
  bool change = true;
  int num_iter = 0;
  while (change) {
    change = false;
    num_iter++;
    #pragma omp parallel for
    for (NodeT_ u=0; u < g.num_nodes(); u++) {
      for (NodeT_ v : g.out_neigh(u)) {
        NodeT_ comp_u = comp[u];
        NodeT_ comp_v = comp[v];
        if (comp_u == comp_v) continue;
        // Hooking condition so lower component ID wins independent of direction
        NodeT_ high_comp = comp_u > comp_v ? comp_u : comp_v;
        NodeT_ low_comp = comp_u + (comp_v - high_comp);
        if (high_comp == comp[high_comp]) {
          change = true;
          comp[high_comp] = low_comp;
        }
      }
    }
    #pragma omp parallel for
    for (NodeT_ n=0; n < g.num_nodes(); n++) {
      while (comp[n] != comp[comp[n]]) {
        comp[n] = comp[comp[n]];
      }
    }
  }
  std::cout << "Shiloach-Vishkin took " << num_iter << " iterations" << std::endl;
  return comp;
}
}

#endif