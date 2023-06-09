#ifndef GLOBAL_H
#define GLOBAL_H

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <set>
#include <chrono>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <climits>
#include <deque>
#include "gapbs-agile/src/benchmark.h"


extern std::map<int, int>iters_per_k;

extern double hooking_time;
extern double compression_time;
extern double sp_edge_time;
extern double summary_merge_time;
extern double parallel_summary_graph_time;
extern double gapbs_sv_conn_time;
extern double gapbs_afforest_conn_time;
extern double serial_sorting_time;
extern double parallel_sorting_time;
extern double serial_splitter_time;
extern double parallel_splitter_time;
extern double sp_edge_whole_time;
extern double intersectlist_time;
extern double truss_read_time;
extern double graph_read_time;
extern double edge2adj_time;
extern double sorting_adj_time;
extern double initializeP_time;
extern double reMapping_compID_time;


#endif