#ifndef EGRAPH_H
#define EGRAPH_H

#include "global.h"
#include "GraphManip.h"

class EGraph
{
	gapbs::WGraph& support;
	GraphManip& gp;
	std::map<int, std::vector<MyEdge>>& trussgroups;
	MyEdgeList& edges;

	public:
	
	EGraph(gapbs::WGraph& support_, GraphManip& gp_, std::map<int, std::vector<MyEdge>>& trussgroups_, MyEdgeList& edges_);
	void initializeParentComp(gapbs::pvector<int>& comp);
	gapbs::pvector<int> conn_comp();
	std::vector<int> out_neigh(int& e_index, int round = 0);
	gapbs::pvector<int> conn_comp_afforest();
	bool directed();
	std::vector<int> in_neigh(int& e_index);
};


#endif