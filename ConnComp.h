#ifndef CONNCOMP_H
#define CONNCOMP_H

#include "global.h"
#include "GraphManip.h"

typedef std::vector<std::vector<std::pair<int, int>>> sm_graph_per_thread;

class ConnComp
{
	bool hooking;
	
	public:
	void edge_hooking(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, bool& hooking, GraphManip& gp);
	void edge_hooking(gapbs::WGraph& support, const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, bool& hooking, GraphManip& gp);
	void edge_compression(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, GraphManip& gp);
	void super_edge_creation(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, GraphManip& gp, std::vector<std::set<std::pair<int, int>>>& super_edges);
	void super_edge_creation(gapbs::WGraph& support, const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, GraphManip& gp, std::vector<std::set<std::pair<int, int>>>& super_edges);
	std::vector<std::set<std::pair<int, int>>> conn_comp_edge(const MyEdgeList& edges, std::map<int, std::vector<MyEdge>>& trussgroups, GraphManip& gp);
	std::vector<std::set<std::pair<int, int>>> conn_comp_edge(gapbs::WGraph& support, const MyEdgeList& edges, std::map<int, std::vector<MyEdge>>& trussgroups, GraphManip& gp);
	std::set<std::pair<int, int>> createSummaryGraph(std::vector<std::set<std::pair<int, int>>> super_edges);
	std::vector<std::pair<int, int>> mergeSummaryGraph(std::vector<std::set<std::pair<int, int>>>& super_edges);
};

#endif
