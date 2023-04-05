#include "EGraph.h"
#include "gapbs-agile/src/cc_sv.h"
#include "gapbs-agile/src/afforest.h"

EGraph::EGraph(gapbs::WGraph& support_, GraphManip& gp_, std::map<int, std::vector<MyEdge>>& trussgroups_, MyEdgeList& edges_): support(support_), gp(gp_), trussgroups(trussgroups_), edges(edges_) {}


std::vector<int> EGraph::out_neigh(int& e_index, int round)
{
	std::vector<int>indices;
		
	MyEdge& e = edges[e_index];
	int u = e.first;
	int v = e.second;
	int k = support.get_unsafe(u,v).w;

	std::vector<int>& temp = gp.intersectlist[e_index];
		
	for(int j = 0; j < temp.size(); j++)
	{
		int w = temp[j];

		MyEdge e1{std::min(u,w), std::max(u,w)};
		MyEdge e2{std::min(v,w), std::max(v,w)};

		int k1 = support.get_unsafe(u,w).w;
		int k2 = support.get_unsafe(v,w).w;

			
		if (k1 >= k && k2 >= k)	//ensuring k-triangle
		{
			if(k1 == k)
			{
				indices.push_back(gp.edge2index.at(e1));
			}
			if(k2 == k)
			{
				indices.push_back(gp.edge2index.at(e2));
			}
		}
	}
	return indices;
}

std::vector<int> EGraph::in_neigh(int& e_index)
{		
	return std::vector<int>();
}

bool EGraph::directed()
{
	return false;
}

void EGraph::initializeParentComp(gapbs::pvector<int>& comp)
{
	#pragma omp parallel for
	for(int i = 0; i < comp.size(); i++)
	{
		comp[i] = i;
	}
}

gapbs::pvector<int> EGraph::conn_comp()
{
	auto conn_start = std::chrono::high_resolution_clock::now();

	gapbs::pvector<int>comp(edges.size());

	initializeParentComp(comp);
	
	for (int k = gp.kmin; k <= gp.kmax; k++)
	{
		if(trussgroups.find(k) != trussgroups.end())
		{
			std::vector<MyEdge>& phi_k = trussgroups.at(k);
			//prepare the list of indices going to send to gapbs' ShiloachVishkin
			std::vector<int>indices(phi_k.size());
			#pragma omp parallel for
			for(int i = 0; i < phi_k.size(); i++)
			{
				auto& e = phi_k[i];
				int& index = gp.edge2index[e];
				indices[i] = index;
			}
			gapbs::ShiloachVishkin(*this, indices, comp);
		}
	}

	auto conn_end = std::chrono::high_resolution_clock::now();

	gapbs_sv_conn_time += std::chrono::duration_cast<std::chrono::nanoseconds>(conn_end - conn_start).count();

	return comp;
}


gapbs::pvector<int> EGraph::conn_comp_afforest()
{
	auto conn_start = std::chrono::high_resolution_clock::now();

	gapbs::pvector<int>comp_afforest(edges.size());

	initializeParentComp(comp_afforest);
	
	for (int k = gp.kmin; k <= gp.kmax; k++)
	{
		if(trussgroups.find(k) != trussgroups.end())
		{
			std::vector<MyEdge>& phi_k = trussgroups.at(k);
			//prepare the list of indices going to send to gapbs' ShiloachVishkin
			std::vector<int>indices(phi_k.size());
			#pragma omp parallel for
			for(int i = 0; i < phi_k.size(); i++)
			{
				auto& e = phi_k[i];
				int& index = gp.edge2index[e];
				indices[i] = index;
			}
			gapbs::Afforest(*this, indices, comp_afforest);
		}
	}

	auto conn_end = std::chrono::high_resolution_clock::now();

	gapbs_afforest_conn_time += std::chrono::duration_cast<std::chrono::nanoseconds>(conn_end - conn_start).count();

	return comp_afforest;
}

