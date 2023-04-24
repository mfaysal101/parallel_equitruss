#include "ConnComp.h"
#include "MergeSort.h"


void ConnComp::edge_hooking(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, bool& hooking, GraphManip& gp)
{
	#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		auto e = phi_k[i];		//got the edge in the current k

		int index = gp.edge2index[e];

		std::vector<int>& temp = gp.intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e

		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			std::pair<int, int> e1{std::min(u,w), std::max(u,w)};
			std::pair<int, int> e2{std::min(v,w), std::max(v,w)};
			int k1, k2;

			k1 = gp.edge2k[e1];
			k2 = gp.edge2k[e2];
			
			if (k1 >= k && k2 >= k)	//ensuring k-triangle
			{
				int comp_u = gp.edge2P[e];
				int comp_v = gp.edge2P[e1];
				
				int high_comp = comp_u > comp_v ? comp_u : comp_v;
				int low_comp = comp_u + (comp_v - high_comp);

				if(low_comp < high_comp && high_comp == gp.edge2P[edges[high_comp]] && k == k1)
				{
					gp.edge2P[edges[high_comp]] = low_comp;
					hooking = true;
				}
			
				int comp_w = gp.edge2P[e2];

				high_comp = comp_u > comp_w ? comp_u : comp_w;
				low_comp = comp_u + (comp_w - high_comp);
 				
				if(low_comp < high_comp && high_comp == gp.edge2P[edges[high_comp]] && k == k2)
				{
					gp.edge2P[edges[high_comp]] = low_comp;
					hooking = true;
				}
			}
		}
	}
}


///This is edge_hooking function overload for gapbs integration

void ConnComp::edge_hooking(gapbs::WGraph& support, const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, bool& hooking, GraphManip& gp)
{
	#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		auto e = phi_k[i];		//got the edge in the current k

		int index = gp.edge2index[e];

		std::vector<int>& temp = gp.intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e

		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			std::pair<int, int> e1{std::min(u,w), std::max(u,w)};
			std::pair<int, int> e2{std::min(v,w), std::max(v,w)};

			int k1 = support.get_unsafe(u,w).w;
			int k2 = support.get_unsafe(v,w).w;
			
			if (k1 >= k && k2 >= k)	//ensuring k-triangle
			{
				int comp_u = gp.edge2P[e];
				int comp_v = gp.edge2P[e1];
				
				int high_comp = comp_u > comp_v ? comp_u : comp_v;
				int low_comp = comp_u + (comp_v - high_comp);

				if(low_comp < high_comp && high_comp == gp.edge2P[edges[high_comp]] && k == k1)
				{
					gp.edge2P[edges[high_comp]] = low_comp;
					hooking = true;
				}
			
				int comp_w = gp.edge2P[e2];

				high_comp = comp_u > comp_w ? comp_u : comp_w;
				low_comp = comp_u + (comp_w - high_comp);
 				
				if(low_comp < high_comp && high_comp == gp.edge2P[edges[high_comp]] && k == k2)
				{
					gp.edge2P[edges[high_comp]] = low_comp;
					hooking = true;
				}
			}
		}
	}
}


void ConnComp::edge_compression(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, GraphManip& gp)
{
	#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		auto e = phi_k[i];
		
		int p_id = gp.edge2P[e];				// gp.edge2P[e] will return me an index of an edge, using edges[p[e]] will return me the parent edge, applying p[edges[p[e]]] will return me the index of the parent edge
		int pp_id = gp.edge2P[edges[p_id]];

		while (pp_id != p_id)    
		{
			gp.edge2P[e] = pp_id;

			p_id = pp_id;
			pp_id = gp.edge2P[edges[p_id]];

		}
	}
}



void ConnComp::super_edge_creation(const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, GraphManip& gp, std::vector<std::set<std::pair<int, int>>>& super_edges)
{
	#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		size_t tid = omp_get_thread_num();

		auto e = phi_k[i];		//got the edge in the current k

		int index = gp.edge2index[e];

		std::vector<int>& temp = gp.intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e

		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			std::pair<int, int> e1{std::min(u,w), std::max(u,w)};
			std::pair<int, int> e2{std::min(v,w), std::max(v,w)};
			int k1, k2;

			k1 = gp.edge2k[e1];
			k2 = gp.edge2k[e2];

			int lowest_k = std::min(k1, k2);
			lowest_k = std::min(k, lowest_k);
				
			if(k > lowest_k && lowest_k == k1)
			{	
				int par1 = gp.edge2P[e1];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});								
			}

			if(k > lowest_k && lowest_k == k2)
			{
				int par1 = gp.edge2P[e2];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});
			}
		}
	}
}


///This is super_edge_creation function overload for gapbs integration

void ConnComp::super_edge_creation(gapbs::WGraph& support, const MyEdgeList& edges, std::vector<MyEdge>& phi_k, int& k, GraphManip& gp, std::vector<std::set<std::pair<int, int>>>& super_edges)
{
	#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		size_t tid = omp_get_thread_num();

		auto e = phi_k[i];		//got the edge in the current k

		int index = gp.edge2index[e];

		std::vector<int>& temp = gp.intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e

		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			std::pair<int, int> e1{std::min(u,w), std::max(u,w)};
			std::pair<int, int> e2{std::min(v,w), std::max(v,w)};

			int k1 = support.get_unsafe(u,w).w;
			int k2 = support.get_unsafe(v,w).w;

			int lowest_k = std::min(k1, k2);
			lowest_k = std::min(k, lowest_k);
				
			if(k > lowest_k && lowest_k == k1)
			{	
				int par1 = gp.edge2P[e1];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});								
			}

			if(k > lowest_k && lowest_k == k2)
			{
				int par1 = gp.edge2P[e2];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});
			}
		}
	}
}

///This is super_edge_creation function overload for gapbs integration with entire edge set instead of phi_k set only

void ConnComp::super_edge_creation(gapbs::WGraph& support, const MyEdgeList& edges, GraphManip& gp, std::vector<std::set<std::pair<int, int>>>& super_edges)
{

#pragma omp parallel for
	for (size_t i = 0; i < edges.size(); i++)
	{
		size_t tid = omp_get_thread_num();

		auto e = edges[i];

		int index = gp.edge2index[e];

		std::vector<int>& temp = gp.intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e
		
		for (int j = 0; j < temp.size(); j++)
		{
			
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			std::pair<int, int> e1{std::min(u,w), std::max(u,w)};
			std::pair<int, int> e2{std::min(v,w), std::max(v,w)};

			int k = support.get_unsafe(u,v).w;
			int k1 = support.get_unsafe(u,w).w;
			int k2 = support.get_unsafe(v,w).w;

			int lowest_k = std::min(k1, k2);
			lowest_k = std::min(k, lowest_k);

			if(k > lowest_k && lowest_k == k1)
			{	
				int par1 = gp.edge2P[e1];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});								
			}

			if(k > lowest_k && lowest_k == k2)
			{
				int par1 = gp.edge2P[e2];
				int par2 = gp.edge2P[e];
				if(par1 > par2)
				{
					int temp = par1;
					par1 = par2;
					par2 = temp;
				}
				super_edges[tid].insert({par1, par2});
			}
		}
	}
}


std::vector<std::set<std::pair<int, int>>> ConnComp::conn_comp_edge(const MyEdgeList& edges, std::map<int, std::vector<MyEdge>>& trussgroups, GraphManip& gp)
{
	size_t numThreads = 1;
	
	std::vector<std::set<std::pair<int, int>>> super_edges;
	
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
		#pragma omp single
		{
			super_edges.resize(numThreads);
		}
	}

	int it_count = 0;
	
	for (int i = gp.kmin; i <= gp.kmax; i++)
	{
		if(trussgroups.find(i) != trussgroups.end())
		{
			int it_count = 0;
			hooking = true;
			while (hooking)
			{
				it_count++;
				hooking = false;
				auto hook_start = std::chrono::high_resolution_clock::now();
				edge_hooking(edges, trussgroups[i], i, hooking, gp);
				auto hook_end = std::chrono::high_resolution_clock::now();
				hooking_time += std::chrono::duration_cast<std::chrono::nanoseconds>(hook_end - hook_start).count();

				auto comp_start = std::chrono::high_resolution_clock::now();
				edge_compression(edges, trussgroups[i], gp);
				auto comp_end = std::chrono::high_resolution_clock::now();
				compression_time += std::chrono::duration_cast<std::chrono::nanoseconds>(comp_end - comp_start).count();
			}
			it_count++;	//this extra increment is to ensure the one extra iteration for super-edge creation (algorithm 5.5) is counted
			iters_per_k.insert({i, it_count});
			auto se_start = std::chrono::high_resolution_clock::now();
			super_edge_creation(edges, trussgroups[i], i, gp, super_edges);
			auto se_end = std::chrono::high_resolution_clock::now();
			sp_edge_time += std::chrono::duration_cast<std::chrono::nanoseconds>(se_end - se_start).count();	
		}
	}
	
	return super_edges;
}


///this is conn_comp_edge overloaded function with gapbs integration

std::vector<std::set<std::pair<int, int>>> ConnComp::conn_comp_edge(gapbs::WGraph& support, const MyEdgeList& edges, std::map<int, std::vector<MyEdge>>& trussgroups, GraphManip& gp)
{
	size_t numThreads = 1;
	
	std::vector<std::set<std::pair<int, int>>> super_edges;
	std::vector<std::set<std::pair<int, int>>> super_edges_whole;
	
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
		#pragma omp single
		{
			super_edges.resize(numThreads);
			super_edges_whole.resize(numThreads);
		}
	}

	int it_count = 0;
	
	for (int i = gp.kmin; i <= gp.kmax; i++)
	{
		if(trussgroups.find(i) != trussgroups.end())
		{
			int it_count = 0;
			hooking = true;
			while (hooking)
			{
				it_count++;
				hooking = false;
				auto hook_start = std::chrono::high_resolution_clock::now();
				edge_hooking(support, edges, trussgroups[i], i, hooking, gp);
				auto hook_end = std::chrono::high_resolution_clock::now();
				hooking_time += std::chrono::duration_cast<std::chrono::nanoseconds>(hook_end - hook_start).count();

				auto comp_start = std::chrono::high_resolution_clock::now();
				edge_compression(edges, trussgroups[i], gp);
				auto comp_end = std::chrono::high_resolution_clock::now();
				compression_time += std::chrono::duration_cast<std::chrono::nanoseconds>(comp_end - comp_start).count();
			}
			it_count++;	//this extra increment is to ensure the one extra iteration for super-edge creation (algorithm 5.5) is counted
			iters_per_k.insert({i, it_count});
			auto se_start = std::chrono::high_resolution_clock::now();
			super_edge_creation(support, edges, trussgroups[i], i, gp, super_edges);
			auto se_end = std::chrono::high_resolution_clock::now();
			sp_edge_time += std::chrono::duration_cast<std::chrono::nanoseconds>(se_end - se_start).count();	
		}
	}
	
	auto se_whole_start = std::chrono::high_resolution_clock::now();
	super_edge_creation(support, edges, gp, super_edges_whole);
	auto se_whole_end = std::chrono::high_resolution_clock::now();
	sp_edge_whole_time = std::chrono::duration_cast<std::chrono::nanoseconds>(se_whole_end - se_whole_start).count();
	return super_edges;

}


//this function builds the summary graph by using multiple threads in parallel
std::vector<std::pair<int, int>> ConnComp::mergeSummaryGraph(std::vector<std::set<std::pair<int, int>>>& super_edges)
{
	auto parallel_summary_start = std::chrono::high_resolution_clock::now();
	
	double t1, t2, t3, t4;
	t1 = 0.0;
	t2 = 0.0;
	t3 = 0.0;
	t4 = 0.0;
	int hash_scale = 107;
	std::vector<std::pair<int, int>> final_sp_graph;
	std::vector<sm_graph_per_thread> sm_graph;
	size_t numThreads = 1;
	size_t total_num_sp_edges = 0;
	std::vector<size_t> num_thread_local_SpEdges;
	std::vector<size_t> displacments;
	sm_graph_per_thread combined_thread_local_sm_g;
	
	auto start_t1 = std::chrono::high_resolution_clock::now();
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
		#pragma omp single
		{
			sm_graph.resize(numThreads);
			num_thread_local_SpEdges.resize(numThreads);
			displacments.resize(numThreads);
			combined_thread_local_sm_g.resize(numThreads, std::vector<std::pair<int, int>>());
		}

		int tid = omp_get_thread_num();
		
		sm_graph_per_thread		sm_graph_t(numThreads, std::vector<std::pair<int, int>>());		//allocating rows for total number of threads
		
		for(auto it = super_edges[tid].begin(); it != super_edges[tid].end(); it++)
		{
			int x = (*it).first;
			int y = (*it).second;
			int dest_t = (((x + y) * (x + y + 1) / 2) + y) & (numThreads - 1);
			sm_graph_t[dest_t].push_back((*it));
		}
		sm_graph[tid] = sm_graph_t;
		
		std::set<std::pair<int, int>>().swap(super_edges[tid]);
	}
	
	std::vector<std::set<std::pair<int, int>>>().swap(super_edges);
	
	auto end_t1 = std::chrono::high_resolution_clock::now();
	t1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t1 - start_t1).count();
	
	auto start_t2 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel reduction(+:total_num_sp_edges)
	{
		int tid = omp_get_thread_num();
		
		for(auto& sm_t: sm_graph)
		{
			combined_thread_local_sm_g[tid].insert(end(combined_thread_local_sm_g[tid]), sm_t[tid].begin(), sm_t[tid].end());
			std::vector<std::pair<int,int>>().swap(sm_t[tid]);
		}
		
		sort(combined_thread_local_sm_g[tid].begin(), combined_thread_local_sm_g[tid].end(),
		[](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});
		
		combined_thread_local_sm_g[tid].erase(unique(combined_thread_local_sm_g[tid].begin(), combined_thread_local_sm_g[tid].end()), combined_thread_local_sm_g[tid].end());
		
		num_thread_local_SpEdges[tid] = combined_thread_local_sm_g[tid].size();
		
		total_num_sp_edges += num_thread_local_SpEdges[tid];
		
	}
	
	//printf("Total num of sp edges:%lu\n", total_num_sp_edges);
	auto end_t2 = std::chrono::high_resolution_clock::now();
	t2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t2 - start_t2).count();
	
	displacments[0] = 0;
		
	for(size_t i = 1; i < numThreads; i++)
	{
		displacments[i] = displacments[i - 1] + num_thread_local_SpEdges[i - 1];
	}
	
	final_sp_graph.resize(total_num_sp_edges);

	auto start_t3 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		
		//printf("threads id:%d, number of local sp edge:%lu\n", tid, combined_thread_local_sm_g[tid].size());
		
		for(int i = 0; i < combined_thread_local_sm_g[tid].size(); i++)
		{
			final_sp_graph[displacments[tid] + i] = combined_thread_local_sm_g[tid][i];
		}
		std::vector<std::pair<int, int>>().swap(combined_thread_local_sm_g[tid]);
	}
	
	sm_graph_per_thread().swap(combined_thread_local_sm_g);
	
	auto end_t3 = std::chrono::high_resolution_clock::now();
	t3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t3 - start_t3).count();
	
	auto start_t4 = std::chrono::high_resolution_clock::now();
	
	/*
	sort(final_sp_graph.begin(), final_sp_graph.end(), [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});
	*/
	MergeSort(final_sp_graph.begin(), final_sp_graph.end(), [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});

	auto end_t4 = std::chrono::high_resolution_clock::now();
	t4 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t4 - start_t4).count();
	
	auto parallel_summary_end = std::chrono::high_resolution_clock::now();
	parallel_summary_graph_time = std::chrono::duration_cast<std::chrono::nanoseconds>(parallel_summary_end - parallel_summary_start).count();
	
	printf("t1:%0.9f, t2:%0.9f, t3:%0.9f, t4:%0.9f\n", t1*(1e-9), t2*(1e-9), t3*(1e-9), t4*(1e-9));
	
	return final_sp_graph;	
}


std::set<std::pair<int, int>> ConnComp::createSummaryGraph(std::vector<std::set<std::pair<int, int>>> super_edges)
{
	
	auto summary_start = std::chrono::high_resolution_clock::now();
	
	std::set<std::pair<int, int>> summary_graph;
	int total = 0;
	int validate = 0;

	for (int i = 0; i < super_edges.size(); i++)
	{
		std::set<std::pair<int, int>> temp = super_edges[i];
		total += temp.size();
		for (auto it = temp.begin(); it != temp.end(); it++)
		{
			summary_graph.insert((*it));
		}
	}
	auto summary_end = std::chrono::high_resolution_clock::now();
	summary_merge_time = std::chrono::duration_cast<std::chrono::nanoseconds>(summary_end - summary_start).count();
	return summary_graph;
}