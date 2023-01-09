#include <stdio.h>
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
#include<climits>

using namespace std;




typedef std::vector<std::pair<int, int> > EdgeList;
typedef std::vector<std::vector<int> > AdjList;


AdjList adjlist;

std::vector<int>parent;
std::map<std::pair<int, int>, int> p;
int numEdges;
int numVertices;
bool hooking;
size_t maxk;
std::vector<std::vector<int>> intersectlist;
std::map<pair<int, int>, int>trussk;
std::map<int, vector<pair<int, int>>>trussgroups;
std::map<pair<int, int>, int> edge2index;
//std::vector<std::set<pair<pair<int, int>, pair<int, int>>>> super_edges;
std::vector<std::set<pair<int, int>>> super_edges;
std::vector<std::pair<int, int>> summary_graph;
int kmin, kmax;
double hooking_time, compression_time, sp_edge_time, summary_merge_time;
std::map<int, int>iters_per_k;

namespace
{
	template <class InputIterator1, class InputIterator2>
	vector<int> IntersectionSize(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, int n)
	{
		vector<int> result;
		while (first1 != last1 && first2 != last2)
		{
			if (*first1 >= n || *first2 >= n)
			{
				break;
			}
			if (*first1 < *first2)
			{
				++first1;
			}
			else if (*first1 > *first2)
			{
				++first2;
			}
			else
			{
				result.push_back(*first1);
				++first1;
				++first2;
			}
		}
		return result;
	}
}


size_t getNumVertices(const EdgeList& edges)
{
	int num = 0;

#pragma omp parallel for reduction (max:num)
	for (size_t i = 0; i < edges.size(); i++)
	{
		num = max(num, 1 + max(edges[i].first, edges[i].second));
	}

	numVertices = num;

	cout << "Total number of vertices:" << numVertices << endl;

	return numVertices;
}

EdgeList ReadEdgeListFromFile(const char* filename)
{
	EdgeList edgelist;

	cout << "Reading network " << filename << " file\n";

	ifstream infile(filename);

	if (infile)
	{
		printf("File open successful\n");
	}
	else
	{
		printf("Failed to read file\n");
		exit(1);
	}

	string line = "";
	numEdges = 0;

	int index = 0;

	while (getline(infile, line))
	{
		istringstream iss(line);
		int src, dst;
		if ((iss >> src >> dst))
		{
			if(src == dst)
			{
				continue;
			}
			edgelist.push_back(make_pair(src, dst));
		}
	}

	numEdges = edgelist.size();

	cout << "Total number of edges:" << numEdges << endl;

	sort(edgelist.begin(), edgelist.end(), [](const pair<int, int>& edge1, const pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});

	for(auto edge:edgelist)
	{
		edge2index[edge] = index++;
	}
	return edgelist;

}

void EdgeToAdjList(const EdgeList& edges)
{
	adjlist.resize(getNumVertices(edges));

	for (auto edge : edges)
	{
		adjlist[edge.first].push_back(edge.second);
		adjlist[edge.second].push_back(edge.first);
	}

}

void printParent(ofstream& out)
{
	for (int i = 0; i < numVertices; i++)
	{
		out << i << "\t" << parent[i] << endl;
	}
	out.close();
}

void initializeParent()
{
	parent.resize(numVertices);

	for (int i = 0; i < numVertices; i++)
	{
		parent[i] = i;
	}
}

void parallel_hooking()
{
#pragma omp parallel for
	for (int u = 0; u < numVertices; u++)
	{
		vector<int>& temp = adjlist[u];

#pragma omp parallel for
		for (int j = 0; j < temp.size(); j++)
		{
			int v = temp[j];

			if (parent[u] < parent[v] && parent[v] == parent[parent[v]])
			{
				parent[parent[v]] = parent[u];
				hooking = true;
			}
		}
	}
}

void serial_hooking()
{
	for (int u = 0; u < numVertices; u++)
	{
		vector<int>& temp = adjlist[u];

		for (int j = 0; j < temp.size(); j++)
		{
			int v = temp[j];

			if (parent[u] < parent[v] && parent[v] == parent[parent[v]])
			{
				parent[parent[v]] = parent[u];
				hooking = true;
			}
		}
	}
}

void parallel_compression()
{
#pragma omp parallel for
	for (int v = 0; v < numVertices; v++)
	{
		while (parent[parent[v]] != parent[v])
		{
			parent[v] = parent[parent[v]];
		}
	}
}

void serial_compression()
{
	for (int v = 0; v < numVertices; v++)
	{
		while (parent[parent[v]] != parent[v])
		{
			parent[v] = parent[parent[v]];
		}
	}
}

void shiloach_vishkin_parallel()
{
	hooking = true;

	while (hooking)
	{
		hooking = false;
		parallel_hooking();
		parallel_compression();
	}
}

void shiloach_vishkin_serial()
{
	hooking = true;

	while (hooking)
	{
		hooking = false;
		serial_hooking();
		serial_compression();
	}
}

void populateIntersectList(const EdgeList& edges)
{
	const int n = adjlist.size();

	maxk = 0;

	for (int i = 0; i < edges.size(); i++)
	{
		int u = edges[i].first;
		int v = edges[i].second;

		intersectlist[i] = IntersectionSize(adjlist[u].begin(), adjlist[u].end(), adjlist[v].begin(), adjlist[v].end(), n);

		maxk = max(maxk, intersectlist[i].size());
	}
}

void bucketSortedEdgelist(int kmax, const EdgeList& edges, vector<vector<int>>& sp, vector<pair<int, int>>& sorted_elbys, map<int, int>& svp, map<pair<int, int>, int>& sorted_ep)
{
	vector<int> bucket((kmax + 1), 0);

	for (auto it = sp.begin(); it != sp.end(); it++)
	{
		bucket[(*it).size()]++;
	}

	int temp;

	int pt = 0;

	for (int i = 0; i < kmax + 1; i++)
	{
		temp = bucket[i];
		bucket[i] = pt;
		pt = pt + temp;
	}

	for (int i = 0; i < sp.size(); i++)
	{
		sorted_elbys[bucket[sp[i].size()]] =  edges[i];
		sorted_ep.insert(make_pair(edges[i], bucket[sp[i].size()]));
		if (svp.find(sp[i].size()) == svp.end())
		{
			svp.insert(make_pair(sp[i].size(), bucket[sp[i].size()]));
		}
		bucket[sp[i].size()] = bucket[sp[i].size()] + 1;
	}

}


void computeTruss(const EdgeList& edges)
{
	map<int, pair<int, int>> klistdict;
	vector<pair<int, int>> sorted_elbys(p.size());
	map<pair<int, int>, int> sorted_ep;
	map<int, int> svp;

	bucketSortedEdgelist(maxk, edges, intersectlist, sorted_elbys, svp, sorted_ep);

}


void edge_hooking(const EdgeList& edges, vector<pair<int, int>>& phi_k, int& k)
{

#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{

		auto e = phi_k[i];		//got the edge in the current k

		int index = edge2index[e];

		//printf("k:%d, (u,v):(%d,%d)\n", k, e.first, e.second);

		vector<int>& temp = intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e

//#pragma omp parallel for
		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			pair<int, int> e1;
			pair<int, int> e2;
			int k1, k2;

			auto it1 = trussk.find({u, w});
			if(it1 != trussk.end())
			{
				e1.first = u;
				e1.second = w;
				k1 = it1->second;
			}
			else
			{
				e1.first = w;
				e1.second = u;
				k1 = trussk[e1];
			}
			
			auto it2 = trussk.find({v, w});
			if(it2 != trussk.end())
			{
				e2.first = v;
				e2.second = w;
				k2 = it2->second;
			}
			else
			{
				e2.first = w;
				e2.second = v;
				k2 = trussk[e2];
			}

			/*
			if (k1 >= k && k2 >= k)	//ensuring k-triangle
			{
				if (p[e] < p[e1] && p[e1] == p[edges[p[e1]]] && k == k1)
				{
					p[edges[p[e1]]] = p[e];
					hooking = true;
				}

				if (p[e] < p[e2] && p[e2] == p[edges[p[e2]]] && k == k2)
				{
					p[edges[p[e2]]] = p[e];
					hooking = true;
				}
			}
			*/

			if (k1 >= k && k2 >= k)	//ensuring k-triangle
			{
				int comp_u = p[e];
				int comp_v = p[e1];
				
				int high_comp = comp_u > comp_v ? comp_u : comp_v;
				int low_comp = comp_u + (comp_v - high_comp);

				if(low_comp < high_comp && high_comp == p[edges[high_comp]] && k == k1)
				{
					p[edges[high_comp]] = low_comp;
					hooking = true;
				}
			
				int comp_w = p[e2];

				high_comp = comp_u > comp_w ? comp_u : comp_w;
				low_comp = comp_u + (comp_w - high_comp);
 				
				if(low_comp < high_comp && high_comp == p[edges[high_comp]] && k == k2)
				{
					p[edges[high_comp]] = low_comp;
					hooking = true;
				}
			}
		}
	}
}


void edge_compression(const EdgeList& edges, vector<pair<int, int>>& phi_k)
{
#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		auto e = phi_k[i];

		//printf("(u,v):(%d,%d), p[(u,v)]:%d, p[edges[p[e]]]:%d)\n", e.first, e.second, p[e], p[edges[p[e]]]);
		
		int p_id = p[e];
		int pp_id = p[edges[p_id]];

		while (pp_id != p_id)    // p[e] will return me an index of an edge, using edges[p[e]] will return me the parent edge, applying p[edges[p[e]]] will return me the index of the parent edge
		{
			p[e] = pp_id;

			p_id = pp_id;
			pp_id = p[edges[p_id]];

		}
	}
}


void super_edge_creation(const EdgeList& edges, vector<pair<int, int>>& phi_k, int& k)
{

#pragma omp parallel for
	for (int i = 0; i < phi_k.size(); i++)
	{
		size_t tid = omp_get_thread_num();

		auto e = phi_k[i];		//got the edge in the current k

		int index = edge2index[e];

		vector<int>& temp = intersectlist[index];		//got the corresponding intersecting edges that makes triangle with e
		
		//printf("Shiloach-Vishkin parallel inside sp_edge_creation:(%d,%d), numtriangle:%d\n", e.first, e.second, temp.size());

//#pragma omp parallel for
		for (int j = 0; j < temp.size(); j++)
		{
			int w = temp[j];
			int u = e.first;
			int v = e.second;

			pair<int, int> e1;
			pair<int, int> e2;
			int k1, k2;
			
			auto it1 = trussk.find({u, w});
			if(it1 != trussk.end())
			{
				e1.first = u;
				e1.second = w;
				k1 = it1->second;
			}
			else
			{
				e1.first = w;
				e1.second = u;
				k1 = trussk[e1];
			}
			
			auto it2 = trussk.find({v, w});
			if(it2 != trussk.end())
			{
				e2.first = v;
				e2.second = w;
				k2 = it2->second;
			}
			else
			{
				e2.first = w;
				e2.second = v;
				k2 = trussk[e2];
			}

			int lowest_k = min(k1, k2);
			lowest_k = min(k, lowest_k);
				
			if(k > lowest_k && lowest_k == k1)
			{
				super_edges[tid].insert({p[e1], p[e]});				
			}

			if(k > lowest_k && lowest_k == k2)
			{
				super_edges[tid].insert({p[e2], p[e]});
			}
		}
	}
}


void conn_comp_edge(const EdgeList& edges)
{

	size_t numThreads = 1;
	#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
		#pragma omp single
		{
			super_edges.resize(numThreads);
		}
	}

	int it_count = 0;

	for (int i = kmin; i <= kmax; i++)
	{
		int it_count = 0;
		hooking = true;
		//printf("Shiloach-Vishkin parallel, k:%d\n", i);
		while (hooking)
		{
			it_count++;
			hooking = false;
			auto hook_start = std::chrono::high_resolution_clock::now();
			edge_hooking(edges, trussgroups[i], i);
			auto hook_end = std::chrono::high_resolution_clock::now();
			hooking_time += std::chrono::duration_cast<std::chrono::nanoseconds>(hook_end - hook_start).count();

			auto comp_start = std::chrono::high_resolution_clock::now();
			edge_compression(edges, trussgroups[i]);
			auto comp_end = std::chrono::high_resolution_clock::now();
			compression_time += std::chrono::duration_cast<std::chrono::nanoseconds>(comp_end - comp_start).count();
		}
		it_count++;	//this extra increment is to ensure the one extra iteration for super-edge creation (algorithm 5.5) is counted
		iters_per_k.insert({i, it_count});
		auto se_start = std::chrono::high_resolution_clock::now();
		super_edge_creation(edges, trussgroups[i], i);
		auto se_end = std::chrono::high_resolution_clock::now();
		sp_edge_time += std::chrono::duration_cast<std::chrono::nanoseconds>(se_end - se_start).count();
		
	}
}

void readTruss(ifstream& in)
{
	string line = "";

	if (in.is_open())
	{
		while (getline(in, line))
		{
			istringstream iss(line);
			string token;
			vector<int>tuple;

			while (getline(iss, token, ','))
			{
				tuple.push_back(stoi(token));
			}

			pair<int, int> temp;
			temp.first = tuple[0];
			temp.second = tuple[1];
			kmin = min(kmin, tuple[2]);
			kmax = max(kmax, tuple[2]);
			trussk.insert(make_pair(temp, tuple[2]));
			trussgroups[tuple[2]].push_back(temp);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		trussgroups.erase(i);
	}
	kmin = max(kmin, 3);

	in.close();
}


void sortEachAdjList()
{
	int len = adjlist.size();

#pragma omp parallel for
	for (int i = 0; i < len; i++)
	{
		sort(adjlist[i].begin(), adjlist[i].end());
	}
}

void printEdgeParent(ofstream& out)
{
	for (auto it = p.begin(); it != p.end(); it++)
	{
		out << it->first.first << "\t" << it->first.second << "\t" << it->second << endl;
	}
	out.close();
}

void printSummaryGraph(ofstream& out)
{
	for (auto sp_edge : summary_graph)
	{
		out << sp_edge.first << "\t" << sp_edge.second << endl;
	}
	out.close();
}

void initializeEdgeParent(const EdgeList& edges)
{
	for (int i = 0; i < edges.size(); i++)
	{
		p[edges[i]] = i;
	}
}

void createSummaryGraph(const EdgeList& edges)
{
	
	auto summary_start = std::chrono::high_resolution_clock::now();
	set<pair<int, int>> visited;

	int total = 0;

	for (int i = 0; i < super_edges.size(); i++)
	{
		set<pair<int, int>> temp = super_edges[i];
		printf("superedge vector:%d, size:%lu\n", i, temp.size());
		total += temp.size();
		for (auto it = temp.begin(); it != temp.end(); it++)
		{

			int p1 = p[edges[(*it).first]];
			int p2 = p[edges[(*it).second]];
			
			//printf("%d\t%d\t%d\n", p1, p2, i);

			//printf("thread id:%d, (c1, c2):(%d, %d)\t(p1, p2):(%d, %d)\n", i, (*it).first, (*it).second, p1, p2);	

			if (!visited.count({ p1, p2 }))
			{
				summary_graph.push_back({ p1, p2 });
				visited.insert({ p1, p2 });
			}
		}
	}
	printf("total pair of edges:%d\n", total);
	auto summary_end = std::chrono::high_resolution_clock::now();
	summary_merge_time += std::chrono::duration_cast<std::chrono::nanoseconds>(summary_end - summary_start).count();
}

int main(int argc, char* argv[])
{
	string networkfile;
	string trussfile;

	networkfile = argv[1];
	trussfile = argv[2];

	kmin = INT_MAX;
	kmax = INT_MIN;

	EdgeList edgelist = ReadEdgeListFromFile(networkfile.c_str());

	//ofstream out1("out1.txt"), out2("out2.txt");
	
	ofstream out3("out3.txt"), out4("out4.txt");

	ifstream trussinput(trussfile.c_str());

	EdgeToAdjList(edgelist);

	intersectlist.resize(edgelist.size());

	/*initializeParent();

	auto s_start = std::chrono::high_resolution_clock::now();

	shiloach_vishkin_serial();

	auto s_end = std::chrono::high_resolution_clock::now();

	auto s_time = std::chrono::duration_cast<std::chrono::nanoseconds>(s_end - s_start).count();

	printParent(out1);

	cout << "\nnow the parallel\n\n";

	initializeParent();

	auto p_start = std::chrono::high_resolution_clock::now();

	shiloach_vishkin_parallel();

	auto p_end = std::chrono::high_resolution_clock::now();

	auto p_time = std::chrono::duration_cast<std::chrono::nanoseconds>(p_end - p_start).count();

	printParent(out2);

	printf("========serial_time:%0.9f===========\n", s_time * (1e-9));
	printf("========parallel_time:%0.9f===========\n", p_time * (1e-9));
	
	*/
	
	// This function initialize the parent component ID of each edge to itself
	initializeEdgeParent(edgelist);

	// Sorting the adjacency list for faster edge intersection
	sortEachAdjList();

	populateIntersectList(edgelist);

	// This function currently read the trussness of all the edges from
	// a file generated by serial code k-truss decomposition
	readTruss(trussinput);

	hooking_time = 0.0;
	compression_time = 0.0;
	sp_edge_time = 0.0;
	summary_merge_time = 0.0;

	auto start = std::chrono::high_resolution_clock::now();

	// This function is the Shiloach-Vishkin connected component over edges
	conn_comp_edge(edgelist);

	// This function is right now serial, need to design parallel implementation
	createSummaryGraph(edgelist);

	auto end = std::chrono::high_resolution_clock::now();

	double totalExecutionTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	printf("========totalExecutionTime:%0.9f===========\n", totalExecutionTime*(1e-9));
	printf("========Time for hooking phase:%0.9f===========\n", hooking_time*(1e-9));
	printf("========Time for compression phase:%0.9f===========\n", compression_time*(1e-9));
	printf("========Time for super edges creation phase:%0.9f===========\n", sp_edge_time*(1e-9));
	printf("========Time for merging to summary graph creation:%0.9f===========\n", summary_merge_time*(1e-9));

	int total_iters = 0;

	for(auto it = iters_per_k.begin(); it != iters_per_k.end(); it++)
	{
		total_iters += it->second;
		printf("k:%d, iterations:%d\n", it->first, it->second);
	}

	printf("total_iterations_for_hooking:%d\n", total_iters);

	printEdgeParent(out3);

	printSummaryGraph(out4);

	return 0;
}

class Edge
{
public:

	int s;
	int t;


	Edge()
	{

	}

	Edge(int source, int target)
	{
		s = source;
		t = target;
	}

	// "<" operator overloading required by c++ map for custom object type
	bool operator<(const Edge& ob) const
	{
		return s < ob.s || (s == ob.s && t < ob.t);
	}

};
