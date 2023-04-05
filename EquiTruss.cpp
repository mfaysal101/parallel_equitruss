#include "GraphManip.h"
#include "ConnComp.h"
#include "EGraph.h"
#include "MergeSort.h"

double hooking_time;
double compression_time;
double sp_edge_time;
double summary_merge_time;
double parallel_summary_graph_time;
double gapbs_sv_conn_time;
double gapbs_afforest_conn_time;
double serial_sorting_time;
double parallel_sorting_time;
std::map<int, int>iters_per_k;


class CLEquiTruss : public gapbs::CLApp 
{
    std::string _trussness_filename;
public:
    CLEquiTruss(int argc, char** argv, std::string name)
        : CLApp(argc, argv, name) {
        get_args_ += "t:";
        AddHelpLine('t', "file", "Input graph with edge trussness");
    }

    void HandleArg(signed char opt, char* opt_arg) override {
        switch (opt) {
        case 't':  _trussness_filename = std::string(opt_arg); break;
        default: CLApp::HandleArg(opt, opt_arg);
        }
    }

    std::string trussness_filename() const { return _trussness_filename; }
};

void printComp(std::ofstream& out, const MyEdgeList& edges, gapbs::pvector<int>& comp)
{
	for(int i = 0; i < comp.size(); i++)
	{
		out << edges[i].first << "\t" << edges[i].second << "\t" << comp[i] << std::endl;
	}
}

void printComp(std::ofstream& out, const MyEdgeList& edges, gapbs::pvector<int>& comp, gapbs::WGraph& support_)
{
	for(int i = 0; i < comp.size(); i++)
	{
		if(support_.get_unsafe(edges[i].first, edges[i].second).w > 2)
		{
			out << edges[i].first << "\t" << edges[i].second << "\t" << comp[i] << std::endl;	
		}
	}
}


void printEdgeParent(std::ofstream& out, GraphManip& gp)
{
	for (auto it = gp.edge2P.begin(); it != gp.edge2P.end(); it++)
	{
		out << it->first.first << "\t" << it->first.second << "\t" << it->second << std::endl;
	}
	out.close();
}

void printSummaryGraph(std::ofstream& out, std::set<std::pair<int, int>>& summary_graph)
{
	for (auto sp_edge : summary_graph)
	{
		out << sp_edge.first << "\t" << sp_edge.second << std::endl;
	}
	out.close();
}


void printSummaryGraph(std::ofstream& out, std::vector<std::pair<int, int>>& smG)
{
	for (auto sp_edge : smG)
	{
		out << sp_edge.first << "\t" << sp_edge.second << std::endl;
	}
	out.close();
}

std::vector<int> sortingSpNodeParallel(gapbs::pvector<int>& comp_afforest)
{
	auto start = std::chrono::high_resolution_clock::now();
	int N = comp_afforest.size();
	std::vector<int>SNData(N, 0);

	#pragma omp parallel for
	for(int i = 0; i < N; i++)
	{
		SNData[i] = i;
	}

	omp_set_dynamic(0);

	MergeSort(SNData.begin(), SNData.end(), [&](const int edge_id1, const int edge_id2){
		return (comp_afforest[edge_id1] < comp_afforest[edge_id2]);
	});
	auto end = std::chrono::high_resolution_clock::now();
	parallel_sorting_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	return SNData;   
}

std::vector<int> sortingSpNode(gapbs::pvector<int>& comp_afforest)
{	
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<int>SNData(comp_afforest.size(), 0);
	for(int i = 0; i < comp_afforest.size(); i++)
	{
		SNData[i] = i;
	}

	std::sort(SNData.begin(), SNData.end(), [&](const int edge_id1, const int edge_id2){
		return (comp_afforest[edge_id1] < comp_afforest[edge_id2]); 
	});
	auto end = std::chrono::high_resolution_clock::now();
	serial_sorting_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	return SNData;
}

void printSpNdClusters(std::ofstream& out, const std::vector<int>& sortedSNData, const gapbs::pvector<int>& comp_afforest)
{
	//out<<"length of the totaledges:"<<sortedSNData.size()<<std::endl;

	std::map<int, std::set<int>>sorted_indices;

	for(int i = 0; i < sortedSNData.size(); i++)
	{
		int index = sortedSNData[i];
		sorted_indices[comp_afforest[index]].insert(index);
		//out<<index<<"\t"<<comp_afforest[index]<<std::endl;
	}

	for(auto it = sorted_indices.begin(); it != sorted_indices.end(); it++)
	{
		auto comp_id = it->first;
		auto comp_edges = it->second;
		for(auto it2 = comp_edges.begin(); it2 != comp_edges.end(); it2++)
		{
			out<<comp_id<<"\t"<<(*it2)<<std::endl;
		}
	}
	out.close();
}

//This function renumbers the super node ids starting from 0. It converts the super node ids to splitter ids
std::map<size_t, size_t> RenameSNIDs(const std::vector<int>& sortedSNData, const gapbs::pvector<int>& comp_afforest)
{
	std::map<size_t, size_t>reMap_SNIDs;

	size_t splitter_id = 0;

	for(int i = 0; i < sortedSNData.size(); i++)
	{
		size_t index = sortedSNData[i];
		size_t snID = comp_afforest[index];
		if(reMap_SNIDs.find(snID) == reMap_SNIDs.end())
		{
			reMap_SNIDs[snID] = splitter_id++;
		}	
	}
	
	return reMap_SNIDs;
}

void printRenamedSNIDs(std::ofstream& out, std::map<size_t, size_t> reMap_SNIDs)
{
	for(auto it = reMap_SNIDs.begin(); it != reMap_SNIDs.end(); it++)
	{
		out<<it->first<<"\t"<<it->second<<std::endl;
	}
	out.close();
}

int main(int argc, char* argv[])
{
	std::string networkfile;
	std::string trussfile;
	

	CLEquiTruss cli(argc, argv, "Equitruss-only");
       	if (!cli.ParseArgs())
	{
		return -1;
	}
                
       	gapbs::WeightedBuilder b(cli);
       	gapbs::WGraph support = b.MakeGraph();
       	if ( support.directed() ) 
	{
            std::cout << "Input graph is directed but equitruss requires undirected" << std::endl;
       	}
       	support.PrintStats();


	networkfile = argv[1];
	trussfile = argv[2];

	std::string ktruss_input = cli.trussness_filename();

	GraphManip gp;

	//This is the old way (before gapbs) of reading edgelist from file, I am commenting out this one.
	//MyEdgeList edgelist = gp.ReadEdgeListFromFile(networkfile.c_str());
	//This is the new way of reading edgelist after gapbs integration
	MyEdgeList edgelist = gp.BuildEdgeListFrom(support);
	
	std::ofstream out3("spNodes.txt"), out4("spEdges.txt"), out5("out5.txt"), out6("gapSV.txt"), out7("gapAfforest.txt"), out8("svNdTrimmed.txt"), out9("spNdClusters.txt"), out10("spNdClustersPar.txt"), out11("renamed_snIDs.txt");

	ListOfList adjlist = gp.EdgeToAdjList(edgelist);
	
	// Sorting the adjacency list in-place for faster edge intersection
	gp.sortInPlaceAdjList(adjlist);
	
	// declaring parent component ID for each edge
	//map<Edge, int> p;
	
	// This function initialize the parent component ID of each edge to itself
	gp.initializeEdgeParent(edgelist);

	gp.populateIntersectList(edgelist, adjlist);

	std::ifstream trussinput(ktruss_input.c_str());
	// This function currently read the trussness of all the edges from
	// a file generated by serial code k-truss decomposition
	std::map<int, std::vector<MyEdge>> trussgroups = gp.readTruss(trussinput, support);

	hooking_time = 0.0;
	compression_time = 0.0;
	sp_edge_time = 0.0;
	summary_merge_time = 0.0;
	parallel_summary_graph_time = 0.0;
	gapbs_sv_conn_time = 0.0;
	gapbs_afforest_conn_time = 0.0;
	serial_sorting_time = 0.0;
	parallel_sorting_time = 0.0;

	auto start = std::chrono::high_resolution_clock::now();


	EGraph egraph(support, gp, trussgroups, edgelist);

	gapbs::pvector<int> comp = egraph.conn_comp();

	gapbs::pvector<int> comp_afforest = egraph.conn_comp_afforest();


	ConnComp sv;
	// This function is the old way (before gapbs integration) Shiloach-Vishkin connected component over edges, commenting out this one
	//std::vector<std::set<std::pair<int, int>>> super_edges = sv.conn_comp_edge(edgelist, trussgroups, gp);
	// This is the new way of calling sv.conn_comp_edge() after gapbs integration
	std::vector<std::set<std::pair<int, int>>> super_edges = sv.conn_comp_edge(support, edgelist, trussgroups, gp);

	// This function is right now serial, need to design parallel implementation
	//std::set<std::pair<int, int>> summary_graph = sv.createSummaryGraph(super_edges);

	std::vector<std::pair<int,int>> summaryG = sv.mergeSummaryGraph(super_edges);

	auto end = std::chrono::high_resolution_clock::now();
	
	double totalExecutionTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	
	
	int total_iters = 0;

	for(auto it = iters_per_k.begin(); it != iters_per_k.end(); it++)
	{
		total_iters += it->second;
		printf("k:%d, iterations:%d\n", it->first, it->second);
	}

	printf("total_iterations_for_hooking:%d\n", total_iters);

	//printEdgeParent(out3, gp);

	//printSummaryGraph(out4, summary_graph);
	
	//printSummaryGraph(out5, summaryG);

	std::vector<int> sortedSNData = sortingSpNode(comp_afforest);

	std::vector<int> par_sortedSNData = sortingSpNodeParallel(comp_afforest);

	std::map<size_t, size_t> reMap_SNIDs = RenameSNIDs(par_sortedSNData, comp_afforest);

	printf("========totalExecutionTime:%0.9f===========\n", totalExecutionTime*(1e-9));
	printf("========Time for hooking phase:%0.9f===========\n", hooking_time*(1e-9));
	printf("========Time for compression phase:%0.9f===========\n", compression_time*(1e-9));
	printf("========Time for super edges creation phase:%0.9f===========\n", sp_edge_time*(1e-9));
	printf("========Time for merging to summary graph creation:%0.9f===========\n", summary_merge_time*(1e-9));
	printf("========Time for parallel summary graph creation:%0.9f===========\n", parallel_summary_graph_time*(1e-9));
	printf("========Time for gapbs ShiloachVishkin connected components:%0.9f===========\n", gapbs_sv_conn_time*(1e-9));
	printf("========Time for gapbs Afforest connected components:%0.9f===========\n", gapbs_afforest_conn_time*(1e-9));
	printf("========Time for sequential sorting connected components:%0.9f===========\n", serial_sorting_time*(1e-9));
	printf("========Time for parallel sorting connected components:%0.9f===========\n", parallel_sorting_time*(1e-9));

	
	printComp(out6, edgelist, comp);

	printComp(out7, edgelist, comp_afforest);

	//printComp(out8, edgelist, comp, support);

	printSpNdClusters(out9, sortedSNData, comp_afforest);

	printSpNdClusters(out10, par_sortedSNData, comp_afforest);

	printRenamedSNIDs(out11, reMap_SNIDs);


	return 0;
}
