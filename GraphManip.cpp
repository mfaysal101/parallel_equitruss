#include "GraphManip.h"

namespace
{
	template <class InputIterator1, class InputIterator2>
	std::vector<int> IntersectionSize(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, int n)
	{
		std::vector<int> result;
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

GraphManip::GraphManip()
{
	numEdges = 0;
	numVertices = 0;
	kmin = std::numeric_limits<gapbs::WeightT>::max();
	kmax = std::numeric_limits<gapbs::WeightT>::min();
}

long long GraphManip::getNumVertices(const MyEdgeList& edges)
{
	int num = 0;

#pragma omp parallel for reduction (max:num)
	for (size_t i = 0; i < edges.size(); i++)
	{
		num = std::max(num, 1 + std::max(edges[i].first, edges[i].second));
	}

	numVertices = num;

	std::cout << "Total number of vertices:" << numVertices << std::endl;

	return numVertices;
}

MyEdgeList GraphManip::ReadEdgeListFromFile(const char* filename)
{
	MyEdgeList edgelist;

	std::cout << "Reading network " << filename << " file\n";

	std::ifstream infile(filename);

	if (infile)
	{
		printf("File open successful\n");
	}
	else
	{
		printf("Failed to read file\n");
		exit(1);
	}

	std::string line = "";
	numEdges = 0;

	int index = 0;

	while (getline(infile, line))
	{
		std::istringstream iss(line);
		int src, dst;
		if ((iss >> src >> dst))
		{
			if(src == dst)
			{
				continue;
			}
			edgelist.push_back({src, dst});
		}
	}

	numEdges = edgelist.size();

	std::cout << "Total number of edges:" << numEdges << std::endl;

	sort(edgelist.begin(), edgelist.end(), [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});

	for(auto edge:edgelist)
	{
		edge2index[edge] = index++;
	}
	
	return edgelist;
}


MyEdgeList GraphManip::BuildEdgeListFrom(const gapbs::WGraph& support)
{
	MyEdgeList edgelist;
	numEdges = support.num_edges();
	edgelist.reserve(support.num_edges());
	for(int u = 0; u < support.num_nodes(); ++u) 
	{
    		for(auto& dest : support.out_neigh(u)) 
		{
			auto v = dest.v;
			if ( u < v )
			{
				edgelist.push_back({u,v});
			}
		}
	}

	std::cout << "Total number of edges:" << numEdges << std::endl;

	sort(edgelist.begin(), edgelist.end(), [](const std::pair<int, int>& edge1, const std::pair<int, int>& edge2) {
		return (edge1.first < edge2.first) || (edge1.first == edge2.first && edge1.second < edge2.second);
		});

	int index{0};
	for(auto edge:edgelist)
	{
		edge2index[edge] = index++;
	}
	
	return edgelist;
}



ListOfList GraphManip::EdgeToAdjList(const MyEdgeList& edges)
{
	ListOfList adjlist(getNumVertices(edges));

	for (auto edge : edges)
	{
		adjlist[edge.first].push_back(edge.second);
		adjlist[edge.second].push_back(edge.first);
	}
	
	return adjlist;
}



void GraphManip::sortInPlaceAdjList(ListOfList& adjlist)
{
	int len = adjlist.size();

#pragma omp parallel for
	for (int i = 0; i < len; i++)
	{
		sort(adjlist[i].begin(), adjlist[i].end());
	}
}


void GraphManip::populateIntersectList(const MyEdgeList& edges, ListOfList& adjlist)
{
	const int n = adjlist.size();
	
	intersectlist.resize(edges.size());

	for (int i = 0; i < edges.size(); i++)
	{
		int u = edges[i].first;
		int v = edges[i].second;
		intersectlist[i] = IntersectionSize(adjlist[u].begin(), adjlist[u].end(), adjlist[v].begin(), adjlist[v].end(), n);
	}
}


std::map<int, std::vector<MyEdge>> GraphManip::readTruss(std::ifstream& in)
{
	std::string line = "";

	std::map<int, std::vector<MyEdge>> trussgroups;
	
	if (in.is_open())
	{
		while (getline(in, line))
		{
			std::istringstream iss(line);
			std::string token;
			std::vector<int>tuple;

			while (getline(iss, token, ','))
			{
				tuple.push_back(std::stoi(token));
			}

			std::pair<int, int> temp;
			temp.first = tuple[0];
			temp.second = tuple[1];
			kmin = std::min(kmin, tuple[2]);
			kmax = std::max(kmax, tuple[2]);
			edge2k.insert(make_pair(temp, tuple[2]));
			trussgroups[tuple[2]].push_back(temp);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		trussgroups.erase(i);
	}
	kmin = std::max(kmin, 3);

	in.close();
	
	return trussgroups;
}

std::map<int, std::vector<MyEdge>> GraphManip::readTruss(std::ifstream& in, gapbs::WGraph& support)
{
	
	std::map<int, std::vector<MyEdge>> trussgroups;
	std::string line = "";
	std::vector<int>tuple;

	auto setTrussness = [&support](int u, int v, int k) 
	{
		auto& dest_ = support.get_unsafe(u,v);
		//operator() for NodeId returns the node id
		if (dest_.v != v) 
		{
			std::cerr << "Error: unable to find vertex {"<<u<<","
			 << v<<"} in graph" << std::endl;
		}
		dest_.w = k;
	};

	
	if (in.is_open())
	{
		while (getline(in, line))
		{
			tuple.clear();
			std::istringstream iss(line);
			std::string token;

			while (getline(iss, token, ','))
			{
				tuple.push_back(std::stoi(token));
			}

			std::pair<int, int> temp;
			temp.first = tuple[0];
			temp.second = tuple[1];
			kmin = std::min(kmin, tuple[2]);
			kmax = std::max(kmax, tuple[2]);
			//duplicated so that we store full trussness
                        	setTrussness(temp.first, temp.second, tuple[2]);
                        	setTrussness(temp.second, temp.first, tuple[2]);
			trussgroups[tuple[2]].push_back(temp);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		trussgroups.erase(i);
	}
	kmin = std::max(kmin, 3);

	in.close();
	
	return trussgroups;
}


void GraphManip::initializeEdgeParent(const MyEdgeList& edges)
{
	for (int i = 0; i < edges.size(); i++)
	{
		edge2P[edges[i]] = i;
	}
}
