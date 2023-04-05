#ifndef GRAPHMANIP_H
#define GRAPHMANIP_H

#include "global.h"

typedef std::pair<int, int> MyEdge;
typedef std::vector<MyEdge> MyEdgeList;
typedef std::vector<std::vector<int>> ListOfList;
typedef std::map<MyEdge, int> Edge2Item; 

class GraphManip
{	
	public:
	
	long long numEdges;
	long long numVertices;
	int kmin;
	int kmax;
	ListOfList intersectlist;
	Edge2Item edge2k;
	Edge2Item edge2index;
	Edge2Item edge2P;
	
	GraphManip();
	long long getNumVertices(const MyEdgeList& edges);
	MyEdgeList ReadEdgeListFromFile(const char* filename);
	MyEdgeList BuildEdgeListFrom(const gapbs::WGraph& support);
	ListOfList EdgeToAdjList(const MyEdgeList& edges);
	void sortInPlaceAdjList(ListOfList& adjlist);
	void initializeEdgeParent(const MyEdgeList& edges);
	void populateIntersectList(const MyEdgeList& edges, ListOfList& adjlist);
	std::map<int, std::vector<MyEdge>> readTruss(std::ifstream& in);
	std::map<int, std::vector<MyEdge>> readTruss(std::ifstream& in, gapbs::WGraph& support);
}; 

#endif
