#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <limits.h>

using namespace std;

typedef unsigned long long ul_long;
typedef pair<ul_long, ul_long> edge;

int main(int argc, char* argv[])
{
	string inputfile = argv[1];
	string outputfile = argv[2];
	
	cout << "Reading network " << inputfile << " file\n";
	ifstream infile(inputfile);
	
	if (infile)
	{
		printf("File open successful\n");
	}
	else
	{
		printf("Failed to read file\n");
		exit(1);
	}

	auto cmp = [](const edge& e1, const edge& e2)
	{
		if(e1.first == e2.first)
		{
			return (e1.second < e2.second);
		}
		return (e1.first < e2.first);
	};
	
	set<edge, decltype(cmp)> edgelist(cmp);
	set<ul_long> vertices;
	map<ul_long, ul_long> ordering;
	ul_long i = 0;
	
	ofstream outfile(outputfile);
	string line = "";

	while (getline(infile, line))
	{
		istringstream iss(line);
		ul_long src, dst;
		if ((iss >> src >> dst))
		{
			if(src == dst)
			{
				continue;
			}
			if(src > dst)
			{
				ul_long temp = src;
				src = dst;
				dst = temp;
			}
			vertices.insert(src);
			vertices.insert(dst);
			edgelist.insert({src, dst});
			//outfile << src << "\t" << dst << endl;
		}
	}
	
	for(auto it = vertices.begin(); it != vertices.end(); it++)
	{
		ordering[(*it)] = i;
		++i;
	}
	printf("Number of vertices:%lu\n", vertices.size());
	printf("Number of edges:%lu\n", edgelist.size());
	
	for(auto it = edgelist.begin(); it != edgelist.end(); it++)
	{
		ul_long u = ordering[(*it).first];
		ul_long v = ordering[(*it).second];
		outfile << u << "\t" << v << endl;
	}
	
	
	infile.close();
	outfile.close();
	
	return 0;

}