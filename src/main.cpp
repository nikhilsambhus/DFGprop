#include "DFGPart.h"
#include "DFGAnaly.h"
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>

int main(int argc, char **argv) {
	string fname = argv[1];
	DAG graph(fname);

	DFGAnaly dfgA(graph);
	DFGPart dfgP(graph);
	double par = dfgA.getParallelism();
	cout << "Parallelism in graph is " << par << endl;
	uint32_t clen = dfgA.criticalPathLen();
	cout << "Critical path length is " << clen << endl;
	
	dfgA.getBasicProps();

	//partitionDFG(graph);
	int cgra_size = atoi(argv[2]);
	int routing_size = atoi(argv[3]);
	int map_size = cgra_size - routing_size;
	dfgP.partitionDFGVar(map_size);
	return 0;
}
