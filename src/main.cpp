#include "DFGPart.h"
#include "DFGAnaly.h"
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>

int main(int argc, char **argv) {
	string fname = argv[1];
	try {
		DAG graph(fname);
		cout << "Created a DAG of fname " << fname << endl; 
//		return 0;

		DFGAnaly dfgA(graph);
		DFGPart dfgP(graph);
		double par = dfgA.getParallelism();
		cout << "Parallelism in graph is " << par << endl;
		uint32_t clen = dfgA.criticalPathLen();
		cout << "Critical path length is " << clen << endl;

		dfgA.getBasicProps();

		int cgra_size = atoi(argv[2]);
		int percentage = atoi(argv[3]);
		float routing_size = percentage/100 * cgra_size;

		int map_size = cgra_size - routing_size;
		dfgP.partitionDFGVar(map_size);
	}catch(string er) {
		cout << er << endl;
	}
	return 0;
}
