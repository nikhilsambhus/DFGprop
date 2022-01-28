#include "DFGUtils.h"
#include "Node.h"
#include "Graph.h"
#include "Edge.h"
#include "GraphUtils.h"
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>

DAG normalizeDAG(DAG &dfg) {
	map<uint32_t, uint32_t> NodeMap;
	DAG newDfg;
	int count = 0;
	for(list<Node>::iterator it = dfg.nodeBegin(); count < dfg.getNumNodes(); it++, count++) {
		//cout << it->getID() << it->getLabel() << endl;
		NodeMap[it->getID()] = count;
		newDfg.addNode(count, it->getLabel());

	}
	
	count = 0;
	for(list<Edge>::iterator it = dfg.edgeBegin(); count < dfg.getNumEdges(); it++, count++) {
		cout << it->getLabel() << "\n";
		newDfg.addEdge(count, NodeMap[it->getSrcNodeID()], NodeMap[it->getDestNodeID()], it->getLabel());
	}
	
	return newDfg;
	
}
int main(int argc, char **argv) {
	string fname = argv[1];
	try {
		DAG graph(fname);
		DAG newDfg = normalizeDAG(graph);
		string tfname = std::string("norm_") + fname;
		toDOT(tfname, newDfg);
	} catch (string rx) {
		cout << rx << endl;
	}
	return 0;
}
