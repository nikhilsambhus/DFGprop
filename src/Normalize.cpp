#include "DFGUtils.h"
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>

void normalizeDAG(DAG &dfg) {
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
		newDfg.addEdge(count, NodeMap[it->getSrcNodeID()], NodeMap[it->getDestNodeID()], it->getLabel());
	}
	
	string tfname = "roialign_norm.dot";
	toDOT(tfname, newDfg);
	
}
int main(int argc, char **argv) {
	string fname = argv[1];
	DAG graph(fname);
	normalizeDAG(graph);
	return 0;
}
