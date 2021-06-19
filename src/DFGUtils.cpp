#include "DFGUtils.h"
#include <vector>
#include <stack>
#include <list>

void topoSortHelper(DAG &gp, uint32_t node, vector<bool> &visited, stack<uint32_t> &st) {
	visited[node] = true;
	list<Node> succs;
	gp.getSuccessors(node, succs);
	for(Node nd : succs) {
		uint32_t ndId = nd.getID();
		if(visited[ndId] == false) {
			topoSortHelper(gp, ndId, visited, st);
		}
	}
	st.push(node);
}

vector<uint32_t> topoSort(DAG &gp) {
	vector<uint32_t> order;
	stack<uint32_t> stNodes;

	vector<bool> visited;
	for(uint32_t i = 0; i < gp.getNumNodes(); i++) {
		visited.push_back(false);
	}

	for(uint32_t i = 0; i < gp.getNumNodes(); i++) {
		if(visited[i] == false) {
			topoSortHelper(gp, i, visited, stNodes);
		}
	}

	while(!stNodes.empty()) {
		order.push_back(stNodes.top());
		stNodes.pop();
	}

	return order;
}
int main(int argc, char **argv) {
	string fname = argv[1];
	cout << fname << "\n";
	DAG graph(fname);
	vector<uint32_t> topOrder = topoSort(graph);

	cout << "Topological sort order is : ";
	for(uint32_t n : topOrder) {
		cout << graph.findNode(n)->getLabel() << " ";
	}
	cout << endl;

	return 0;
}
