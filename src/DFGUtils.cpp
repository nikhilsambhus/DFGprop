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
uint32_t assignTime(DAG &gp, vector<uint32_t> &topOrder, vector<uint32_t> &timeSt) {
	for(uint32_t i = 0; i < gp.getNumNodes(); i++) {
		timeSt.push_back(-1);
	}

	uint32_t timeMax = 0;
	for(uint32_t n : topOrder) {
		list<Node> preds;
		gp.getPredecessors(n, preds);
		uint32_t max = 0;
		for(Node prNode : preds) {
			uint32_t id = prNode.getID();
			if(timeSt[id] == -1) {
				cout << "Error, predecessor time cannot be -1\n";
				return -1;
			}
			else if(max < timeSt[id]) {
				max = timeSt[id];
			}
		}
		timeSt[n] = max + 1;
		if(timeMax < timeSt[n]) {
			timeMax = timeSt[n];
		}
	}

	return timeMax;
}
uint32_t getParallelism(DAG &gp) {
	vector<uint32_t> topOrder = topoSort(gp);
	vector<uint32_t> timeSt; 
	uint32_t timeMax = assignTime(gp, topOrder, timeSt);
	cout << "Time Max is " << timeMax << endl;
	cout << "Topological sort order is\n";
	for(uint32_t n : topOrder) {
		cout << gp.findNode(n)->getLabel() << " has time " << timeSt[n] << endl;
	}
	
	cout << endl;
	return 0;
}
int main(int argc, char **argv) {
	string fname = argv[1];
	DAG graph(fname);

	getParallelism(graph);
	return 0;
}
