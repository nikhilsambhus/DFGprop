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

vector<uint32_t> getAllPreds(DAG &dfg, uint32_t ldv) {
	vector<uint32_t> retIds; //return all pred Ids
	queue<uint32_t> predQ; //queue of predecessors
	predQ.push(ldv);//push first the load
	while(!predQ.empty()) {
		uint32_t vert = predQ.front(); //get front element 
		predQ.pop();//pop out / delete the front element
		
		//add to ret list of preds only if this is not the initial load
		if(vert != ldv) {
			retIds.push_back(vert);
		}

		//get all preds of this vert and add them to the queue
		list<Node> predList;
		dfg.getPredecessors(vert, predList);
		for(auto nd : predList) {
			predQ.push(nd.getID());
		}
	}

	return retIds;
}

DAG loadSanDAG(DAG &dfg) {
	map<uint32_t, bool> deleteNodes; //ids of Nodes before loads to be deleted
	
	//iterate through each node and find predec of load nodes; store nodes to be deleted
	for(auto it = dfg.nodeBegin(); it != dfg.nodeEnd(); ++it) {
		//check if load node
		string opName = it->getLabel();
		if(opName.find("LOD") != string::npos) {
			vector<uint32_t> allPreds = getAllPreds(dfg, it->getID());
			//go through all preds and add them in deleteNodes map..these are the ones to be deleted and map takes care it is accounted for once only
			for(uint32_t id : allPreds) {
				deleteNodes[id] = true;
			}
		}
	}

	map<uint32_t, uint32_t> oldNew; //mapping olf old values to new normalized values

	int norm = 0;//normalized vertex values start from 0

	DAG newDfg;
	
	//scan nodes and skip over deleted nodes
	for(auto it = dfg.nodeBegin(); it != dfg.nodeEnd(); ++it) {
		//check if node id has to be deleted
		uint32_t id = it->getID();
		if(deleteNodes.find(id) != deleteNodes.end()) {
			continue; //if yes then skip this node
		}
		else { //not to be deleted, store node with new id
			oldNew[id] = norm;
			newDfg.addNode(norm, it->getLabel());
			norm = norm + 1; 
		}
	}

	//scan the edges and skip if any of source or dest in deleted nodes
	//store the remaining edges using source and dest ids from oldNew map
	for(auto it = dfg.edgeBegin(); it != dfg.edgeEnd(); it++) {
		uint32_t src = it->getSrcNodeID();
		uint32_t dest = it->getDestNodeID();
		//add edge only if none of src and dest in delete nodes
		if((deleteNodes.find(src) == deleteNodes.end()) && (deleteNodes.find(dest) == deleteNodes.end())) {
			newDfg.addEdge(it->getID(), oldNew[src], oldNew[dest], it->getLabel());
		}
	}

	//print stats of orig and new dfg
	cout << "Original DFG node count " << dfg.getNumNodes() << " ";
	cout << "Edge count " << dfg.getNumEdges() << endl;
	cout << "New DFG node count " << newDfg.getNumNodes() << " ";
	cout << "Edge count " << newDfg.getNumEdges() << endl;
	return newDfg;
}
/* args : InputDotFileName, Outputfilename*/
int main(int argc, char **argv) {
	string fname = argv[1];
	try {
		DAG graph(fname);
		DAG newDfg = loadSanDAG(graph);
		string tfname = argv[2];
		toDOT(tfname, newDfg);
	} catch (string rx) {
		cout << rx << endl;
	}
	return 0;
}
