// --------------------------------------------------------------------------------
//
// partition.cpp - Use the node and branch callbacks for optimizing a 01ILP problem
//
// Use:
//     partition  data.lp
//

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <cstdint>
#include "Edge.h"
#include "Graph.h"
#include <vector>
ILOSTLBEGIN
using namespace std;
class PartitionILP {
	private:
	IloEnv env;
	IloCplex *cplexPtr;
	IloModel *modelPtr;
	IloNumVarArray *varPtr;
	DAG graph;
	int RSize;
	int TSize;
	int numVertices;
	int numEdges;
	int numParts;
	map<pair<int, int>, IloBoolVar> ijMap; //map of ij values
	vector<map<pair<int, int>, IloBoolVar>> klMapVec;//map of kl values for cross partition edges
	public:
	PartitionILP(DAG gp, int rsize, int tsize) {
		modelPtr = new IloModel(env);
		cplexPtr = new IloCplex(env);
		varPtr = new IloNumVarArray(env);
		graph = gp;
		RSize = rsize; //partition size
		TSize = tsize; //transaction limit size
		numVertices = gp.getNumNodes();
		numEdges = gp.getNumEdges();
		numParts = ceil(float(numVertices) / float(RSize)); //set initial partition size to total vertices divided by partition size
		cout << "Num Parts trying with " << numParts << endl;
	}

	void addColVars() {
		//add all vertices-parts mapping as cols

		this->ijMap.clear();
		//define Xij's
		//add Xij's to objective function
		IloObjective objective = IloMinimize(env);
		int count = 0;
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				ijMap[{i, j}] = IloBoolVar(env);
				varPtr->add(ijMap[{i, j}]);
				objective.setLinearCoef(ijMap[{i, j}], 1);
				count++;
			}
		}
		cout << "Xij variable added count = " << count << endl;

		this->klMapVec.clear();
		/*count = 0;
		//add columns for each edge (parts * parts) for communication objective function
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			map<pair<int, int>, IloBoolVar> klMap;
			for(int k = 0; k < numParts; k++) {
				for(int l = k; l < numParts; l++) {
					klMap[{k ,l}] = IloBoolVar(env);
					varPtr->add(klMap[{k, l}]);
					if(k != l) {
						objective.setLinearCoef(klMap[{k, l}], 1);
					}
					else 
						objective.setLinearCoef(klMap[{k, l}], 0);
				}
				count++;
			}
			this->klMapVec.push_back(klMap);
		}
		
		cout << "Xij^kl variable added count = " << count << endl;
		*/
		modelPtr->add(objective);
	}

	//add uniqueness constraint of mapping n vertices to p partitions
	void addUniqueCons() {
		int nCons = 0;
		for(int i = 0; i < numVertices; i++) {
			IloRange range = IloRange(env, 1, 1);
			for(int j = 0; j < numParts; j++) {
				range.setLinearCoef(ijMap[{i,j}], 1);
			}
			modelPtr->add(range);
			nCons++;
		}
		cout << "Number of uniquness constraints rows added " << nCons << endl;
	}
	
	void addSizeCons() {
		int nCons = 0;
		for(int i = 0; i < numParts; i++) {
			IloRange range = IloRange(env, 0, RSize);
			for(int j = 0; j < numVertices; j++) {
				range.setLinearCoef(ijMap[{j, i}], 1);
			}
			modelPtr->add(range);
		}
	}

	bool solve() {
		cplexPtr->extract(*modelPtr);
		if(!cplexPtr->solve()) {
			cout << "Failed to optimize" << endl;
			return false;	
		}

		IloNumArray vals(env);
		
		cplexPtr->getValues(vals, *varPtr);
		cout << "Status value = " << cplexPtr->getStatus() << endl;
		cout << "Objective function value = " << cplexPtr->getObjValue() << endl;
		cout << "Solution vector = " << vals << endl; 
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				if(cplexPtr->getValue(ijMap[{i, j}])) {
					cout << i << j << endl;
				}
			}
		}

	}
};
static void usage (const char *progname);

int main (int argc, char **argv)
{
	if(argc != 4) {
		cout << "Too few arguments, 3 expected" << endl;
		return -1;
	}
	DAG gp;

	char *inpName = argv[1];
	try {
		gp = DAG(inpName);
	} catch(string ex) {
		cout << ex << endl;
	}

	int size = atoi(argv[2]);
	int trans_limit = atoi(argv[3]);
	int iterations = 10;
	PartitionILP *gp1 = new PartitionILP(gp, size, trans_limit);
	for(int i = 1; i <= 10; i++) {
		gp1->addColVars();//set objective function; define all vars
		gp1->addUniqueCons();
		gp1->addSizeCons();
		gp1->solve();
		break;
	}

	return 0;
}
