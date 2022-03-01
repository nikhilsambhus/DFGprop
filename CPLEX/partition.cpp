// --------------------------------------------------------------------------------
//
// partition.cpp - Use the node and branch callbacks for optimizing a 01ILP problem
//
// Use:
//     partition  data.lp
//

#include <ilcplex/ilocplex.h>
#include <chrono>
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
	int loadWeight = 10;
	//for each pair, there are {k,l} pair mappings
	map<pair<int, int>, IloBoolVar> ijMap; //map of ij values
	map<pair<int, int>, IloBoolVar> WipMap; //map of wip values
	map<pair<int, int>, IloBoolVar> YipMap; //map of yip values
	vector<map<pair<int, int>, IloBoolVar>> klMapVec;//map of kl values for cross partition edges
	public:
	PartitionILP(DAG gp, int rsize, int tsize, int nPts) {
		modelPtr = new IloModel(env);
		cplexPtr = new IloCplex(env);
		varPtr = new IloNumVarArray(env);
		graph = gp;
		RSize = rsize; //partition size
		TSize = tsize; //transaction limit size
		numVertices = gp.getNumNodes();
		numEdges = gp.getNumEdges();
		numParts = nPts;
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
				string name = "x";
				name = name + to_string(i);
				name = name + ",";
				name = name + to_string(j);
				ijMap[{i, j}] = IloBoolVar(env, name.c_str());
				varPtr->add(ijMap[{i, j}]);
				objective.setLinearCoef(ijMap[{i, j}], 0);
				count++;
			}
		}
		cout << "Xij variable added count = " << count << endl;
		
		count = 0;
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts - 1; j++) {
				string yname = "y";
				yname = yname + to_string(i);
				yname = yname + ",";
				yname = yname + to_string(j);

				YipMap[{i, j}] = IloBoolVar(env, yname.c_str());

				string wname = "w";
				wname = wname + to_string(i);
				wname = wname + ",";
				wname = wname + to_string(j);

				WipMap[{i, j}] = IloBoolVar(env, wname.c_str());

				varPtr->add(YipMap[{i, j}]);
				varPtr->add(WipMap[{i, j}]);
				count += 2;
				objective.setLinearCoef(WipMap[{i, j}], 1);
			}
		}
		cout << "Y and W variables added count = " << count << endl;

		/*this->klMapVec.clear();
		count = 0;
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
			nCons++;
			modelPtr->add(range);
		}
		cout << "Number of size constraint rows " << nCons << endl;
	}

	void addEdgePrec() {
		int nCons = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src_i = it->getSrcNodeID();
			uint32_t dest_j = it->getDestNodeID();
			IloRange range = IloRange(env, -IloInfinity, 0);

			//sum of sources
			for(int i = 0; i < numParts; i++) {
				range.setLinearCoef(ijMap[{src_i, i}], i);			
			}

			//sum of dests
			for(int j = 0; j < numParts; j++) {
				range.setLinearCoef(ijMap[{dest_j, j}], -j);
			}

			modelPtr->add(range);
			nCons++;

		}
		cout << "Number of edge precedence constraints " << nCons << endl;
	}

	void addWriteCons() {
		int nCons = 0;
		for(list<Node>::iterator it = graph.nodeBegin(); it != graph.nodeEnd(); it++) {
			uint32_t i_v = it->getID();
			for(int p = 0; p < numParts - 1; p++) {
				IloRange range1 = IloRange(env, -IloInfinity, 0);
				IloRange range2 = IloRange(env, -IloInfinity, 0);
				//sum over all successors of this node
				list<Node> succ;
				graph.getSuccessors(i_v, succ);
				if(succ.size() == 0) {
					//just set yip and wip to 0 and continue
					IloRange y1 = IloRange(env, 0, 0);
					IloRange w1 = IloRange(env, 0, 0);
					y1.setLinearCoef(YipMap[{i_v, p}], 1);
					w1.setLinearCoef(WipMap[{i_v, p}], 1);
					modelPtr->add(y1);
					modelPtr->add(w1);
					nCons += 2;
					continue;
				}
				for(auto nd : succ) {
					uint32_t j_v = nd.getID();
					for(int l = p + 1; l < numParts; l++) { 
						range1.setLinearCoef(ijMap[{j_v, l}], -1);
						range2.setLinearCoef(ijMap[{j_v, l}], 1);
					}
				}
				range1.setLinearCoef(YipMap[{i_v, p}], 1);
				range2.setLinearCoef(YipMap[{i_v, p}], -1 * (int)succ.size());

				modelPtr->add(range1);
				modelPtr->add(range2);
				
				IloRange range3 = IloRange(env, -IloInfinity, 1);
				IloRange range4 = IloRange(env, -IloInfinity, 0);

				range3.setLinearCoef(ijMap[{i_v, p}], 1);
				range3.setLinearCoef(YipMap[{i_v, p}], 1);
				range3.setLinearCoef(WipMap[{i_v, p}], -1);
				
				range4.setLinearCoef(ijMap[{i_v, p}], -1);
				range4.setLinearCoef(YipMap[{i_v, p}], -1);
				range4.setLinearCoef(WipMap[{i_v, p}], 2);

				modelPtr->add(range3);
				modelPtr->add(range4);

				nCons += 4;
			}
		}

		cout << "Total set of write constraints " << nCons << endl;

	}

	void addCommCons() {
		int i = 0;
		
		//Constraints modelling edges as communication which contain both intrapartition and same partition edges
		int nCons = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src_i = it->getSrcNodeID();
			uint32_t dest_j = it->getDestNodeID();

			//first constraint sum (p < l) Xi_j^p_l = Xj_l
			for(int l = 0; l < numParts; l++) {
				IloRange range = IloRange(env, 0, 0);
				for(int p = 0; p <= l; p++) {
					range.setLinearCoef(klMapVec[i][{p, l}], 1);
				}
				range.setLinearCoef(ijMap[{dest_j, l}], -1);
				modelPtr->add(range);
				nCons++;
			}

			//second constraint sum (p > k) Xi_j^k_p = Xi_k
			for(int k = 0; k < numParts; k++) {
				IloRange range = IloRange(env, 0, 0);
				for(int p = k; p < numParts; p++) {
					range.setLinearCoef(klMapVec[i][{k, p}], 1);
				}
				range.setLinearCoef(ijMap[{src_i, k}], -1);
				modelPtr->add(range);
				nCons++;
			}

			i++;
		}

		cout << "Number of communication constraint rows added " << nCons << endl;
	}
	
	void addTransWriteCons() {
		int nCons = 0;

		for(int p = 0; p < numParts - 1; p++) {
			IloRange range = IloRange(env, 0, TSize);
			for(int v = 0; v < numVertices; v++) {
				range.setLinearCoef(WipMap[{v, p}], 1);
			}
			nCons++;
			modelPtr->add(range);
		}
		cout << "Transaction write constraints added = " << nCons << endl;
	}
	void addTransCons() {

		int nCons = 0;
		//(l > k) X^kl_ij < T
		for(int k = 0; k < numParts - 1; k++) {
			IloRange range = IloRange(env, 0, TSize);
			for(int l = k + 1; l < numParts; l++) {
				for(auto klMap: this->klMapVec) {
					range.setLinearCoef(klMap[{k, l}], 1);
				}
			}
			modelPtr->add(range);
			nCons++;
		}

		//(k < l) X^kl_ij < T
		for(int l = 1; l < numParts; l++) {
			IloRange range = IloRange(env, 0, TSize);
			for(int k = 0; k < l; k++) {
				for(auto klMap: this->klMapVec) {
					range.setLinearCoef(klMap[{k, l}], 1);
				}
			}
			modelPtr->add(range);
			nCons++;
		}

		cout << "Number of transaction constraint rows added " << nCons << endl;

	}

	bool solve() {
		int countMaps = 0;
		try {
			cplexPtr->extract(*modelPtr);
			cplexPtr->exportModel("test.lp");
			if(!cplexPtr->solve()) {
				cout << "Failed to optimize" << endl;
				return false;	
			}

			cout << "Status value = " << cplexPtr->getStatus() << endl;
			cout << "Objective function value = " << cplexPtr->getObjValue() << endl;
			cout << "Number of rows = " << cplexPtr->getNrows() << endl;
			cout << "Number of cols = " << cplexPtr->getNcols() << endl;

			IloNumArray vals(env);
			cplexPtr->getValues(vals, *varPtr);
			//cout << "Solution vector = " << vals << endl; 
			/*for(int i = 0; i < numVertices; i++) {
				for(int j = 0; j < numParts; j++) {
					if(cplexPtr->getValue(ijMap[{i, j}])) {
						countMaps++;
						cout << "Node " << i << " is mapped to " << j << endl;
					}
				}
			}*/
		}
		catch (IloException ex) {
			cout << ex << endl;
		}
		//return (countMaps == numVertices); //return true if 1-1 mapping done
		return true;
	}

	
	bool compareEqual(double val1, double val2) {
		if(fabs(val1 - val2) < 1e-5) {
			return true;
		}
		return false;
	}
	void ValidateSoln() {
		cout << "Asserting uniqueness constraints" << endl;
		//Check if vertex mapped to only one partition
		for(int i = 0; i < numVertices; i++) {
			int count = 0;
			for(int j = 0; j < numParts; j++) {
				double val = cplexPtr->getValue(ijMap[{i, j}]);
				if(compareEqual(val, 1) == true) {
					count++;
				}
			}
			/*if(count != 1) {
				cout << "Vertex " << i << " is mapped to these many parts " << count << endl;
				for(int j = 0; j < numParts; j++) {
					cout << i << " " << j << " pos is " << cplexPtr->getValue(ijMap[{i, j}]) << " ";
				}
			}*/
			assert(count == 1); //only one vexter needs to be mapped to some partition
		}

		cout << "Asserting size constraints" << endl;
		for(int i = 0; i < numParts; i++) {
			int count = 0;
			for(uint32_t j = 0; j < graph.getNumNodes(); j++) {
				double val = cplexPtr->getValue(ijMap[{j, i}]);
				if(compareEqual(val, 1)) {
					count++; //add if vertex present in this partition
				}
			}

			assert(count <= RSize);
		}

		cout << "Asserting edge precedences" << endl;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src = it->getSrcNodeID();
			uint32_t dest = it->getDestNodeID();
			//find source partition
			int srcPart = -1;

			for(int j = 0; j < numParts; j++) {
				double val = cplexPtr->getValue(ijMap[{src, j}]);
				if(compareEqual(val, 1)) {
					srcPart = j;
					break;
				}
			}
			assert(srcPart != -1); //mapping should be found

			//find destination partition
			int destPart = -1;

			for(int j = 0; j < numParts; j++) {
				double val = cplexPtr->getValue(ijMap[{dest, j}]);
				if(compareEqual(val, 1)) {
					destPart = j;
					break;
				}
			}
			assert(srcPart != -1); //mapping should be found
			assert(srcPart <= destPart); //partition of source should be less than destination


		}

		cout << "Asserting Wips and trans limits " << endl;
		
		map<int, int> writeCount; //for each partition
		for(int i = 0; i < numParts; i++) {
			writeCount[i] = 0;
		}
		for(int v = 0; v < numVertices; v++) {
			for(int p = 0; p < numParts - 1; p++) {
				//check if v is in partition p
				double val = cplexPtr->getValue(ijMap[{v, p}]);
				if(compareEqual(val, 1) == true) {
					list<Node> succ;
					graph.getSuccessors(v, succ);
					//check if any successor mapped to partition greater than p
					//check successors of v where they are mapped
					bool someNodeMap = false;
					for(auto nd : succ) {
						int j = nd.getID();
						bool isNodeMap = false;
						for(int l = p + 1; l < numParts; l++) {
							double val = cplexPtr->getValue(ijMap[{j, l}]);
							if(compareEqual(val, 1)) {
								isNodeMap = true;
								break;
							}
						}
						if(isNodeMap == true) {
							someNodeMap = true;
							break;
						}
					}
					
					//if some succ node is mapped to next partition 
					if(someNodeMap == true) {
						double val = cplexPtr->getValue(WipMap[{v, p}]);
						assert(compareEqual(val, 1) == true);
						//inc count of base node writes
						writeCount[p] = writeCount[p] + 1;
					}else {
						double val = cplexPtr->getValue(WipMap[{v, p}]);
						assert(compareEqual(val, 0) == true);
					}
				}
				else { //if v is not in p them Wvp must be 0
					double val = cplexPtr->getValue(WipMap[{v, p}]);
					assert(compareEqual(val, 0) == true);
				}
			}
		}


		cout << "Write counts ";
		for(int i = 0; i < numParts; i++) {
			//write trans count
			cout << writeCount[i] << " ";
			assert(writeCount[i] <= TSize);
		}
		cout << endl;
	}
};

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
	int iterations = 100;

	auto start = chrono::high_resolution_clock::now();
	int numParts = ceil(float(gp.getNumNodes()) / float(size)); //set initial partition size to total vertices divided by partition size
	for(int i = 1; i <= iterations; i++) {
		PartitionILP *gp1 = new PartitionILP(gp, size, trans_limit, numParts);
		gp1->addColVars();//set objective function; define all vars
		gp1->addUniqueCons();
		gp1->addSizeCons();
		gp1->addEdgePrec();
		gp1->addWriteCons();
		gp1->addTransWriteCons();
		//gp1->addCommCons();
		//gp1->addTransCons();
		if(gp1->solve() == true) {
			gp1->ValidateSoln();
			cout << "Solution found in iteration number " << i << " with partitions " << numParts << endl;
			break;
		}
		numParts++;
		delete gp1;
	}
	auto stop = chrono::high_resolution_clock::now();

	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
	
	cout << "Time taken in seconds " << duration/1000.0 << endl;
	
	return 0;
}
