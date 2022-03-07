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

	int nUniqCons = 0; //uniqness constraints rows count
	int nCapCons = 0; //capacity contsraints rows count
	int nPrecCons = 0; //precedence constraints rows count
	int nIntParCons = 0; //inter partition constraints rows count
	int nLdStCons = 0; //load/store constraints count
	int nTransCons = 0; //Transaction constraints count

	int nXij = 0; //number of assignment variables
	int nLp = 0; //number of loadgroup variables
	int nInPa = 0; //number of Xikl and Yikl variables

	int RSize;
	int TSize;
	int numVertices;
	int numEdges;
	int numParts;
	int loadWeight = 1; //weight for load from memory
	int writeWeight = 1; //weight for intermediate writes
	map<int, vector<int>> loadGroups;
	//for each pair, there are {k,l} pair mappings
	map<pair<int, int>, IloBoolVar> ijMap; //map of Xij vars
	vector<map<pair<int, int>, IloBoolVar>> XiklMap; //map of Xijk vars for inter communication..for each vertex i there are k and l partition number mappings
	vector<map<pair<int, int>, IloBoolVar>> YiklMap; //for each vertex i, there are kl pairs belonging to k and l partition number each
	
	map<int, vector<IloBoolVar>> lpMap; //loadgroup_id->lpvector..one lp vector for each load group..lp vector has each element for one partition

	vector<map<pair<int, int>, IloBoolVar>> klMapVec;//map of kl values for cross partition edges
	public:
	PartitionILP(DAG gp, int rsize, int tsize, int nPts, int loadWt) {
		modelPtr = new IloModel(env);
		cplexPtr = new IloCplex(env);
		varPtr = new IloNumVarArray(env);
		graph = gp;
		RSize = rsize; //partition size
		TSize = tsize; //transaction limit size
		loadWeight = loadWt; //weight of load
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
		this->nXij = count; //append count

		//iterate through graph nodes and store load ids in corresponding group vector
		for(list<Node>::iterator it = graph.nodeBegin(); it != graph.nodeEnd(); it++) {
			string op = it->getLabel();
			if(op.find("load") != string::npos || op.find("LOD") != string::npos) { ///two possible vals for describing load nodes
				//find group id which is placed after load;"
				int pos  = op.find(";");
				int group_id = stoi(op.substr(pos + 1, string::npos)); //from : + 1 till end of string
				
				if(loadGroups.find(group_id) != loadGroups.end()) {
					vector<int> &loadsV = loadGroups[group_id]; //group id data present..append to the vector
					loadsV.push_back(it->getID());
				} else {
					vector<int> loadsV;
					loadsV.push_back(it->getID());
					loadGroups[group_id] = loadsV;
				}
			}
		}

		count = 0;
		//add group number of lp vectors
		for(auto elem : loadGroups) {
			//add partition number of load boolean variable..each represent if there is atleast one load of the group in given partition
			vector<IloBoolVar> lpV;
			for(int i = 0; i < numParts; i++) {
				string name = "lp" + to_string(elem.first) + "," + to_string(i);
				lpV.push_back(IloBoolVar(env, name.c_str()));
				count++; //increment count of load group variables
				objective.setLinearCoef(lpV[i], loadWeight); //set objective function's coefficient to loadweight
			}
			lpMap[elem.first] = lpV;
		}
		this->nLp = count; //assign count to load group variables
		
		count = 0;
		//define kl variables for each i for inter communication..define for both y and x variables
		for(int i = 0; i < numVertices; i++) {
			map<pair<int, int>, IloBoolVar> xdMap; //define kl variables for xikl
			map<pair<int, int>, IloBoolVar> ydMap; //define kl variables for yikl
			for(int k = 0; k < numParts - 1; k++) {
				for(int l = k + 1; l < numParts; l++) {
					count += 2; //increment twice
					xdMap[{k, l}] = IloBoolVar(env, ("x_" + to_string(i) + "_" + to_string(k) + "_" + to_string(l)).c_str());
					varPtr->add(xdMap[{k, l}]); //add variable in the model
					ydMap[{k, l}] = IloBoolVar(env, ("y_" + to_string(i) + "_" + to_string(k) + "_" + to_string(l)).c_str());
					varPtr->add(ydMap[{k, l}]); //add variable in the model
				}
			}
			XiklMap.push_back(xdMap);//ap:pend the kl map for this particular vertex
			YiklMap.push_back(ydMap);//ap:pend the kl map for this particular vertex
		}

		this->nInPa = count;

		//read map to check which xikl should be made 1
		//this is done because set cofficient 1 simply sets it but if reads and writes are both present then we want set coefficient to be 2
		vector<map<pair<int, int>, bool>> RiklMap;
		//sum objective function over all each vertex i..for accounting reads based on Xikl
		for(int i = 0; i < numVertices; i++) {
			map<pair<int, int>, bool> tmp; 
			RiklMap.push_back(tmp);//add into read temp map for accounting later
			//sum of all the partitions starting from 2 to end (numbering starts from 0)
			for(int l = 1; l < numParts; l++) {
				//sum over all xikl variables where p < l
				for(int p = 0; p < l; p++) {
					//objective.setLinearCoef(XiklMap[i][{p, l}], 1);
					RiklMap[i][{p, l}] = true;//set pl to 1 and count it later while accounting objective for reads
				}
			}
		}

		//sum objective function over all each vertex i..for accounting writes based on Xikl
		for(int i = 0; i < numVertices; i++) {
			//sum of all partitions start from 1 to k - 1 (numbering starts from 0)
			for(int k = 0; k < numParts - 1; k++) {
				//sum over xikl variables for p > k
				for(int p = k + 1; p < numParts; p++) {
					if(RiklMap[i].find({k, p}) == RiklMap[i].end()) {//none was added while accounting for reads
						objective.setLinearCoef(XiklMap[i][{k, p}], 1);
					}
					else {
						//added while accounting for reads hence take the coefficient as 2 instead
						objective.setLinearCoef(XiklMap[i][{k,p}], 2);
					}
				}
			}
		}
		

		//summation in the objective function
		//print out number of loads present in the graph
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
		this->nUniqCons = nCons;
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
		this->nCapCons = nCons;
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
		nPrecCons = nCons;
		cout << "Number of edge precedence constraints " << nCons << endl;
	}

	//add constraints for inter partition
	void addInterPartCons() {
		int nCons = 0;
		for(list<Node>::iterator it = graph.nodeBegin(); it != graph.nodeEnd(); it++) {
			int i = it->getID(); //for each vertex i
			list<Node> succ;
			graph.getSuccessors(i, succ);
			if(succ.size() == 0) {
				//set all Xikls and Yikls to 0 for this particular source node
				for(int k = 0; k < numParts - 1; k++) {
					for(int l = k + 1; l < numParts; l++) {
						IloRange xikl = IloRange(env, 0, 0);
						IloRange yikl = IloRange(env, 0, 0);
						xikl.setLinearCoef(XiklMap[i][{k, l}], 1);
						yikl.setLinearCoef(YiklMap[i][{k, l}], 1);
						modelPtr->add(xikl);
						modelPtr->add(yikl);
						nCons += 2;
					}
				}

				continue;
			}
			
			//for each partition pair kl : sum Xjl (where j is sucessor) >= Yikl
			//for each pariition pair kl : sum Xjl (where j is successor) <= Size of succ * Yikl
			for(int k = 0; k < numParts - 1; k++) {
				for(int l = k + 1; l < numParts; l++) {
					IloRange Xy1 = IloRange(env, -IloInfinity, 0);
					IloRange Xy2 = IloRange(env, -IloInfinity, 0);
					for(auto nd: succ) {
						int j = nd.getID();
						Xy1.setLinearCoef(ijMap[{j, l}], -1);
						Xy2.setLinearCoef(ijMap[{j, l}], 1);
					}
					Xy1.setLinearCoef(YiklMap[i][{k, l}], 1);
					Xy2.setLinearCoef(YiklMap[i][{k, l}], -1 * (int)succ.size());

					nCons += 2;
					modelPtr->add(Xy1);
					modelPtr->add(Xy2);
				}
			}

			//for each kl pair of r define the following equations
			//Xik + Yikl <= 1 + Xikl
			//Xik + Yikl >= 2.Xikl
			for(int k = 0; k < numParts - 1; k++) {
				for(int l = k + 1; l < numParts; l++) {
					IloRange Xy1 = IloRange(env, -IloInfinity, 1);
					IloRange Xy2 = IloRange(env, -IloInfinity, 0);
					
					//Xik + Yikl -Xikl <= 1
					Xy1.setLinearCoef(ijMap[{i, k}], 1);
					Xy1.setLinearCoef(YiklMap[i][{k, l}], 1);
					Xy1.setLinearCoef(XiklMap[i][{k, l}], -1);
					
					//-Xik - Yikl +  2.Xikl <= 0
					Xy2.setLinearCoef(ijMap[{i, k}], -1);
					Xy2.setLinearCoef(YiklMap[i][{k, l}], -1);
					Xy2.setLinearCoef(XiklMap[i][{k, l}], 2);
					
					nCons += 2;
					modelPtr->add(Xy1);
					modelPtr->add(Xy2);
				}
			}

		}

		this->nIntParCons = nCons;
		cout << "Total inter partition constraints " << nCons << endl;
	}
	/*void addWriteCons() {
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
	}*/
	
	void addTransCons() {

		int nCons = 0;
		
		//for reads transaction constraints..for reads k starts from 1 but 0 should be considered for loads
		for(int k = 0; k < numParts; k++) {
			IloRange rdLd = IloRange(env, 0, TSize);
			//sum over all vertices
			for(int i = 0; i < numVertices; i++) {
				//for reads.. p < k sum all Xipk
				for(int p = 0; p < k; p++) {
					rdLd.setLinearCoef(XiklMap[i][{p, k}], 1);
				}
			}

			//for loads consider for this partition k
			for(auto elem : lpMap) {
				auto ldV = elem.second; //vector of load variables for each partition
				rdLd.setLinearCoef(ldV[k], 1);
			}
			
			nCons++;
			modelPtr->add(rdLd);
		}

		//for now considering only writes not stores, hence going from 0 to numparts - 1
		////// later when consider store groups might simply let go k upto numparts to account for stores
		for(int k = 0; k < numParts - 1; k++) {
			IloRange wrSt = IloRange(env, 0, TSize);
			
			//sum over all vertices
			for(int i = 0; i < numVertices; i++) {
				//for writes : p > k
				for(int p = k + 1; p < numParts; p++) {
					wrSt.setLinearCoef(XiklMap[i][{k, p}], 1);	
				}
			}

			/////skip for now but to add constraints for stores
			nCons++;
			modelPtr->add(wrSt);
		}
		
		this->nTransCons = nCons;
		cout << "Number of transaction constraint rows added " << nCons << endl;

	}
	
	void addLoadReuse() {
		int nCons = 0;
		for(auto elem : loadGroups) {
			int group_id = elem.first;
			//two equations for each partition for a given load group
			for(int p = 0; p < numParts; p++) {
				IloRange range1 = IloRange(env, -IloInfinity, 0);
				IloRange range2 = IloRange(env, -IloInfinity, 0);

				//negative sum all Xlp where ld is load
				for(int ld: elem.second) {
					range1.setLinearCoef(ijMap[{ld, p}], -1);
				}
				//plus Lp
				range1.setLinearCoef(lpMap[group_id][p], 1);

				//positive sum all Xlp where ld is load
				for(int ld : elem.second) {
					range2.setLinearCoef(ijMap[{ld, p}], 1);
				}

				//- num of loads * lp
				range2.setLinearCoef(lpMap[group_id][p], -1 * (int)elem.second.size());
				
				nCons += 2;
				modelPtr->add(range1);
				modelPtr->add(range2);
			}

		}
		this->nLdStCons = nCons;
	}

	//function to print all variables and row constraints
	void printVarCons() {
		cout << graph.getNumNodes() << " ";
		cout << graph.getNumEdges() << " ";

		int vals = 0; // summation of cols
		cout << nXij << " ";
		vals += nXij;
		cout << nLp <<  " ";
		vals += nLp;
		cout << nInPa << " ";
		vals += nInPa;
		cout << vals << " ";
		
		vals = 0; //summation of rows
		cout << nUniqCons << " ";
		vals += nUniqCons;
		cout << nCapCons << " ";
		vals += nCapCons;
		cout << nPrecCons << " ";
		vals += nPrecCons;
		cout << nIntParCons << " ";
		vals += nIntParCons;
		cout << nLdStCons << " ";
		vals += nLdStCons;
		cout << nTransCons << " ";
		vals += nTransCons;
		cout << vals << " ";

		cout << endl;

	}
	//print stats after iteration
	void printStats(double time, int iteration) {
		cout << time << " ";
		cout << numParts << " ";
		cout << iteration << endl;
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

			//get count of total load nodes by iterating through loadgroup map
			int loadCount = 0;
			for(auto elem : loadGroups) {
				loadCount = loadCount + elem.second.size();
				cout << "Group " << elem.first << " has " << elem.second.size() << " elements\n";
			}
			cout << "Number of load nodes " << loadCount << endl;


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

	//find partition to which this vertex is mapped to
	int getMapPart(int v) {
		for(int p = 0; p < numParts; p++) {
			double val = cplexPtr->getValue(ijMap[{v, p}]);
			if(compareEqual(val, 1) == true) {
				return p;
			}
		}
		return -1;//return -1 if no partition found, ideally should not happen
	}
	void ValidateSoln() {
		cout << "Asserting uniqueness constraints" << endl;
		//Check if vertex mapped to only one partition
		for(int i = 0; i < numVertices; i++) {
			int count = 0;
			for(int j = 0; j < numParts; j++) {
				double val = cplexPtr->getValue(ijMap[{i, j}]);
				if (compareEqual(val, 1) == true) {
					count++;
				}
			}
			assert(count == 1); //only one vertex needs to be mapped to some partition
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

		cout << "Asserting Xikls " << endl;
		
		map<int, int> writeCount; //for each partition
		map<int, int> readCount; //reads required by a partiion key
		map<int, int> outEdgesCount;//number of out edges to some subsequent partition emerging from key partition
		map<int, int> inEdgesCount; //number of incoming edges onto this key partition from some previous partition
		for(int i = 0; i < numParts; i++) {
			writeCount[i] = 0;
			readCount[i] = 0;
			outEdgesCount[i] = 0;
			inEdgesCount[i] = 0;
		}

		
		//find write, out Edges
		//also counting number of reads as it is easy to do from successors logic given below
		for(int v = 0; v < numVertices; v++) {
			int k = getMapPart(v);
			list<Node> succ;
			graph.getSuccessors(v, succ);
			map<int, bool> uniqDest; //map of unique subsequent partitions to which an out edge goes
			//uniq dest partitions because they will cause only one read on the destination partition
			bool isSomeSucc = false; //is there some successor in subsequent partition to which vertex v's output goes
			for(auto nd : succ) {
				int s_id = nd.getID(); //get successor id
				int l = getMapPart(s_id);
				if(l > k) { //this successor node is mapped to some subsequent partition
					outEdgesCount[k]++;
					isSomeSucc = true;
					uniqDest[l] = true;
					//assert that X (write) for vertex v starting at partition k and landing in partition l is true
					double val = cplexPtr->getValue(XiklMap[v][{k, l}]);
					assert(compareEqual(val, 1) == true);
				}
			}
			
			//increment respective counts of destinations having unique partitions 
			for(auto uq : uniqDest) {
				readCount[uq.first] += 1;
			}

			if(isSomeSucc) { //atleast one succ mapped in subsequent partitions
				writeCount[k]++;
			}
		}

		//in Edges
		for(int v = 0; v < numVertices; v++) {
			int k = getMapPart(v);
			list<Node> preds;
			graph.getPredecessors(v, preds);
			for(auto nd : preds) {
				int p_id = nd.getID(); //get predecessor id;
				int l = getMapPart(p_id);
				if(l < k) {// predecessor mapped to earlier partition
					inEdgesCount[k]++; //increment incoming edges count to this partition
				}
			}
		}
			

		//count distinct cluster of loads
		int loadTrans = 0;
		for(auto elem : loadGroups) {
			//count consider each group nodes
			vector<int> loadsV = elem.second;
			map<int, bool> partMapd; //maintain partition to which a load node in this group is mapped
			for(int ld : loadsV) {
				int p = getMapPart(ld);
				assert(p != -1);
				partMapd[p] = true;
			}
			
			loadTrans = loadTrans + partMapd.size();
			//size of partMapd tells how many different partitions loads are present..or how many load transactions need to be issued
		}
	
		cout << "Load transactions = " << loadTrans << endl;
		

		cout << "Write counts|Out edges counts ";
		for(int i = 0; i < numParts; i++) {
			//write trans count
			cout << writeCount[i] << "|";
			cout << outEdgesCount[i] << ", ";
			assert(writeCount[i] <= TSize);
		}
		cout << endl;

		//print out edges count
		cout << "Reads counts|In edges count ";
		for(int i = 0; i < numParts; i++) {
			cout << readCount[i] << "|";
			cout << inEdgesCount[i] << ", ";
			assert(readCount[i] <= TSize);
		}

		cout << endl;
	}
};


/*args required: dotfilename, size of partition, transaction limit, load weight*/
int main (int argc, char **argv)
{
	if(argc != 5) {
		cout << "Too few arguments, 4 expected" << endl;
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
	int loadWt = atoi(argv[4]);
	int iterations = 100;

	std::cout << std::fixed << std::setprecision(2); //set precision to 2 decimal places

	auto start = chrono::high_resolution_clock::now();
	int numParts = ceil(float(gp.getNumNodes()) / float(size)); //set initial partition size to total vertices divided by partition size
	for(int i = 1; i <= iterations; i++) {
		PartitionILP *gp1 = new PartitionILP(gp, size, trans_limit, numParts, loadWt);
		gp1->addColVars();//set objective function; define all vars
		gp1->addUniqueCons(); //add uniqness constraints
		gp1->addSizeCons(); //add size constraints
		gp1->addEdgePrec(); //add edge precedence constraints

		gp1->addInterPartCons(); //add constraints w.r.t inter partition communication
		gp1->addLoadReuse(); //add load reuse constraints
		gp1->addTransCons();
		gp1->printVarCons();
		if(gp1->solve() == true) {
			gp1->ValidateSoln();
			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
			double secs = duration/1000.0;
			gp1->printStats(secs, i);
			cout << "Solution found in iteration number " << i << " with partitions " << numParts << endl;
			break;
		}
		numParts++;
		delete gp1;
	}
	
	return 0;
}
