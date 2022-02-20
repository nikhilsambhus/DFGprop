#include <iostream>
#include <string>
#include "glpk.h"
#include "Edge.h"
#include "Graph.h"
#include <cmath>
#include <bits/stdc++.h>
using namespace std;

class GraphILP {
	private:
		glp_prob *lp; //lp object
		int numEdges, numVertices, numParts;
		int RSize; //size of partition
		int TSize;//size of transaction
		DAG graph;//input graph
		vector<map<pair<int, int>, int>> klMapVec;//map of kl values for cross partition edges
		map<pair<int, int>, int> ijMap; //map of ij values

	public:

	GraphILP(string name, DAG gp, int rsize, int tsize) {
		this->graph = gp;
		this->RSize = rsize;
		this->TSize = tsize;
		this->numVertices = gp.getNumNodes();
		this->numEdges = gp.getNumEdges();
		cout << "Num Edges " << numEdges << " Num vertices " << numVertices << endl;

		this->numParts = ceil(float(numVertices) / float(RSize)); //set initial partition size to total vertices divided by partition size
		//numParts = 1;

		this->lp = glp_create_prob();
		glp_set_prob_name(this->lp, name.c_str()); //create lp object with name in constructor
	}

	~GraphILP() {
		glp_delete_prob(lp);
	}

	void eraseProb() {
		glp_erase_prob(lp);
	}
	//increment number of partitions
	void incParts() {
		numParts++;
	}

	int getNumParts() {
		return numParts;
	}


	void addColVars() {
		//add all vertices-parts mapping as cols

		this->ijMap.clear();
		int cls = glp_add_cols(lp, numVertices * numParts);
		int count = 0;
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				ijMap[{i, j}] = cls;
				glp_set_obj_coef(lp, cls, 0.0);
				glp_set_col_bnds(lp, cls, GLP_DB, 0.0, 1.0);
				glp_set_col_kind(lp, cls, GLP_BV);
				cls++;
				count++;
			} 
		}

		cout << "Xij variable added count = " << count << endl;

		this->klMapVec.clear();
		count = 0;
		//add columns for each edge (parts * parts) for communication objective function
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			map<pair<int, int>, int> klMap;
			for(int k = 0; k < numParts; k++) {
				for(int l = k; l < numParts; l++) {
					cls = glp_add_cols(lp, 1);
					klMap[{k ,l}] = cls;
					if(k != l) {// do not add edges in same partition
						glp_set_obj_coef(lp, cls, 1.0);
					}
					glp_set_col_bnds(lp, cls, GLP_DB, 0.0, 1.0);
					glp_set_col_kind(lp, cls, GLP_BV);
					count++;
				}
			}
			this->klMapVec.push_back(klMap);
		}
		cout << "Xij^kl variable added count = " << count << endl;


	}

	void addTransCons() {
		
		int nCons = 0;
		//(l > k) X^kl_ij < T
		for(int k = 0; k < numParts - 1; k++) {
			vector<double> arr;
			vector<int> ja;
			arr.push_back(0);
			ja.push_back(0);//0th element not used
			for(int l = k + 1; l < numParts; l++) {
				for(auto klMap: this->klMapVec) {
					ja.push_back(klMap[{k, l}]);
					arr.push_back(1);
				}
			}

			int nr = glp_add_rows(lp, 1); //add one for k's outgoing edges
			glp_set_mat_row(lp, nr, arr.size() - 1, ja.data(), arr.data());
			glp_set_row_bnds(lp, nr, GLP_UP, 0.0, TSize);
			nCons++;
		}
		
		//(k < l) X^kl_ij < T
		for(int l = 1; l < numParts; l++) {
			vector<double> arr;
			vector<int> ja;
			arr.push_back(0);
			ja.push_back(0);//0th element not used
			for(int k = 0; k < l; k++) {
				for(auto klMap: this->klMapVec) {
					ja.push_back(klMap[{k, l}]);
					arr.push_back(1);
				}
			}

			int nr = glp_add_rows(lp, 1); //add one for k's outgoing edges
			glp_set_mat_row(lp, nr, arr.size() - 1, ja.data(), arr.data());
			glp_set_row_bnds(lp, nr, GLP_UP, 0.0, TSize);
			nCons++;

		}

		cout << "Number of transaction constraint rows added " << nCons << endl;

	}

	//add uniqueness constraint of mapping n vertices to p partitions
	void addUniqueCons() {
		//numVertices * numParts coef matrix for constraints
		int *ja = new int [1 + numParts];
		double *arr = new double [1 + numParts];
		int nCons = 0;
		
		int nr = glp_add_rows(lp, numVertices); //numVertices rows constraints
		
			//set rows bound sum euqal to 1.0
		for(int j = 0; j < numVertices; j++) {
			glp_set_row_bnds(lp, nr + j, GLP_FX, 1.0, 1.0);
		}

		//set constraints matrix to 1 for a given vertex in all partitions as only vertex needs to get mapped
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				//shift by numParts * row_id + index
				ja[j + 1] = ijMap[{i, j}];
				arr[j + 1] = 1.0;
			}
			glp_set_mat_row(lp, nr + i, numParts, ja, arr);
			nCons++;
		}

		cout << "Number of uniquness constraints rows added " << nCons << endl;

	}

	//validate uniqueness constraints
	void ValidateUniq() {
		cout << "Asserting uniqueness constraints " << endl;
		for(uint32_t i = 0; i < this->graph.getNumNodes(); i++) {
			int count = 0;
			for(int j = 0; j < numParts; j++) {
				if(glp_mip_col_val(lp, ijMap[{i, j}]) == 1) {
					count++;
				}
			}
			assert(count == 1);
		}
	}
	void addSizeCons() {
		int *ja = new int [1 + numVertices];
		double *arr = new double [1 + numVertices];

		int nr = glp_add_rows(lp, numParts); //numParts rows constraints
		for(int j = 0; j < numParts; j++) {
			glp_set_row_bnds(lp, nr + j, GLP_UP, 0.0, RSize);
		}

		//set constraints to less than Rsize for all vertices of a partition
		for(int i = 0; i < numParts; i++) {
			for(int j = 0; j < numVertices; j++) {
				//row + shift by i + partition number id
				ja[j + 1] = i + j*numParts + 1;
				arr[j + 1] = 1.0;
			}
			glp_set_mat_row(lp, nr + i, numVertices, ja, arr);
		}




	}

	void ValidateSize() {
		cout << "Asserting size constraints" << endl;
		for(int i = 0; i < numParts; i++) {
			int count = 0;
			for(uint32_t j = 0; j < graph.getNumNodes(); j++) {
				if(glp_mip_col_val(lp, ijMap[{j, i}]) == 1) {
					count++; //add if vertex present in this partition
				}
			}

			assert(count <= RSize);
		}
	}

	void addCommCons() {

		int *ja = new int[1 + 3]; //3 terms in the each equation
		double *arr = new double[1 + 3]; //3 terms in the each equation
		int i = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src = it->getSrcNodeID();
			uint32_t dest = it->getDestNodeID();
			
			for(int k = 0; k < numParts - 1; k++) {
				for(int l = k + 1; l < numParts; l++) {
					int nr = glp_add_rows(lp, 2);

					//first equation
					//Xi_k
					ja[1] = (src * numParts) + k + 1;
					//Xj_l
					ja[2] = (dest * numParts) + l + 1;
					//Xi_j^k_l
					ja[3] = this->klMapVec[i][{k ,l}];
					//ja[3] = (numVertices * numParts) + (i * numParts * numParts) +  (k * numParts) + l + 1;
					arr[1] = arr[2] = 1.0;
					arr[3] = -1.0;

					glp_set_mat_row(lp, nr, 3, ja, arr);
					glp_set_row_bnds(lp, nr, GLP_UP, 0.0, 1.0);


					//second equation
					//format same, coefs different
					arr[1] = arr[2] = -1.0;
					arr[3] = 2.0;

					glp_set_mat_row(lp, nr + 1, 3, ja, arr);
					glp_set_row_bnds(lp, nr + 1, GLP_UP, 0.0, 0.0);

				}

			}
			i++;
		}



	}

	void addCommCons2() {
		int i = 0;

		int nCons = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src_i = it->getSrcNodeID();
			uint32_t dest_j = it->getDestNodeID();

			//first constraint sum (p < l) Xi_j^p_l = Xj_l
			for(int l = 0; l < numParts; l++) {
				vector<int> ja;
				vector<double> arr;
				ja.push_back(0);
				arr.push_back(0); //0th not used
				for(int p = 0; p <= l; p++) {
					ja.push_back(klMapVec[i][{p, l}]);
					arr.push_back(1);
				}
				ja.push_back(ijMap[{dest_j, l}]);
				arr.push_back(-1);
				int nr = glp_add_rows(lp, 1);
				nCons++;
				glp_set_mat_row(lp, nr, arr.size() - 1, ja.data(), arr.data());
				glp_set_row_bnds(lp, nr, GLP_FX, 0.0, 0.0);
			}

			//second constraint sum (p > k) Xi_j^k_p = Xi_k
			for(int k = 0; k < numParts; k++) {
				vector<int> ja;
				vector<double> arr;
				ja.push_back(0);
				arr.push_back(0); //0th not used
				for(int p = k; p < numParts; p++) {
					ja.push_back(klMapVec[i][{k, p}]);
					arr.push_back(1);
				}
				ja.push_back(ijMap[{src_i, k}]);
				arr.push_back(-1);
				int nr = glp_add_rows(lp, 1);
				nCons++;
				glp_set_mat_row(lp, nr, arr.size() - 1, ja.data(), arr.data());
				glp_set_row_bnds(lp, nr, GLP_UP, 0.0, 0.0);
			}

			i++;
		}

		cout << "Comm 2 Constraints added = " << nCons << endl;
	}

	void addEdgePrec() {

		int *ja = new int [1 + 2*numParts];
		double *arr = new double [1 + 2*numVertices];
		int nr = glp_add_rows(lp, numEdges); //numEdges rows constraints

		//add numEdges rows with limit to <= 0
		int i = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			glp_set_row_bnds(lp, nr + i, GLP_UP, 0.0, 0.0);
			i++;
		}


		i = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			//for src of edge do summation

			uint32_t src = it->getSrcNodeID();
			uint32_t dest = it->getDestNodeID();

			//partition number * Xvertex_partition
			for(int i = 0; i < numParts; i++) {
				ja[i + 1] = (src*numParts) + i + 1;
				arr[i + 1] = (i + 1) * 1;
			}

			//-partition number * Xvertex_partition
			for(int i = 0; i < numParts; i++) {
				ja[numParts + i + 1] = (dest*numParts) + i + 1;
				arr[numParts + i + 1] = -(i + 1) * 1;
			}

			glp_set_mat_row(lp, nr + i, 2*numParts, ja, arr);
			i++;
		}
		
	}

	void ValidatePrecTrans() {
		cout << "Asserting edge precedence; transaction limits" << endl;
		map<int, int> inPartCounts;//store incoming edges onto this partition 
		map<int, int> outPartCounts; //store outgoing edges from this partition
		for(int i = 0; i < numParts; i++) {
			inPartCounts[i] = 0;
			outPartCounts[i] = 0;
		}
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src = it->getSrcNodeID();
			uint32_t dest = it->getDestNodeID();
			//find source partition
			int srcPart = -1;
			
			for(int j = 0; j < numParts; j++) {
				if(glp_mip_col_val(lp, ijMap[{src, j}]) == 1) {
					srcPart = j;
					break;
				}
			}
			assert(srcPart != -1); //mapping should be found

			//find destination partition
			int destPart = -1;
			
			for(int j = 0; j < numParts; j++) {
				if(glp_mip_col_val(lp, ijMap[{dest, j}]) == 1) {
					destPart = j;
					break;
				}
			}
			assert(srcPart != -1); //mapping should be found
			assert(srcPart <= destPart); //partition of source should be less than destination

			if(srcPart < destPart) { //for edge not present in the same partition
				outPartCounts[srcPart]++; //increment appropriate transaction counts
				inPartCounts[destPart]++;
			}

		}

		//Assert cross partition counts are meet
		for(int i = 0; i < numParts; i++) {
			assert(inPartCounts[i] <= TSize);
			assert(outPartCounts[i] <= TSize);
		}
	}
	
	void printProb() {
		int ja[1024];
		double arr[1024];
		int nrows = glp_get_num_rows(lp);
		for(int i = 0; i < nrows; i++) {
			int nelems = glp_get_mat_row(lp, i + 1, ja, arr);
			for(int j = 0; j < nelems; j++) {
				cout << arr[j + 1] << "*" << ja[j + 1] << " "; 
			}
			double ub = glp_get_row_ub(lp, i + 1);
			cout << "\t" << ub;
			cout << endl;
		}

		int ncols = glp_get_num_cols(lp);
		for(int  i = 0; i < ncols; i++) {
			double coef = glp_get_obj_coef(lp, i + 1);
			cout << coef << " ";
		}
		cout << endl;
	}

	bool solve() {
		//solve equations
		glp_set_obj_dir(lp, GLP_MIN);
		glp_iocp parm;
		glp_init_iocp(&parm);
		parm.presolve = GLP_ON;
		glp_simplex(lp, NULL);
		if(glp_intopt(lp, &parm) != 0) {
			return false;
		}
		double z = glp_mip_obj_val(lp);
		cout << "Objective function output "<< z << endl;

		//print values
		int countVP = 0;
		cout << "Mappings ";
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				//cout << glp_get_col_prim(lp, i * numParts + j + 1) << endl; 
				if(glp_mip_col_val(lp, ijMap[{i, j}])) {
					//cout << "Vertex id " << i  << " is mapped to " << j + 1 <<  endl;
					cout << i <<  "|" << j << " ";
					countVP++;
				}
			}
		}
		
		cout << endl;

		/*int ncols = glp_get_num_cols(lp);
		cout << "Output vals of kl for each edge on newline in form part1:part2=bool present value" << endl;
		for(auto klMap : this->klMapVec) {
			for(int l = 0; l < numParts; l++) {
				for(int k = l + 1; k < numParts; k++) {
					cout << l << ":" << k << "=" << glp_mip_col_val(lp, klMap[{l, k}]) << " ";
				}

			}
			cout << endl;
		}*/
		
		return (countVP == numVertices); //return success if all vertices mapped to some partition

	}

};

/*Arguments required 
	graph file name
	Rsize
*/
int main(int argc, char **argv) {
	
	if(argc != 4) {
		cout << "Too few arguments, 3 expected" << endl;
		return -1;
	}
	DAG gp;
	try {
		gp = DAG(argv[1]);
	} catch(string ex) {
		cout << ex << endl;
	}
	
	GraphILP *gp1 = new GraphILP("basic", gp, atoi(argv[2]), atoi(argv[3]));

	
	int iterations = 10;
	while(iterations) {
		cout << "Trying next with number of partitions " << gp1->getNumParts() << endl;
		gp1->eraseProb();
		gp1->addColVars();
		gp1->addUniqueCons();
		gp1->addSizeCons();
		gp1->addEdgePrec();
		gp1->addCommCons2();
		gp1->addTransCons();
		//gp1->printProb();
		if(gp1->solve() == true) {
			cout << "Converged at total number of partitions equal to " << gp1->getNumParts() << endl;
			gp1->ValidateUniq();
			gp1->ValidateSize();
			gp1->ValidatePrecTrans();
			break;
		}
		gp1->incParts();
		iterations--;
	}
	return 0;
}
