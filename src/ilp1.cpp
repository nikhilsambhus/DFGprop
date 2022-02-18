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
		DAG graph;
		vector<map<pair<int, int>, int>> klMapVec;

	public:

	GraphILP(string name, DAG gp, int rsize) {
		this->graph = gp;
		this->RSize = rsize;
		this->numVertices = gp.getNumNodes();
		this->numEdges = gp.getNumEdges();

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
		int cls = glp_add_cols(lp, numVertices * numParts);
		for(int i = 0; i < (numVertices * numParts); i++) {
			glp_set_obj_coef(lp, cls + i, 0.0);
			glp_set_col_bnds(lp, cls + i, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, cls + i, GLP_BV);
		} 
		

		//add columns for each edge (parts * parts) for communication objective function
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			cls = glp_add_cols(lp, numParts * numParts);
			map<pair<int, int>, int> klMap;
			for(int k = 0; k < numParts; k++) {
				for(int l = k + 1; l < numParts; l++) {
					klMap[{k ,l}] = cls;
					glp_set_obj_coef(lp, cls, 1.0);
					glp_set_col_bnds(lp, cls, GLP_DB, 0.0, 1.0);
					glp_set_col_kind(lp, cls, GLP_BV);
					cls++;
				}
				cout << endl;
			}
			this->klMapVec.push_back(klMap);
		}


	}
	//add uniqueness constraint of mapping n vertices to p partitions
	void addUniqueCons() {
		//numVertices * numParts coef matrix for constraints
		int *ja = new int [1 + numParts];
		double *arr = new double [1 + numParts];
		
		int nr = glp_add_rows(lp, numVertices); //numVertices rows constraints
		
			//set rows bound sum euqal to 1.0
		for(int j = 0; j < numVertices; j++) {
			glp_set_row_bnds(lp, nr + j, GLP_FX, 1.0, 1.0);
		}

		//set constraints matrix to 1 for a given vertex in all partitions as only vertex needs to get mapped
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				//shift by numParts * row_id + index
				ja[j + 1] = i * numParts + j + 1;
				arr[j + 1] = 1.0;
			}
			glp_set_mat_row(lp, nr + i, numParts, ja, arr);
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

	void addCommCons() {

		int *ja = new int[1 + 3]; //3 terms in the each equation
		double *arr = new double[1 + 3]; //3 terms in the each equation
		int i = 0;
		for(list<Edge>::iterator it = graph.edgeBegin(); it != graph.edgeEnd(); it++) {
			uint32_t src = it->getSrcNodeID();
			uint32_t dest = it->getDestNodeID();
			
			for(int k = 0; k < numParts; k++) {
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
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				//cout << glp_get_col_prim(lp, i * numParts + j + 1) << endl; 
				if(glp_mip_col_val(lp, i * numParts + j + 1)) {
					cout << "Vertex id " << i  << " is mapped to " << j + 1 <<  endl;
				}
			}
		}


		/*int ncols = glp_get_num_cols(lp);
		for(int  i = 0; i < ncols; i++) {
			double ans = glp_mip_col_val(lp, i + 1);
			cout << ans << " ";
		}
		cout << endl;
		*/
		return true;

	}

};

/*Arguments required 
	graph file name
	Rsize
*/
int main(int argc, char **argv) {
	
	if(argc != 3) {
		cout << "Too few arguments, 2 expected" << endl;
		return -1;
	}
	DAG gp;
	try {
		gp = DAG(argv[1]);
	} catch(string ex) {
		cout << ex << endl;
	}
	
	GraphILP *gp1 = new GraphILP("basic", gp, atoi(argv[2]));

	
	int iterations = 10;
	while(iterations) {
		gp1->eraseProb();
		gp1->addColVars();
		gp1->addUniqueCons();
		gp1->addSizeCons();
		gp1->addEdgePrec();
		gp1->addCommCons();
		cout << "Trying next with number of partitions " << gp1->getNumParts() << endl;
		//gp1->printProb();
		if(gp1->solve() == true) {
			cout << "Converged at total number of partitions equal to " << gp1->getNumParts() << endl;
			break;
		}
		gp1->incParts();
		iterations--;
	}
	return 0;
}
