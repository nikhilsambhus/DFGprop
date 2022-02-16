#include <iostream>
#include <string>
#include "glpk.h"
#include "Edge.h"
#include "Graph.h"
using namespace std;

class GraphILP {
	private:
		glp_prob *lp; //lp object
		int numEdges, numVertices, numParts;
		int RSize; //size of partition
		DAG graph;
	public:

	GraphILP(string name, DAG gp, int rsize) {
		this->graph = gp;
		this->RSize = rsize;
		this->numVertices = gp.getNumNodes();
		this->numEdges = gp.getNumEdges();

		this->numParts = numVertices / RSize + 1; //set initial partition size to total vertices divided by partition size
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
		glp_add_cols(lp, numVertices * numParts);
		for(int i = 0; i < (numVertices * numParts); i++) {
			glp_set_obj_coef(lp, i + 1, 1.0);
			glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, i + 1, GLP_BV);
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

	bool solve() {
		//solve equations
		glp_set_obj_dir(lp, GLP_MAX);
		glp_iocp parm;
		glp_init_iocp(&parm);
		parm.presolve = GLP_ON;
		glp_simplex(lp, NULL);
		if(glp_intopt(lp, &parm) != 0) {
			return false;
		}
		double z = glp_get_obj_val(lp);
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
		cout << "Trying next with number of partitions " << gp1->getNumParts() << endl;
		if(gp1->solve() == true) {
			cout << "Converged at total number of partitions equal to " << gp1->getNumParts() << endl;
			break;
		}
		gp1->incParts();
		iterations--;
	}
	return 0;
}
