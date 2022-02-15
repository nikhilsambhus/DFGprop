#include <iostream>
#include <string>
#include "glpk.h"
using namespace std;

class GraphPart1 {
	private:
		glp_prob *lp; //lp object
		int numEdges, numVertices, numParts;
		int RSize; //size of partition
	public:

	GraphPart1(string name) {
		lp = glp_create_prob();
		glp_set_prob_name(lp, name.c_str()); //create lp object with name in constructor
	}

	~GraphPart1() {
		glp_delete_prob(lp);
	}

	void eraseProb() {
		glp_erase_prob(lp);
	}
	//set properties of graph and partitions
	void setGraphProp(int n, int e, int p, int size) {
		numEdges = e;
		numVertices = n;
		numParts = p;
		RSize = size;
	}

	//add uniqueness constraint of mapping n vertices to p partitions
	void addUniqueCons() {
		//numVertices * numParts coef matrix for constraints
		int *ja = new int [1 + numParts];
		double *arr = new double [1 + numParts];
		
		int nr = glp_add_rows(lp, numVertices); //numVertices rows constraints
		
		//add all vertices-parts mapping as cols
		glp_add_cols(lp, numVertices * numParts);
		for(int i = 0; i < (numVertices * numParts); i++) {
			glp_set_obj_coef(lp, i + 1, 1.0);
			glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, i + 1, GLP_BV);
		}

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
			//	arr[i*(numVertices * numParts) + i + j*numParts + 1] = 1.0;
			}
			glp_set_mat_row(lp, nr + i, numVertices, ja, arr);
			cout << endl;
		}




	}

	void addEdgePrec() {

		int numEdges = 2;
		int *ja = new int [1 + 2*numParts];
		double *arr = new double [1 + 2*numVertices];
		int nr = glp_add_rows(lp, numEdges); //numEdges rows constraints

		//add numEdges rows with limit to <= 0
		for(int j = 0; j < numEdges; j++) {
			glp_set_row_bnds(lp, nr + j, GLP_UP, 0.0, 0.0);
		}
	
		for(int i = 0; i < numEdges; i++) {
			//for src of edge do summation
			
			int src = 4;
			int dest = 0;
			
			if(i == 1) { 
				src = 3;
				dest = 0;
			}
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
					cout << "Vertex " << i + 1 << " is mapped to " << j + 1 <<  endl;
				}
			}
		}
		
		return true;

	}

};
int main() {
	
	int numVer = 5;
	int numEdges = 5;
	int numParts = 0;
	int Rsize = 2;

	GraphPart1 *gp1 = new GraphPart1("basic");

	do {
		gp1->eraseProb();
		numParts = numParts + 1;
		gp1->setGraphProp(numVer, numEdges, numParts, Rsize);
		gp1->addUniqueCons();
		gp1->addSizeCons();
		gp1->addEdgePrec();
		cout << "Trying with number of partitions " << numParts << endl;
	} while(gp1->solve() == false);

	return 0;
}
