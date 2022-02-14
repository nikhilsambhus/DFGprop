#include <iostream>
#include <string>
#include "glpk.h"
using namespace std;

class GraphPart1 {
	private:
		glp_prob *lp; //lp object
		int numEdges, numVertices, numParts;
		int RSize = 2; //size of partition
	public:

	GraphPart1(string name) {
		lp = glp_create_prob();
		
		glp_set_prob_name(lp, name.c_str()); //create lp object with name in constructor
	}

	//set properties of graph and partitions
	void setGraphProp(int n, int e, int p) {
		numEdges = e;
		numVertices = n;
		numParts = p;
	}

	//add uniqueness constraint of mapping n vertices to p partitions
	void addUniqueCons1() {
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

	void addSizeCons2() {
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

	void solve() {

		//solve equations
		glp_set_obj_dir(lp, GLP_MAX);
		glp_iocp parm;
		glp_init_iocp(&parm);
		parm.presolve = GLP_ON;
		glp_simplex(lp, NULL);
		glp_intopt(lp, &parm);
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
		

	}

};
int main() {
	
	GraphPart1 gp1("basic");
	gp1.setGraphProp(5, 5, 3);
	gp1.addUniqueCons1();
	gp1.addSizeCons2();
	gp1.solve();
	
	return 0;
}
