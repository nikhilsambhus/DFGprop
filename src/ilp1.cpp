#include <iostream>
#include <string>
#include "glpk.h"
using namespace std;

class GraphPart1 {
	private:
		glp_prob *lp; //lp object
		int numEdges, numVertices, numParts;
		int RSize = 3; //size of partition
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
		int *ia = new int [1 + numVertices * numVertices * numParts];
		int *ja = new int [1 + numVertices * numVertices * numParts];
		double *arr = new double [1 + numVertices * numVertices * numParts];
		
		glp_add_rows(lp, numVertices); //numVertices rows constraints
		for(int j = 0; j < numVertices; j++) {
			glp_set_row_bnds(lp, j + 1, GLP_FX, 1.0, 1.0);
		}
	
		//initialize to 0 all coefs
		int count = 1;
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numVertices * numParts; j++) {
				ia[count] = i + 1;
				ja[count] = j + 1;
				arr[count] = 0.0;
				count++;
			}
		}


		//set constraints to 1
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numParts; j++) {
				//row + shift + index
				arr[i*(numVertices * numParts) + i*numParts + j + 1] = 1.0;
			}
		}

		//debug
		for(int i = 0; i < numVertices; i++) {
			for(int j = 0; j < numVertices * numParts; j++) {
				cout << arr[i*(numVertices * numParts) + j + 1] << " ";
			}
			cout << endl;
		}

		//add all vertices-parts mapping as cols
		glp_add_cols(lp, numVertices * numParts);
		for(int i = 0; i < (numVertices * numParts); i++) {
			glp_set_obj_coef(lp, i + 1, 1.0);
			glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
			glp_set_col_kind(lp, i + 1, GLP_BV);
		}

		glp_load_matrix(lp, numVertices * numVertices * numParts, ia, ja, arr);
	}

	void addSizeCons2() {
		int *ia = new int [1 + numParts * numVertices * numParts];
		int *ja = new int [1 + numParts * numVertices * numParts];
		double *arr = new double [1 + numParts * numVertices * numParts];

		glp_add_rows(lp, numParts); //numParts rows constraints
		for(int j = 0; j < numParts; j++) {
			glp_set_row_bnds(lp, numVertices + j + 1, GLP_UP, 0.0, RSize);
		}

		//initialize to 0 all coefs
		int count = 1;
		for(int i = 0; i < numParts; i++) {
			for(int j = 0; j < numVertices * numParts; j++) {
				ia[count] = numVertices + i + 1;
				ja[count] = j + 1;
				arr[count] = 0.0;
				count++;
			}
		}
		
		//set constraints to less than Rsize
		for(int i = 0; i < numParts; i++) {
			for(int j = 0; j < numVertices; j++) {
				//row + shift by i + partition number id
				arr[i*(numVertices * numParts) + i + j*numParts + 1] = 1.0;
			}
		}

		//debug
		for(int i = 0; i < numParts; i++) {
			for(int j = 0; j < numVertices * numParts; j++) {
				cout << arr[i*(numVertices * numParts) + j + 1] << " ";
			}
			cout << endl;
		}

		glp_load_matrix(lp, numParts * numVertices * numParts, ia, ja, arr);


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
		cout << z << endl;

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
