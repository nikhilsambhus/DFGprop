#ifndef DFGANALY_H
#define DFGANALY_H
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>
#include "Graph.h"
#include "GraphUtils.h"
#define STRIDE_MIN 8
using namespace std;
class DFGAnaly {
	private:
		DAG gp;
		map<string, int> nodeWts = {
			{"load", 1},
			{"store", 1},
			{"fcmp", 1},
			{"br", 1},
			{"phi", 1},
			{"sext", 1},
			{"getelementptr", 1},
			{"fneg", 1},
			{"call exp", 1},
			{"fadd", 1},
			{"fsub", 1},
			{"add", 1},
			{"mul", 1},
			{"fmul", 1},
			{"fdiv", 1},
			{"call_max", 1},
			{"call_min", 1}
		};


	public:
		DFGAnaly(DAG grph);
		DFGAnaly();
		void topoSortHelper(uint32_t node, vector<bool> &visited, stack<uint32_t> &st);
		vector<uint32_t> topoSort();
		uint32_t assignTime(vector<uint32_t> &topOrder, vector<int32_t> &timeSt);
		double getParallelism();
		uint32_t criticalPathLen();
		void getBasicProps();
};
#endif
