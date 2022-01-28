#ifndef DFGPART_H
#define DFGPART_H
#include <vector>
#include <stack>
#include <list>
#include <bits/stdc++.h>
#include <string>
#include "Graph.h"
typedef struct partData {
	int total;
	int out;
}partData;

typedef struct partDef {
	int start, end;
	partData pData;
}partDef;


class DFGPart {
	private: 
		DAG gp;

	public:
	DFGPart(DAG grph);
	vector<vector<int>> getCombs(int timeMax, int k);
	partData getInterNds(int start, int end, vector<int32_t> &timeSt);
	int partitionDFGnP(int npart, int map_size);
	void partitionDFGVar(int map_size);
	void getBasicProfs();
	int printParts(vector<partDef>& AllParts);
};
#endif 
