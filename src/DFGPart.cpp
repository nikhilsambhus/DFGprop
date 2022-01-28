#include "DFGPart.h"
#include "DFGAnaly.h"
DFGPart::DFGPart(DAG grph) {
	gp = grph;
}

vector<vector<int>> DFGPart::getCombs(int timeMax, int k) {
	int n = timeMax;
	std::string bitmask(k, 1); // K leading 1's
	bitmask.resize(n, 0); // N-K trailing 0's

	vector<vector<int>> combs;
	do {
		vector<int> splits;
		for (int i = 0; i < n; ++i) 
		{
			if (bitmask[i]) { 
				//std::cout << " " << i + 1;
				splits.push_back(i+1);
			}
		}
		combs.push_back(splits);
		//std::cout << std::endl;
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	
	return combs;
}

partData DFGPart::getInterNds(int start, int end, vector<int> &timeSt) {
	uint32_t total = 0, out = 0;
	partData pData;
	for(uint32_t nd = 0; nd < timeSt.size(); nd++) {
		if(timeSt[nd] >= start && timeSt[nd] <= end) {
			list<Node> succs;
			string label = gp.findNode(nd)->getLabel();
			total++;
			gp.getSuccessors(nd, succs);
			bool outYes = false;
			for(Node nd : succs) {
				if(timeSt[nd.getID()] > end) {
					outYes = true;
					break;
				}
			}
			if(outYes == true) {
				out++;
			}
		}
	}
	pData.total = total;
	pData.out = out;
	return pData;
}

int DFGPart::printParts(vector<partDef>& AllParts) {
	int pno = 1;
	int totalCoast = 0; //for total intermediate outputs
	for(auto pd: AllParts) {
		totalCoast += pd.pData.out;
		cout << " Pno. " << pno << "(" << pd.start << "-" << pd.end << ") Total nodes: " << pd.pData.total << " Intermediate outputs: " << pd.pData.out; 
		pno++;
	}
	cout << " Total coast = " << totalCoast << endl;
	return totalCoast;
}
int DFGPart::partitionDFGnP(int npart, int map_size) {
	DFGAnaly danl = DFGAnaly(gp);
	vector<uint32_t> topOrder = danl.topoSort();
	vector<int32_t> timeSt; 
	int applicableParts = 0;
	int32_t timeMax = danl.assignTime(topOrder, timeSt);
	vector<vector<int>> combs = getCombs(timeMax, npart - 1);
	vector<partDef> selectedMin;
	int minCoast = INT_MAX;
	for(vector<int> splits: combs) {
		int prev = 1;
		int totalCoast = 0; //for total intermediate outputs
		vector<partDef> AllParts;
		splits.push_back(timeMax);
		for(int split : splits) {
			partData pData = getInterNds(prev, split, timeSt);
			if(pData.total > map_size) {
				break;
			}
			partDef pd;
			pd.start = prev;
			pd.end = split;
			pd.pData = pData;
			prev = split + 1;
			AllParts.push_back(pd);
		}
		if(AllParts.size() == npart) {
			totalCoast = printParts(AllParts);
			applicableParts++;
			if(totalCoast < minCoast) {
				selectedMin = AllParts;
				minCoast = totalCoast;
			}
		}
	}

	if(applicableParts) {
		cout << "Choosen min total coast partition ";
		printParts(selectedMin);
	}
	return applicableParts;
}

void DFGPart::partitionDFGVar(int map_size) {
	for(int nparts = 2; nparts <= 4; nparts++) {
		cout << "Trying with " << nparts << endl;
		if(partitionDFGnP(nparts, map_size)) {
			cout << "Success\n" << endl;
			break;
		}
	}
}


