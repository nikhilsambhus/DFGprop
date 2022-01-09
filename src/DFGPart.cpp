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
				splits.push_back(i);
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

int DFGPart::partitionDFGnP(int npart, int map_size) {
	DFGAnaly danl = DFGAnaly(gp);
	vector<uint32_t> topOrder = danl.topoSort();
	vector<int32_t> timeSt; 
	int applicableParts = 0;
	int32_t timeMax = danl.assignTime(topOrder, timeSt);
	vector<vector<int>> combs = getCombs(timeMax, npart - 1);
	for(vector<int> splits: combs) {
		int prev = 0;
		vector<partDef> AllParts;
		splits.push_back(timeMax);
		for(int split : splits) {
			partData pData = getInterNds(prev, split, timeSt);
			if(pData.total > map_size) {
				break;
			}
			//cout << " Pno. " << pno << "(" << prev << "-" << split << ") Total nodes: " << pData.total << " Intermediate outputs: " << pData.out; 
			partDef pd;
			pd.start = prev;
			pd.end = split;
			pd.pData = pData;
			prev = split + 1;
			AllParts.push_back(pd);
		}
		//partData pData = getInterNds(prev, timeMax, gp, timeSt);
		//cout << " Pno. " << pno << "(" << prev << "-" << timeMax << ") Total nodes: " << pData.total << " Intermediate outputs: " << pData.out; 
		if(AllParts.size() == npart) {
			int pno = 1;
			for(auto pd: AllParts) {
				cout << " Pno. " << pno << "(" << pd.start << "-" << pd.end << ") Total nodes: " << pd.pData.total << " Intermediate outputs: " << pd.pData.out; 
				pno++;
			}
			applicableParts++;
			cout << endl;
		}
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


