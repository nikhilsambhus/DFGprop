#include <stdint.h>
#include <iostream>
#include <tuple>
#include <map>
using namespace std;
#define STRIDE_MIN 8
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
typedef struct partData {
	int total;
	int out;
}partData;

typedef struct partDef {
	int start, end;
	partData pData;
}partDef;
#include "Graph.h"
#include "GraphUtils.h"

