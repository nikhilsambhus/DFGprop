#include <stdint.h>
#include <iostream>
#include <tuple>
#include <map>
using namespace std;
#define STRIDE_MIN 8
map<string, int> nodeWts = {
{"load", 1},
{"store", 1},
{"fadd", 1},
{"fmul", 1},
{"fdiv", 1},
{"call max", 1},
{"call min", 1}
};

#include "Graph.h"
#include "GraphUtils.h"

