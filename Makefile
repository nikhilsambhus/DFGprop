CC=g++
SRC=./src
INC=./include
GR_LIB=/home/sambhusn/llvm-project/llvm/lib/Transforms/LLVMAssngs/cgramap/Graph/lib/libgraph.a

DFGUtils.out : ${SRC}/DFGUtils.cpp ${INC}/DFGUtils.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/DFGUtils.cpp -I ${INC} ${GR_LIB} -o DFGUtils.out

Normalize.out : ${SRC}/Normalize.cpp ${INC}/DFGUtils.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/Normalize.cpp -I ${INC} ${GR_LIB} -o Normalize.out
