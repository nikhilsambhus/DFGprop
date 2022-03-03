CC=g++ -Wall -g
SRC=./src
INC=./include
GR_LIB=/home/sambhusn/llvm-project/llvm/lib/Transforms/LLVMAssngs/cgramap/Graph/lib/libgraph.a -lm -lglpk

DFGUtils.out : ${SRC}/DFGUtils.cpp ${INC}/DFGUtils.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/DFGUtils.cpp -I ${INC} ${GR_LIB} -o DFGUtils.out

DFGAnaly.o: ${SRC}/DFGAnaly.cpp ${INC}/DFGAnaly.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/DFGAnaly.cpp -I ${INC} ${GR_LIB} -c 

DFGPart.o : ${SRC}/DFGPart.cpp ${INC}/DFGPart.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/DFGPart.cpp -I ${INC} ${GR_LIB} -c

main.o : ${SRC}/main.cpp ${INC}/* ${GR_LIB} DFGPart.o DFGAnaly.o
	${CC} -std=c++11 ${SRC}/main.cpp -I ${INC} ${GR_LIB} DFGPart.o DFGAnaly.o -o main.o

Normalize.out : ${SRC}/Normalize.cpp ${INC}/DFGUtils.h ${GR_LIB}
	${CC} -std=c++11 ${SRC}/Normalize.cpp -I ${INC} ${GR_LIB} -o Normalize.out

ConvLoadSan.out : ${SRC}/ConvLoadSan.cpp ${INC}/* ${GR_LIB}
	${CC} -std=c++11 ${SRC}/ConvLoadSan.cpp -I ${INC} ${GR_LIB} -o ConvLoadSan.out

ilp1.o : ${SRC}/ilp1.cpp ${INC}/* ${GR_LIB}
	${CC} -std=c++11 ${SRC}/ilp1.cpp -I ${INC} ${GR_LIB} -o ilp1.o

dotconv1.o : dotconv1.cpp ${INC}/* ${GR_LIB}
	${CC} -std=c++11 dotconv1.cpp -I ${INC} ${GR_LIB} -o dotconv1.o
