SRC_DIR = ${CURDIR}/src
OBJ_DIR = ${CURDIR}/objs
EXP_DIR = ${CURDIR}/examples

FN_PRECS = ${OBJ_DIR}/precs.txt
FN_COMP = ${OBJ_DIR}/comp.txt

lab: ${OBJ_DIR}/lab.o
	${OBJ_DIR}/lab.o

research: ${SRC_DIR}/graph.py ${FN_PRECS}
	python3 ${SRC_DIR}/graph.py ${FN_PRECS} ${FN_COMP}
	
${FN_PRECS}: ${OBJ_DIR}/research.o
	${OBJ_DIR}/research.o ${FN_PRECS} ${FN_COMP}

${OBJ_DIR}:
	mkdir ${OBJ_DIR}

${OBJ_DIR}/lab.o: ${OBJ_DIR} ${SRC_DIR}/lab.cpp
	@g++ -g3 ${SRC_DIR}/lab.cpp -o ${OBJ_DIR}/lab.o

${OBJ_DIR}/research.o: ${OBJ_DIR} ${SRC_DIR}/research.cpp
	@g++ -g3 ${SRC_DIR}/research.cpp -o ${OBJ_DIR}/research.o

clear:
	rm ${OBJ_DIR}/*
