SRC_DIR = ${CURDIR}/src
OBJ_DIR = ${CURDIR}/objs
EXP_DIR = ${CURDIR}/examples

FN 		= ${OBJ_DIR}/check.txt

lab: ${OBJ_DIR}/lab.o
	${OBJ_DIR}/lab.o

research: ${SRC_DIR}/graph.py
	python3 ${SRC_DIR}/graph.py ${FN}
	
${FN}: ${OBJ_DIR}/research.o
	${OBJ_DIR}/research.o ${FN}

${OBJ_DIR}:
	mkdir ${OBJ_DIR}

${OBJ_DIR}/lab.o: ${OBJ_DIR} ${SRC_DIR}/lab.c
	@gcc -g3 ${SRC_DIR}/lab.c -o ${OBJ_DIR}/lab.o -Wformat=0

${OBJ_DIR}/research.o: ${OBJ_DIR} ${SRC_DIR}/research.c
	@gcc -g3 ${SRC_DIR}/research.c -o ${OBJ_DIR}/research.o -Wformat=0

clear:
	rm ${OBJ_DIR}/*

help:
	@echo "Available options:"
	@echo "	make [lab] -- solve matrix via Gauss"
	@echo " make speed-comp -- compare speed with np.linalg"
	@echo "	make input_errors -- compare precision with np.linalg"
