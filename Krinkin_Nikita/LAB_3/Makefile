SRC_DIR = ${CURDIR}/src
WORK_DIR = ${SRC_DIR}/workdir
MAT_DIR = ${SRC_DIR}/Materials


all: ${WORK_DIR}/precision.txt ${WORK_DIR}/graph.py
	python3 ${WORK_DIR}/graph.py ${WORK_DIR}/precision.txt

input_errors: ${WORK_DIR}/input_errors.o
	${WORK_DIR}/input_errors.o

${WORK_DIR}/precision.txt: ${WORK_DIR}/precision.o
	${WORK_DIR}/precision.o ${WORK_DIR}/precision.txt

${WORK_DIR}/precision.o: ${WORK_DIR}/precision.cpp ${MAT_DIR}
	g++ -g -I${MAT_DIR} ${WORK_DIR}/precision.cpp -o ${WORK_DIR}/precision.o

${WORK_DIR}/input_errors.o: ${WORK_DIR}/input_errors.cpp ${MAT_DIR}
	g++ -g -I${MAT_DIR} ${WORK_DIR}/input_errors.cpp -o ${WORK_DIR}/input_errors.o

help:
	@echo "Available options:"
	@echo "	make [all] -- calculate f(x) = 0 and show results"
	@echo "	make input_errors -- show input data sensitivity"