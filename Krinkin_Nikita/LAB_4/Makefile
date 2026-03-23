SRC_DIR = ${CURDIR}/src
WORK_DIR = ${SRC_DIR}/workdir
MAT_DIR = ${SRC_DIR}/Materials

lab: ${WORK_DIR}/lab
	${WORK_DIR}/lab

precision: ${WORK_DIR}/precision.txt ${WORK_DIR}/graph.py
	python3 ${WORK_DIR}/graph.py ${WORK_DIR}/precision.txt

input_errors: ${WORK_DIR}/input_errors.o
	${WORK_DIR}/input_errors.o

${WORK_DIR}/precision.txt: ${WORK_DIR}/precision.o
	${WORK_DIR}/precision.o ${WORK_DIR}/precision.txt


${WORK_DIR}/lab: ${WORK_DIR}/lab.cpp ${MAT_DIR}
	g++ -g -I${MAT_DIR} ${WORK_DIR}/lab.cpp -o ${WORK_DIR}/lab

${WORK_DIR}/precision.o: ${WORK_DIR}/precision.cpp ${MAT_DIR}
	g++ -g -I${MAT_DIR} ${WORK_DIR}/precision.cpp -o ${WORK_DIR}/precision.o

${WORK_DIR}/input_errors.o: ${WORK_DIR}/input_errors.cpp ${MAT_DIR}
	g++ -g -I${MAT_DIR} ${WORK_DIR}/input_errors.cpp -o ${WORK_DIR}/input_errors.o

help:
	@echo "Available options:"
	@echo "	make [lab] -- calculate f(x) = 0 and show results"
	@echo " make precision -- calculate f(x) for different epsilon for given delta"
	@echo "	make input_errors -- show input data sensitivity"