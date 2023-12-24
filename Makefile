CFLAGS = -O3 -Wall -Wextra
CC = g++
OUTPUT_FILE = MAtrixByMatrixMult
SOURCE_DIRECTORY = ./src
MATRIX_GENERATION_DIRECTORY = ${SOURCE_DIRECTORY}/mtx

.PHONY:
	all
	sequential
	clean

all: sequential


sequential:init random_matrix_generation
	@echo
	@echo "[INFO] Starting compilation for sequential code..." 
	${CC} ${CFLAGS} init.o random_matrix_generation.o -o ${OUTPUT_FILE}
	rm init.o
	rm random_matrix_generation.o 
	@echo "[INFO] Computation completed succesfully!"
	@echo

init:
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/init.c -o init.o 

random_matrix_generation:
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/random_generator.cpp	-o random_matrix_generation.o



clean:
	@echo "[INFO] Cleaning up..."
	rm ${OUTPUT_FILE}
	rm *.o

