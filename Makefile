CFLAGS = -O3 -Wall -Wextra
CC = g++
OUTPUT_FILE = MAtrixByMatrixMult
SOURCE_DIRECTORY = ./src
MATRIX_GENERATION_DIRECTORY = ${SOURCE_DIRECTORY}/mtx
UTILS_DIRECTORY = ${SOURCE_DIRECTORY}/utils

.PHONY:
	all
	sequential
	generate_square_matrix
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
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/random_generator.cpp	-o random_generator.o
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/generator.c -o generator.o
	${CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/file_io.c -o file_io.o
	${CC} ${CFLAGS} generator.o random_generator.o file_io.o -o generate
	@mv generate src/utils/
	@rm *.o

generate_square_matrix: random_matrix_generation
	@./src/utils/generate 16 16
	@./src/utils/generate 32 32
	@./src/utils/generate 64 64
	@./src/utils/generate 128 128
	@./src/utils/generate 256 256
	@./src/utils/generate 512 512
	@./src/utils/generate 1024 1024
	@./src/utils/generate 2048 2048
	@./src/utils/generate 4096 4096

	@rm src/utils/generate

clean:
	@echo "[INFO] Cleaning up..."
	rm ${OUTPUT_FILE}
	rm *.o

