CFLAGS = -O3 -Wall -Wextra
CC = g++
OUTPUT_FILE = MAtrixByMatrixMult
SOURCE_DIRECTORY = ./src
MATRIX_GENERATION_DIRECTORY = ${SOURCE_DIRECTORY}/mtx
UTILS_DIRECTORY = ${SOURCE_DIRECTORY}/utils
LOGGER_DIRECTORY = ${SOURCE_DIRECTORY}/logger


.PHONY:
	all
	sequential
	generate_square_matrix
	clean

all: sequential


file_io:
	${CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/file_io.c -o file_io.o

logger:
	${CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o

sequential_computation:
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/matrix_computation.c -o matrix_computation.o

time_formatter:
	${CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/time_formatter.cpp -o time_formatter.o



sequential_compile: file_io logger sequential_computation time_formatter
	@echo
	@echo "[INFO] Starting compilation for sequential code..." 
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/sequential.c -o sequential.o 
	${CC} ${CFLAGS} file_io.o logger.o  sequential.o matrix_computation.o time_formatter.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

sequential: sequential_compile
	@./${OUTPUT_FILE} 16 16
	@./${OUTPUT_FILE} 32 32
	@./${OUTPUT_FILE} 64 64
	@./${OUTPUT_FILE} 128 128
	@./${OUTPUT_FILE} 256 256
	@./${OUTPUT_FILE} 512 512
	@./${OUTPUT_FILE} 1024 1024
	@./${OUTPUT_FILE} 2048 2048
	@./${OUTPUT_FILE} 4096 4096

	@rm ${OUTPUT_FILE}


init:
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/init.c -o init.o 




random_matrix_generation: file_io logger
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/random_generator.cpp	-o random_generator.o
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/generator.c -o generator.o
	${CC} ${CFLAGS} generator.o random_generator.o file_io.o logger.o -o generate
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

