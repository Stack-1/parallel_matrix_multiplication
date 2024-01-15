CFLAGS = -O3 -Wall -Wextra
CC = g++
MPI_CC = mpic++
OUTPUT_FILE = MarixByMatrixMult
SOURCE_DIRECTORY = ./src
MATRIX_GENERATION_DIRECTORY = ${SOURCE_DIRECTORY}/mtx
UTILS_DIRECTORY = ${SOURCE_DIRECTORY}/utils
LOGGER_DIRECTORY = ${SOURCE_DIRECTORY}/logger
NUMBER_OF_PROCESS = 1

.PHONY:
	all
	sequential
	generate_square_matrix
	clean

all: sequential


--load_modules:
	@sh ./src/load_modules.sh


--file_io:
	${CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/file_io.c -o file_io.o

--logger:
	${CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o

--time_formatter:
	${CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/time_formatter.cpp -o time_formatter.o


--file_io_mpi:
	${MPI_CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/file_io.c -o file_io.o

--logger_mpi:
	${MPI_CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o

--time_formatter_mpi:
	${MPI_CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/time_formatter.cpp -o time_formatter.o




--sequential_compile: --load_modules --file_io --logger --time_formatter
	@echo
	@echo "[INFO] Starting compilation for sequential code..." 
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/sequential.c -o sequential.o 
	${CC} ${CFLAGS} file_io.o logger.o  sequential.o time_formatter.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

sequential: --sequential_compile
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


--mpi: --load_modules --logger_mpi --time_formatter_mpi
	@echo
	@echo "[INFO] Starting compilation for process parallel (MPI) code..."
	$(MPI_CC) $(CFLAGS) -c ${SOURCE_DIRECTORY}/mpi_computation.c -o mpi_computation.o 
	${MPI_CC} ${CFLAGS} logger.o mpi_computation.o time_formatter.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

process-parallel: --mpi
	@mpirun -np ${NUMBER_OF_PROCESS} ./${OUTPUT_FILE} 16 16

	rm ${OUTPUT_FILE}


--init:
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/init.c -o init.o 




--random_matrix_generation: --file_io --logger
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/random_generator.cpp	-o random_generator.o
	${CC} ${CFLAGS} -c ${MATRIX_GENERATION_DIRECTORY}/generator.c -o generator.o
	${CC} ${CFLAGS} generator.o random_generator.o file_io.o logger.o -o generate
	@mv generate src/utils/
	@rm *.o

generate_square_matrix: --random_matrix_generation
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

