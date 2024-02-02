CFLAGS = -O3 -Wall -Wextra
CC = g++
MPI_CC = mpic++
OUTPUT_FILE = MarixByMatrixMult
MATRIX_GENERATION_DIRECTORY = ../mtx_generation
UTILS_DIRECTORY = ../utils
LOGGER_DIRECTORY = ${UTILS_DIRECTORY}/logger
PROCS = 1

.PHONY:
	all
	process-parallel
	clean

all: process-parallel


--load_modules:
	@sh ./load_modules_mpi.sh


--file_io_mpi:
	${MPI_CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/file_io/file_io.c -o file_io.o

--logger_mpi:
	${MPI_CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o

--time_formatter_mpi:
	${MPI_CC} ${CFLAGS} -c ${UTILS_DIRECTORY}/string_formatter/time_formatter.cpp -o time_formatter.o



--mpi: --load_modules --logger_mpi --time_formatter_mpi --file_io_mpi
	@echo
	@echo "[INFO] Starting compilation for process parallel (MPI) code..."
	$(MPI_CC) $(CFLAGS) -c c/mpi_computation.c -o mpi_computation.o 
	${MPI_CC} ${CFLAGS} logger.o mpi_computation.o time_formatter.o file_io.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

process-parallel: --mpi
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 16 16 16

	rm ${OUTPUT_FILE}


--init:
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/init.c -o init.o 


clean:
	@echo "[INFO] Cleaning up..."
	rm ${OUTPUT_FILE}
	rm *.o

