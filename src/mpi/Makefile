CFLAGS = -O3 -Wall -Wextra
CC = g++
MPI_CC = mpic++
OUTPUT_FILE = MarixByMatrixMult
MATRIX_GENERATION_DIRECTORY = ../mtx_generation
UTILS_DIRECTORY = ../utils
LOGGER_DIRECTORY = ${UTILS_DIRECTORY}/logger
PROCS = 20
PROCS_GRID_ROWS = 5
PROCS_GRID_COLS = 4

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
	$(MPI_CC) $(CFLAGS) -c c/mpi_matrix_distriibution.c -o mpi_matrix_distriibution.o 
	$(MPI_CC) $(CFLAGS) -c c/mpi_computation.c -o mpi_computation.o 
	${MPI_CC} ${CFLAGS} logger.o mpi_computation.o mpi_matrix_distriibution.o time_formatter.o file_io.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

process-parallel: --mpi
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 128 128 128 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 256 256 256 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 500 500 500 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 512 512 512 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1000 1000 1000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1024 1024 1024 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 2048 2048 2048 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 3000 3000 3000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 4096 4096 4096 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 5000 5000 5000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 8192 8192 8192 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 10000 10000 10000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}

	@echo "[INFO] Square matrices mpi computation ended"

	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 256 32 256 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 500 32 500 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 512 32 512 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1000 32 1000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1024 32 1024 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 2048 32 2048 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 3000 32 3000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 4096 32 4096 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 5000 32 5000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 8192 32 8192 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 10000 32 10000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}

	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 256 64 256 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 500 64 500 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 512 64 512 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1000 64 1000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1024 64 1024 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 2048 64 2048 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 3000 64 3000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 4096 64 4096 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 5000 64 5000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 8192 64 8192 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 10000 64 10000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}

	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 256 128 256 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 500 128 500 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 512 128 512 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1000 128 1000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1024 128 1024 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 2048 128 2048 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 3000 128 3000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 4096 128 4096 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 5000 128 5000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 8192 128 8192 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 10000 128 10000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}

	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 256 156 256 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 500 156 500 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 512 156 512 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1000 156 1000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 1024 128 1024 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 2048 156 2048 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 3000 128 3000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 4096 156 4096 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 5000 156 5000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 8192 156 8192 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}
	@mpirun -np ${PROCS} ./${OUTPUT_FILE} 10000 156 10000 ${PROCS_GRID_ROWS} ${PROCS_GRID_COLS}

	@echo "[INFO] Rectangular matrices mpi computation ended"

	rm ${OUTPUT_FILE}

clean:
	@echo "[INFO] Cleaning up..."
	rm ${OUTPUT_FILE}
	rm *.o

