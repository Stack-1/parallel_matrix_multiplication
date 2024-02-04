CFLAGS = -O3 -Wall -Wextra
CC = g++
OUTPUT_FILE = MarixByMatrixMult
SOURCE_DIRECTORY = ./c
UTILS_DIRECTORY = ../utils
LOGGER_DIRECTORY = ${UTILS_DIRECTORY}/logger
FILE_IO_DIRECTORY = ${UTILS_DIRECTORY}/file_io
FORMATTER_DIRECTORY = ${UTILS_DIRECTORY}/string_formatter
DATA_DIRECTORY = ../data/

--load_modules:
	@sh ./load_modules_sequential.sh


--file_io:
	${CC} ${CFLAGS} -c ${FILE_IO_DIRECTORY}/file_io.c -o file_io.o 

--logger:
	${CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o

--time_formatter:
	${CC} ${CFLAGS} -c ${FORMATTER_DIRECTORY}/time_formatter.cpp -o time_formatter.o


--sequential_compile: --load_modules --file_io --logger --time_formatter
	@echo
	@echo "[INFO] Starting compilation for sequential code..." 
	${CC} ${CFLAGS} -c ${SOURCE_DIRECTORY}/sequential.c -o sequential.o 
	${CC} ${CFLAGS} file_io.o logger.o  sequential.o time_formatter.o -o ${OUTPUT_FILE}
	@rm *.o
	@echo "[INFO] Compilation completed succesfully!"
	@echo

sequential: --sequential_compile
	@./${OUTPUT_FILE} 100 100 100
	@./${OUTPUT_FILE} 128 128 128
	@./${OUTPUT_FILE} 256 256 256
	@./${OUTPUT_FILE} 500 500 500
	@./${OUTPUT_FILE} 512 512 512
	@./${OUTPUT_FILE} 1000 1000 1000
	@./${OUTPUT_FILE} 1024 1024 1024
	@./${OUTPUT_FILE} 2048 2048 2048
	@./${OUTPUT_FILE} 3000 3000 3000
	@./${OUTPUT_FILE} 4096 4096 4096
	@./${OUTPUT_FILE} 5000 5000 5000
	@./${OUTPUT_FILE} 8192 8192 8192
	@./${OUTPUT_FILE} 10000 10000 10000

	@echo "[INFO] Square matrices sequential computation ended"

	@./${OUTPUT_FILE} 256 32 256
	@./${OUTPUT_FILE} 500 32 500
	@./${OUTPUT_FILE} 512 32 512
	@./${OUTPUT_FILE} 1000 32 1000
	@./${OUTPUT_FILE} 1024 32 1024
	@./${OUTPUT_FILE} 2048 32 2048
	@./${OUTPUT_FILE} 3000 32 3000
	@./${OUTPUT_FILE} 4096 32 4096
	@./${OUTPUT_FILE} 5000 32 5000
	@./${OUTPUT_FILE} 8192 32 8192
	@./${OUTPUT_FILE} 10000 32 10000

	@./${OUTPUT_FILE} 256 64 256
	@./${OUTPUT_FILE} 500 64 500
	@./${OUTPUT_FILE} 512 64 512
	@./${OUTPUT_FILE} 1000 64 1000
	@./${OUTPUT_FILE} 1024 64 1024
	@./${OUTPUT_FILE} 2048 64 2048
	@./${OUTPUT_FILE} 3000 64 3000
	@./${OUTPUT_FILE} 4096 64 4096
	@./${OUTPUT_FILE} 5000 64 5000
	@./${OUTPUT_FILE} 8192 64 8192
	@./${OUTPUT_FILE} 10000 64 10000

	@./${OUTPUT_FILE} 256 128 256
	@./${OUTPUT_FILE} 500 128 500
	@./${OUTPUT_FILE} 512 128 512
	@./${OUTPUT_FILE} 1000 128 1000
	@./${OUTPUT_FILE} 1024 128 1024
	@./${OUTPUT_FILE} 2048 128 2048
	@./${OUTPUT_FILE} 3000 128 3000
	@./${OUTPUT_FILE} 4096 128 4096
	@./${OUTPUT_FILE} 5000 128 5000
	@./${OUTPUT_FILE} 8192 128 8192
	@./${OUTPUT_FILE} 10000 128 10000

	@./${OUTPUT_FILE} 256 156 256
	@./${OUTPUT_FILE} 500 156 500
	@./${OUTPUT_FILE} 512 156 512
	@./${OUTPUT_FILE} 1000 156 1000
	@./${OUTPUT_FILE} 1024 128 1024
	@./${OUTPUT_FILE} 2048 156 2048
	@./${OUTPUT_FILE} 3000 128 3000
	@./${OUTPUT_FILE} 4096 156 4096
	@./${OUTPUT_FILE} 5000 156 5000
	@./${OUTPUT_FILE} 8192 156 8192
	@./${OUTPUT_FILE} 10000 156 10000

	@echo "[INFO] Rectangular matrices sequential computation ended"

	@rm ${OUTPUT_FILE}
