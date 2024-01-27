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
	@./${OUTPUT_FILE} 2 2 2 
	@./${OUTPUT_FILE} 16 16 16 
	@./${OUTPUT_FILE} 32 32 32 
	@./${OUTPUT_FILE} 64 64 64 
	@./${OUTPUT_FILE} 128 128 128 
	@./${OUTPUT_FILE} 256 256 256 
	@./${OUTPUT_FILE} 512 512 512 
	@./${OUTPUT_FILE} 1024 1024 1024 
	@./${OUTPUT_FILE} 2048 2048 2048 
	@./${OUTPUT_FILE} 4096 4096 4096 
	@./${OUTPUT_FILE} 5000 5000 5000 


	@rm ${OUTPUT_FILE}
