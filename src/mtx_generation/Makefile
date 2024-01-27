CFLAGS = -O3 -Wall -Wextra
CC = g++
OUTPUT_FILE = Generator
UTILS_DIRECTORY = ../utils
LOGGER_DIRECTORY = ${UTILS_DIRECTORY}/logger
FILE_IO_DIRECTORY = ${UTILS_DIRECTORY}/file_io
FORMATTER_DIRECTORY = ${UTILS_DIRECTORY}/string_formatter


--file_io:
	${CC} ${CFLAGS} -c ${FILE_IO_DIRECTORY}/file_io.c -o file_io.o 

--logger:
	${CC} ${CFLAGS} -c ${LOGGER_DIRECTORY}/logger.c -o logger.o


--random_matrix_generation: --file_io --logger
	${CC} ${CFLAGS} -c random_generator.cpp	-o random_generator.o
	${CC} ${CFLAGS} -c generator.c -o generator.o
	${CC} ${CFLAGS} generator.o random_generator.o file_io.o logger.o -o ${OUTPUT_FILE}
	@rm *.o

generate_square_matrix: --random_matrix_generation
	@./${OUTPUT_FILE} 2 2
	@./${OUTPUT_FILE} 16 16
	@./${OUTPUT_FILE} 32 32
	@./${OUTPUT_FILE} 64 64
	@./${OUTPUT_FILE} 128 128
	@./${OUTPUT_FILE} 256 256
	@./${OUTPUT_FILE} 512 512
	@./${OUTPUT_FILE} 1024 1024
	@./${OUTPUT_FILE} 2048 2048
	@./${OUTPUT_FILE} 4096 4096
	@./${OUTPUT_FILE} 5000 5000

	@rm ${OUTPUT_FILE}