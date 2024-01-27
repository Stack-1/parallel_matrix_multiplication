#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rndgen.h"
#include "../utils/file_io/stackio.h"
#include "../utils/logger/logger.h"

int main(int argc, char *argv[]){
    int rows = 0;
    int cols = 0;
    char *matrix_A_name = (char *)"A";
    char *matrix_B_name = (char *)"B";
    char *matrix_C_name = (char *)"C";

    if(argc != 3){
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%d", &rows) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%d", &cols) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

	puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");



    float *matrix;

    // Allocate memory
    matrix = (float *)malloc(sizeof(float)*rows*cols);
    


    puts("MATRIX A:");
    // Generate matrix
    generate_real_value_matrix(matrix,rows,cols);
    write_matrix_to_file(matrix,rows,cols,matrix_A_name);


    puts("MATRIX B:");
    generate_real_value_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix,rows,cols,matrix_B_name);



    puts("MATRIX C:");
    generate_real_value_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix,rows,cols,matrix_C_name);
    //test();

    puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");
}