#include <stdio.h>
#include <stdlib.h>
#include "rndgen.h"
#include "../utils/stackio.h"


int main(int argc, char *argv[]){
    size_t rows = 0;
    size_t cols = 0;

    if(argc != 3){
        printf("[ERROR] Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%ld", &rows) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        printf("[ERROR] Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%ld", &cols) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        printf("[ERROR] Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    float *matrix[rows];

    // Allocate memory
    for(size_t i = 0;i<rows;i++){
        matrix[i] = (float *)malloc(sizeof(float)*cols);
    }


    // Generate matrix
    generate_real_value_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix,rows,cols);

    //test();
}