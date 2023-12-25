#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "mtx/rndgen.h"


void compute_sequential_matrix_by_matrix_multiplication(float **matrix_A,float **matrix_B,float **matrix_C,size_t rows,size_t cols){
	for (size_t i=0;i<rows;i++) 
        for(size_t j=0;j<cols;j++) {
            for(size_t k=0;k<cols;k++){    
                matrix_C[i][j]=matrix_C[i][j]+matrix_A[i][k]*matrix_B[k][j];
			}
    }
}



int main(int argc, char *argv[]){
	size_t rows = 0;
    size_t cols = 0;
	float temp = 0;     
    clock_t start, end;
	double cpu_time_used;

    
	if(argc != 3){
        printf("[ERROR] Program should be called like ./<elf file name> <# rows> <# cols>\n");
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

	float *matrix_A[rows];
	float *matrix_B[rows];
	float *matrix_C[rows];
	// TODO: See if it is more efficent to keep matrix as 1D dimencional array


    puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");

	for(int i = 0;i<cols;i++){
		matrix_A[i] = (float *)malloc(cols*sizeof(float));
		matrix_B[i] = (float *)malloc(cols*sizeof(float));
		matrix_C[i] = (float *)malloc(cols*sizeof(float));
	}

	puts("[INFO] Matrix memory correctly allocated on HOST\n");


	//TODO: Get data from file
	generate_real_value_matrix(matrix_A,rows,cols);
	generate_real_value_matrix(matrix_B,rows,cols);
	generate_real_value_matrix(matrix_C,rows,cols);

	puts("[INFO] Matrix values generation correctly finished on HOST\n");

	if(cols <= 16 && rows <= 16){
		puts("Matrix A:");
		print_matrix(matrix_A,rows,cols);
	
		puts("Matrix B:");
		print_matrix(matrix_B,rows,cols);
	
		puts("Matrix C:");
		print_matrix(matrix_C,rows,cols);
	}

	start = clock();
	compute_sequential_matrix_by_matrix_multiplication(matrix_A,matrix_B,matrix_C,rows,cols);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("[INFO] Sequential computation ended in %f s\n",cpu_time_used);
	
	if(cols <= 16 && rows <= 16){
		puts("New matrix C:");
		print_matrix(matrix_C,rows,cols);
	}

	puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");

	return 0;
}