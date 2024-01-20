#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "mtx/rndgen.h"
#include "utils/stackio.h"
#include "logger/logger.h"

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
    clock_t read_start, read_end, computation_start, computation_end;
	double cpu_time_used;
	double total_time = 0.0;
	char logger_message[LOG_MESSAGE_SIZE];
	char formatted_string[64];

    
	if(argc != 3){
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%ld", &rows) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%ld", &cols) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols>");
        exit(EXIT_FAILURE);
    }

	float *matrix_A[rows];
	float *matrix_B[rows];
	float *matrix_C[rows];
	// TODO: See if it is more efficent to keep matrix as 1D dimencional array


    puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");


	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Starting sequential computation for matrix %ld x %ld",rows,cols);
	logger_info(logger_message);

	for(size_t i = 0;i<cols;i++){
		matrix_A[i] = (float *)malloc(cols*sizeof(float));
		matrix_B[i] = (float *)malloc(cols*sizeof(float));
		matrix_C[i] = (float *)malloc(cols*sizeof(float));
	}

	logger_info("Matrix memory correctly allocated on HOST");


	// Getting data from previously generated file
	read_start = clock();
	read_matrix_from_file(matrix_A,rows,cols,(char *)"A");
	read_matrix_from_file(matrix_B,rows,cols,(char *)"B");
	read_matrix_from_file(matrix_C,rows,cols,(char *)"C");
	read_end = clock();
	cpu_time_used = ((double) (read_end - read_start)) / CLOCKS_PER_SEC;
	total_time += cpu_time_used;

	getFormattedTime(cpu_time_used,(char *)formatted_string);
	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Reading matrix A,B and C from memory ended in %s",formatted_string);
	logger_info(logger_message);



	logger_info("Matrix values acquired correctly on HOST");

	if(cols <= 16 && rows <= 16){
		puts("Matrix A:");
		print_matrix(matrix_A,rows,cols);
	
		puts("Matrix B:");
		print_matrix(matrix_B,rows,cols);
	
		puts("Matrix C:");
		print_matrix(matrix_C,rows,cols);
	}

	
	computation_start = clock();
	compute_sequential_matrix_by_matrix_multiplication(matrix_A,matrix_B,matrix_C,rows,cols);
	computation_end = clock();
	cpu_time_used = ((double) (computation_end - computation_start)) / CLOCKS_PER_SEC;
	total_time += cpu_time_used;

	memset(formatted_string,0,64);
	getFormattedTime(cpu_time_used,(char *)formatted_string);

	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Sequential computation ended in %s",formatted_string);
	logger_info(logger_message);
	
	if(cols <= 16 && rows <= 16){
		puts("New matrix C:");
		print_matrix(matrix_C,rows,cols);
	}


	memset(logger_message,0,LOG_MESSAGE_SIZE);
	sprintf(logger_message,"Computation + I/O time: %f\n",total_time);
	logger_info(logger_message);
	puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");

	return 0;
}