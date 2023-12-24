#include <stdio.h>
#include <malloc.h>
#include "mtx/rndgen.h"

#define ROWS 16
#define COLS 16

int main(int argc, char *argv[]){
	float *matrix_A[ROWS];
	float *matrix_B[ROWS];
	float *matrix_C[ROWS];
	// TODO: See if it is more efficent to keep matrix as 1D dimencional array

	puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");

	for(int i = 0;i<COLS;i++){
		matrix_A[i] = (float *)malloc(COLS*sizeof(float));
		matrix_B[i] = (float *)malloc(COLS*sizeof(float));
		matrix_C[i] = (float *)malloc(COLS*sizeof(float));
	}

	puts("[INFO] Matrix memory correctly allocated on HOST");

	generate_real_value_matrix(matrix_A,ROWS,COLS);
	generate_real_value_matrix(matrix_B,ROWS,COLS);
	generate_zero_value_matrix(matrix_C,ROWS,COLS);

	puts("[INFO] Matrix values generation correctly finished on HOST");

	puts("Matrix A:");
	print_matrix(matrix_A,ROWS,COLS);
	
	puts("Matrix B:");
	print_matrix(matrix_B,ROWS,COLS);
	
	puts("Matrix C:");
	print_matrix(matrix_C,ROWS,COLS);


	puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");

	return 0;
}
