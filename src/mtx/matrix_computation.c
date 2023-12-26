#include <aio.h>
#include "rndgen.h"

void compute_sequential_matrix_by_matrix_multiplication(float **matrix_A,float **matrix_B,float **matrix_C,size_t rows,size_t cols){
	for (size_t i=0;i<rows;i++) 
        for(size_t j=0;j<cols;j++) {
            for(size_t k=0;k<cols;k++){    
                matrix_C[i][j]=matrix_C[i][j]+matrix_A[i][k]*matrix_B[k][j];
			}
    }
}