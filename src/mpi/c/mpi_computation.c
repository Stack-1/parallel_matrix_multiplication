#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "utils/stackio.h"
#include "logger/logger.h"

#define BLOCK_ROWS 3
#define BLOCK_COLS 3
#define AUDIT if(my_rank == 0)


/*

/// @brief This function computes the matrix multiplication between A and B, and saves the result in the array C_temp.
float *multiplyMatrices(float *A, int A_rows, int A_cols, float *B, int B_rows, int B_cols, float *C_temp, int rank) {
    
    if (A_cols != B_rows) {
        fprintf(stderr, "The number of columns in A must be equal to the number of rows in B.\n");
        return NULL;
    }



    C_temp = (float *)calloc(A_rows * B_cols, sizeof(float));
    if (C_temp == NULL) {
        fprintf(stderr, "C_temp allocation failed");
        return NULL;
    }

    for (int i = 0; i < A_rows; ++i) {
        for (int j = 0; j < B_cols; ++j) {
            for (int k = 0; k < A_cols; ++k) {
                C_temp[i * B_cols + j] += A[i * A_cols + k] * B[k * B_cols + j];
            }
        }
    }

    return C_temp;
}


void *send_matrix_to_root(float *matrix, int num_rows, int num_cols, int rank, int proc_rows, int proc_cols, MPI_Comm comm) {
    
    //
    //   E.g.: 
    //        Process grid:
    //        0 1 2 3 4
    //        5 6 7 8 9
//
  //          the processes 0, 1, 2, 3, 4 will send their C_temp to process 0 (1/5 = 2/5 = 3/5 = 4/5 = 0)
    //        the processes 5, 6, 7, 8, 9 will send their C_temp to process 1 (6/5 = 7/5 = 8/5 = 9/5 = 1)
      //      the processes 0 and 5 are the receivers
    

    int row_index = rank / proc_cols;           // this is the index of the row of the process grid in which the current process is located
    int root_process = row_index * proc_cols;   // this is the index of the first process in the same row in which the current process is located

    // There will be a group for each row in the process grid (so, proc_rows groups).
    //  We're gonna create a new communicator for each group, having color equals to the row index.
    
    MPI_Comm group_comm;
    MPI_Comm_split(comm, row_index, rank, &group_comm);    // the processes on the same row will be assigned to the same group

    // this array will contain the element-wise sum of all the C_temp matrices belonging to the processes in the same group
    float *result = (float *)calloc(num_rows * num_cols, sizeof(float));  
    if (result == NULL) {
        fprintf(stderr, "result allocation failed");
        return NULL;
    }

    MPI_Reduce(matrix, result, num_rows*num_cols, MPI_FLOAT, MPI_SUM, 0, group_comm);

    MPI_Barrier(group_comm);    // sync barrier

    #ifdef DEBUG
    int group_rank;
    MPI_Comm_rank(group_comm, &group_rank);

    if (group_rank == 0)  {
        printf("\n[Process %d] Printing the resulting matrix C_temp (%dx%d)\n", rank, num_rows, num_cols);
        for (int i = 0; i < num_rows; i++) {
            puts("");
            for (int j = 0; j < num_cols; j++) {
                printf("%d ", (int)result[i*num_cols + j]);
            }
        }
        puts("");
        fflush(stdout);
    }
    #endif

    return result;
}


void write_result_to_file(char *filename, float *result, int result_len, int rank, int matrix_rows, int matrix_cols, int block_rows, int proc_rows, int proc_cols, MPI_Comm comm) {
    
    printf("[Process %d] Trying to write result (length %d) on the resulting matrix (%dx%d) file\n", rank, result_len, matrix_rows, matrix_cols);

    MPI_Datatype filetype;
    MPI_File file;

    MPI_Status status;

    int rank_norm = rank / proc_cols;
    printf("[Process %d] rank_norm: %d\n", rank, rank_norm);

    int gsizes[2], distribs[2], dargs[2], psizes[2];

    gsizes[0] = matrix_rows; // rows of the original matrix 
    gsizes[1] = matrix_cols; // columns of the original matrix

    // block cyclic distribution 
    distribs[0] = MPI_DISTRIBUTE_CYCLIC;
    distribs[1] = MPI_DISTRIBUTE_CYCLIC;

    dargs[0] = block_rows; // rows of the block 
    dargs[1] = matrix_cols; // columns of the block 

    psizes[0] = proc_rows; // no. of processes in vertical dimension of process grid 
    psizes[1] = 1; // no. of processes in horizontal dimension of process grid (1, because every process in the same row gets the same rows of the original matrix) 

    MPI_Type_create_darray(proc_rows, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    
    MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file.\n");
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(file, result, result_len, MPI_FLOAT, &status);

    MPI_File_close(&file);
}
*/



/**
 * @brief This function computes how many values should be assigned to every process:
 *        it determines the max number of blocks that can be assigned to a process, by computing how many blocks
 *        there are on the horizontal and vertical "axis" (i.e., rows and cols dimensions) and how these must be
 *        distributed in the process grid, according to a block cyclic distribution
 *        Note that the output of this function is an upper limit: most of the processes need less space
 * @param matrix_rows Number of rows of the matrix to be distributed
 * @param matrix_cols Number of columns of the matrix to be distributed
 * @param block_rows Number of rows of the block
 * @param block_cols Number of columns of the block
 * @param proc_rows Number of rows of the process grid
 * @param proc_cols Number of columns of the process grid
 * @return Number of floats the process should allocate to receive its portion of the matrix
*/
int compute_submatrix_dimension(
        size_t matrix_rows, 
        size_t matrix_cols, 
        int block_rows, 
        int block_cols, 
        int proc_rows, 
        int proc_cols, 
        int my_rank) {

    int my_rank_coord_x, my_rank_coord_y;

    my_rank_coord_x = my_rank / proc_cols;
    my_rank_coord_y = my_rank % proc_cols; 

    //printf("[Process %d] my coordinates are (%d, %d)\n", my_rank, my_rank_coord_x, my_rank_coord_y);
    

    /* COLS */
    int tot_num_blocks = ceil((float)matrix_cols / (float)block_cols); // 12 / 3 = 4

    int tot_blocks_assigned = tot_num_blocks / proc_cols; // 4 / 2 = 2

    if (tot_num_blocks % proc_cols != 0 && my_rank_coord_y < tot_num_blocks % proc_cols) {
            tot_blocks_assigned++;
    } 

    /* assign last block portion */

    int tot_nums_assigned_x = tot_blocks_assigned * block_cols; // 2 * 3 = 6

    int last_process_x = tot_num_blocks % proc_cols == 0 ? proc_cols -1 : (tot_num_blocks % proc_cols) -1;
    if (my_rank_coord_y == last_process_x && matrix_cols % block_cols != 0) {
        tot_nums_assigned_x -= (block_cols - (matrix_cols % block_cols));
    }

    /* ROWS */ 
    // my_rank_coord_x = 1 
    tot_num_blocks = ceil((float)matrix_rows / (float)block_rows);    // 25 / 3 = 9

    tot_blocks_assigned = tot_num_blocks / proc_rows;   // 9 / 3 = 3

    if (tot_num_blocks % proc_rows != 0 && my_rank_coord_x < tot_num_blocks % proc_rows) {
            tot_blocks_assigned++;
    } 

    /* assign last block portion */

    int tot_nums_assigned_y = tot_blocks_assigned * block_rows; // 3 * 3 = 9

    int last_process_y = ((tot_num_blocks % proc_rows) == 0) ? (proc_rows -1) : ((tot_num_blocks % proc_rows) -1);  // 2
    if (my_rank_coord_x == last_process_y && matrix_rows % block_rows != 0) {
        tot_nums_assigned_y -= (block_rows - (matrix_rows % block_rows));
    }

    printf("[Process %d] Assigned %d rows\n", my_rank, tot_nums_assigned_y);
    printf("[Process %d] Assigned %d cols\n", my_rank, tot_nums_assigned_x);
    printf("[Process %d] Assigned %d elements\n", my_rank, tot_nums_assigned_y * tot_nums_assigned_x);
    puts("");

    return tot_nums_assigned_x * tot_nums_assigned_y;

}

/**
 * @brief This function divides the matrices A, B and C between n_proc processes, according to a block-cyclic distribution.
 * @param matrix_file_name Name of the file containing A matrix in binary format
 * @param number_of_processes number of processes to divide the matrices between
 * @param matrix_rows matrix dimension on y axis
 * @param matrix_cols matrix dimension on x axis
 * @param my_rank ID of the current process (PID)
 * @param comm MPI private comunicator used in thi kernel of computation
*/
void generate_block_cyclic_distribution_matrix_A(
        char *matrix_file_name, 
        size_t matrix_rows, 
        size_t matrix_cols, 
        int block_rows, 
        int block_cols, 
        int process_grid_y, 
        int process_grid_x, 
        int my_rank,
        MPI_Comm comm) {

    int matrix_size;
    float *matrix;
    char logger_message[LOG_MESSAGE_SIZE];
    
    MPI_Datatype filetype;
    MPI_File file;

    MPI_Status status;

    /**
     * Reading from file.
     */
    int gsizes[2], distribs[2], dargs[2], psizes[2];

    gsizes[0] = matrix_rows; /* x dimension of the original matrix */
    gsizes[1] = matrix_cols; /* y dimension of the original matrix */

    /* block cyclic distribution */
    distribs[0] = MPI_DISTRIBUTE_CYCLIC;
    distribs[1] = MPI_DISTRIBUTE_CYCLIC;

    dargs[0] = block_rows; /* x dimension of the block */
    dargs[1] = block_cols; /* y dimension of the block */

    psizes[0] = process_grid_y; /* no. of processes in vertical dimension of process grid */
    psizes[1] = process_grid_x; /* no. of processes in horizontal dimension of process grid */

    MPI_Type_create_darray(process_grid_y * process_grid_x, my_rank, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
	    memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Error opening file %s\n",matrix_file_name);
	    logger_info(logger_message);
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);

    
    matrix_size = compute_submatrix_dimension(matrix_rows, matrix_cols, block_rows, block_cols, process_grid_y, process_grid_x, my_rank);

    matrix = (float*) malloc(matrix_size * sizeof(float));
    MPI_File_read_all(file, matrix, matrix_size,
            MPI_FLOAT, &status);

    MPI_File_close(&file);

    printf("[Process %d] received : ", my_rank);
    for (int i = 0; i < matrix_size; i++) {
        printf("%f ", matrix[i]);
    }
    puts("\n");
    fflush(stdout);

}

/**
 * @brief This function has the responsability to devide the matrix B in rows and give each group of rows
 * to the correct process in order to succesively compute the molutiplication between A block and B rows
 * @param matrix_file_name Name of the file containing A matrix in binary format
 * @param matrix_rows Number of rows of the matrix to be distributed
 * @param matrix_cols Number of columns of the matrix to be distributed
 * @param block_rows Number of rows of the block
 * @param proc_rows Number of rows of the process grid
 * @param my_rank ID of the current process (PID)
 * @param comm MPI private comunicator used in thi kernel of computation
*/
float *generate_rows_cyclyc_distribution(
        char *matrix_file_name, 
        size_t matrix_rows, 
        size_t matrix_cols, 
        int block_rows, 
        int proc_rows, 
        int my_rank, 
        MPI_Comm comm) {

    int matrix_size;
    int num_rows = (int)matrix_rows / (block_rows*proc_rows);
    float *matrix;
    
    MPI_Datatype filetype;
    MPI_File file;

    MPI_Status status;

    int rank_norm = my_rank % proc_rows;

    int gsizes[2], distribs[2], dargs[2], psizes[2];

    gsizes[0] = matrix_rows; // rows of the original matrix 
    gsizes[1] = matrix_cols; // columns of the original matrix 

    // block cyclic distribution 
    distribs[0] = MPI_DISTRIBUTE_CYCLIC;
    distribs[1] = MPI_DISTRIBUTE_CYCLIC;

    dargs[0] = block_rows; // rows of the block 
    dargs[1] = matrix_cols; // columns of the block 

    psizes[0] = proc_rows; // no. of processes in vertical dimension of process grid 
    psizes[1] = 1; // no. of processes in horizontal dimension of process grid (1, because every process in the same row gets the same rows of the original matrix) 

    MPI_Type_create_darray(proc_rows, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file.\n");
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    matrix_size = num_rows*matrix_cols;

    matrix = (float*) malloc(sizeof(float));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed");
        return NULL;
    }
    MPI_File_read_all(file, matrix, matrix_size,
            MPI_FLOAT, &status);

    MPI_File_close(&file);

    printf("\n======\ni = %d\nBlock received : ", my_rank);
    for (int i = 0; i < matrix_size; i++) {
        printf("%d ", (int)matrix[i]);
    }
    puts("\n");
    fflush(stdout);

}



/**
 * @brief Startup of the parallel computation using MPI
 * @param matrix_A_file_name Name of the file in which is stored the matrix A
 * @param matrix_B_file_name Name of the file in which is stored the matrix B
 * @param matrix_C_file_name Name of the file in which is stored the matrix C
 * @param rows Number of rows of the matrix
 * @param cols Number of cols of the matrix
 * @param number_of_processes Number of processes launched by mpirun, know at runtime
 * @param my_rank ID of the current process (PID)
 * @param comm MPI private comunicator used in thi kernel of computation
*/
void compute_process_parallel_matrix_by_matrix_multiplication(
        char *matrix_A_file_name,
        char *matrix_B_file_name,
        char *matrix_C_file_name,
        size_t rows,
        size_t cols,
        int number_of_processes,
        int my_rank,
        MPI_Comm comm){

    int process_grid_rows, process_grid_cols;
    //float *matrix_data_process_subarray = (float *)malloc(sizeof());

	process_grid_cols = ceil(sqrt(number_of_processes)); //Prefer getting more process on column to scale better
	process_grid_rows = (int)(sqrt(number_of_processes));

    printf("%d\n",compute_submatrix_dimension(rows,cols,BLOCK_ROWS,BLOCK_COLS,process_grid_rows,process_grid_cols,my_rank));



    generate_block_cyclic_distribution_matrix_A(matrix_A_file_name, rows, cols, BLOCK_ROWS, BLOCK_COLS,process_grid_rows,process_grid_cols,my_rank,comm);
    
    //generate_rows_cyclyc_distribution(matrix_B_file_name, rows, cols, BLOCK_ROWS, process_grid_rows, my_rank, comm);    
    

}



int main(int argc, char *argv[]){
	size_t rows = 0;
    size_t cols = 0;
    clock_t read_start, read_end, computation_start, computation_end;
	double cpu_time_used;
	double total_time = 0.0;
	char logger_message[LOG_MESSAGE_SIZE];
	char formatted_string[64];
    char matrix_A_file_name[64];
    char matrix_B_file_name[64];
    char matrix_C_file_name[64];
    int number_of_processes;
    int my_rank;
    MPI_Comm comm;


    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm); //Get private comunicator

    MPI_Comm_size(comm, &number_of_processes); //Initialize number of process
    MPI_Comm_rank(comm, &my_rank); //Initialize rank (process_id - PID)


	if(argc != 3){
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols> ");
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

    memset(matrix_A_file_name,0,64);
    memset(matrix_B_file_name,0,64);
    memset(matrix_C_file_name,0,64);
    sprintf(matrix_A_file_name,"src/data/matrix%lux%lu/matrix_A_%lux%lu.bin", rows, cols, rows, cols);
    sprintf(matrix_B_file_name,"src/data/matrix%lux%lu/matrix_B_%lux%lu.bin",rows, cols, rows, cols);
    sprintf(matrix_C_file_name,"src/data/matrix%lux%lu/matrix_C_%lux%lu.bin",rows, cols, rows, cols);

    AUDIT{
        puts("");
	    puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");
    	memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Starting process parallel computation (MPI) for matrix %ld x %ld",rows,cols);
	    logger_info(logger_message);
    }


    // Starting actual parallel matrix multiplication
	computation_start = clock();
	compute_process_parallel_matrix_by_matrix_multiplication(
        matrix_A_file_name,
        matrix_B_file_name,
        matrix_C_file_name,
        rows,
        cols,
        number_of_processes,
        my_rank,
        comm
        );

	computation_end = clock();
	cpu_time_used = ((double) (computation_end - computation_start)) / CLOCKS_PER_SEC;
	total_time += cpu_time_used;

	memset(formatted_string,0,64);
	getFormattedTime(cpu_time_used,(char *)formatted_string);


    AUDIT{
        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Process parallel (MPI) computation ended in %s",formatted_string);
	    logger_info(logger_message);

        memset(logger_message,0,LOG_MESSAGE_SIZE);
        sprintf(logger_message,"Computation + I/O time: %f\n",total_time);
        logger_info(logger_message);
        puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
        puts("");
    }

    MPI_Finalize();

	return 0;
}