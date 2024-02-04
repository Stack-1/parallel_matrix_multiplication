#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "../../utils/logger/logger.h"
#include "../../utils/string_formatter/formatter.h"
#include "../../utils/file_io/stackio.h"


#define BLOCK_ROWS 3
#define BLOCK_COLS 3
#define AUDIT if(my_rank == 0)

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
        __typeof__ (b) _b = (b); \
        _a > _b ? _a : _b; })

typedef struct computed_dimensions{
    int rows;
    int cols;
}computed_dimensions;

#ifndef FORMATTED_STRING_SIZE
#define FORMATTED_STRING_SIZE 64
#endif


void write_result_to_file(char *filename, float *result, int result_len, int rank, int matrix_rows, int matrix_cols, int block_rows, int proc_rows, int proc_cols, MPI_Comm comm) {
    
    MPI_Datatype filetype;
    MPI_File file;

    MPI_Status status;

    int rank_norm = rank / proc_cols;
    int gsizes[2], distribs[2], dargs[2], psizes[2];

    gsizes[0] = matrix_rows; /* rows of the original matrix */
    gsizes[1] = matrix_cols; /* columns of the original matrix */

    /* block cyclic distribution */
    distribs[0] = MPI_DISTRIBUTE_CYCLIC;
    distribs[1] = MPI_DISTRIBUTE_CYCLIC;

    dargs[0] = block_rows; /* rows of the block */
    dargs[1] = matrix_cols; /* columns of the block */

    psizes[0] = proc_rows; /* no. of processes in vertical dimension of process grid */
    psizes[1] = 1; /* no. of processes in horizontal dimension of process grid (1, because every process in the same row gets the same rows of the original matrix) */

    MPI_Type_create_darray(proc_rows, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    
    MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file %s to write.\n",filename);
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(file, result, result_len, MPI_FLOAT, &status);

    MPI_File_close(&file);
}


void send_matrix_to_root(float *matrix, float **return_matrix, int num_rows, int num_cols, int my_rank, int proc_rows, int proc_cols, MPI_Comm comm) {
    
    //
    //   E.g.: 
    //        Process grid:
    //        0 1 2 3 4
    //        5 6 7 8 9
//
  //          the processes 0, 1, 2, 3, 4 will send their C_temp to process 0 (1/5 = 2/5 = 3/5 = 4/5 = 0)
    //        the processes 5, 6, 7, 8, 9 will send their C_temp to process 1 (6/5 = 7/5 = 8/5 = 9/5 = 1)
      //      the processes 0 and 5 are the receivers
    

    int row_index = my_rank / proc_cols;           // this is the index of the row of the process grid in which the current process is located

    // There will be a group for each row in the process grid (so, proc_rows groups).
    //  We're gonna create a new communicator for each group, having color equals to the row index.
    
    MPI_Comm group_comm;
    MPI_Comm_split(comm, row_index, my_rank, &group_comm);    // the processes on the same row will be assigned to the same group



    MPI_Reduce(matrix, *return_matrix, num_rows*num_cols, MPI_FLOAT, MPI_SUM, 0, group_comm);

    MPI_Barrier(group_comm);    // sync barrier

    #ifdef DEBUG
    int group_rank;
    MPI_Comm_rank(group_comm, &group_rank);

    if (group_rank == 0)  {
        printf("\n[Process %d] Printing the resulting matrix C_temp (%dx%d)\n", my_rank, num_rows, num_cols);
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

}




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
computed_dimensions *compute_submatrix_dimension(
        int matrix_rows, 
        int matrix_cols, 
        int block_rows, 
        int block_cols, 
        int proc_rows, 
        int proc_cols, 
        int my_rank) {

    int my_rank_coord_x, my_rank_coord_y;
    computed_dimensions *dimensions = (computed_dimensions *)malloc(sizeof(computed_dimensions));

    my_rank_coord_x = my_rank / proc_cols;
    my_rank_coord_y = my_rank % proc_cols; 

#ifdef DEBUG
    AUDIT
        printf("[Process %d] my coordinates are (%d, %d)\n", my_rank, my_rank_coord_x, my_rank_coord_y);
#endif

    /* COLS */
    int tot_num_blocks = ceil( (float)matrix_cols / (float)block_cols); // 12 / 3 = 4

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
#ifdef DEBUG
    printf("[Process %d] Assigned %d rows\n", my_rank, tot_nums_assigned_y);
    printf("[Process %d] Assigned %d cols\n", my_rank, tot_nums_assigned_x);
    printf("[Process %d] Assigned %d elements\n", my_rank, tot_nums_assigned_y * tot_nums_assigned_x);
    puts("");
#endif
    dimensions->rows = tot_nums_assigned_y;
    dimensions->cols = tot_nums_assigned_x;
    return dimensions;

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
        float **matrix_A,
        int N, 
        int K, 
        int N_subarray,
        int K_subarray,
        int block_rows, 
        int block_cols, 
        int process_grid_y, 
        int process_grid_x, 
        int my_rank,
        MPI_Comm comm) {

    char logger_message[LOG_MESSAGE_SIZE];
    
    MPI_Datatype filetype;
    MPI_File file;
    MPI_Status status;

    /**
     * Reading from file.
     */
    int dims[2] = {N,K};
    int distribs[2] =  {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows,block_cols};  // Arguments to pass to the MPI_Datatype constructor
    int proc_dims[2] = {process_grid_y,process_grid_x};
    int num_procs = process_grid_x * process_grid_y;

    MPI_Type_create_darray(num_procs, my_rank, 2, dims, distribs, dargs, proc_dims, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
	    memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Error opening file %s\n",matrix_file_name);
	    logger_info(logger_message);
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 2*sizeof(int), MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_read_all(file, (*matrix_A), N_subarray * K_subarray, MPI_FLOAT, &status);

    MPI_File_close(&file);

    /*printf("[Process %d] received : ", my_rank);
    for (int i = 0; i < matrix_size; i++) {
        printf("%f ", matrix[i]);
    }
    puts("\n");
    fflush(stdout);*/
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
void generate_rows_distribution(
        char *matrix_file_name, 
        float **matrix,
        int K, 
        int M, 
        int K_subarray,
        int M_subarray,
        int block_rows, 
        int proc_cols, 
        int my_rank, 
        MPI_Comm comm) {

    MPI_Datatype filetype;
    MPI_File file;
    MPI_Status status;

    int rank_norm = my_rank % proc_cols;
    int gsizes[2] = {K, M};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows, M_subarray}; // I want to have all the cols of the matrix, having only the correct rows 
    int psizes[2] = {proc_cols,1}; // We want to divide the matrix only in one dimension, by rows
    

    MPI_Type_create_darray(proc_cols, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file.\n");
        MPI_Abort(comm, 1);
    }


    MPI_File_set_view(file, 2*sizeof(int), MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_read_all(file, (*matrix), K_subarray*M_subarray, MPI_FLOAT, &status);


    MPI_File_close(&file);

    /*printf("\n======\ni = %d\nBlock received : ", my_rank);
    for (int i = 0; i < matrix_size; i++) {
        printf("%f ", matrix[i]);
    }
    puts("\n");
    fflush(stdout);*/
}


void generate_rows_distribution_C(
        char *matrix_file_name, 
        float **matrix,
        int N, 
        int M, 
        int N_subarray,
        int M_subarray,
        int block_rows, 
        int proc_rows,
        int proc_cols, 
        int my_rank, 
        MPI_Comm comm) {

    MPI_Datatype filetype;
    MPI_File file;
    MPI_Status status;

    int rank_norm = my_rank / proc_cols;
    int gsizes[2] = {N, M};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows, M_subarray}; // I want to have all the cols of the matrix, having only the correct rows 
    int psizes[2] = {proc_rows,1}; // We want to divide the matrix only in one dimension, by rows

    MPI_Type_create_darray(proc_rows, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file %s.\n",matrix_file_name);
        MPI_Abort(comm, 1);
    }


    MPI_File_set_view(file, 2*sizeof(int), MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_read_all(file, (*matrix), N_subarray*M_subarray, MPI_FLOAT, &status);


    MPI_File_close(&file);

    /*printf("\n======\ni = %d\nBlock received : ", my_rank);
    for (int i = 0; i < matrix_size; i++) {
        printf("%f ", matrix[i]);
    }
    puts("\n");
    fflush(stdout);*/
}


/** @brief This function computes the matrix multiplication between A and B, and saves the result in the array C_temp.
 * 
*/
void multiplyMatrices(float *matrix_A, float *matrix_B, float **matrix_C, int N, int K, int M, int my_rank) {
#ifdef DEBUG
    printf("[Process %d] Computation of matrix %dx%d\n", my_rank, N, M);
#endif
	for (int i=0; i<N; ++i){ 
        for(int k=0; k<K; ++k) {
            for(int j=0; j<M; ++j){    
                (*matrix_C)[i * M + j] += matrix_A[i * K + k] * matrix_B[k * M + j];
			}
    	}
	}
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
        int N,
        int K,
        int M,
        int process_grid_rows,
        int process_grid_cols,
        int my_rank,
        MPI_Comm comm){

    computed_dimensions *matrix_sizes;

    float *matrix_A_subarray;
    float *matrix_B_subarray;
    float *matrix_C_subarray;
    float *matrix_C_subarray_rows;
    float *matrix_C_subarray_from_file;

    int N_subarray;
    int K_subarray;
    int M_subarray;

    AUDIT{ 
        printf("[INFO] Process grid is: %d %d\n",process_grid_rows,process_grid_cols);
        puts("----------------------------------------------------------------");
        for(int i = 0;i<process_grid_rows;++i){
            printf("|");
            for(int j = 0;j<process_grid_cols;++j){
                printf("%d ",i*process_grid_cols+j);
            }
            puts("|");
        }
        puts("----------------------------------------------------------------");
        fflush(stdout);
    }

    matrix_sizes = compute_submatrix_dimension(N, K, BLOCK_ROWS, BLOCK_COLS, process_grid_rows, process_grid_cols, my_rank);

    N_subarray = matrix_sizes->rows;
    K_subarray = matrix_sizes->cols;
    M_subarray = M;

    
    matrix_A_subarray = (float *)malloc(N_subarray * K_subarray * sizeof(float)); 
    matrix_B_subarray = (float *)malloc(K_subarray * M_subarray * sizeof(float));
    matrix_C_subarray = (float *)malloc(N_subarray * M_subarray * sizeof(float));

    matrix_C_subarray_rows = (float *)malloc(N_subarray * M_subarray * sizeof(float));

    memset(matrix_C_subarray,0,N_subarray * M_subarray * sizeof(float));

    generate_block_cyclic_distribution_matrix_A(matrix_A_file_name,&matrix_A_subarray ,N, K, N_subarray, K_subarray,BLOCK_ROWS, BLOCK_COLS,process_grid_rows,process_grid_cols,my_rank,comm);

    generate_rows_distribution(matrix_B_file_name,&matrix_B_subarray, K, M, K_subarray , M_subarray, BLOCK_ROWS, process_grid_cols, my_rank,comm);


    multiplyMatrices(matrix_A_subarray,matrix_B_subarray,&matrix_C_subarray, N_subarray, K_subarray ,M_subarray,my_rank);

    send_matrix_to_root(matrix_C_subarray, &matrix_C_subarray_rows,N_subarray, M_subarray, my_rank, process_grid_rows, process_grid_cols, comm);

    MPI_Barrier(comm);

    matrix_C_subarray_from_file = (float *)malloc(N_subarray * M_subarray * sizeof(float));
            if(my_rank == 2){
            printf("%d %d %d %d %d %d %d %d\n",N,M,N_subarray,M_subarray,BLOCK_ROWS,process_grid_rows,process_grid_cols,my_rank/process_grid_cols);
                        fflush(stdout);

            }
    
    generate_rows_distribution_C(matrix_C_file_name,&matrix_C_subarray_from_file, N, M, N_subarray , M_subarray, BLOCK_ROWS, process_grid_rows, process_grid_cols,my_rank,comm);

    MPI_Barrier(comm);


    if (my_rank % process_grid_cols == 0) {
        for (int i = 0; i < N_subarray * M_subarray; ++i) {
            matrix_C_subarray_rows[i] += matrix_C_subarray_from_file[i];
        }
    }

    MPI_Barrier(comm);

    MPI_Comm first_element_comm;
    MPI_Comm_split(comm, my_rank % process_grid_cols, my_rank, &first_element_comm);

    int first_element_rank;
    MPI_Comm_rank(first_element_comm, &first_element_rank);

     /* Write result on file */
    if (my_rank % process_grid_cols == 0) {
        char file_name[128];
        if(N==K && K==M){
            sprintf(file_name, "../../data/square/matrix%dx%d/matrix_C_mpi_%dx%d.bin",N,M,N,M);
        }else{
            sprintf(file_name, "../../data/rectangular/matrix%dx%d/matrix_C_mpi_%dx%d.bin",N,M,N,M);  
        }
        
        write_result_to_file(file_name, matrix_C_subarray_rows, N_subarray * M_subarray, my_rank, N, M, BLOCK_ROWS, process_grid_rows, process_grid_cols, first_element_comm);
    }

}



int main(int argc, char *argv[]){
	int N = 0;
    int K = 0;
    int M = 0;

	char logger_message[LOG_MESSAGE_SIZE];
	char formatted_string[FORMATTED_STRING_SIZE];

    char matrix_A_file_name[FORMATTED_STRING_SIZE];
    char matrix_B_file_name[FORMATTED_STRING_SIZE];
    char matrix_C_file_name[FORMATTED_STRING_SIZE];

    int number_of_processes;
    int process_grid_rows,process_grid_cols;

    double total_time = 0.0f;
    double start, end;

    bool is_square = false;
        
    int my_rank;
    MPI_Comm comm;


    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm); //Get private communicator

    MPI_Comm_size(comm, &number_of_processes); //Initialize number of process
    MPI_Comm_rank(comm, &my_rank); //Initialize rank (process_id - PID)


	if(argc != 6){
        logger_error("Program should be called like ./<elf file name> <# rows> <# cols> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%d", &N) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%d", &K) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[3], "%d", &M) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[4], "%d", &process_grid_rows) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[5], "%d", &process_grid_cols) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M> <# process grid rows> <# process grid cols>");
        exit(EXIT_FAILURE);
    }

    // As first thing I want to check if the process grid is coherent with the number of processes
    if(process_grid_rows*process_grid_cols != number_of_processes){
        AUDIT
            logger_error("Invalid process grid! The program will not start.");
        exit(EXIT_FAILURE);
    }
    
    // Check if matrix is a square matrix
    is_square = (N==K && K==M) ? true : false;

    memset(matrix_A_file_name,0,FORMATTED_STRING_SIZE);
    memset(matrix_B_file_name,0,FORMATTED_STRING_SIZE);
    memset(matrix_C_file_name,0,FORMATTED_STRING_SIZE);
    if(is_square){
        sprintf(matrix_A_file_name,"../../data/square/matrix%dx%d/matrix_A_%dx%d.bin", N, K, N, K);
        sprintf(matrix_B_file_name,"../../data/square/matrix%dx%d/matrix_B_%dx%d.bin",K, M, K, M);
        sprintf(matrix_C_file_name,"../../data/square/matrix%dx%d/matrix_C_%dx%d.bin",N, M, N, M);
    }else{
        sprintf(matrix_A_file_name,"../../data/rectangular/matrix%dx%d/matrix_A_%dx%d.bin", N, K, N, K);
        sprintf(matrix_B_file_name,"../../data/rectangular/matrix%dx%d/matrix_B_%dx%d.bin",K, M, K, M);
        sprintf(matrix_C_file_name,"../../data/rectangular/matrix%dx%d/matrix_C_%dx%d.bin",N, M, N, M); 
    }
    AUDIT{
        puts("");
	    puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");
    	memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Starting process parallel computation (MPI) for matrix %d x %d",N,M);
	    logger_info(logger_message);
        fflush(stdout);
    }


    // Starting actual parallel matrix multiplication
    start = MPI_Wtime();
	compute_process_parallel_matrix_by_matrix_multiplication(
        matrix_A_file_name,
        matrix_B_file_name,
        matrix_C_file_name,
        N,
        K,
        M,
        process_grid_rows,
        process_grid_cols,
        my_rank,
        comm
        );

	end = MPI_Wtime();
	total_time = (double) (end - start);;

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(total_time,(char *)formatted_string);


    AUDIT{
        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Process parallel (MPI) computation ended in %s",formatted_string);
	    logger_info(logger_message);

        memset(logger_message,0,LOG_MESSAGE_SIZE);
        sprintf(logger_message,"Computation + I/O time: %f\n",total_time);
        logger_info(logger_message);
        puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
        puts("");
        fflush(stdout);
    }


    AUDIT{
        float *matrix_c_sequential;
        float *mpi;
        float diff = 0.0f;
        float max_diff = 0.0f;
        float max_rel_diff = 0.0f;


        matrix_c_sequential = (float *)malloc(sizeof(float)*N*M);
        mpi = (float *)malloc(sizeof(float)*N*M);


        read_matrix_from_file(matrix_c_sequential,N,M,(char *)"sequential_C",is_square);
        read_matrix_from_file_mpi(mpi,N,M,(char *)"C_mpi",is_square);
        
#ifdef DEBUG
        fflush(stdout);
        print_matrix(matrix_c_sequential,N,M);
        printf("\n\n\n\n");
        fflush(stdout);
        print_matrix(mpi,N,M);
        fflush(stdout);
#endif

        for(int i=0; i<N; ++i){
            for(int j=0; j<M; ++j){
                diff += abs(mpi[i*M+j] - matrix_c_sequential[i*M+j]);

                max_diff = max( abs(mpi[i*M+j] - matrix_c_sequential[i*M+j]), max_diff);
                max_rel_diff = max(  
                    abs(mpi[i*M+j] - matrix_c_sequential[i*M+j]) /  max(mpi[i*M+j], matrix_c_sequential[i*M+j]), 
                    max_rel_diff
                    );
            }
        }
        fflush(stdout);

        puts("");
        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Diff is: %f\n",diff);
	    logger_info(logger_message);
        
        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Max diff is: %f\n",max_diff);
	    logger_info(logger_message);

        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Max rel diff is: %f\n",max_rel_diff);
	    logger_info(logger_message);


    }


    MPI_Finalize();

	return 0;
}