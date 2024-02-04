#include "mpi_computation.h"


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


void reduce_partial_results(float *matrix, float **return_matrix, int num_rows, int num_cols, int my_rank, int proc_rows, int proc_cols, MPI_Comm comm) {
    
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


    int group_rank;
    MPI_Comm_rank(group_comm, &group_rank);
#ifdef DEBUG
    if (group_rank == 0)  {
        printf("\n[Process %d] Printing the resulting matrix C_temp (%dx%d)\n", my_rank, num_rows, num_cols);
        for (int i = 0; i < num_rows; i++) {
            puts("");
            for (int j = 0; j < num_cols; j++) {
                printf("%f ", (*return_matrix)[i*num_cols + j]);
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

    my_rank_coord_x = my_rank / proc_cols; // 0
    my_rank_coord_y = my_rank % proc_cols; // 0

#ifdef DEBUG
    AUDIT
        printf("[Process %d] my coordinates are (%d, %d)\n", my_rank, my_rank_coord_x, my_rank_coord_y);
#endif

    /* COLS */
    int tot_num_blocks = ceil( (float)matrix_cols / (float)block_cols); 

    int tot_blocks_assigned = tot_num_blocks / proc_cols; 

    if (tot_num_blocks % proc_cols != 0 && my_rank_coord_y < tot_num_blocks % proc_cols) {
            tot_blocks_assigned++; 
    } 

    /* assign last block portion */

    int tot_nums_assigned_x = tot_blocks_assigned * block_cols; 

    int last_process_x = tot_num_blocks % proc_cols == 0 ? proc_cols -1 : (tot_num_blocks % proc_cols) -1; // 1
    if (my_rank_coord_y == last_process_x && matrix_cols % block_cols != 0) {
        tot_nums_assigned_x -= (block_cols - (matrix_cols % block_cols));
    }

    /* ROWS */ 
    tot_num_blocks = ceil((float)matrix_rows / (float)block_rows);    

    tot_blocks_assigned = tot_num_blocks / proc_rows;   
    if (tot_num_blocks % proc_rows != 0 && my_rank_coord_x < tot_num_blocks % proc_rows) {
            tot_blocks_assigned++;
    } 

    /* assign last block portion */

    int tot_nums_assigned_y = tot_blocks_assigned * block_rows; // 42 * 3 = 

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
    int dargs[2] = {block_rows, N_subarray}; // I want to have all the cols of the matrix, having only the correct rows 
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