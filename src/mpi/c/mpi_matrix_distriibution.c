#include "mpi_computation.h"

/**
 * @brief Function that has the responsability to write the result of the MPI computation on file.
 * To be sure that each process writes the right part of the file, we define a new datatype, specifying
 * MPI_DISTRIBUTE_CYCLIC to be coherent with the block cyclic distribution scheme used so far.
 * @param matrix_file_name Name of the file in which we want to write the final computatation of C
 * @param matrix The sub-matrix of C containing the values we want to save
 * @param sub_matrix_size The size of the sub-matrix of C we want to save on file
 * @param matrix_rows matrix dimension on y axis
 * @param matrix_cols matrix dimension on x axis
 * @param block_rows Dimension of block on y axis
 * @param process_grid_rows Dimension of the process grid on the y axis
 * @param process_grid_cols Dimension of the process grid on the x axis
 * @param my_rank ID of the current process (PID)
 * @param comm MPI private comunicator used in thi kernel of computation
*/
void write_result_to_file(
        char *matrix_file_name, 
        float *matrix, 
        int sub_matrix_size,  
        int matrix_rows, 
        int matrix_cols, 
        int block_rows, 
        int process_grid_rows, 
        int process_grid_cols, 
        int my_rank,
        MPI_Comm comm) {
    
    MPI_Datatype filetype;
    MPI_File file;

    MPI_Status status;

    int raws_index = my_rank / process_grid_cols;

    int gsizes[2] = {matrix_rows, matrix_cols};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows, matrix_cols}; // I want to have all the cols of the matrix, having only the correct rows 
    int psizes[2] = {process_grid_rows,1}; // We want to divide the matrix only in one dimension, by rows
    

    MPI_Type_create_darray(process_grid_rows, raws_index, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    
    MPI_File_open(comm, matrix_file_name, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file %s to write.\n",matrix_file_name);
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(file, matrix, sub_matrix_size, MPI_FLOAT, &status);

    MPI_File_close(&file);
}

/**
 * @brief Function that has the responsability of putting toghether all the partial computation obtained 
 * by each process on the same row of the process grid.
 * @param matrix The sub-matrix containing the partial result of the computation
 * @param return_matrix The address of the memory in which the reduction will be saved
 * @param rows Number of rows of the sub-matrix
 * @param cols Number of columns of the sub-matrix
 * @param process_grid_cols Number of columnas of the process grid
 * @param my_rank PID od the process running
 * @param comm Communicator in use
*/
void reduce_partial_results(
        float *matrix, 
        float **return_matrix, 
        int rows, 
        int cols, 
        int process_grid_cols,
        int my_rank, 
        MPI_Comm comm) {
    

    int row_index = my_rank / process_grid_cols;        
    int my_group_rank;

    MPI_Comm group_comm;
    MPI_Comm_split(comm, row_index, my_rank, &group_comm); // Create communicator for the group

    MPI_Reduce(matrix, *return_matrix, rows*cols, MPI_FLOAT, MPI_SUM, 0, group_comm); // Actual reduction

    MPI_Barrier(group_comm);   

    MPI_Comm_rank(group_comm, &my_group_rank);
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
 * @param process_gird_rows Number of rows of the process grid
 * @param process_gird_cols Number of columns of the process grid
 * @return Number of floats the process should allocate to receive its portion of the matrix
*/
computed_dimensions *compute_submatrix_dimension(
        int matrix_rows, 
        int matrix_cols, 
        int block_rows, 
        int block_cols, 
        int process_grid_rows, 
        int process_grid_cols, 
        int my_rank) {

    int my_rank_coord_cols, my_rank_coord_rows;
    computed_dimensions *dimensions = (computed_dimensions *)malloc(sizeof(computed_dimensions));

    my_rank_coord_cols = my_rank / process_grid_cols; 
    my_rank_coord_rows = my_rank % process_grid_cols; 

#ifdef DEBUG
    AUDIT
        printf("[Process %d] my coordinates are (%d, %d)\n", my_rank, my_rank_coord_cols, my_rank_coord_rows);
#endif

    /* Compute for columns */
    int blocks_on_dimension = ceil( (float)matrix_cols / (float)block_cols); 

    int block_per_process = blocks_on_dimension / process_grid_cols; 

    if (blocks_on_dimension % process_grid_cols != 0 && my_rank_coord_rows < blocks_on_dimension % process_grid_cols) {
            block_per_process++; 
    } 

    /* Assign last block */
    int block_per_process_on_cols = block_per_process * block_cols; 

    int last_process_x = blocks_on_dimension % process_grid_cols == 0 ? process_grid_cols -1 : (blocks_on_dimension % process_grid_cols) -1; // 1
    if (my_rank_coord_rows == last_process_x && matrix_cols % block_cols != 0) {
        block_per_process_on_cols -= (block_cols - (matrix_cols % block_cols));
    }

    /* Compute for rows */ 
    blocks_on_dimension = ceil((float)matrix_rows / (float)block_rows);    

    block_per_process = blocks_on_dimension / process_grid_rows;   
    if (blocks_on_dimension % process_grid_rows != 0 && my_rank_coord_cols < blocks_on_dimension % process_grid_rows) {
            block_per_process++;
    } 

    /* Assign last block */
    int block_per_process_on_rows = block_per_process * block_rows; 

    int last_process_y = ((blocks_on_dimension % process_grid_rows) == 0) ? (process_grid_rows -1) : ((blocks_on_dimension % process_grid_rows) -1);  // 2
    if (my_rank_coord_cols == last_process_y && matrix_rows % block_rows != 0) {
        block_per_process_on_rows -= (block_rows - (matrix_rows % block_rows));
    }


    dimensions->rows = block_per_process_on_rows;
    dimensions->cols = block_per_process_on_cols;


#ifdef DEBUG
    printf("[Process %d] %d rows\n", my_rank, tot_nums_assigned_y);
    printf("[Process %d] %d cols\n", my_rank, tot_nums_assigned_x);
    printf("[Process %d] %d values\n", my_rank, tot_nums_assigned_y * tot_nums_assigned_x);
    puts("");
#endif


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
}

/**
 * @brief This function has the responsability to devide the matrix B in rows and give each group of rows
 * to the correct process in order to succesively compute the molutiplication between A block and B rows
 * @param matrix_file_name Name of the file containing A matrix in binary format
 * @param matrix_rows Number of rows of the matrix to be distributed
 * @param matrix_cols Number of columns of the matrix to be distributed
 * @param block_rows Number of rows of the block
 * @param process_gird_rows Number of rows of the process grid
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
        int process_gird_cols, 
        int my_rank, 
        MPI_Comm comm) {

    MPI_Datatype filetype;
    MPI_File file;
    MPI_Status status;

    int rank_norm = my_rank % process_gird_cols;
    int gsizes[2] = {K, M};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows, M_subarray}; // I want to have all the cols of the matrix, having only the correct rows 
    int psizes[2] = {process_gird_cols,1}; // We want to divide the matrix only in one dimension, by rows
    

    MPI_Type_create_darray(process_gird_cols, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file.\n");
        MPI_Abort(comm, 1);
    }


    MPI_File_set_view(file, 2*sizeof(int), MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_read_all(file, (*matrix), K_subarray*M_subarray, MPI_FLOAT, &status);


    MPI_File_close(&file);
}


void generate_rows_distribution_C(
        char *matrix_file_name, 
        float **matrix,
        int N, 
        int M, 
        int N_subarray,
        int M_subarray,
        int block_rows, 
        int process_gird_rows,
        int process_gird_cols, 
        int my_rank, 
        MPI_Comm comm) {

    MPI_Datatype filetype;
    MPI_File file;
    MPI_Status status;

    int rank_norm = my_rank / process_gird_cols;

    int gsizes[2] = {N, M};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {block_rows, N_subarray}; // I want to have all the cols of the matrix, having only the correct rows 
    int psizes[2] = {process_gird_rows,1}; // We want to divide the matrix only in one dimension, by rows

    MPI_Type_create_darray(process_gird_rows, rank_norm, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);

                AUDIT{
            puts("HERE");
            fflush(stdout);
        }

    MPI_File_open(comm, matrix_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (file == MPI_FILE_NULL) {
        fprintf(stderr, "Error in opening file %s.\n",matrix_file_name);
        MPI_Abort(comm, 1);
    }

    MPI_File_set_view(file, 2*sizeof(int), MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
    
    MPI_File_read_all(file, (*matrix), N_subarray*M_subarray, MPI_FLOAT, &status);


    MPI_File_close(&file);
}