#include "mpi_computation.h"


/** @brief This function computes the matrix multiplication between A and B, and saves the result in the array C_temp.
 * 
*/
void mpi_matrix_multiplication(float *matrix_A, float *matrix_B, float **matrix_C, int N, int K, int M, int my_rank) {
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
 * @param N Number of rows of the matrix A
 * @param K Number of cols of the matrix A
 * @param M Number of cols of the matrix B
 * @param process_grid_rows Dimension of processes grid rows 
 * @param process_grid_cols Dimension of processes grid columns 
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

#ifdef DEBUG
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
#endif

    // First of all we have to compute sub-matrix sizes for each process according to block cyclic distribution and save them
    matrix_sizes = compute_submatrix_dimension(N, K, BLOCK_ROWS, BLOCK_COLS, process_grid_rows, process_grid_cols, my_rank);

    /* Save the results in the struct defined in the header file*/
    N_subarray = matrix_sizes->rows;
    K_subarray = matrix_sizes->cols;
    M_subarray = M;

    /* Allocate memory for each sub-matrix of the process */
    matrix_A_subarray = (float *)malloc(N_subarray * K_subarray * sizeof(float)); 
    matrix_B_subarray = (float *)malloc(K_subarray * M_subarray * sizeof(float));
    matrix_C_subarray = (float *)malloc(N_subarray * M_subarray * sizeof(float));
    matrix_C_subarray_rows = (float *)malloc(N_subarray * M_subarray * sizeof(float)); // Matrix needed to reducce the partial results obtained by each process
    memset(matrix_C_subarray,0,N_subarray * M_subarray * sizeof(float));
    matrix_C_subarray_from_file = (float *)malloc(N_subarray * M_subarray * sizeof(float)); // Matrix needed to store the sub-matrix of C read from file

#ifdef DEBUG
    AUDIT
        logger_info("Memory correctly allocatedand initialized!");
#endif

    /* Get sub-matrix of A according to block cyclic distribution */
    generate_block_cyclic_distribution_matrix_A(matrix_A_file_name,&matrix_A_subarray ,N, K, N_subarray, K_subarray,BLOCK_ROWS, BLOCK_COLS,process_grid_rows,process_grid_cols,my_rank,comm);

    /* Get sub-matrix of B according to rows cyclic distribution */
    generate_rows_distribution(matrix_B_file_name,&matrix_B_subarray, K, M, K_subarray , M_subarray, BLOCK_ROWS, process_grid_cols, my_rank,comm);

    /* Compute the actual multiplication, giving each process a partial result of the final computation and a partial view of C by rows*/
    mpi_matrix_multiplication(matrix_A_subarray,matrix_B_subarray,&matrix_C_subarray, N_subarray, K_subarray ,M_subarray,my_rank);

    /* Give the first process of each row in the process grid all the results of the processes in the row of the process grid*/
    reduce_partial_results(matrix_C_subarray, &matrix_C_subarray_rows,N_subarray, M_subarray, process_grid_cols, my_rank, comm);

    MPI_Barrier(comm); // Wait for all process to finish

    
    /* Read C values from file to add them to the computation */
    generate_rows_distribution_C(matrix_C_file_name,&matrix_C_subarray_from_file, N, M, N_subarray , M_subarray, BLOCK_ROWS, process_grid_rows, process_grid_cols,my_rank,comm);
        AUDIT{
            puts("HERE");
            fflush(stdout);
        }
    MPI_Barrier(comm);

    if (my_rank % process_grid_cols == 0) {
        for (int i = 0; i < N_subarray * M_subarray; ++i) {
            matrix_C_subarray_rows[i] += matrix_C_subarray_from_file[i];
        }
    }

    MPI_Barrier(comm);

    /* Get a new communicator to write on file using only the first processes of the param grid row */
    MPI_Comm first_element_comm;
    MPI_Comm_split(comm, my_rank % process_grid_cols, my_rank, &first_element_comm);

    int first_element_rank;
    MPI_Comm_rank(first_element_comm, &first_element_rank);

     /* Write C on file */
    if (my_rank % process_grid_cols == 0) {
        char file_name[128];
        if(N==K && K==M){
            sprintf(file_name, "../../data/square/matrix%dx%dx%d/matrix_C_mpi_%dx%d.bin",N,K,M,N,M);
        }else{
            sprintf(file_name, "../../data/rectangular/matrix%dx%dx%d/matrix_C_mpi_%dx%d.bin",N,K,M,N,M);  
        }

        write_result_to_file(file_name, matrix_C_subarray_rows, N_subarray * M_subarray, N, M, BLOCK_ROWS, process_grid_rows, process_grid_cols, my_rank,first_element_comm);

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
        sprintf(matrix_A_file_name,"../../data/square/matrix%dx%dx%d/matrix_A_%dx%d.bin", N, K, M, N, K);
        sprintf(matrix_B_file_name,"../../data/square/matrix%dx%dx%d/matrix_B_%dx%d.bin",N,K, M, K, M);
        sprintf(matrix_C_file_name,"../../data/square/matrix%dx%dx%d/matrix_C_%dx%d.bin",N,K, M, N, M);
    }else{
        sprintf(matrix_A_file_name,"../../data/rectangular/matrix%dx%dx%d/matrix_A_%dx%d.bin", N, K, M, N, K);
        sprintf(matrix_B_file_name,"../../data/rectangular/matrix%dx%dx%d/matrix_B_%dx%d.bin",N, K, M, K, M);
        sprintf(matrix_C_file_name,"../../data/rectangular/matrix%dx%dx%d/matrix_C_%dx%d.bin",N, K, M, N, M); 
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
        // Compare sequential and mpi result
        float *matrix_c_sequential;
        float *mpi;
        float max_diff = 0.0f;
        float max_rel_diff = 0.0f;


        matrix_c_sequential = (float *)malloc(sizeof(float)*N*M);
        mpi = (float *)malloc(sizeof(float)*N*M);

        read_matrix_from_file(matrix_c_sequential,N,K,M,N,M,(char *)"sequential_C",is_square);
        read_matrix_from_file_mpi(mpi,N,K,M,N,M,(char *)"C_mpi",is_square);
        
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

                max_diff = max( abs(mpi[i*M+j] - matrix_c_sequential[i*M+j]), max_diff);
                max_rel_diff = max(  
                    (float) abs(mpi[i*M+j] - matrix_c_sequential[i*M+j]) /  (float)(max(mpi[i*M+j], matrix_c_sequential[i*M+j]) ), 
                    max_rel_diff
                    );
            }
        }
        fflush(stdout);

        puts("");
        
        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Max diff is: %f\n",max_diff);
	    logger_info(logger_message);

        memset(logger_message,0,LOG_MESSAGE_SIZE);
	    sprintf(logger_message,"Max rel diff is: %f\n",max_rel_diff);
	    logger_info(logger_message);

        // Stats computation and saving
        float gflops = 0.0f;
        float partial = (total_time != 0.0f) ? (float)(2.0*N*K*M)/total_time : (float)(2.0*N*K*M)/0.0000000001;
	    gflops = (float)( partial / 1000000000);

        write_mpi_stats(N,K,M,total_time,max_diff,max_rel_diff,gflops);
        
        
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

    

    MPI_Finalize();

	return 0;
}