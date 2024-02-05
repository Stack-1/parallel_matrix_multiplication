#include "stackio.h"

#ifndef DATA_DIR
#define DATA_DIR "../../data/"
#endif


/**
 * @brief Utilit function to create a directory
 * @return Return code of mkdir call
*/
static int create_dir(char *rel_path,mode_t mode) {
  return mkdir(rel_path, mode);
}

/**
 * @brief Utilit function to print matrix for debugging
*/
void print_matrix(float *matrix,int rows, int cols){
    for(int i = 0;i<rows;i++){
        for(int j = 0;j<cols;j++){
            printf("%f ",matrix[i*cols + j]);
        }
        puts("");
    }
}



/**
 * @brief Utitliy function used to write generated matrices to file
 * @param matrix The matrix to be written to file
 * @param rows Number of rows of the matrix
 * @param cols Number of columns of the matrix
 * @param matrix_name The name of the matrix we want to save in file (A/B/C)
 * @param is_squared Treu if the matrix we are saving is part of the square matrices computation
 * @return Number of bytes written in file
*/
int write_matrix_to_file(float *matrix, int N, int K, int M,int rows, int cols,char *matrix_name,bool is_square){
    int ret = 0;
    struct stat st;
    FILE *matrix_file;
    char file_name[64];
    char dir_name[64]; 
    char log_string[LOG_MESSAGE_SIZE];

    sprintf(dir_name,DATA_DIR);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    if(is_square){
        sprintf(&(dir_name)[strlen(dir_name)],"square/");
    }else{
        sprintf(&(dir_name)[strlen(dir_name)],"rectangular/");
    }

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(&(dir_name)[strlen(dir_name)],"matrix%dx%dx%d/",N,K,M);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    memcpy(file_name,dir_name,strlen(dir_name));

    sprintf(&(file_name)[strlen(dir_name)],"matrix_%s_%dx%d.bin",matrix_name,rows,cols);

    // Open or create file
	if((matrix_file=fopen(file_name, "w+"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}

    // Write rows number as first element of the binary file
    ret = fwrite(&rows,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }
    
    // Write cols number as second element of the binary file
    ret = fwrite(&cols,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }

    ret = fwrite(matrix,sizeof(float),rows*cols,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }
    

    fclose(matrix_file);


    memset(log_string,0,LOG_MESSAGE_SIZE);
    sprintf(log_string,"File %s populated correctly!\n",file_name);
	logger_info(log_string);

    return sizeof(float)*cols*rows + sizeof(int) + sizeof(int);
}
/**
 * @brief Utitliy function used to read generated matrices from file
 * @param matrix The matrix to be read from file
 * @param rows_expected Number of rows of the matrix
 * @param cols_expected Number of columns of the matrix
 * @param matrix_name The name of the matrix we want to save in file (A/B/C)
 * @param is_squared Treu if the matrix we are saving is part of the square matrices computation
*/
void read_matrix_from_file(float *matrix,int N, int K, int M,int rows_expected, int cols_expected, char *matrix_name,bool is_square){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    int rows = 0;
    int cols = 0;
    int ret = 0;

    sprintf(file_name,DATA_DIR);

    if(is_square){
        sprintf(&(file_name)[strlen(file_name)],"square/");
    }else{
        sprintf(&(file_name)[strlen(file_name)],"rectangular/");
    }

    sprintf(&(file_name)[strlen(file_name)],"matrix%dx%dx%d/matrix_%s_%dx%d.bin",N,K,M,matrix_name,rows_expected,cols_expected);
    printf("%s\n",file_name);

    // Open or create file
	if((matrix_file=fopen(file_name, "r"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_error(log_string);
		exit(EXIT_FAILURE);
	}

    // Read rows number as first element of the binary file
    ret = fread(&rows,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(rows != rows_expected){    
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }

    // Read cols number as first element of the binary file
    ret = fread(&cols,sizeof(int),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(cols != cols_expected){
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }


    ret = fread(matrix,sizeof(float),rows*cols,matrix_file);
    if(ret == 0 && errno == EOF){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error reading bytes from file %s to matrix! Error returned: %s\n",file_name, strerror(errno));
        logger_error(log_string); 
        exit(EXIT_FAILURE);
    }
    


    fclose(matrix_file);
}


/**
 * @brief Utitliy function used to read generated matrices from file for MPI_darray data type. 
 * The difference with the previous functions is that this one consider the files without the 
 * number of rows and columns saved at the beginning of it.
 * @param matrix The matrix to be read from file
 * @param rows_expected Number of rows of the matrix
 * @param cols_expected Number of columns of the matrix
 * @param matrix_name The name of the matrix we want to save in file (A/B/C)
 * @param is_squared Treu if the matrix we are saving is part of the square matrices computation
*/
void read_matrix_from_file_mpi(float *matrix,int N, int K, int M,int rows_expected, int cols_expected, char *matrix_name,bool is_square){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    int rows = rows_expected;
    int cols = cols_expected;
    int ret = 0;

    sprintf(file_name,DATA_DIR);

    if(is_square){
        sprintf(&(file_name)[strlen(file_name)],"square/");
    }else{
        sprintf(&(file_name)[strlen(file_name)],"rectangular/");
    }

    sprintf(&(file_name)[strlen(file_name)],"matrix%dx%dx%d/matrix_%s_%dx%d.bin",N,K,M,matrix_name,rows_expected,cols_expected);

    // Open or create file
	if((matrix_file=fopen(file_name, "r"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_error(log_string);
		exit(EXIT_FAILURE);
	}

    ret = fread(matrix,sizeof(float),rows*cols,matrix_file);
    if(ret == 0 && errno == EOF){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error reading bytes from file %s to matrix! Error returned: %s\n",file_name, strerror(errno));
        logger_error(log_string); 
        exit(EXIT_FAILURE);
    }
    


    fclose(matrix_file);
}

/**
 * @brief Utitlity function to save times obtained in matrix generation to a txt file
 * @param matrix_name The name of the matrix we are considering (A/B/C/Total time)
 * @param rows Number of rows of the matrix
 * @param cols Number of columns of the matrix
 * @param time Time of generating the matrix in seconds
*/
void write_times_to_txt_file(char *matrix_name, int rows, int cols, double time){
    int ret = 0;
    struct stat st;
    FILE *matrix_file;
    char *file_name = (char *)"data/times.txt";
    char *dir_name = (char *)"data/"; 
    char log_string[LOG_MESSAGE_SIZE];

    // Create directory if needed
    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }


    // Open or create file
	if((matrix_file=fopen(file_name, "a"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}


    fprintf(matrix_file,"Matrix %s, %dx%d %f sec\n",matrix_name,rows,cols,time);

    fclose(matrix_file);
}


void write_sequential_computation_csv(int N,int K,int M,double read_time, double comptation_time,double write_time,double gflops){
    int ret = 0;
    struct stat st;
    FILE *matrix_file;
    char *file_name = (char *)"data/stats.csv";
    char *dir_name = (char *)"data/"; 
    char log_string[LOG_MESSAGE_SIZE];

    // Create directory if needed
    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }


    // Open or create file
	if((matrix_file=fopen(file_name, "a"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}

    fprintf(matrix_file,"%d,%d,%d,%f,%f,%f,%f\n",N,K,M,read_time,write_time,comptation_time,gflops);

    fclose(matrix_file);
}


void write_mpi_stats(int N, int K, int M,double total_time,float max_diff,float max_rel_diff,float gflops){
    int ret = 0;
    struct stat st;
    FILE *matrix_file;
    char *file_name = (char *)"data/stats.csv";
    char *dir_name = (char *)"data/"; 
    char log_string[LOG_MESSAGE_SIZE];

    // Create directory if needed
    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); 
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }


    // Open or create file
	if((matrix_file=fopen(file_name, "a"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}

    fprintf(matrix_file,"%d,%d,%d,%f,%f,%f,%f\n",N,K,M,total_time,max_diff,max_rel_diff,gflops);

    fclose(matrix_file);
}






void test(){
    int rows = 16;
    int cols = 16;
    float *matrix;
    float *matrix_recv;

    matrix = (float *)malloc(rows*cols*sizeof(float));
    matrix_recv = (float *)malloc(rows*cols*sizeof(float));
    

    
    for(int i = 0;i<rows;i++){
        for(int j = 0;j<cols;j++){
            matrix[i*cols + j] = i*cols + j;
        }
    }

#ifdef DEBUG    
    puts("");
    print_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix,rows,rows,cols, rows, cols,"A",true);
#endif

    read_matrix_from_file(matrix_recv,rows, rows, cols,rows,cols,(char *)"A",true);

    print_matrix(matrix_recv,rows,cols);
}