#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "stackio.h"
#include "../logger/logger.h"



#ifndef DATA_DIR
#define DATA_DIR "../../data/"
#endif





static int create_dir(char *rel_path,mode_t mode) {
  return mkdir(rel_path, mode);
}


void print_matrix(double *matrix,int rows, int cols){
    for(int i = 0;i<rows;i++){
        for(int j = 0;j<cols;j++){
            printf("%f ",matrix[i*cols + j]);
        }
        puts("");
    }
}




// Return bytes written in file
int write_matrix_to_file(double *matrix, int rows, int cols, char *matrix_name){
    int ret = 0;
    struct stat st;
    //mode_t mode = st.st_mode & (S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
    FILE *matrix_file;
    char file_name[63];
    char dir_name[32]; 
    char log_string[LOG_MESSAGE_SIZE];

    sprintf(dir_name,DATA_DIR);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(&(dir_name)[strlen(DATA_DIR)],"matrix%dx%d/",rows,cols);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
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

    ret = fwrite(matrix,sizeof(double),rows*cols,matrix_file);
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

    return sizeof(double)*cols*rows + sizeof(int) + sizeof(int);
}


void read_matrix_from_file(double *matrix,int rows_expected, int cols_expected, char *matrix_name){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    int rows = 0;
    int cols = 0;
    int ret = 0;

    sprintf(file_name,DATA_DIR);
    sprintf(&(file_name)[strlen(DATA_DIR)],"matrix%dx%d/matrix_%s_%dx%d.bin",rows_expected,cols_expected,matrix_name,rows_expected,cols_expected);

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


    ret = fread(matrix,sizeof(double),rows*cols,matrix_file);
    if(ret == 0 && errno == EOF){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error reading bytes from file %s to matrix! Error returned: %s\n",file_name, strerror(errno));
        logger_error(log_string); 
        exit(EXIT_FAILURE);
    }
    


    fclose(matrix_file);
}

void read_matrix_from_file_mpi(double *matrix,int rows_expected, int cols_expected, char *matrix_name){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    int rows = rows_expected;
    int cols = cols_expected;
    int ret = 0;

    sprintf(file_name,DATA_DIR);
    sprintf(&(file_name)[strlen(DATA_DIR)],"matrix%dx%d/matrix_%s_%dx%d.bin",rows_expected,cols_expected,matrix_name,rows_expected,cols_expected);

    // Open or create file
	if((matrix_file=fopen(file_name, "r"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_error(log_string);
		exit(EXIT_FAILURE);
	}

    ret = fread(matrix,sizeof(double),rows*cols,matrix_file);
    if(ret == 0 && errno == EOF){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error reading bytes from file %s to matrix! Error returned: %s\n",file_name, strerror(errno));
        logger_error(log_string); 
        exit(EXIT_FAILURE);
    }
    


    fclose(matrix_file);
}

void test(){
    int rows = 16;
    int cols = 16;
    double *matrix;
    double *matrix_recv;

    matrix = (double *)malloc(rows*cols*sizeof(double));
    matrix_recv = (double *)malloc(rows*cols*sizeof(double));
    

    
    for(int i = 0;i<rows;i++){
        for(int j = 0;j<cols;j++){
            matrix[i*cols + j] = i*cols + j;
        }
    }

#ifdef DEBUG    
    puts("");
    print_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix, rows, cols,"A");
#endif

    read_matrix_from_file(matrix_recv,rows,cols,(char *)"A");

    print_matrix(matrix_recv,rows,cols);
}