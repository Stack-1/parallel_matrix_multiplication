#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "stackio.h"
#include "../logger/logger.h"


static int create_dir(char *rel_path,mode_t mode) {
  return mkdir(rel_path, mode);
}


void print_matrix(float **matrix,size_t rows, size_t cols){
    for(size_t i = 0;i<rows;i++){
        for(size_t j = 0;j<cols;j++){
            printf("%f ",matrix[i][j]);
        }
        puts("");
    }
}




// Return bytes written in file
size_t write_matrix_to_file(float **matrix, size_t rows, size_t cols, char *matrix_name){
    int ret = 0;
    struct stat st;
    //mode_t mode = st.st_mode & (S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP);
    FILE *matrix_file;
    char file_name[63];
    char dir_name[32]; 
    char log_string[LOG_MESSAGE_SIZE];
    size_t r[1]; 
    size_t c[1];
    r[0] = rows;
    c[0] = cols;


    sprintf(dir_name,"src/data/");

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(dir_name,"src/data/matrix%ldx%ld/",rows,cols);

    if ((ret = stat(dir_name, &st)) == -1) { // If directory not found
        create_dir(dir_name,0770); // TODO: Check for the right parametric bitmask
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Created data directory, with permissions: -rw-r--r-- %s\n",dir_name);
        logger_info(log_string);
    }

    sprintf(file_name,"src/data/matrix%ldx%ld/matrix_%s_%ldx%ld.bin",rows,cols,matrix_name,rows,cols);

    // Open or create file
	if((matrix_file=fopen(file_name, "w+"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
		exit(EXIT_FAILURE);
	}

    // Write rows number as first element of the binary file
    ret = fwrite(r,sizeof(size_t),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }
    
    // Write cols number as second element of the binary file
    ret = fwrite(c,sizeof(size_t),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_info(log_string);
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0;i<rows;i++){
        ret = fwrite(matrix[i],sizeof(float),cols,matrix_file);
        if(ret == 0){
            memset(log_string,0,LOG_MESSAGE_SIZE);
            sprintf(log_string,"Error opening file %s\n",file_name);
		    logger_info(log_string);
            exit(EXIT_FAILURE);
        }
    }

    fclose(matrix_file);


    memset(log_string,0,LOG_MESSAGE_SIZE);
    sprintf(log_string,"File %s populated correctly!\n",file_name);
	logger_info(log_string);

    return sizeof(float)*cols*rows + sizeof(size_t) + sizeof(size_t);
}


void read_matrix_from_file(float **matrix,size_t rows_expected, size_t cols_expected, char *matrix_name){
    FILE *matrix_file;
    char file_name[64];
    char log_string[LOG_MESSAGE_SIZE];
    size_t rows[1] = {0};
    size_t cols[1] = {0};
    int ret = 0;

    sprintf(file_name,"src/data/matrix%ldx%ld/matrix_%s_%ldx%ld.bin",rows_expected,cols_expected,matrix_name,rows_expected,cols_expected);

    // Open or create file
	if((matrix_file=fopen(file_name, "r"))==NULL) {
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error opening file %s\n",file_name);
		logger_error(log_string);
		exit(EXIT_FAILURE);
	}


    // Read rows number as first element of the binary file
    ret = fread(rows,sizeof(size_t),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(rows[0] != rows_expected){    
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }

    // Read cols number as first element of the binary file
    ret = fread(cols,sizeof(size_t),1,matrix_file);
    if(ret == 0){
        memset(log_string,0,LOG_MESSAGE_SIZE);
        sprintf(log_string,"Error writing bytes on file %s\n",file_name);
        logger_error(log_string); //TODO: Check with errno
        exit(EXIT_FAILURE);
    }else if(cols[0] != cols_expected){
        logger_error("Matrix passed to the function must be choerent with the size of the matrix in the file!\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0;i<rows[0];i++){
        ret = fread(matrix[i],sizeof(float),cols[0],matrix_file);
        if(ret == 0){
            memset(log_string,0,LOG_MESSAGE_SIZE);
            sprintf(log_string,"Error reading bytes from file %s to matrix!\n",file_name);
            logger_error(log_string); //TODO: Check with errno
            exit(EXIT_FAILURE);
        }
    }


    fclose(matrix_file);
}



void test(){
    size_t rows = 16;
    size_t cols = 16;
    float *matrix[rows];
    float *matrix_recv[rows];

    for(size_t i = 0;i<rows;i++){
        matrix[i] = (float *)malloc(cols*sizeof(float));
        matrix_recv[i] = (float *)malloc(cols*sizeof(float));
    }

    
    for(size_t i = 0;i<rows;i++){
        for(size_t j = 0;j<cols;j++){
            matrix[i][j] = i*cols + j;
        }
    }

        
    /*puts("");
    print_matrix(matrix,rows,cols);

    write_matrix_to_file(matrix, rows, cols,"A");*/

    read_matrix_from_file(matrix_recv,rows,cols,(char *)"A");

    print_matrix(matrix_recv,rows,cols);
}
