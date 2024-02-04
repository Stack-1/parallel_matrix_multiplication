#include "rndgen.h"


#ifndef LOGGER_MESSAGE_SIZE
#define LOGGER_MESSAGE_SIZE 256
#endif

#ifndef FORMATTED_STRING_SIZE
#define FORMATTED_STRING_SIZE 64
#endif

int main(int argc, char *argv[]){
    int N = 0;
    int K = 0;
    int M = 0;
    char *matrix_A_name = (char *)"A";
    char *matrix_B_name = (char *)"B";
    char *matrix_C_name = (char *)"C";
    double cpu_time_used_A = 0.0;
    double cpu_time_used_B = 0.0;
    double cpu_time_used_C = 0.0;
	double total_time = 0.0;
    struct timeval stop, start;
    char logger_buffer[LOGGER_MESSAGE_SIZE];
    char formatted_string[FORMATTED_STRING_SIZE];
    bool is_square;

    if(argc != 4){
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[1], "%d", &N) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[2], "%d", &K) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    if (sscanf (argv[3], "%d", &M) != 1) {
        fprintf(stderr, "[ERROR] - not an integer");
        logger_error("Program should be called like ./<elf file name> <N> <K> <M>");
        exit(EXIT_FAILURE);
    }

    is_square = (N==K && K==M) ? true : false;

	puts("");
	puts("---------------------------------------------------------------------------[START]---------------------------------------------------------------------------");



    float *matrix_A;
    float *matrix_B;
    float *matrix_C;

    // Allocate memory
    matrix_A = (float *)malloc(sizeof(float)*N*K);
    matrix_B = (float *)malloc(sizeof(float)*K*M);
    matrix_C = (float *)malloc(sizeof(float)*N*M);

    


    puts("MATRIX A:");
    // Generate matrix
    gettimeofday(&start, NULL);
    generate_real_value_matrix(matrix_A,N,K);
    write_matrix_to_file(matrix_A,N,K,M,N,K,matrix_A_name,is_square);
	gettimeofday(&stop, NULL);
	cpu_time_used_A = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += cpu_time_used_A;

    puts("MATRIX B:");

    gettimeofday(&start, NULL);
    generate_real_value_matrix(matrix_B,K,M);
    write_matrix_to_file(matrix_B,N,K,M,K,M,matrix_B_name,is_square);
	gettimeofday(&stop, NULL);
	cpu_time_used_B = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += cpu_time_used_B;


    puts("MATRIX C:");

    gettimeofday(&start, NULL);
    generate_real_value_matrix(matrix_C,N,M);
    write_matrix_to_file(matrix_C,N,K,M,N,M,matrix_C_name,is_square);
	gettimeofday(&stop, NULL);
	cpu_time_used_C = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
	total_time += cpu_time_used_C;

#ifdef DEBUG
    test();
#endif

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(cpu_time_used_A,(char *)formatted_string);
    memset(logger_buffer,0,LOG_MESSAGE_SIZE);
    sprintf(logger_buffer,"Matrix A generation time: %s",formatted_string);
    logger_info(logger_buffer);

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(cpu_time_used_B,(char *)formatted_string);
    memset(logger_buffer,0,LOG_MESSAGE_SIZE);
    sprintf(logger_buffer,"Matrix B generation time: %s",formatted_string);
    logger_info(logger_buffer);

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(cpu_time_used_C,(char *)formatted_string);
    memset(logger_buffer,0,LOG_MESSAGE_SIZE);
    sprintf(logger_buffer,"Matrix C generation time: %s",formatted_string);
    logger_info(logger_buffer);

	memset(formatted_string,0,FORMATTED_STRING_SIZE);
	getFormattedTime(total_time,(char *)formatted_string);
    memset(logger_buffer,0,LOG_MESSAGE_SIZE);
    sprintf(logger_buffer,"Total generation time: %s",formatted_string);
    logger_info(logger_buffer);


    write_times_to_txt_file((char *)"A",N,K,cpu_time_used_A);
    write_times_to_txt_file((char *)"B",K,M,cpu_time_used_B);
    write_times_to_txt_file((char *)"C",N,M,cpu_time_used_C);
    write_times_to_txt_file((char *)"Total time",N,M,total_time);


    puts("----------------------------------------------------------------------------[END]----------------------------------------------------------------------------");
	puts("");
}