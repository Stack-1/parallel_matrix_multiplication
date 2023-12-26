#include "logger.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void logger(const char* tag, const char* message) {
   time_t now;
   time(&now);
   printf("%s [%s]: %s\n", ctime(&now), tag, message);
}


void logger_info(const char* message) {
    logger((char *)"INFO",message);
}


void logger_error(const char* message) {
    logger((char *)"ERROR",message);
}


void logger_debug(const char* message) {
    logger((char *)"DEBUG",message);
}