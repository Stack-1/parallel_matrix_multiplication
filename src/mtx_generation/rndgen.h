#ifndef __MTX__
#define __MTX__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include "../utils/file_io/stackio.h"
#include "../utils/logger/logger.h"
#include "../utils/string_formatter/formatter.h"


void generate_real_value_matrix(float *matrix,int rows, int cols);
void generate_zero_value_matrix(float *matrix,int rows, int cols);


#endif

