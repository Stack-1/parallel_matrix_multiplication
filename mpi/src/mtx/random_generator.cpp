#include <random>
#include <iostream>
#include "rndgen.h"

using namespace std;

constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 100;

void generate_real_value_matrix(float **matrix,size_t rows, size_t cols){
    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<float> dist(FLOAT_MIN,FLOAT_MAX); // distribute results between FLOAT_MIN and FLOAT_MAX inclusive.
    
    
    for(size_t i = 0;i<rows;i++){
        for(size_t j = 0;j<cols;j++){
            matrix[i][j] = dist(gen);
        }
    }
}

void generate_zero_value_matrix(float **matrix,size_t rows, size_t cols){
    for(size_t i = 0;i<rows;i++){
        for(size_t j = 0;j<cols;j++){
            matrix[i][j] = 0;
        }
    }
}


