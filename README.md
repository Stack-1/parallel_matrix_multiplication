# Prallel Matrix Multiplication
This repository contains the code implemented for the SCPA(Sistemi di Calcolo Parallelo e Applicazioni - PCSA Parallel Computing Systems and Applications) class. The goal of this project is to use MPI, OpnMP and CUDA to measure the performances of a matrix by matrix multiplication. 

To excecute the code:
- Download the source code
- Move to src/mtx_generation/ and excecute make generate-random to store the matrix generations on file
- Move to src/sequential/ and excecute make sequential for the sequential computation
- Move to src/mpi/, excecute module load mpi and excecute make process-parallel for the mpi computation
- Move to src/cuda/ , excecute cmake . to setup the enivronment, excecute make MatrixMultiplication to compile and run ./cuda_launcher.sh to excecute CUDA computation.

