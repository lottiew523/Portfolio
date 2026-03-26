#!/bin/bash
#PBS -P ud25
#PBS -N BoNorm
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=8000MB
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ud25
#PBS -l wd
#PBS -j oe

module purge

module load intel-mpi
export OMP_NUM_THREADS=1

mpirun -np 4  ./SWE.exe > output.dat

## Before running the program do compilation by
## mpif90 -c -real-size 64 -double-size 64 mpi_parameters*.f90 && mpif90 -c -real-size 64 -double-size 64 MY_MPI.f90 && mpif90 -real-size 64 -double-size 64 -O2 -o SWE.exe mpi_main.f90 *.o
