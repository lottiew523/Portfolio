rm -r *.bin
rm -r *.gif
rm -r *.mod
rm -r *.o
gfortran -c mpi_parameters.f90
mpif90 -c MY_MPI.f90
mpif90 -o SWE.exe mpi_main.f90 *.o
mpirun -np 4 SWE.exe
echo Simulation Complete!
python3 view.py --dt=0.02