Shallow Water Equations Solver

This project was developed as part of academic coursework and/or collaborative work. 
Certain components may reflect group contributions. The code is provided for 
demonstration purposes only and is not licensed for reuse or distribution.
		
#----------------------------------------------------------------------------

Modified 26/03/2026 by Charlotte Woods - view.py, bug fixes to I/O, runHPC.sh, run.sh

## My Contributions

- Parallelised serial Lax-Wendroff Solver (finite volume shallow water equations, indexing conventions)
- Wrote post-processing and visualisation scripts
- Debugging and verification (stability, CFL condition testing, speedup)

#----------------------------------------------------------------------------


# Files in this package
	mpi_parameters.f90
	MY_MPI.f90
	mpi_main.f90
	run.sh
	view.py
	runHPC.sh
	
This package contains code for solving the Shallow Water Equations using a Lax-Wendroff solver.
The package solves the transient state of a 2D area of size Lx by Ly with no external forces applies.

# Instructions on Running

Method 1 (recommended)
	
	All physical conditions of the code can be changed in mpi_parameters.f90 .
	
	
	
	LOCAL: Use 
	
	./run.sh 
	
	to run suite
	
	
	
	HPC: Use 
	
	mpif90 -c -real-size 64 -double-size 64 mpi_parameters*.f90 && mpif90 -c -real-size 64 -double-size 64 MY_MPI.f90 && mpif90 -real-size 64 -double-size 64 -O2 -o SWE.exe mpi_main.f90 *.o
	
	Followed by 
	
	qsub runHPC.sh
	
	Note that for HPC all files must exist in same WD, and ncpus must be pinned matching -np flag. e.g.
	!				#PBS -l ncpus=4
	!				mpirun -np 4  ./SWE.exe > output.dat
	

Method 2: manual compile (for debugging or partial runs)

	
	# Parameter set
	All physical conditions of the code can be changed in parameters, including start and end time.

	To compile run the following in a folder with all .f90 files:
		gfortran -c mpi_parameters.f90
	!		mpif90 -c MY_MPI.f90
	!		mpif90 -o SWE.exe mpi_main.f90 *.o
	To run the code use the following command in the same folder (-np 4 can be changed to the required number of processors):
	!		mpirun -np 4 SWE.exe

	To create gif:
	!		python3 view.py

