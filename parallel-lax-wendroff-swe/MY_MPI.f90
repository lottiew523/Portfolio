MODULE MY_MPI
   
	USE mpi_parameters

   IMPLICIT NONE
    INCLUDE 'mpif.h'
	
	
	! initialize variables
	INTEGER :: pid, nprocs, ierr
    INTEGER :: local_nx, local_ny
    INTEGER :: CART_COMM
	INTEGER :: cart_rank, coord(2)
	INTEGER :: ndims, dims(2), shift, dimx, dimy
	INTEGER	:: left, right, top, bottom
	INTEGER :: ROW_TYPE
	LOGICAL :: periodic(2), reorder
	
	CONTAINS
	
!-------------------------- init_mpi subroutine -------------------------------
    SUBROUTINE init_mpi()
		! This subroutine sets up a 2D Cartesian topology for domain decomposition and enables periodic boundaries in both directions. 
		! The subroutine also creates a custom MPI datatype (ROW_TYPE) for efficient row communication since fortran is column major so
		! the row communication is non-contiguous and therefore need MPI_Type_vector to describe the correct pattern.
		
		
		! Initializes MPI and gets rank and size.
        CALL MPI_Init(ierr)
        CALL MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
        CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
		
		! Initialize Cartesian Mesh/ Grid Topology
		ndims = 2							! Number of dimensions
		dims = (/0, 0/)  					! Dimension vector. let MPI calculate the correct dimension
		periodic = (/ .true., .true. /)		! Boundaries in x- and y-direction are periodic
		reorder = .true.					! MPI can change the rank numbering
		CALL MPI_Dims_create(nprocs, ndims, dims, ierr)
		!coord(1) = 0


		CALL MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, CART_COMM, ierr)
		CALL MPI_Comm_rank(CART_COMM, cart_rank, ierr)	!Obtain new rank in new communicator
		CALL MPI_Cart_coords(CART_COMM, cart_rank, ndims, coord, ierr)
		
	
		! Obtain rank of processors to the different directions
		dimx  = 0 														! shift along the x-axis to find left/right neighbors.
		dimy  = 1 														! shift along the y-axis to find top/bottom neighbors
		shift = 1 														! Specifies how far to move along specific dimension. 1 means one step (one neighbor away). If 2 = neighbors two blocks away.
		CALL MPI_Cart_shift(CART_COMM, dimx, shift, left, right, ierr)  ! x-direction
		CALL MPI_Cart_shift(CART_COMM, dimy, shift, bottom, top, ierr)	! y-direction
	
		
		! Example local grid size (can be adjusted or passed in)

		
		! Create row type for column communication
 		! CALL MPI_Type_free(ROW_TYPE, ierr)
!		CALL MPI_Type_vector(local_ny, 1, local_nx+2, MPI_DOUBLE_PRECISION, ROW_TYPE, ierr)
!		CALL MPI_Type_commit(ROW_TYPE, ierr)
		
		
    END SUBROUTINE init_mpi


!-------------------------- exhange_ghost_points SUBROUTINE -------------------------------	
		
	SUBROUTINE variable_communication(var)
	
		REAL(KIND=8), INTENT(INOUT) :: var(0:local_nx+1, 0:local_ny+1)
		INTEGER :: req(8), reqs(2), stat(MPI_STATUS_SIZE,8), ReqNr
		INTEGER :: tagL, tagR, tagB, TagT
		INTEGER :: ierr
		
		! Initialize tag (this ensures the correct send/receive pair matches.)
		tagL = 100
		tagR = 101
		TagT = 102
		TagB = 103
	
		ReqNr = 8

	
		!-------- Non-blocking communication for exchanging ghost points -----------
		! TODO: Chck if the communication can work for several variables in one or if it has to be done seperately.
		
		! Top/Bottom communication
        CALL MPI_Irecv(var(1, 			 0), local_nx, MPI_DOUBLE_PRECISION, top, tagT, CART_COMM, req(1), ierr) 		      	! neighbor above will send its bottom interior row to fill ghost row j=0
        CALL MPI_Isend(var(1, 			 1), local_nx, MPI_DOUBLE_PRECISION, top, tagB, CART_COMM, req(2), ierr)				! send first interior row (j=1) to top neighbor.
        CALL MPI_Irecv(var(1, local_ny + 1), local_nx, MPI_DOUBLE_PRECISION, bottom, tagB, CART_COMM, req(3), ierr) ! Receive into top ghost row from bottom neighbor. Neighbor sends its top interior row.
        CALL MPI_Isend(var(1, 	  local_ny), local_nx, MPI_DOUBLE_PRECISION, bottom, tagT, CART_COMM, req(4), ierr) 	! Send top interior row to bottom neighbor, who will store it in its bottom ghost (j=0).

        ! Left/Right communication
        CALL MPI_Irecv(var(0, 		 	 1), 1, ROW_TYPE, left, tagL, CART_COMM, req(5), ierr) 			! receive into left ghost column (i=0). data from left neighbor sends its rightmost interior column
        CALL MPI_Isend(var(local_nx, 	 1), 1, ROW_TYPE, right, tagL, CART_COMM, req(6), ierr) 	! Send rightmost interior column (i=local_nx) to right neighbor.
        CALL MPI_Irecv(var(local_nx + 1, 1), 1, ROW_TYPE, right, tagR, CART_COMM, req(7), ierr)	! receive into right ghost column from right neighbor’s left interior
        CALL MPI_Isend(var(1, 			 1), 1, ROW_TYPE, left, tagR, CART_COMM, req(8), ierr) 			! send left interior column (i=1) to left neighbor

		
		! Wait for the nonblocking to complete before using the ghost points
        CALL MPI_WAITALL(8, req, stat, ierr)
		
    END SUBROUTINE variable_communication



	! ------------------- exchange ghost point for each variable --------------------

	SUBROUTINE exchange_ghost_points(eta, u, v)
		
		REAL(KIND=8), INTENT(INOUT) :: eta(0:local_nx+1, 0:local_ny+1)
		REAL(KIND=8), INTENT(INOUT) :: u  (0:local_nx+1, 0:local_ny+1)
		REAL(KIND=8), INTENT(INOUT) :: v  (0:local_nx+1, 0:local_ny+1)

		
		CALL variable_communication(eta)
		CALL variable_communication(u)
		CALL variable_communication(v)
		
	END SUBROUTINE exchange_ghost_points



	!-------------------------------- Parallel I/O ----------------------------------
	 
	SUBROUTINE write_var_parallel_io(filename, local_var)
		! Sources: https://www.cscs.ch/fileadmin/user_upload/contents_publications/tutorials/fast_parallel_IO/MPI-IO_NS.pdf 
		!		   https://wgropp.cs.illinois.edu/courses/cs598-s16/lectures/lecture32.pdf
		
		! Collective parallel I/O approach
			!Independent: MPI_File_write_at (each rank writes alone)
			!Collective: MPI_File_write_at_all (We use this)
		
		
		CHARACTER(*), INTENT(IN) :: filename
		REAL(KIND = 8), INTENT(IN) :: local_var(1:local_nx, 1:local_ny) ! Local var without ghost points 
																	   
		INTEGER :: fh, filetype, ierr
		INTEGER :: gsizes(2), lsizes(2), starts(2)
		INTEGER(KIND = MPI_OFFSET_KIND) :: disp

		! Define global and local sizes for the subarray                                             
		gsizes = (/ nx, ny /)								! Global grid dimensions
		lsizes = (/ local_nx, local_ny /)					! Local block dimensions
		starts = (/ coord(1)*local_nx, coord(2)*local_ny /)	! Starting indices for this rank's block
	  
		! Create a derived datatype representing this rank's subarray in the global file
		CALL MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
		CALL MPI_Type_commit(filetype, ierr)

		! Open the shared file for writing (collectively)
		disp = 0_MPI_OFFSET_KIND
		CALL MPI_File_open(CART_COMM, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
		
		! Set the file view so each rank writes to its correct position
		CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, filetype, 'native', MPI_INFO_NULL, ierr)

		! Collective write: all ranks participate
		CALL MPI_File_write_all(fh, local_var, local_nx*local_ny, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

		! Close file and free datatype
		CALL MPI_File_close(fh, ierr)
		CALL MPI_Type_free(filetype, ierr)

		! Print confirmation from rank 0
		IF (pid == 0) WRITE(*,*) 'Parallel-IO: wrote file ', TRIM(filename)
	END SUBROUTINE write_var_parallel_io



    !-------------------------- finalize_mpi subroutine -------------------------------
	 

END MODULE MY_MPI
