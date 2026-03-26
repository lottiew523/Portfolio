PROGRAM mpi_main

	use MY_MPI
	use mpi_parameters
	!use mpi
	
    IMPLICIT NONE

	! ---------------------- Declare Variables ----------------------- 
	integer :: iteration
	integer :: iters
	integer :: steps
	character(len=100) :: filename
	REAL(KIND = 8) :: local_min, local_max, global_min, global_max
	
	iteration = 0
	!step = 100
	steps = tmax / dt / step
	
	! -------------------- Precompute Parameters ---------------------
	! Starting and ending index for x- and y-direction
	iStart 	= 1
	iEnd 	= nx
	jStart 	= 1
	jEnd	= ny
	
	
	! Delta spacesteps
	dx = Lx/(nx)
	dy = Ly/(ny)

	! Precompute - delta spacesteps / delta timestep
	dt_dx = dt/dx
	dt_dy = dt/dy
	
	! Defining the points we want to print out
	pointx = nx/2	! point in x
	pointy = ny/2	! point in y

	

	! ------------------- init_mpi subroutine -----------------------
	! This subroutine sets up 2D Cartesian topology for domain decomposition and enables periodic boundaries in both directions. 
	! It also creates a custom MPI datatype (ROW_TYPE) for efficient row communication since fortran is column major.
	call init_mpi()
	
	if (pid == 0) then
		write(*,*) "Local Index Bounds:"
		write(*,*) "       PID:   Location:                Indices:"
		
		open(99, file='run_info.txt', status='replace')
	
		write(99,'(A,1X,F10.6)') 'dt', dt
		write(99,'(A,1X,I0)')     'nx', nx
		write(99,'(A,1X,I0)')     'ny', ny
		write(99,'(A,1X,F10.6)') 'tmax', tmax
		write(99,'(A,1X,I0)')     'step', step
	
		close(99)
	end if
	
	
	! ----- Define local subdomain indices for each MPI process -----
	! Example local grid size (can be adjusted or passed in)
	nLocX = nx / dims(1) ! nx ny num cells
	nLocY = ny / dims(2)
	
	if ((coord(1)) == 0) then
		iStart = 1          
		iEnd = nLocX  +  MOD(nx, dims(1))                   ! rank 0 block absorbs extra
	else
		iStart = coord(1) * nLocX +  MOD(nx, dims(1))+1	    ! first owned column
		iEnd = (coord(1)+1) * nLocX +  MOD(nx, dims(1))	    ! last owned column
	end if
	

	
	if ((coord(2)) == 0) then
		jStart = 1          
		jEnd = nLocY  +  MOD(ny, dims(2))                   ! rank 0 block absorbs extra
	else
		jStart = (coord(2)) * nLocY +  MOD(ny, dims(2))+1	! first owned column
		jEnd = (coord(2)+1) * nLocY +  MOD(ny, dims(2))	    ! last owned column
	end if
	
	local_ny = jEnd + 1 - jStart
	local_nx = iEnd + 1 - iStart


		


	! Index for ghost cell
	iGW	= iStart- 1
	iGE	= iEnd	+ 1
	jGS	= jStart- 1
	jGN	= jEnd	+ 1
	
	
	CALL MPI_Type_vector(local_ny, 1, local_nx+2, MPI_DOUBLE_PRECISION, ROW_TYPE, ierr)
	CALL MPI_Type_commit(ROW_TYPE, ierr)


		
	! ---------------------- Allocate Arrays -----------------------
	
	ALLOCATE(x(iStart:iEnd), y(jStart:jEnd))
	
	! Allocate height and velocity variables (with ghost cells)
	ALLOCATE(eta(iGW:iGE, jGS:jGN))
	ALLOCATE(u	(iGW:iGE, jGS:jGN))
	ALLOCATE(v	(iGW:iGE, jGS:jGN))
	
	! Allocate for updated height and velocity variables
	ALLOCATE(eta_new(iGW:iGE, jGS:jGN))
	ALLOCATE(u_new	(iGW:iGE, jGS:jGN))
	ALLOCATE(v_new	(iGW:iGE, jGS:jGN))
	
	! Allocate height flux arrays
	Allocate(ueta(iGW:iGE,JGS:JGN))
	allocate(veta(iGW:iGE,JGS:JGN))
	
	! Allocate momentum flux arrays
	Allocate(Ux(iGW:iGE,JGS:JGN))
	allocate(Uy(iGW:iGE,JGS:JGN))
	Allocate(Vx(iGW:iGE,JGS:JGN))
	allocate(Vy(iGW:iGE,JGS:JGN))

	! Allocate staggered (midpoint) variables
	allocate(eta_mid_x(igw:iend, jstart:jend))
	allocate(eta_mid_y(istart:iend, jgs:jend))
	allocate(ueta_mid_x(igw:iend, jstart:jend))
    allocate(veta_mid_x(igw:iend, jstart:jend))
	allocate(ueta_mid_y(istart:iend, jgs:jend))
	allocate(veta_mid_y(istart:iend, jgs:jend))
	
	! Allocate midpoint momentum fluxes 
	allocate(Ux_mid(igw:iend, jstart:jend))
	allocate(Uy_mid(istart:iend, jgs:jend))
	allocate(Vx_mid(igw:iend, jstart:jend))
	allocate(Vy_mid(istart:iend, jgs:jend))
	
	
	! Allocate new timestep fluxes
	allocate(veta_new(iGW:iGE,JGS:JGN))
	allocate(ueta_new(iGW:iGE,JGS:JGN))
	
	
	
	! ---------------------- Initialise Arrays -----------------------
	eta = 0.0
	u	= 0.0
	v	= 0.0
	eta_new = 0.0
	u_new	= 0.0
	v_new	= 0.0
	
	! Position vectors
	DO j = jStart,jEnd
		DO i = iStart,iEnd
			x(i) = dx*(i-1)+dx/2
			y(j) = dy*(j-1)+dy/2
		END DO
	END DO

	
	! ---------------------- Apply intitial conditions -----------------------
	! Initial condition
	DO j = jGS,jGN
		DO i = iGW,iGE
			eta(i,j)= eta0 + amp * exp(-((x(i)-x0*Lx)**2 + (y(j)-y0*Ly)**2)/(2.0*(sigma*Lx)**2))! 0.1*sin(2*PI*x(i)/(Lx)) + eta0
			u(i,j)	= 0.0
			v(i,j)	= 0.0
		END DO
	END DO


	
	! ---------------------- Apply Boundary Conditions -----------------------
		! ---------------- Exchange ghost cells for eta, u, v -----------------
	CALL MPI_Barrier(CART_COMM, ierr)
		
		
		CALL variable_communication(eta)
		CALL variable_communication(u)
		CALL variable_communication(v)	

	! Debug: Print processor grid coordinates and local domain extents
	write(*,*) pid, coord(1), coord(2), iStart, iEnd, jstart, jend
	
	
	! ================================= START ITERATION LOOP ==================================
	! =========================================================================================
	DO WHILE (t < tmax)
	
		! Update time counter
		t = t + dt
		
	
		! ========================== Lax-Wendroff solver ============================
		! centre computed nodes
		ueta = u*eta	! x-dir height flux
		veta = v*eta	! y-dir height flux
		!
		Ux = ueta*u+0.5*g*eta**2		! EW-dir Flux: Advection (momentum flux + hydrostatic pressure)
		Uy = ueta*v						! Cross Flux: (Uy = eta*u*v)
		!! NS-dir Momentum flux
		Vx = Uy							! Cross Flux: (Vx = eta*v*u = eta*u*v = Uy)	{NOTE: same node so Uy = Vx }
		Vy = veta*v+0.5*g*eta**2; 		! NS-dir Flux: Advection (momentum flux + hydrostatic pressure)	
		
		
		! Calculate staggered nodes
		do j = jstart,jend
			DO i = iGW,iEnd
				eta_mid_x(i,j) = 0.5*((eta(i,j) + eta(i+1,j)) - dt_dx*(ueta(i+1,j) - ueta(i,j))) ! Height in staggered node: ew direction
				ueta_mid_x(i,j) = 0.5*((ueta(i,j) + ueta(i+1,j)) - dt_dx*(Ux(i+1,j) - Ux(i,j)))  ! ew-dir momentum in staggered node: ew grid direction
				veta_mid_x(i,j) = 0.5*((veta(i,j) + veta(i+1,j)) - dt_dx*(Vx(i+1,j) - Vx(i,j)))  ! ns-dir momentum in staggered node: ew grid direction
			END DO
		END DO

		
		do j = jgs,jend
			DO i = istart,iEnd
				eta_mid_y(i,j) = 0.5*((eta(i,j) + eta(i,j+1)) - dt_dy*(veta(i,j+1) - veta(i,j))) ! Height in staggered node: ns direction
				ueta_mid_y(i,j)  = 0.5*((ueta(i,j) + ueta(i,j+1)) - dt_dy*(Uy(i,j+1) - Uy(i,j))) ! ew-dir momentum in staggered node: ns grid direction
				veta_mid_y(i,j) = 0.5*((veta(i,j) + veta(i,j+1)) - dt_dy*(Vy(i,j+1) - Vy(i,j)))  ! ns-dir momentum in staggered node: ns grid direction
			END DO
		END DO

		
		! Use staggered node values to find new centred grid values
		! ew-dir flux
		Ux_mid = ueta_mid_x * ueta_mid_x/eta_mid_x + 0.5*g*eta_mid_x**2	! ew-dir Flux: Advection (momentum flux + hydrostatic pressure)
		Uy_mid = ueta_mid_y * veta_mid_y/eta_mid_y						! Cross Flux: (Uy = eta*u*v = ueta*veta/eta)
		
		! ns-dir flux
		Vx_mid = veta_mid_x * ueta_mid_x/eta_mid_x						! Cross Flux: (Vx = eta*u*v = ueta*veta/eta) {NOTE: staggered nodes so Uy(i,j) /= Vx(i,j)}
		Vy_mid = veta_mid_y * veta_mid_y/eta_mid_y + 0.5*g*eta_mid_y**2	! ns-dir Flux: Advection (momentum flux + hydrostatic pressure)



		! Calculate height and momentum in centred grid nodes
		Do j = jstart,jend
			do i = istart,iend
				eta_new(i,j) = eta(i,j) - dt_dx*(ueta_mid_x(i,j) - ueta_mid_x(i-1,j)) - dt_dy*(veta_mid_y(i,j) - veta_mid_y(i,j-1)) ! New Height
				ueta_new(i,j) = ueta(i,j) - dt_dx*(Ux_mid(i,j) - Ux_mid(i-1,j)) - dt_dy*(Uy_mid(i,j) - Uy_mid(i,j-1)) 				! New x-dir momentum
				veta_new(i,j)= veta(i,j) - dt_dx*(Vx_mid(i,j) - Vx_mid(i-1,j)) - dt_dy*(Vy_mid(i,j) - Vy_mid(i,j-1)) 				! New y-dir momentum
			end do
		end do

		! ======================== Lax-Wendroff solver ends =========================
		
	
	
		u_new = ueta_new/eta_new		! x-dir velocity (divide momentum by height (height analogous to mass)
		v_new = veta_new/eta_new		! y-dir velocity (divide momentum by height (height analogous to mass)

		! Update values 
		eta = eta_new
		u   = u_new
		v   = v_new

			

		! ---------------- Exchange ghost cells for eta, u, v -----------------
		!CALL MPI_WAITALL(8, req, stat, ierr)
		CALL MPI_Barrier(CART_COMM, ierr)

		!CALL exchange_ghost_points(eta, u, v)
		CALL variable_communication(eta)
		CALL variable_communication(u)
		CALL variable_communication(v)

		if (mod(iteration, steps) == 0) then
			write(filename, '(A,I5.5,A)') 'SWE_eta_', iteration, '.bin'
			CALL write_var_parallel_io(filename, eta(istart:iEnd,jstart:jend))
		end if
		
		iteration = iteration + 1
	
	END DO  ! end iterations	
	! ================================ END OF ITERATION LOOP ==================================
	! =========================================================================================
	
	
	! Compute local min and max for the height
	local_min = minval(eta(iStart:iEnd, jStart:jEnd))
	local_max = maxval(eta(iStart:iEnd, jStart:jEnd))

	! Reduce across all MPI processes
	call MPI_Allreduce(local_min, global_min, 1, MPI_REAL8, MPI_MIN, CART_COMM, ierr)
	call MPI_Allreduce(local_max, global_max, 1, MPI_REAL8, MPI_MAX, CART_COMM, ierr)

	! Write out global min/max from rank 0
	if (coord(1) == 0 .and. coord(2) == 0) then
		write(*,*) 'Global Min eta = ', global_min
		write(*,*) 'Global Max eta = ', global_max
	end if
	
	

	! Finalize MPI
	CALL MPI_Type_free(ROW_TYPE, ierr)
	CALL MPI_Finalize(ierr)
	

END PROGRAM mpi_main

