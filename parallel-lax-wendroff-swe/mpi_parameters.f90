MODULE mpi_parameters
    IMPLICIT NONE

	! ------------- Declare parameters -------------
	INTEGER, PARAMETER 			::	nx = 100					! Number of grid points x
	INTEGER, PARAMETER 			::	ny = 100					! Number of grid points y
	INTEGER, PARAMETER			::	step = 100					! TARGET number of export frames, exports every floor(tmax/dt/steps)
	
	integer:: nLocX
	integer:: nLocY
	
	REAL(KIND = 8), PARAMETER 	::	Lx = 50						! Length of domain x (m)
	REAL(KIND = 8), PARAMETER 	::	Ly = 50						! Length of domain y (m)
	REAL(KIND = 8), PARAMETER	::	PI = 4.*ATAN(1d0)			! PI	
	REAL(KIND = 8), PARAMETER	::	eta0 = 1.0					! Mean depth (m)
	REAL(KIND = 8), PARAMETER	::	amp = 0.5 					! Initial bump amplitude (m)
	REAL(KIND = 8), PARAMETER	::	sigma = 0.1					! Relative width (m)
	REAL(KIND = 8), PARAMETER	::	x0 = 0.5					! Midpoint of the domain normalized x
	REAL(KIND = 8), PARAMETER	::	y0 = 0.5 					! Midpoint of the domain normalized y
	REAL(KIND = 8), PARAMETER	::	g = 9.81					! Gravitational Acceleration (ms^-2)
	REAL(KIND = 8), PARAMETER	::	CFL = 0.01					! CFL
	REAL(KIND = 8), PARAMETER	::	eps_eta = 1E-10				! Division Safeguard
	REAL(KIND = 8), PARAMETER	::	tmax = 50.					! Maximum time 
	REAL(KIND = 8), PARAMETER	::	dt = 0.02					! time step
	

	! Declare variables
	REAL(KIND = 8)	::	t = 0.0	! Initial time
	
	! ------------- Declare allocatable -------------
	INTEGER 					:: i,j,iStart, iEnd, jStart, jEnd,iGW, iGE, jGS, jGN ! Indeces
	INTEGER			 			:: pointx, pointy									 ! Specific grid points
	REAL(KIND = 8) 				:: dx, dy, dt_dx, dt_dy								 ! Grid spacing and time-stepping
	REAL(KIND = 8), ALLOCATABLE :: x(:), y(:)										 ! Coordinate vectors for grid positions
	
	REAL(KIND = 8), ALLOCATABLE :: eta(:,:)		, u(:,:)	, v(:,:)					! Height elevation and velocity fields
	REAL(KIND = 8), ALLOCATABLE :: eta_new(:,:)	, u_new(:,:), v_new(:,:), veta_new(:,:) ! Updated Height elevation and velocity fields
	
	! Allocatable arrays for fluxes and staggered grid variables
	REAL(KIND = 8), ALLOCATABLE :: ueta(:,:), Ux(:,:), ueta_mid_x(:,:), veta_mid_x(:,:)
	REAL(KIND = 8), ALLOCATABLE :: Ux_mid(:,:), Vx_mid(:,:), ueta_new(:,:), Uy(:,:), veta(:,:), Vx(:,:), Vy(:,:)
	REAL(KIND = 8), ALLOCATABLE :: eta_mid_x(:,:), eta_mid_y(:,:)
	REAL(KIND = 8), ALLOCATABLE :: ueta_mid_y(:,:), veta_mid_y(:,:)
	REAL(KIND = 8), ALLOCATABLE :: Uy_mid(:,:), Vy_mid(:,:)
END MODULE mpi_parameters
