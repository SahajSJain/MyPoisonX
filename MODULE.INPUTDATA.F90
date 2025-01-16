MODULE INPUTDATA
  use iso_fortran_env, only: real32,real64,real128
  IMPLICIT NONE
  INTEGER, PARAMETER :: CGR = real64
  INTEGER, PARAMETER :: ifINPUTDAT =69
  INTEGER :: Npx, Npy, Npz !! No. of procs
  INTEGER :: NGL !! No. of Ghost layers
  INTEGER :: Nx_glbl, Ny_glbl, Nz_glbl !! No. of face centers global
  INTEGER :: Nx, Ny, Nz !! No. of grid points
  INTEGER :: Nxc_glbl, Nyc_glbl, Nzc_glbl !! No. of cell centers global
  INTEGER :: Nxc, Nyc, Nzc !! No. of cell centers
  INTEGER :: Nx_f, Ny_f, Nz_f !!No. of cell faces
  INTEGER :: x_unif, y_unif, z_unif
  REAL(KIND=CGR) :: x_start, x_end, y_start, y_end, z_start, z_end
  INTEGER :: iter
  INTEGER :: DEBUGGER, DEBUGGEROUT, EXPORTFILES
  INTEGER :: i_s, i_e, j_s, j_e, k_s, k_e
  INTEGER :: SOLVEMODE, SRJTYPE, NORMTYPE
  INTEGER, PARAMETER :: JACOBIMODE = 1, SORMODE =2, SRJMODE =3
  LOGICAL :: JACOBION, SORON, SRJON
  REAL(KIND=CGR) :: RELAXFAC, TOLERANCE
  INTEGER :: PINITFLAG, SOURCEFLAG,COMM_MODE_FLAG
  REAL(KIND=CGR) :: PINITVAL
  INTEGER :: MAXITER, INNERITER, SYNCITER
  INTEGER :: BC_EAST, BC_WEST, BC_NORTH, BC_SOUTH , BC_FRONT ,   BC_BACK
  INTEGER :: DIRBC=1, NEUMBC=2, CUSTOMBC=3
  REAL(KIND=CGR) :: TIME_MAIN_S, TIME_SOLVE_S, TIME_COMM_S
  REAL(KIND=CGR) :: TIME_MAIN_E, TIME_SOLVE_E, TIME_COMM_E
  REAL(KIND=CGR), PARAMETER ::   PIE = 4.0_CGR*ATAN(1.0_CGR)
  REAL ::   VAL_EAST,    VAL_WEST ,   VAL_NORTH  , VAL_SOUTH  , VAL_FRONT  , VAL_BACK
  REAL :: BCT_EAST, BCT_WEST, BCT_NORTH, BCT_SOUTH , BCT_FRONT ,   BCT_BACK
end module INPUTDATA
