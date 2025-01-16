subroutine startup
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
!! open files
  OPEN(ifINPUTDAT, File="input.dat", ACTION='READ', STATUS='OLD')
  READ(ifINPUTDAT,*)!======================== MyPoisonX ==============================
  READ(ifINPUTDAT,*)!Npx, Npy, Npz
  READ(ifINPUTDAT,*) Npx, Npy, Npz
  READ(ifINPUTDAT,*)!Ngl
  READ(ifINPUTDAT,*) Ngl
  READ(ifINPUTDAT,*)!====================== Grid Data =========================
  READ(ifINPUTDAT,*)!Nx_GLBL     Ny_GLBL     Nz_GLBL
  READ(ifINPUTDAT,*) Nx_GLBL,     Ny_GLBL,     Nz_GLBL
  READ(ifINPUTDAT,*)!x_unif      x_start     x_end
  READ(ifINPUTDAT,*) x_unif  ,    x_start  ,   x_end
  READ(ifINPUTDAT,*)!y_unif      y_start     y_end
  READ(ifINPUTDAT,*) y_unif  ,    y_start  ,   y_end
  READ(ifINPUTDAT,*)!z_unif      z_start     z_end
  READ(ifINPUTDAT,*) z_unif  ,    z_start  ,   z_end
  READ(ifINPUTDAT,*)!====================SOLVER DETAILS=======================
  READ(ifINPUTDAT,*)!DEBUGGER    DEBUGGEROUT     EXPORTFILES
  READ(ifINPUTDAT,*) DEBUGGER   , DEBUGGEROUT   ,  EXPORTFILES
  READ(ifINPUTDAT,*)!SOLVEMODE   RELAXFAC    SRJTYPE    COMM_MODE_FLAG
  READ(ifINPUTDAT,*) SOLVEMODE ,  RELAXFAC  ,  SRJTYPE, COMM_MODE_FLAG
  READ(ifINPUTDAT,*)!PINITFLAG   PINITVAL    SOURCEFLAG
  READ(ifINPUTDAT,*) PINITFLAG  , PINITVAL  ,  SOURCEFLAG
  READ(ifINPUTDAT,*) !MAXITER     INNERITER   TOLERANCE   NORMTYPE
  READ(ifINPUTDAT,*) MAXITER   , INNERITER  , TOLERANCE ,  NORMTYPE
  READ(ifINPUTDAT,*)  !BC_EAST     BC_WEST     BC_NORTH    BC_SOUTH    BC_FRONT    BC_BACK
  READ(ifINPUTDAT,*) BC_EAST   ,  BC_WEST  ,   BC_NORTH ,   BC_SOUTH   , BC_FRONT ,   BC_BACK
  READ(ifINPUTDAT,*)  !VAL_EAST    VAL_WEST    VAL_NORTH   VAL_SOUTH   VAL_FRONT   VAL_BACK
  READ(ifINPUTDAT,*) VAL_EAST,    VAL_WEST ,   VAL_NORTH  , VAL_SOUTH  , VAL_FRONT  , VAL_BACK

  Nx=(Nx_glbl-1)/Npx+1;
  Ny=(Ny_glbl-1)/Npy+1;
  Nz=(Nz_glbl-1)/Npz+1;
  Nxc=Nx-1;
  Nyc=Ny-1;
  Nzc=Nz-1;

! WRITE(*,*) "Task ", MY_RANK," out of ", COMM_SZ," tasks sees ", Npx* Npy* Npz

end subroutine startup
