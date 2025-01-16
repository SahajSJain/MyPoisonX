program mypoisonx
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
#ifdef MPI
  use mpi
#endif
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER IERROR
  call MPI_INIT(ierror)
  PROGRAM_TIME_START=MPI_WTIME();
  call startup
  call MPISTARTUP
  CALL MPI_Barrier(CommCart,ierror)
  CALL MEMORY_ALLOCATE

  CALL initializeMyPoisonX
  CALL MPI_Barrier(CommCart,ierror)
  IF(SOLVEMODE>0) CALL write_dump
  CALL MPI_Barrier(CommCart,ierror)
  IF(ImtheBOSS) WRITE(*,*) "MyPoisonX DIED BECAUSE SIMULATION IS OVER"
  CALL MPIKILL()
end program mypoisonx
