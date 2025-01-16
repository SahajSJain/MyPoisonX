subroutine MPISTARTUP
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER IERROR
  logical, parameter ::  reorder=.true.
  call MPI_COMM_SIZE(MPI_COMM_WORLD, COMM_SZ, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, ierror)
  PROGRAM_TIME_START=MPI_WTIME();
  IF(MY_RANK .eq. 0) THEN
    ImtheBOSS=.true.
  ELSE
    ImtheBOSS=.false.
  ENDIF
  Nptot=Npx*Npy*Npz
  ! print *, 'Hello World from process: ', MY_RANK, 'of ', COMM_SZ
  IF(Nptot .ne. COMM_SZ) THEN
    IF(ImtheBOSS) WRITE(*,*) "MyPoisonX DIED BECAUSE Nptot<COMM_SZ"
    CALL MPIKILL()
  ENDIF
  IF(BC_EAST .eq. 4) THEN
    IF(BC_WEST .eq. 4) THEN
      PeriodicArr(1)=1;
    ELSE
      IF(ImtheBOSS) WRITE(*,*) "MyPoisonX DIED BECAUSE INVALID BCS IN EAST_WEST"
      CALL MPIKILL()
    ENDIF
  ELSE
    PeriodicArr(1)=0;
  ENDIF
  IF(BC_NORTH .eq. 4) THEN
    IF(BC_SOUTH .eq. 4) THEN
      PeriodicArr(2)=1;
    ELSE
      IF(ImtheBOSS) WRITE(*,*) "MyPoisonX DIED BECAUSE INVALID BCS IN NORTH_SOUTH"
      CALL MPIKILL()
    ENDIF
  ELSE
    PeriodicArr(2)=0;
  ENDIF
  IF(BC_FRONT .eq. 4) THEN
    IF(BC_BACK .eq. 4) THEN
      PeriodicArr(3)=1;
    ELSE
      IF(ImtheBOSS) WRITE(*,*) "MyPoisonX DIED BECAUSE INVALID BCS IN FRONT_BACK"
      CALL MPIKILL()
    ENDIF
  ELSE
    PeriodicArr(3)=0;
  ENDIF
  NumProcArr=(/Npx, Npy, Npz /);
  !! Do Cartesian decomposition
  call MPI_Cart_create(MPI_COMM_WORLD, 3, NumProcArr,PeriodicArr, reorder, CommCart, ierror)
  call  MPI_Comm_rank(CommCart,myCartRank ,ierror);
  call   MPI_Cart_coords(CommCart, myCartRank, 3, myCOORDS);
  WRITE(*,*) "I am MPI Process" , myCartRank," out of ", COMM_SZ,", and I am located at",myCOORDS
  !!! GET NEIGHBOURS
  CALL MPI_CART_SHIFT(CommCart, 0, 1, myWest, myEast, ierror)
  ! CALL MPI_CART_SHIFT(CommCart, 0, 1, myCartRank, myEast, ierror)
  ! CALL MPI_CART_SHIFT(CommCart, 1, -1, myCartRank, mySouth, ierror)
  CALL MPI_CART_SHIFT(CommCart, 1, 1, mySouth, myNorth, ierror)
  ! CALL MPI_CART_SHIFT(CommCart, 2, -1, myCartRank, myBack, ierror)
  CALL MPI_CART_SHIFT(CommCart, 2, 1, myBack, myFront, ierror)

  CALL MPI_Barrier(CommCart,ierror)
  WRITE(*,*) "I am MPI Process" , myCartRank," out of ", COMM_SZ,", and my neighbours are"&
    ,myWest, myEast, mySouth, myNorth, myBack,myFront
  CALL MPI_Barrier(CommCart,ierror)
  myIs=myCOORDS(1);
  myJs=myCOORDS(2);
  myKs=myCOORDS(3);
end subroutine MPISTARTUP

subroutine MPIKILL()
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER IERROR
  PROGRAM_TIME_END=MPI_WTIME();
  CALL MPI_Barrier(CommCart,ierror)
  CALL MPI_FINALIZE(IERROR)
  STOP
end subroutine MPIKILL

SUBROUTINE MPI_SETUP_DATATYPES( Nxc1, Nyc1, Nzc1, NGL1)
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1, NGL1

  INTEGER :: ierr
  integer, dimension(3) :: array_sizes
  integer, dimension(3) :: array_subsizes
  integer, dimension(3) :: array_starts
  array_sizes = (/ Nxc1 + 2 * NGL1, Nyc1 + 2 * NGL1, Nzc1 + 2 * NGL1 /)
  array_subsizes = (/ NGL1, Nyc1 + 2 * NGL1, Nzc1 + 2 * NGL1 /)
  ! Define the start of the subarray within the full array
  ! Note: Fortran arrays are 1-based, MPI expects 0-based
  array_starts = (/ 0, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_W_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_W_recv, ierr)
  array_starts = (/ NGL1, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_W_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_W_send, ierr)

  array_starts = (/ Nxc1, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_E_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_E_send, ierr)
  array_starts = (/ Nxc1+NGL1, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_E_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_E_recv, ierr)

!! north south arrays

  array_subsizes = (/ Nxc1 + 2 * NGL1, NGL1, Nzc1 + 2 * NGL1 /)
  ! Define the start of the subarray within the full array
  ! Note: Fortran arrays are 1-based, MPI expects 0-based
  array_starts = (/ 0, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_S_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_S_recv, ierr)
  array_starts = (/ 0, NGL1, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_S_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_S_send, ierr)

  array_starts = (/ 0, Nyc1, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_N_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_N_send, ierr)
  array_starts = (/ 0, Nyc1+NGL1, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_N_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_N_recv, ierr)

!! back front arrays

  array_subsizes = (/ Nxc1 + 2 * NGL1, Nyc1 + 2 * NGL1,NGL1 /)
  ! Define the start of the subarray within the full array
  ! Note: Fortran arrays are 1-based, MPI expects 0-based
  array_starts = (/ 0, 0, 0 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_B_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_B_recv, ierr)
  array_starts = (/ 0, 0, NGL1 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_B_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_B_send, ierr)

  array_starts = (/ 0, 0, Nzc1 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_F_send, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_F_send, ierr)
  array_starts = (/ 0, 0, Nzc1+NGL1 /)
  call MPI_Type_create_subarray(3, array_sizes, array_subsizes, array_starts, &
    MPI_ORDER_FORTRAN, MPI_DOUBLE, type_F_recv, ierr) !! SEND TO WEST TYPE
  call MPI_Type_commit(type_F_recv, ierr)

end subroutine MPI_SETUP_DATATYPES

subroutine MPI_Communicate_Block_data(var, TIMEOUT,Nxc1, Nyc1, Nzc1)
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var
  REAL(KIND=CGR), INTENT(out) :: TIMEOUT
  INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1
  interface
    SUBROUTINE MPI_COMM_SEND_RECV(var1, TIMEOUT1,Nxc2, Nyc2, Nzc2)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var1
      REAL(KIND=CGR), INTENT(out) :: TIMEOUT1
      INTEGER, INTENT(in) :: Nxc2, Nyc2, Nzc2
    end subroutine MPI_COMM_SEND_RECV
    SUBROUTINE MPI_COMM_NB(var1, TIMEOUT1,Nxc2, Nyc2, Nzc2)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var1
      REAL(KIND=CGR), INTENT(out) :: TIMEOUT1
      INTEGER, INTENT(in) :: Nxc2, Nyc2, Nzc2
    end subroutine MPI_COMM_NB
  end interface
  SELECT CASE (COMM_MODE_FLAG)
   CASE (1) !! sendrecv
    CALL MPI_COMM_SEND_RECV(var,TIMEOUT,Nxc1, Nyc1, Nzc1)
   CASE (2) !! non blocking
    CALL MPI_COMM_NB(var,TIMEOUT,Nxc1, Nyc1, Nzc1)
   CASE DEFAULT
    CALL MPI_COMM_SEND_RECV(var,TIMEOUT,Nxc1, Nyc1, Nzc1)
  END SELECT
end subroutine MPI_Communicate_Block_data

subroutine MPI_COMM_SEND_RECV(var, TIMEOUT,Nxc1, Nyc1, Nzc1)
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var
  REAL(KIND=CGR), INTENT(out) :: TIMEOUT
  INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1
  INTEGER :: status(MPI_STATUS_SIZE)
  REAL(KIND=CGR) :: startTime, endTime
  INTEGER :: ierror
  TIMEOUT=0.0;
  startTime = MPI_WTIME()

! SEND TO WEST AND RECEIVE FROM EAST
  CALL MPI_SENDRECV(var,1,type_W_send,myWest,tagWs,  &
    var,1,type_E_recv,myEast,tagEr,  &
    CommCart,status,ierror)

!   Opposite east west
!   ------------
  CALL MPI_SENDRECV(var,1,type_E_send,myEast,tagEs,  &
    var,1,type_W_recv,myWest,tagWr,  &
    CommCart,status,ierror)

! SEND TO South AND RECEIVE FROM North
  CALL MPI_SENDRECV(var,1,type_S_send,mySouth,tagSs,  &
    var,1,type_N_recv,myNorth,tagNr,  &
    CommCart,status,ierror)

!   Opposite north south
!   ------------
  CALL MPI_SENDRECV(var,1,type_N_send,myNorth,tagNs,  &
    var,1,type_S_recv,mySouth,tagSr,  &
    CommCart,status,ierror)

! SEND TO Back AND RECEIVE FROM Front
  CALL MPI_SENDRECV(var,1,type_B_send,myBack,tagBs,  &
    var,1,type_F_recv,myFront,tagFr,  &
    CommCart,status,ierror)

!   Opposite front back
!   ------------
  CALL MPI_SENDRECV(var,1,type_F_send,myFront,tagFs,  &
    var,1,type_B_recv,myBack,tagBr,  &
    CommCart,status,ierror)

  endTime = MPI_WTIME()
  TIMEOUT=endTime-startTime
end subroutine MPI_COMM_SEND_RECV

subroutine MPI_COMM_NB(var, TIMEOUT,Nxc1, Nyc1, Nzc1)
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var
  REAL(KIND=CGR), INTENT(out) :: TIMEOUT
  INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1
  INTEGER :: status(MPI_STATUS_SIZE)
  REAL(KIND=CGR) :: startTime, endTime
  INTEGER :: ierror
  TIMEOUT=0.0;
  startTime = MPI_WTIME()

! SEND TO WEST AND RECEIVE FROM EAST
  CALL MPI_SENDRECV(var,1,type_W_send,myWest,tagWs,  &
    var,1,type_E_recv,myEast,tagEr,  &
    CommCart,status,ierror)
!   Opposite east west
!   ------------
  CALL MPI_SENDRECV(var,1,type_E_send,myEast,tagEs,  &
    var,1,type_W_recv,myWest,tagWr,  &
    CommCart,status,ierror)
! SEND TO South AND RECEIVE FROM North
  CALL MPI_SENDRECV(var,1,type_S_send,mySouth,tagSs,  &
    var,1,type_N_recv,myNorth,tagNr,  &
    CommCart,status,ierror)
!   Opposite north south
!   ------------
  CALL MPI_SENDRECV(var,1,type_N_send,myNorth,tagNs,  &
    var,1,type_S_recv,mySouth,tagSr,  &
    CommCart,status,ierror)
! SEND TO Back AND RECEIVE FROM Front
  CALL MPI_SENDRECV(var,1,type_B_send,myBack,tagBs,  &
    var,1,type_F_recv,myFront,tagFr,  &
    CommCart,status,ierror)
!   Opposite front back
!   ------------
  CALL MPI_SENDRECV(var,1,type_F_send,myFront,tagFs,  &
    var,1,type_B_recv,myBack,tagBr,  &
    CommCart,status,ierror)

  endTime = MPI_WTIME()
  TIMEOUT=endTime-startTime;
end subroutine MPI_COMM_NB

subroutine MPI_FREE
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER :: IERROR
  CALL MPI_TYPE_FREE(type_W_recv,IERROR)
  CALL MPI_TYPE_FREE(type_W_send,IERROR)
  CALL MPI_TYPE_FREE(type_E_recv,IERROR)
  CALL MPI_TYPE_FREE(type_E_send,IERROR)
  CALL MPI_TYPE_FREE(type_S_recv,IERROR)
  CALL MPI_TYPE_FREE(type_S_send,IERROR)
  CALL MPI_TYPE_FREE(type_N_recv,IERROR)
  CALL MPI_TYPE_FREE(type_N_send,IERROR)
  CALL MPI_TYPE_FREE(type_B_recv,IERROR)
  CALL MPI_TYPE_FREE(type_B_send,IERROR)
  CALL MPI_TYPE_FREE(type_F_recv,IERROR)
  CALL MPI_TYPE_FREE(type_F_send,IERROR)

end subroutine MPI_FREE

subroutine MPI_SUMMATION(QUANT)
  ! USE INPUTDATA
  ! USE FIELDDATA
  ! USE COEFFDATA
  ! use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  REAL(KIND=CGR), INTENT(inout) :: QUANT
  REAL(KIND=CGR):: SUMOFQUANT
  INTEGER :: IERROR
  CALL MPI_ALLREDUCE(QUANT, SUMOFQUANT, 1, MPI_DOUBLE, MPI_SUM, CommCart, IERROR)
  QUANT=SUMOFQUANT
end subroutine MPI_SUMMATION

subroutine MPI_MAXCALC(QUANT)
  ! USE INPUTDATA
  ! USE FIELDDATA
  ! USE COEFFDATA
  ! use FLOWDATA
  USE mpiDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  REAL(KIND=CGR), INTENT(inout) :: QUANT
  REAL(KIND=CGR):: SUMOFQUANT
  INTEGER :: IERROR
  CALL MPI_ALLREDUCE(QUANT, SUMOFQUANT, 1, MPI_DOUBLE, MPI_MAX, CommCart, IERROR)
  QUANT=SUMOFQUANT
end subroutine MPI_MAXCALC



