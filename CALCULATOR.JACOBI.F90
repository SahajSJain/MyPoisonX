SUBROUTINE JACOBI(kiter1, error1 )
  USE INPUTDATA
  USE COEFFDATA
  USE MPIDATA
  USE FIELDDATA
  USE FLOWDATA

  interface
    SUBROUTINE BOUNDARY_CONDITIONS(var)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var
    end subroutine BOUNDARY_CONDITIONS
    SUBROUTINE MPI_SETUP_DATATYPES( Nxc1, Nyc1, Nzc1, NGL1)
      USE INPUTDATA
      INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1, NGL1
    end subroutine MPI_SETUP_DATATYPES
    SUBROUTINE MPI_Communicate_Block_data(var1, TIMEOUT1,Nxc1, Nyc1, Nzc1)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var1
      REAL(KIND=CGR), INTENT(out) :: TIMEOUT1
      INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1
    end subroutine MPI_Communicate_Block_data
    subroutine MPI_SUMMATION(QUANT)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(inout) :: QUANT
    end subroutine MPI_SUMMATION
    subroutine MPI_MAXCALC(QUANT)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(inout) :: QUANT
    end subroutine MPI_MAXCALC
  end interface

  INTEGER, INTENT(OUT) :: kiter1
  REAL(KIND=CGR), INTENT(OUT) :: error1
  INTEGER :: i, j, k

  REAL(KIND=CGR) :: TSTART1, TSTART2, TSTART3
  REAL(KIND=CGR) :: TEND1, TEND2, TEND3

  REAL(KIND=CGR) :: ITERTIME,innertime, TOTALINNERTIME, TOTALTIME, RESTIME, TOTALRESTIME
  REAL(KIND=CGR) :: COMMTIME, TOTALCOMMTIME, ALLREDTIME, TOTALALLREDTIME, INNERCOMMTIME
  REAL(KIND=CGR) :: LOCAL_RES, GLOBAL_RES
  INTEGER :: INITER
  INTEGER :: IERROR

  INTEGER :: ierr
  REAL(KIND=CGR) :: rhoe, rhow, rhon, rhos, rhob, rhof
  REAL(KIND=CGR) :: SCALEFACTOR, TIMEOUTx
  iter=0;
  TOTALALLREDTIME=0;
  TOTALTIME=0;
  TOTALRESTIME=0;
  TOTALCOMMTIME=0;
  INNERCOMMTIME=0;
  kITER1=0;
  GLOBAL_RES=99999999;
  SCALEFACTOR=(Nxc*Nyc*Nzc*COMM_SZ);
  IF(ImtheBOSS)  WRITE (*,*) "------------------------------------------------------------------------"
  IF(ImtheBOSS)  WRITE (*,*) " iter     local res     global res    time"
  IF(ImtheBOSS)  WRITE (*,*) "------------------------------------------------------------------------"

  do while (kITER1 .lt. MAXITER .and. GLOBAL_RES .gt. TOLERANCE)
    TSTART1 = MPI_WTIME()
    INNERCOMMTIME=0.0;
    do INITER=1, INNERITER
      ! IF(ImtheBOSS) WRITE(*,*) "INNER ITER NUMBER" , INITER
      ! CALL BOUNDARY_CONDITIONS(phi)
      do k=1-NGL, Nzc+Ngl
        do j=1-NGL, Nyc+Ngl
          do i=1-NGL, Nxc+Ngl
            phi0(i,j,k) = phi(i,j,k)
          enddo !i
        ENDDO !j
      ENDDO !k
      DO k=1,Nzc
        DO j=1, Nyc
          DO i=1, Nxc
            phi(i,j,k) = (source(i,j,k) - &
              (Aw(i,j,k)*phi0(i-1,j,k) + Ae(i,j,k)*phi0(i+1,j,k) + &
              As(i,j,k)*phi0(i,j-1,k) + An(i,j,k)*phi0(i,j+1,k) + &
              Ab(i,j,k)*phi0(i,j,k-1) + Af(i,j,k)*phi0(i,j,k+1))) * Apinv(i,j,k)
          enddo !i
        ENDDO !j
      ENDDO !k
      CALL MPI_SETUP_DATATYPES( Nxc, Nyc, Nzc, NGL)
      CALL MPI_Communicate_Block_data(phi, TIMEOUTx,Nxc, Nyc, Nzc)
      CALL MPI_FREE
      TOTALCOMMTIME=TOTALCOMMTIME+COMMTIME;
      INNERCOMMTIME=INNERCOMMTIME+COMMTIME;
    ENDDO
    kITER1=kITER1+INNERITER
    !! CALCULATE RESIDUAL
    residual=0.0;
    LOCAL_RES=0.0;
    GLOBAL_RES=0.0;
    DO k=1,Nzc
      DO j=1, Nyc
        DO i=1, Nxc
          residual(i,j,k)=Aw(i,j,k)*phi(i-1,j,k)+Ae(i,j,k)*phi(i+1,j,k)&
            +As(i,j,k)*phi(i,j-1,k)+An(i,j,k)*phi(i,j+1,k)&
            +Ab(i,j,k)*phi(i,j,k-1)+Af(i,j,k)*phi(i,j,k+1)&
            +phi(i,j,k)*Ap(i,j,k)-source(i,j,k)
        enddo !i
      ENDDO !j
    ENDDO !k
    LOCAL_RES=MAXVAL(abs(residual(1:NXC,1:NYC,1:NZC)))
    GLOBAL_RES=LOCAL_RES
    CALL MPI_MAXCALC(GLOBAL_RES)
    TEND1=MPI_WTIME()
    IF(myCartRank .eq. 0)  WRITE(*,*) kiter1, LOCAL_RES,GLOBAL_RES, TEND1-TSTART1
    CALL MPI_Barrier(CommCart,ierror)
  end do ! while loop
  error1=GLOBAL_RES;
  CALL MPI_Barrier(CommCart,ierror)
  IF(ImtheBOSS)   WRITE (*,*) "------------------------------------------------------------------------"
  WRITE(*,*) "I am MPI Process" , myCartRank," out of ", COMM_SZ,", and my local residual is = ", LOCAL_RES
end subroutine JACOBI
