SUBROUTINE MPITEST(kiter1, error1 )
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
    SUBROUTINE MPI_COMM_SEND_RECV(var1, TIMEOUT1,Nxc1, Nyc1, Nzc1)
      USE INPUTDATA
      REAL(KIND=CGR), INTENT(INOUT),DIMENSION(:,:,:) :: var1
      REAL(KIND=CGR), INTENT(out) :: TIMEOUT1
      INTEGER, INTENT(in) :: Nxc1, Nyc1, Nzc1
    end subroutine MPI_COMM_SEND_RECV
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
  CHARACTER*20 :: fname1, fname2

  do k=1-NGL, Nzc+Ngl
    do j=1-NGL, Nyc+Ngl
      do i=1-NGL, Nxc+Ngl
        phi_ew(i,j,k)=REAL(i);
        phi_ns(i,j,k)=REAL(j);
        phi_fb(i,j,k)=REAL(k);
      enddo
    enddo
  enddo
  uphi_ew=phi_ew;
  uphi_ns=phi_ns;
  uphi_fb=phi_fb;
  CALL MPI_COMM_SEND_RECV(phi_ew, TIMEOUTx,Nxc, Nyc, Nzc)
  CALL MPI_COMM_SEND_RECV(phi_ns, TIMEOUTx,Nxc, Nyc, Nzc)
  CALL MPI_COMM_SEND_RECV(phi_fb, TIMEOUTx,Nxc, Nyc, Nzc)

  WRITE(fname1,"('q.',I3.3,'.',I7.7,'.dat')") MY_RANK,0
  OPEN(71,FILE=fname1)
  WRITE(71,*) " Title = MyPoisonX "
  ! WRITE(71,*) 'VARIABLES = "X","Y", "Z","I", "J", "K", "phi","source","res","error","Aw","Ae","As","An","Ab","Af","Ap"'
  WRITE(71,*) 'VARIABLES = "X","Y", "Z","pEW","pNS","pFB","upEW","upNS","upFB"'
  write(71,*) 'ZONE T="ZONE ',MY_RANK,'", I=',nxc+2*ngl,', J=',nyc+2*ngl,', k =', nzc+2*NGL
  ! WRITE(71,*) 'ZONE  I=',nxc,', J=',nyc,', k =', nzc

  DO k=1-Ngl,nzc+Ngl
    DO j=1-Ngl,nyc+Ngl
      DO i=1-Ngl,nxc+Ngl
        WRITE(71, '(3F12.5, 6F12.5)') xc(i), yc(j), zc(k), &
          phi_ew(i,j,k), phi_ns(i,j,k), phi_fb(i,j,k), &
          uphi_ew(i,j,k), uphi_ns(i,j,k), uphi_fb(i,j,k)
      END DO
    END DO
  END DO

end subroutine MPITEST
