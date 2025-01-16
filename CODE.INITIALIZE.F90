SUBROUTINE initializeMyPoisonX

  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  INTEGER::ierror, i, j, k, kiter
  REAL(KIND=CGR) :: TIMEOUT, jacobierror, phi_exact, xp, yp, zp
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
    SUBROUTINE JACOBI(kiter1, error1 )
      USE INPUTDATA
      INTEGER, INTENT(OUT) :: kiter1
      REAL(KIND=CGR), INTENT(OUT) :: error1
    end subroutine JACOBI
    SUBROUTINE GSSOR(kiter1, error1 )
      USE INPUTDATA
      INTEGER, INTENT(OUT) :: kiter1
      REAL(KIND=CGR), INTENT(OUT) :: error1
    end subroutine GSSOR
    function phicalc(xa, ya, za) result(output)
      USE INPUTDATA
      implicit none
      real(kind=CGR), intent(in) :: xa, ya, za
      real(kind=CGR) :: output
    end function phicalc

  end interface

  CALL setup_grid
  CALL SETUP_COEFFS
  CALL SETUP_SOURCE
  CALL SETUP_BOUNDARYCONDITIONS
  SELECT CASE (PINITFLAG)
   CASE (1) !! initialize as given value
    phi=PINITVAL
   CASE (2) !! initialize as random value
    call random_seed()
    call random_number(phi)
  END SELECT
  CALL BOUNDARY_CONDITIONS(phi)
  phi0=phi;
  phi00=phi;

  CALL MPI_SETUP_DATATYPES( Nxc, Nyc, Nzc, NGL)
  CALL MPI_COMM_SEND_RECV(phi, TIMEOUT,Nxc, Nyc, Nzc)
  ! CALL MPI_FREE

  IF(ImtheBOSS) WRITE(*,*) "INITIALIZATION COMPLETE"
!   IF(ImtheBOSS) WRITE(*,*) phi
  IF(ImtheBOSS) WRITE(*,*) "START JACOBI"
  CALL MPI_Barrier(CommCart,ierror)
  CALL SETUP_COEFFS

  SELECT CASE(SOLVEMODE)
   CASE(0)
    CALL MPITEST(kiter, jacobierror)
   CASE(1)
    CALL JACOBI(kiter, jacobierror );
   CASE(2)
    CALL GSSOR(kiter, jacobierror );
  END SELECT
  IF(SOURCEFLAG .eq. 3) THEN
    do k=1-NGL, Nzc+Ngl
      do j=1-NGL, Nyc+Ngl
        do i=1-NGL, Nxc+Ngl
          xp=xc(i)
          yp=yc(j)
          zp=zc(k)
          phi_exact=phicalc(xp,yp,zp);
          errorphi(i,j,k) = phi_exact - phi(i,j,k)
        enddo !i
      ENDDO !j
    ENDDO !k

  END IF
end subroutine initializeMyPoisonX
