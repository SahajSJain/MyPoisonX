subroutine write_dump
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  USE mpi
  IMPLICIT NONE
  ! include 'mpif.h'
  CHARACTER*20 :: fname1, fname2
  INTEGER:: i,j,k
  WRITE(fname1,"('q.',I3.3,'.',I7.7,'.dat')") MY_RANK,0
  OPEN(71,FILE=fname1)
  WRITE(71,*) "Title = MyPoisonX "
  ! WRITE(71,*) 'VARIABLES = "X","Y", "Z","I", "J", "K", "phi","source","res","error","Aw","Ae","As","An","Ab","Af","Ap"'
  WRITE(71,*) 'VARIABLES = "X","Y", "Z","phi","source","res","error"'

  write(71,'(A,I5,A,I5,A,I5,A,I5)') 'ZONE T="ZONE ',MY_RANK,'", I=',nxc,', J=',nyc,', K=',nzc
  ! WRITE(71,*) 'ZONE  I=',nxc,', J=',nyc,', k =', nzc

  DO k=1,nzc
    DO j=1,nyc
      DO i=1,nxc
        WRITE(71,'(F12.5,1X,F12.5,1X,F12.5,1X,F18.9,1X,F18.9,1X,F18.9,1X,F18.9)') &
        xc(i), yc(j), zc(k), phi(i,j,k), source(i,j,k), residual(i,j,k), errorphi(i,j,k)    
      END DO
    END DO
  END DO
  CLOSE(71)

  WRITE(fname2,"('C.',I3.3,'.',I7.7,'.dat')") MY_RANK,0
  OPEN(72,FILE=fname2)
  WRITE(72,*) "Title = MyPoisonX "
  ! WRITE(71,*) 'VARIABLES = "X","Y", "Z","I", "J", "K", "phi","source","res","error","Aw","Ae","As","An","Ab","Af","Ap"'
  WRITE(72,*) 'VARIABLES = "X","Y", "Z","Aw","Ae","As","An","Ab","Af","Ap","source"'

  write(72,'(A,I5,A,I5,A,I5,A,I5)') 'ZONE T="ZONE ',MY_RANK,'", I=',nxc,', J=',nyc,', K=',nzc
  ! WRITE(71,*) 'ZONE  I=',nxc,', J=',nyc,', k =', nzc

  DO k=1,nzc
    DO j=1,nyc
      DO i=1,nxc
        WRITE(72,'(11F12.4)') xc(i), yc(j), zc(k), Aw(i,j,k), Ae(i,j,k), As(i,j,k), An(i,j,k),&
          Ab(i,j,k), Af(i,j,k), Ap(i,j,k), source(i,j,k)
      END DO
    END DO
  END DO
  CLOSE(72)
end subroutine write_dump
