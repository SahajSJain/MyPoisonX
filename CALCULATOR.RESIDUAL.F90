SUBROUTINE CALCRESIDUAL(var,res )
  USE INPUTDATA
  USE COEFFDATA
  USE MPIDATA
  USE FIELDDATA
  USE FLOWDATA
  ! REAL(KIND=CGR), INTENT(IN), DIMENSION(:,:,:) :: var
  ! REAL(KIND=CGR), INTENT(OUT), DIMENSION(:,:,:) :: res 
  REAL(KIND=CGR), INTENT(IN), DIMENSION((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)) :: var
  REAL(KIND=CGR), INTENT(INOUT), DIMENSION((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)) :: res 

  INTEGER i,j,k 
  res=0;
  DO k=1,Nzc
    DO j=1, Nyc
      DO i=1, Nxc
        res(i,j,k)=Aw(i,j,k)*var(i-1,j,k)+Ae(i,j,k)*var(i+1,j,k)&
          +As(i,j,k)*var(i,j-1,k)+An(i,j,k)*var(i,j+1,k)&
          +Ab(i,j,k)*var(i,j,k-1)+Af(i,j,k)*var(i,j,k+1)&
          +var(i,j,k)*Ap(i,j,k)-source(i,j,k)
      enddo !i
    ENDDO !j
  ENDDO !k

END SUBROUTINE CALCRESIDUAL