subroutine setup_grid
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER i, garbage, iFilein, writeProc
  INTEGER istart, iend, Npoints
  ! INTEGER, DIMENSION(1:(Npx+1)) :: XPOINTS
  ! INTEGER, DIMENSION(1:(Npy+1)) :: YPOINTS
  ! INTEGER, DIMENSION(1:(Npz+1)) :: ZPOINTS
  REAL(KIND=CGR) STEP

  IF(ImtheBOSS)  WRITE(*,*) "START SETUP GRID"
  writeProc=3;
  iFilein=999
  if(x_unif .eq. 1) then
    x_GLBL(1:Nx_GLBL)=LINSPACE(x_start, x_end, Nx_GLBL)
  else
    iFilein=iFilein+1;
    OPEN(iFilein, File="xgrid.dat", ACTION='READ', STATUS='OLD')
    do i=1,Nx_GLBL
      READ(ifINPUTDAT,*) garbage, x_GLBL(i)
    enddo
  endif
  if(y_unif .eq. 1) then
    y_GLBL(1:Ny_GLBL)=LINSPACE(y_start, y_end, Ny_GLBL)
  else
    iFilein=iFilein+1;
    OPEN(iFilein, File="ygrid.dat", ACTION='READ', STATUS='OLD')
    do i=1,Ny_GLBL
      READ(ifINPUTDAT,*) garbage, y_GLBL(i)
    enddo
  endif
  if(z_unif .eq. 1) then
    z_GLBL(1:Nz_GLBL)=LINSPACE(z_start, z_end, Nz_GLBL)
  else
    iFilein=iFilein+1;
    OPEN(iFilein, File="zgrid.dat", ACTION='READ', STATUS='OLD')
    do i=1,Nz_GLBL
      READ(ifINPUTDAT,*) garbage, z_GLBL(i)
    enddo
  endif
  ! IF(ImtheBOSS) WRITE(*,*) z_GLBL


  !! SET GHOST VALUES
  do i=1, NGL
    STEP=x_GLBL(2)-x_GLBL(1);
    x_GLBL(1-i)=x_GLBL(1-i+1)-STEP
    STEP=x_GLBL(Nx_GLBL)-x_GLBL(Nx_GLBL-1);
    x_GLBL(Nx_GLBL+i)=x_GLBL(Nx_GLBL+i-1)+STEP

    STEP=y_GLBL(2)-y_GLBL(1);
    y_GLBL(1-i)=y_GLBL(1-i+1)-STEP
    STEP=y_GLBL(Ny_GLBL)-y_GLBL(Ny_GLBL-1);
    y_GLBL(Ny_GLBL+i)=y_GLBL(Ny_GLBL+i-1)+STEP

    STEP=z_GLBL(2)-z_GLBL(1);
    z_GLBL(1-i)=z_GLBL(1-i+1)-STEP
    STEP=z_GLBL(Nz_GLBL)-z_GLBL(Nz_GLBL-1);
    z_GLBL(Nz_GLBL+i)=z_GLBL(Nz_GLBL+i-1)+STEP
  ENDDO

  XPOINTS(1)=1
  do i=1,Npx
    XPOINTS(i+1)=XPOINTS(i)+Nx-1;
  enddo

  YPOINTS(1)=1
  do i=1,Npy
    YPOINTS(i+1)=YPOINTS(i)+Ny-1;
  enddo

  ZPOINTS(1)=1
  do i=1,Npz
    ZPOINTS(i+1)=ZPOINTS(i)+Nz-1;
  enddo
!   IF(myCartRank .eq. writeProc) WRITE(*,*) XPOINTS

  istart=xpoints(myCOORDS(1)+1)-NGL;
  iend=xpoints(myCOORDS(1)+2)+NGL
!   IF(myCartRank .eq. writeProc) WRITE(*,*) (1-NGL), (Nx+Ngl), istart,iend
  x=x_GLBL(istart:iend )
  istart=ypoints(myCOORDS(2)+1)-NGL;
  iend=ypoints(myCOORDS(2)+2)+NGL
  y=y_GLBL(istart:iend)
  istart=zpoints(myCOORDS(3)+1)-NGL;
  iend=zpoints(myCOORDS(3)+2)+NGL
  z=z_GLBL(istart:iend)
  ! IF(myCartRank .eq. writeProc) WRITE(*,*) x_GLBL
  ! IF(myCartRank .eq. writeProc) WRITE(*,*) x

  do i=1-Ngl, Nxc+NGL
    xc(i)=(x(i)+x(i+1))*0.50d0
    delx(i)=x(i+1)-x(i)
    delxinv(i)=1/delx(i);
  enddo
  do i=1-Ngl, Nyc+NGL
    yc(i)=(y(i)+y(i+1))*0.50d0
    dely(i)=y(i+1)-y(i)
    delyinv(i)=1/dely(i);
  enddo
  do i=1-Ngl, Nzc+NGL
    zc(i)=(z(i)+z(i+1))*0.50d0
    delz(i)=z(i+1)-z(i)
    delzinv(i)=1/delz(i);
  enddo
  do i=1, Nx
    dx(i)=xc(i)-xc(i-1)
    dxinv(i)=1/dx(i)
  enddo
  do i=1, Ny
    dy(i)=yc(i)-yc(i-1)
    dyinv(i)=1/dy(i)
  enddo
  do i=1, Nz
    dz(i)=zc(i)-zc(i-1)
    dzinv(i)=1/dz(i)
  enddo
  if (x_unif .eq. 1) then
    delx = (x_end - x_start) / REAL(Nx_GLBL - 1, KIND=CGR)
    dx = (x_end - x_start) / REAL(Nx_GLBL - 1, KIND=CGR)
    delxinv = REAL(Nx_GLBL - 1, KIND=CGR) / (x_end - x_start)
    dxinv = REAL(Nx_GLBL - 1, KIND=CGR) / (x_end - x_start)
  end if

  if (y_unif .eq. 1) then
    dely = (y_end - y_start) / REAL(Ny_GLBL - 1, KIND=CGR)
    dy = (y_end - y_start) / REAL(Ny_GLBL - 1, KIND=CGR)
    delyinv = REAL(Ny_GLBL - 1, KIND=CGR) / (y_end - y_start)
    dyinv = REAL(Ny_GLBL - 1, KIND=CGR) / (y_end - y_start)
  end if

  if (z_unif .eq. 1) then
    delz = (z_end - z_start) / REAL(Nz_GLBL - 1, KIND=CGR)
    dz = (z_end - z_start) / REAL(Nz_GLBL - 1, KIND=CGR)
    delzinv = REAL(Nz_GLBL - 1, KIND=CGR) / (z_end - z_start)
    dzinv = REAL(Nz_GLBL - 1, KIND=CGR) / (z_end - z_start)
  end if!   IF(myCartRank .eq. writeProc) WRITE(*,*) zc
  ! IF(ImtheBOSS) WRITE(*,*) dxinv
  ! IF(myCartRank .eq. writeProc) WRITE(*,*) dz
CONTAINS
  FUNCTION LINSPACE(FROM, TOf, NUM) RESULT(ARRAY)
    USE INPUTDATA
    IMPLICIT NONE
    REAL(KIND=CGR), INTENT(IN) :: FROM, TOf
    INTEGER, INTENT(IN) :: NUM
    REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE :: ARRAY
    INTEGER :: ILS
    REAL(KIND=CGR) :: STEPLS
    ! Check if NUM is less than or equal to 1
    IF (NUM <= 1) THEN
      ALLOCATE(ARRAY(1))
      ARRAY(1) = FROM
      RETURN
    END IF
    ! Allocate the array for the result
    ALLOCATE(ARRAY(NUM))
    ! Calculate the step size
    STEPLS = (TOf - FROM) / (NUM - 1)
    DO ILS = 1, NUM
      ARRAY(ILS) = FROM + (ILS - 1) * STEPLS
    END DO

  end function LINSPACE
end subroutine setup_grid

SUBROUTINE SETUP_COEFFS
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(KIND=CGR) rhoe, rhow, rhon, rhos, rhob, rhof
  rho=1; !! LATER ADD VARYING RHO FOR MULTIPHASE FLOWS

  DO k=1,Nzc
    DO j=1, Nyc
      DO i=1, Nxc
        rhoe=0.5*(rho(i,j,k)+rho(i+1,j,k))
        rhow=0.5*(rho(i,j,k)+rho(i-1,j,k))
        rhon=0.5*(rho(i,j,k)+rho(i,j+1,k))
        rhos=0.5*(rho(i,j,k)+rho(i,j-1,k))
        rhof=0.5*(rho(i,j,k)+rho(i,j,k+1))
        rhob=0.5*(rho(i,j,k)+rho(i,j,k-1))
        rhoe=1.
        rhow=1.
        rhon=1.
        rhos=1.
        rhof=1.
        rhob=1.
        Aw(i,j,k)=delxinv(i)*dxinv(i)/rhow;
        Ae(i,j,k)=delxinv(i)*dxinv(i+1)/rhoe;
        As(i,j,k)=delyinv(j)*dxinv(j)/rhos;
        An(i,j,k)=delyinv(j)*dxinv(j+1)/rhon;
        Ab(i,j,k)=delzinv(k)*dxinv(k)/rhob;
        Af(i,j,k)=delzinv(k)*dxinv(k+1)/rhof;
        Ap(i,j,k)= -(Aw(i,j,k)+Ae(i,j,k)+As(i,j,k)+An(i,j,k)+Ab(i,j,k)+Af(i,j,k))
        Apinv(i,j,k) = 1/Ap(i,j,k)
      enddo
    ENDDO
  ENDDO
end subroutine SETUP_COEFFS
SUBROUTINE SETUP_SOURCE
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(KIND=CGR) :: xp,yp,zp
  interface
    function rhscalc(xa, ya, za) result(output)
      USE INPUTDATA
      implicit none
      real(kind=CGR), intent(in) :: xa, ya, za
      real(kind=CGR) :: output
    end function rhscalc
  end interface
  SELECT CASE (SOURCEFLAG)
   CASE (1) !! Laplace equation
    source=0.0;
   CASE (2) !! random number
    call random_seed()
    call random_number(source)
   CASE (3) !! Customized
    DO k=1,Nzc
      DO j=1, Nyc
        DO i=1, Nxc
          xp=xc(i)
          yp=yc(j)
          zp=zc(k)
          source(i,j,k)=rhscalc(xp,yp,zp);
        enddo
      ENDDO
    ENDDO
   CASE DEFAULT
    source=0.0
  END SELECT

end subroutine SETUP_SOURCE

SUBROUTINE SETUP_BOUNDARYCONDITIONS
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(KIND=CGR) :: xp, yp, zp
  interface
    function phicalc(xa, ya, za) result(output)
      USE INPUTDATA
      implicit none
      real(kind=CGR), intent(in) :: xa, ya, za
      real(kind=CGR) :: output
    end function phicalc
  end interface



  VAL_w_ARR=VAL_WEST
  VAL_e_ARR=VAL_EAST
  VAL_s_ARR=VAL_SOUTH
  VAL_n_ARR=VAL_NORTH
  VAL_b_ARR=VAL_BACK
  VAL_f_ARR=VAL_FRONT
  !! west
  if(myWest .eq. -1) then
    SELECT CASE (BC_WEST)
     CASE (1) !! DIRICHLET BC
      BCT_WEST=1.0;
      DO j=1, Nyc
        DO k=1, Nzc
          xp=x(1)
          yp=yc(j)
          zp=zc(k)
          VAL_w_ARR(j,k)=VAL_WEST
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_WEST=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_WEST=1.0;
      DO j=1, Nyc
        DO k=1, Nzc
          xp=x(1)
          yp=yc(j)
          zp=zc(k)
          VAL_w_ARR(j,k)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      BC_WEST=1;
      ! WRITE(*,*) "CUSTOM BC WEST"
    END SELECT
  end if
  !! east
  if(myEast .eq. -1) then
    SELECT CASE (BC_EAST)
     CASE (1) !! DIRICHLET BC
      BCT_EAST=1.0;
      DO j=1, Nyc
        DO k=1, Nzc
          xp=x(Nx)
          yp=yc(j)
          zp=zc(k)
          VAL_e_ARR(j,k)=VAL_EAST
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_EAST=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_EAST=1.0;
      DO j=1, Nyc
        DO k=1, Nzc
          xp=x(Nx)
          yp=yc(j)
          zp=zc(k)
          VAL_e_ARR(j,k)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      ! WRITE(*,*) "CUSTOM BC EAST"
      BC_EAST=1;
    END SELECT
  end if
!! South
  if(mySouth .eq. -1) then
    SELECT CASE (BC_SOUTH)
     CASE (1) !! DIRICHLET BC
      BCT_SOUTH=1.0;
      DO i=1, Nx
        DO k=1, Nzc
          xp=xc(i)
          yp=y(1)
          zp=zc(k)
          VAL_s_ARR(i,k)=VAL_SOUTH
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_SOUTH=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_SOUTH=1.0;
      DO i=1, Nx
        DO k=1, Nzc
          xp=xc(i)
          yp=y(1)
          zp=zc(k)
          VAL_s_ARR(i,k)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      BC_SOUTH=1;
    END SELECT
  end if

!! North
  if(myNorth .eq. -1) then
    SELECT CASE (BC_NORTH)
     CASE (1) !! DIRICHLET BC
      BCT_NORTH=1.0;
      DO i=1, Nx
        DO k=1, Nzc
          xp=xc(i)
          yp=y(Ny)
          zp=zc(k)
          VAL_n_ARR(i,k)=VAL_NORTH
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_NORTH=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_NORTH=1.0;
      DO i=1, Nx
        DO k=1, Nzc
          xp=xc(i)
          yp=y(Ny)
          zp=zc(k)
          VAL_n_ARR(i,k)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      BC_NORTH=1;
    END SELECT
  end if
!! Back
  if(myBack .eq. -1) then
    SELECT CASE (BC_BACK)
     CASE (1) !! DIRICHLET BC
      BCT_BACK=1.0;
      DO i=1, Nx
        DO j=1, Nyc
          xp=xc(i)
          yp=yc(j)
          zp=z(1)
          VAL_b_ARR(i,j)=VAL_BACK
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_BACK=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_BACK=1.0;
      DO i=1, Nx
        DO j=1, Nyc
          xp=xc(i)
          yp=yc(j)
          zp=z(1)
          VAL_b_ARR(i,j)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      BC_BACK=1;
    END SELECT
  end if

!! Front
  if(myFront .eq. -1) then
    SELECT CASE (BC_FRONT)
     CASE (1) !! DIRICHLET BC
      BCT_FRONT=1.0;
      DO i=1, Nx
        DO j=1, Nyc
          xp=xc(i)
          yp=yc(j)
          zp=z(Nz)
          VAL_f_ARR(i,j)=VAL_FRONT
        ENDDO
      ENDDO
     CASE (2) !! NEUMANN
      BCT_FRONT=0.0
     CASE (3) !! CUSTOM DIRICHLET
      BCT_FRONT=1.0;
      DO i=1, Nx
        DO j=1, Nyc
          xp=xc(i)
          yp=yc(j)
          zp=z(Nz)
          VAL_f_ARR(i,j)=phicalc(xp,yp,zp);
        ENDDO
      ENDDO
      BC_FRONT=1;
    END SELECT
  end if

  ! West boundary
  if((myWest .eq. -1) .and. (BC_WEST .eq. 1)) then
    DO j=1, Nyc
      DO k=1, Nzc
        source(1,j,k)=source(1,j,k)-2.0*Val_w_arr(j,k)*Aw(1,j,k)
        Ap(1,j,k)=Ap(1,j,k)-Aw(1,j,k)
        Aw(1,j,k)=0.0
      ENDDO
    ENDDO
  end if

! East boundary
  if((myEast .eq. -1) .and. (BC_EAST .eq. 1)) then
    DO j=1, Nyc
      DO k=1, Nzc
        source(Nxc,j,k)=source(Nxc,j,k)-2.0*Val_e_arr(j,k)*Ae(Nxc,j,k)
        Ap(Nxc,j,k)=Ap(Nxc,j,k)-Ae(Nxc,j,k)
        Ae(Nxc,j,k)=0.0
      ENDDO
    ENDDO
  end if

! South boundary
  if((mySouth .eq. -1) .and. (BC_SOUTH .eq. 1)) then
    DO i=1, Nxc
      DO k=1, Nzc
        source(i,1,k)=source(i,1,k)-2.0*Val_s_arr(i,k)*As(i,1,k)
        Ap(i,1,k)=Ap(i,1,k)-As(i,1,k)
        As(i,1,k)=0.0
      ENDDO
    ENDDO
  end if

! North boundary
  if((myNorth .eq. -1) .and. (BC_NORTH .eq. 1)) then
    DO i=1, Nxc
      DO k=1, Nzc
        source(i,Nyc,k)=source(i,Nyc,k)-2.0*Val_n_arr(i,k)*An(i,Nyc,k)
        Ap(i,Nyc,k)=Ap(i,Nyc,k)-An(i,Nyc,k)
        An(i,Nyc,k)=0.0
      ENDDO
    ENDDO
  end if

! Back boundary
  if((myBack .eq. -1) .and. (BC_BACK .eq. 1)) then
    DO i=1, Nxc
      DO j=1, Nyc
        source(i,j,1)=source(i,j,1)-2.0*Val_b_arr(i,j)*Ab(i,j,1)
        Ap(i,j,1)=Ap(i,j,1)-Ab(i,j,1)
        Ab(i,j,1)=0.0
      ENDDO
    ENDDO
  end if

! Front boundary
  if((myFront .eq. -1) .and. (BC_FRONT .eq. 1)) then
    DO i=1, Nxc
      DO j=1, Nyc
        source(i,j,Nzc)=source(i,j,Nzc)-2.0*Val_f_arr(i,j)*Af(i,j,Nzc)
        Ap(i,j,Nzc)=Ap(i,j,Nzc)-Af(i,j,Nzc)
        Af(i,j,Nzc)=0.0
      ENDDO
    ENDDO
  end if

  DO k=1,Nzc
    DO j=1, Nyc
      DO i=1, Nxc
        Apinv(i,j,k) = 1/Ap(i,j,k)
      enddo
    ENDDO
  ENDDO

end subroutine SETUP_BOUNDARYCONDITIONS

function phicalc(xa, ya, za) result(output)
  USE INPUTDATA
  implicit none
  ! Input variables
  real(kind=CGR), intent(in) :: xa, ya, za
  ! Output variable
  real(kind=CGR) :: output
  output = sin(Pie*xa) * sin(2.0d0*Pie*ya) * sin(4.0d0*Pie*za)
end function phicalc

function rhscalc(xa, ya, za) result(output)
  USE INPUTDATA
  implicit none
  ! Input variables
  real(kind=CGR), intent(in) :: xa, ya, za
  ! Output variable
  real(kind=CGR) :: output
  output = -21.0d0 * sin(Pie*xa) * sin(2.0d0*Pie*ya) * sin(4.0d0*Pie*za)
end function rhscalc
