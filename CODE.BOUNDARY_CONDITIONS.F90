SUBROUTINE BOUNDARY_CONDITIONS(var)
  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL :: xp, yp, zp
  REAL, INTENT(INOUT),DIMENSION(:,:,:) :: var
  !! west
!   if(myWest .eq. -1) then
!     SELECT CASE(BC_WEST)
!      CASE (1) !! Dirichlet
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(0,j,k)=(2.*VAL_WEST-var(1,j,k))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(0,j,k)=var(1,j,k)
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(0,j,k)=(2.*VAL_w_ARR(j,k)-var(1,j,k))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if
!   !! east
!   if(myEast .eq. -1) then
!     SELECT CASE(BC_EAST)
!      CASE (1) !! Dirichlet
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(Nxc+1,j,k)=(2.*VAL_EAST-var(Nxc,j,k))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(Nxc+1,j,k)=(var(Nxc,j,k))
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO k=1-Ngl, Nzc+Ngl
!         DO j=1-Ngl, Nyc+Ngl
!           var(Nxc+1,j,k)=(2.*VAL_e_ARR(j,k)-var(Nxc,j,k))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if

! !! south
!   if(mySouth .eq. -1) then
!     SELECT CASE(BC_SOUTH)
!      CASE (1) !! Dirichlet
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,0,k)=(2.*VAL_SOUTH-var(i,1,k))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,0,k)=var(i,1,k)
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,0,k)=(2.*VAL_s_ARR(i,k)-var(i,1,k))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if

! !! north
!   if(myNorth .eq. -1) then
!     SELECT CASE(BC_NORTH)
!      CASE (1) !! Dirichlet
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,Nyc+1,k)=(2.*VAL_NORTH-var(i,Nyc,k))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,Nyc+1,k)=(var(i,Nyc,k))
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO k=1-Ngl, Nzc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,Nyc+1,k)=(2.*VAL_n_ARR(i,k)-var(i,Nyc,k))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if

! !! back
!   if(myBack .eq. -1) then
!     SELECT CASE(BC_BACK)
!      CASE (1) !! Dirichlet
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,0)=(2.*VAL_BACK-var(i,j,1))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,0)=var(i,j,1)
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,0)=(2.*VAL_b_ARR(i,j)-var(i,j,1))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if

! !! front
!   if(myFront .eq. -1) then
!     SELECT CASE(BC_FRONT)
!      CASE (1) !! Dirichlet
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,Nzc+1)=(2.*VAL_FRONT-var(i,j,Nzc))
!         ENDDO
!       ENDDO
!      CASE (2) !! NEUMANN
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,Nzc+1)=(var(i,j,Nzc))
!         ENDDO
!       ENDDO
!      CASE (3) !! CUSTOM DIRICHLET
!       DO j=1-Ngl, Nyc+Ngl
!         DO i=1-Ngl, Nxc+Ngl
!           var(i,j,Nzc+1)=(2.*VAL_f_ARR(i,j)-var(i,j,Nzc))
!         ENDDO
!       ENDDO
!     END SELECT
!   end if
end subroutine BOUNDARY_CONDITIONS
