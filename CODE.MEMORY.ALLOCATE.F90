SUBROUTINE MEMORY_ALLOCATE()

  USE INPUTDATA
  USE FIELDDATA
  USE COEFFDATA
  use FLOWDATA
  USE MPIDATA
  IMPLICIT NONE
  INTEGER :: iAlErr

!! ALLOCATE FIELD DATA
  IF(MY_RANK .eq. 0) WRITE(*,*)  "START MEM ALLOC"
  IF(MY_RANK .eq. 0) WRITE(*,*)  Nxc,Nyc,Nzc,Ngl
  IF(MY_RANK .eq. 0) WRITE(*,*)  Nx,Ny,Nz
  ALLOCATE(XPOINTS((Npx+1)))
  IF(MY_RANK .eq. 0) WRITE(*,*)  "TESTING"
  ALLOCATE(YPOINTS((Npy+1)))
  IF(MY_RANK .eq. 0) WRITE(*,*)  "TESTING2"
  ALLOCATE(ZPOINTS((Npz+1)))
  ALLOCATE(x_GLBL((1-Ngl):(Nx_glbl+Ngl)))
  ALLOCATE(y_GLBL((1-Ngl):(Ny_glbl+Ngl)))
  ALLOCATE(z_GLBL((1-Ngl):(Nz_glbl+Ngl)))
  ALLOCATE(x((1-Ngl):(Nx+Ngl)))
  ALLOCATE(y((1-Ngl):(Ny+Ngl)))
  ALLOCATE(z((1-Ngl):(Nz+Ngl)))
  ALLOCATE(xc((1-Ngl):(Nxc+Ngl)))
  ALLOCATE(yc((1-Ngl):(Nyc+Ngl)))
  ALLOCATE(zc((1-Ngl):(Nzc+Ngl)))
  ALLOCATE(dx((1):(Nx)))
  ALLOCATE(dy((1):(Ny)))
  ALLOCATE(dz((1):(Nz)))
  ALLOCATE(delx((1-Ngl):(Nxc+Ngl)))
  ALLOCATE(dely((1-Ngl):(Nyc+Ngl)))
  ALLOCATE(delz((1-Ngl):(Nzc+Ngl)))
  ALLOCATE(dxinv((1):(Nx)))
  ALLOCATE(dyinv((1):(Ny)))
  ALLOCATE(dzinv((1):(Nz)))
  ALLOCATE(delxinv((1-Ngl):(Nxc+Ngl)))
  ALLOCATE(delyinv((1-Ngl):(Nyc+Ngl)))
  ALLOCATE(delzinv((1-Ngl):(Nzc+Ngl)))
  IF(MY_RANK .eq. 0) WRITE(*,*)  "START MEM ALLOC: Flow data"
!! Allocate memory for flow data
  ALLOCATE(phi((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(rho((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(phi0((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(phi00((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(source((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(errorphi((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(residual((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(phi_ew((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(phi_ns((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(phi_fb((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(uphi_ew((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(uphi_ns((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)
  ALLOCATE(uphi_fb((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)), STAT=iAlErr)

!! Allocate memory for coefficients
  IF(MY_RANK .eq. 0) WRITE(*,*)  "START MEM ALLOC: Coeff data"
  ALLOCATE(Ae((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(Aw((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(An((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(As((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(Af((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(Ab((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(Ap((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  ALLOCATE(Apinv((1):(Nxc),(1):(Nyc),(1):(Nzc)), STAT=iAlErr)
  IF(MY_RANK .eq. 0) WRITE(*,*)  "START MEM ALLOC: BC data"
  ALLOCATE(VAL_w_ARR((1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)))
  ALLOCATE(VAL_e_ARR((1-Ngl):(Nyc+Ngl),(1-Ngl):(Nzc+Ngl)))
  ALLOCATE(VAL_s_ARR((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nzc+Ngl)))
  ALLOCATE(VAL_n_ARR((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nzc+Ngl)))
  ALLOCATE(VAL_b_ARR((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl)))
  ALLOCATE(VAL_f_ARR((1-Ngl):(Nxc+Ngl),(1-Ngl):(Nyc+Ngl)))
  IF(MY_RANK .eq. 0) WRITE(*,*) "END MEM ALLOC"

end subroutine MEMORY_ALLOCATE
