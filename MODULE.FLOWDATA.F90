MODULE FLOWDATA
  USE inputdata

  IMPLICIT NONE
  REAL(KIND=CGR), DIMENSION(:,:,:), ALLOCATABLE  :: phi, phi0, phi00, source, rho, errorphi, residual, PHIL
  REAL(KIND=CGR), DIMENSION(:,:,:), ALLOCATABLE  :: phi_ew, phi_ns, phi_fb
  REAL(KIND=CGR), DIMENSION(:,:,:), ALLOCATABLE  :: uphi_ew, uphi_ns, uphi_fb
end module FLOWDATA
