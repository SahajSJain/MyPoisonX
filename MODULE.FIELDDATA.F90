MODULE FIELDDATA
  USE inputdata

  IMPLICIT NONE
  REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE  :: x, y, z,  xc, yc, zc
  REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE  :: x_GLBL, y_GLBL, z_GLBL,  xc_GLBL, yc_GLBL, zc_GLBL
  REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE  :: dx, dy, dz,  delx, dely, delz !! stored locally
  REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE  :: dxinv, dyinv, dzinv
  REAL(KIND=CGR), DIMENSION(:), ALLOCATABLE  ::  delxinv, delyinv, delzinv
  REAL(KIND=CGR) :: dxmin, dxmax, dymin, dymax, dzmin, dzmax
  REAL(KIND=CGR), DIMENSION(:,:), ALLOCATABLE :: VAL_w_ARR, VAL_e_ARR !j,k
  REAL(KIND=CGR), DIMENSION(:,:), ALLOCATABLE :: VAL_s_ARR, VAL_n_ARR !k,i
  REAL(KIND=CGR), DIMENSION(:,:), ALLOCATABLE :: VAL_f_ARR, VAL_b_ARR !i,j
  INTEGER, DIMENSION(:), ALLOCATABLE  :: XPOINTS , YPOINTS, ZPOINTS
end module FIELDDATA
