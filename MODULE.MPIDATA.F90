MODULE MPIDATA
  USE INPUTDATA
  IMPLICIT NONE
  INTEGER :: COMM_SZ, MY_RANK, Nptot
  REAL(KIND=CGR) :: PROGRAM_TIME_START, PROGRAM_TIME_END
  LOGICAL :: ImtheBOSS
  INTEGER :: myWest, myEast, mySouth, myNorth, myFront, myBack
  INTEGER :: myCartRank
  INTEGER, DIMENSION(3) :: NumProcArr, myCOORDS
  LOGICAL, DIMENSION(3) :: PeriodicArr
  INTEGER :: CommCart, myIs, myJs, myKs
  INTEGER :: type_W_send, type_W_recv
  INTEGER :: type_E_send, type_E_recv
  INTEGER :: type_S_send, type_S_recv
  INTEGER :: type_N_send, type_N_recv
  INTEGER :: type_B_send, type_B_recv
  INTEGER :: type_F_send, type_F_recv

  INTEGER :: tagWs, tagWr
  INTEGER :: tagEs, tagEr
  INTEGER :: tagNs, tagNr
  INTEGER :: tagSs, tagSr
  INTEGER :: tagBs, tagBr
  INTEGER :: tagFs, tagFr

end module MPIDATA
