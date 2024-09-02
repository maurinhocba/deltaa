 SUBROUTINE conini(actio,  maxve)
 !********************************************************************
 !
 !*** contact definition card READ
 !
 !********************************************************************
 !USE c_input
 !USE cont_db
 !USE npoi_db, ONLY : ndime,npoin
 !USE outp_db, ONLY : nprqc, nreqc
 IMPLICIT NONE
 CHARACTER(len=6 ), INTENT(IN OUT) :: actio
 INTEGER(kind=4),INTENT(IN OUT):: maxve
 END SUBROUTINE conini
