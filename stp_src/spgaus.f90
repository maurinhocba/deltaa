SUBROUTINE spgaus()
! Process special Gauss data after smoothing
USE data_db,ONLY: bst_nps
IMPLICIT NONE

  !Local variables

  IF (bst_nps > 0) CALL sgaus12()

RETURN
END SUBROUTINE spgaus
