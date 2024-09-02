SUBROUTINE runend(message)
! This routine stops the run and prints the last read card
IMPLICIT NONE

  !Dummy variables
  CHARACTER(len=*),INTENT(IN):: message
  !Local variables

  IF (LEN_TRIM(message) > 0) WRITE(*,"(//,A)") TRIM(message)
  WRITE(*,"(/,'* * * AN ERROR HAS BEEN DETECTED * * *')")
  STOP

RETURN
END SUBROUTINE runend
