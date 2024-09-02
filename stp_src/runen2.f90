SUBROUTINE runen2(message)
! This routine stops the run when a write error is detected
IMPLICIT NONE

  !Dummy variables
  CHARACTER(len=*), INTENT(IN) :: message

  IF(LEN_TRIM(message) > 0) WRITE(*,"(//,A)") TRIM(message)   !print message to screen
  WRITE(*,1)   !Print "disk full" message
  STOP         !stop execution

RETURN
 1 FORMAT(//,'  AN ERROR HAS BEEN DETECTED:',                             &
           /,'  FAILED WRITTING, THE DISK IS POSSIBLY FULL.',             &
          //,'  PLEASE, CHECK AND FREE DISK SPACE BEFORE CONTINUING.')
END SUBROUTINE runen2
