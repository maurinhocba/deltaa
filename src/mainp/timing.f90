      SUBROUTINE timing (k,ind)

      ! initializes and add CPU times for different tasks

      USE outp_db, ONLY : time
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: k
      LOGICAL,INTENT(IN) :: ind

      REAL (kind=8) :: tcpu

      CALL timuse (tcpu)                          !get clock time
      IF (ind) THEN                            !start task
        time(k+20) = tcpu                         !initial time
      ELSE                                     !end task
        time(k) = time(k) + tcpu - time(k+20)     !adds elapsed time
      END IF
      RETURN

      END SUBROUTINE timing
