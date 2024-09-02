      SUBROUTINE endstr
!********************************************************************
!
!*** END STRATEGY card READ
!
!********************************************************************
      USE c_input, ONLY : listen,exists,backs
      IMPLICIT NONE

      CALL listen('EQNUMS')
      IF( .NOT.exists('ENDSTR') )THEN
        WRITE(*,"(' Warning key-word END_STR not found')")
        backs = .TRUE.
      END IF

      RETURN
      END SUBROUTINE endstr
