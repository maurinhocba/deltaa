SUBROUTINE gausn8( iset )
USE data_db
IMPLICIT NONE
INTEGER, INTENT(IN OUT) :: iset

TYPE (beame), POINTER, SAVE :: eset
INTEGER ::  i,g

iset = iset + 1   !update set number

IF( iset == 1)THEN     !for first set point the head
  eset => beame_head
ELSE                   !point to next set
  eset => eset%next
END IF

DO i=1,eset%nelem
  DO g=1,eset%ngaus
    READ(16)
  END DO
END DO

RETURN
END SUBROUTINE gausn8
