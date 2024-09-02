SUBROUTINE gausn2( iset )
USE data_db
IMPLICIT NONE
INTEGER, INTENT(IN OUT) :: iset

TYPE (truss), POINTER, SAVE :: eset
INTEGER :: i

iset = iset + 1   !update set number

IF( iset == 1)THEN    !for first set point the head
  eset => truss_head
ELSE                  !point to next set
  eset => eset%next
END IF

DO i=1,eset%nelem
  READ(16)
END DO

RETURN
END SUBROUTINE gausn2
