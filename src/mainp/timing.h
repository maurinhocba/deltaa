SUBROUTINE timing (pos,ind)

!routine to increase times devoted to different tasks

IMPLICIT NONE
LOGICAL, INTENT(IN) :: ind    !.TRUE. = START    .FALSE. = END
INTEGER (kind=4),INTENT(IN) :: pos
END SUBROUTINE timing
