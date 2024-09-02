SUBROUTINE cinpu4(maxcc,iwrit,label,npoin,coord)

!.... contact input routine (NEW problem)

!USE cont4_db  !INTENT(OUT)

IMPLICIT NONE
!     arguments
INTEGER (kind=4),INTENT (IN) :: iwrit,npoin,label(:),maxcc
REAL (kind=8), INTENT (IN) :: coord(:,:)
END SUBROUTINE cinpu4
