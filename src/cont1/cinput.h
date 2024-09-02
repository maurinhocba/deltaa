SUBROUTINE cinput(maxve,iwrit,label,npoin,coord)

!.... contact input routine (NEW problem)

USE cont1_db  !INTENT(OUT)

IMPLICIT NONE
!     arguments
INTEGER (kind=4),INTENT (IN) :: iwrit,npoin,label(:),maxve
REAL (kind=8), INTENT (IN) :: coord(:,:)

END SUBROUTINE cinput
