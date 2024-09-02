SUBROUTINE projtr(imn,x,lcseg,xita,eta,y)
!     project over a triangle
IMPLICIT NONE

!     arguments
INTEGER (kind=4), INTENT (IN) :: imn,lcseg(:)
REAL (kind=8), INTENT (IN) :: x(:,:),y(:)
REAL (kind=8), INTENT (OUT) :: xita,eta
END SUBROUTINE projtr
