      SUBROUTINE nearst(xc,nsegm,nearn,xs)

!.... look for the nearest master segment using a GLOBAL searching algorithm
!
!.... input
!....   xc = current spatial coordinates of the segment center
!....   xs = current spatial coordinates of the slave node
!.... output
!....   nearn = number of the nearest master segment

      IMPLICIT NONE
!     arguments
      INTEGER (kind=4), INTENT (IN) :: nsegm
      INTEGER (kind=4), INTENT (OUT) :: nearn
      REAL (kind=8), INTENT (IN) :: xc(:,:),xs(:)
      END SUBROUTINE nearst
