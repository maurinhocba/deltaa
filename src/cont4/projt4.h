SUBROUTINE projt4(xs,nearn,x,lcseg,   &
                  nhseg,issdb,prdat,cutof,gapin,curms,cu)

!.... project the slave node onto a 3-d master triangular surface

IMPLICIT NONE
!     arguments
      LOGICAL, INTENT(IN) :: curms
INTEGER (kind=4), INTENT(IN) :: lcseg(:,:),nhseg(:,:)
INTEGER (kind=4), INTENT(IN OUT) :: nearn,issdb(:)
REAL (kind=8), INTENT (IN) :: x(:,:),cutof,gapin,xs(:),cu(:,:)
REAL (kind=8), INTENT (OUT) :: prdat(:)

END SUBROUTINE projt4
