      SUBROUTINE projt1(xs,nearn,x,isn,lcseg,   &
                        nhseg,issdb,prdat,cutof,gapin,curms,cu)

!.... project the slave node onto a 2-d master segment

      IMPLICIT NONE
!     arguments
      LOGICAL, INTENT(IN) :: curms
      INTEGER (kind=4), INTENT(IN) :: isn,lcseg(:,:),nhseg(:,:)
      INTEGER (kind=4), INTENT(IN OUT) :: nearn,issdb(:)
      REAL (kind=8), INTENT (IN) :: x(:,:),cutof,gapin,xs(:),cu(:)
      REAL (kind=8), INTENT (IN OUT) :: prdat(:)

      END SUBROUTINE projt1
