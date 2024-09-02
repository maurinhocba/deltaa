SUBROUTINE cstie4(rssdb,stiff,cprop,issdb,facto)

!.... compute stiffness matrix for a four-node 3-d contact element

IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: issdb(:)
REAL (kind=8), INTENT(IN) :: rssdb(:),cprop(:),facto
REAL (kind=8), INTENT(OUT) :: stiff(:)
END SUBROUTINE cstie4
