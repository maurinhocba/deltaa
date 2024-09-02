SUBROUTINE cstie1(rssdb,stiff,cprop,issdb,facto,nsymm,astif)

!.... COMPUTE STIFFNESS MATRIX FOR A THREE-NODE 2-D CONTACT ELEMENT

IMPLICIT NONE
!       arguments
INTEGER (kind=4), INTENT(IN) :: nsymm,issdb(:)
REAL (kind=8), INTENT(IN) :: rssdb(:),cprop(:),facto
REAL (kind=8), INTENT(OUT) :: stiff(:),astif(:)
END SUBROUTINE cstie1
