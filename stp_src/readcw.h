SUBROUTINE readcw (areas,coord,ndime,ntype)
!
!   read surface data to compute surface contact wearing
!

IMPLICIT NONE
INTEGER (kind=4), INTENT(IN) :: ndime,ntype
REAL (kind=8), INTENT(IN) :: coord(:,:)
REAL (kind=8), INTENT(OUT) :: areas(:)

END SUBROUTINE readcw
