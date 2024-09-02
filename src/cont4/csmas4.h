SUBROUTINE csmas4(lcseg,x,emass,nsegm,density)
! Compute nodal mass for contact surfaces if required
IMPLICIT NONE
!dummy arguments
INTEGER (kind=4), INTENT(IN) :: nsegm,lcseg(:,:)
REAL (kind=8), INTENT(IN) :: x(:,:),density
REAL (kind=8), INTENT(IN OUT) :: emass(:,:)
END SUBROUTINE csmas4
