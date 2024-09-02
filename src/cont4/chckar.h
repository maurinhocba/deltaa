SUBROUTINE chckar(nsegm,nnseg,lnods,x,label)
!
!     Checks aspect ratio of triangles defining a master surface
!     writes warnings when a specified value is exceeded
!
IMPLICIT NONE
INTEGER (KIND=4), INTENT(IN) :: nsegm,nnseg,lnods(:,:),label(:)
REAL (KIND=8), INTENT(IN) :: x(:,:)

END SUBROUTINE chckar
