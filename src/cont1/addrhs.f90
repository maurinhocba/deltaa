SUBROUTINE addrhs(fcont,eresf,ndime,ien,nen,xfact)

!.... add element right-hand-side to global right-hand-side

IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: ndime,nen,ien(nen)
REAL (kind=8), INTENT(IN) ::  eresf(ndime,nen),xfact
REAL (kind=8), INTENT(IN OUT) ::  fcont(ndime,*)
!     local variables
INTEGER (kind=4) i,n

DO i = 1,nen
  n = ien(i)
  fcont(1:ndime,n) = fcont(1:ndime,n) - eresf(1:ndime,i)*xfact
END DO
RETURN
END SUBROUTINE addrhs
