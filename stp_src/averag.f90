SUBROUTINE averag(npoin,nvarg,vargs,accpn,quads)
!
!  nodal averaging
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npoin,nvarg
REAL(kind=8), INTENT(IN) :: accpn(npoin)
REAL(kind=8), INTENT(IN OUT) :: vargs(nvarg,npoin)
LOGICAL, INTENT(OUT) :: quads

INTEGER :: n

quads = .FALSE.     !initialize, NO missing node
DO n=1,npoin
  IF( accpn(n) > 0) THEN
    vargs(:,n) = vargs(:,n)/accpn(n)
  ELSE
    quads = .TRUE.  !missing node
  END IF
END DO

RETURN
END SUBROUTINE averag
