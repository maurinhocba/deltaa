SUBROUTINE cfrmlm(ien,nen,dtime,xfact,r,indcon)

!.... form lm array for a given element and computes factor
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: nen,ien(nen),indcon
REAL    (kind=8), INTENT(IN) :: dtime,r(*)
REAL    (kind=8), INTENT(OUT) ::  xfact
END SUBROUTINE cfrmlm
