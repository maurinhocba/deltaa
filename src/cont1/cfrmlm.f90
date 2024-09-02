SUBROUTINE cfrmlm(ien,nen,dtime,xfact,r,indcon)

!.... form lm array for a given element and computes factor
USE c_input, ONLY : runend, lures
USE npoi_db  !INTENT(IN) :: ndime,ndofn,emass(:,:),label(:,:)
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: nen,ien(nen),indcon
REAL    (kind=8), INTENT(IN) :: dtime,r(*)
REAL    (kind=8), INTENT(OUT) ::  xfact
!     local variables
INTEGER (kind=4) n,node
REAL    (kind=8) m,maux,sh(5)

!     form shape functions
sh(1) = 1d0
SELECT CASE (nen)
!CASE (2)
!  sh(2) = 1d0
CASE (3)
  sh(2) = 1d0-r(1)
  sh(3) = r(1)
CASE (4)
  sh(2) = 1d0-r(1)-r(2)
  sh(3) = r(1)
  sh(4) = r(2)
!CASE (5)
!  sh(2) = (1d0-r(1))*(1d0-r(2))*0.25d0
!  sh(3) = (1d0+r(1))*(1d0-r(2))*0.25d0
!  sh(4) = (1d0+r(1))*(1d0+r(2))*0.25d0
!  sh(5) = (1d0-r(1))*(1d0+r(2))*0.25d0
END SELECT
!     computes equivalent mass
m = 0d0
DO n = 1, nen
  node = ien(n)
  maux = MAXVAL(emass(1:ndime,node))
  IF(maux > 0) m = m + sh(n)**2/maux
  IF(indcon /= 0)EXIT
END DO
!     computes xfact
IF(m > 0) THEN
  xfact = 2d0/(m*dtime**2)
ELSE
  WRITE(*,*)'  Error in nodal mass ',label(ien(1:nen))
  WRITE(lures,*)'  Error in nodal mass  ',label(ien(1:nen))
  WRITE(lures,*)'  Internal node numbers',ien(1:nen)
  CALL runend('CONTAC: zero masses in contact pair')
END IF
RETURN
END SUBROUTINE cfrmlm
