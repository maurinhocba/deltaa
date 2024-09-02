SUBROUTINE celmn4(ttime,iwric,coord,emass,fcont,coorb,coort)

!.... perform contact searching & computes contact forces

USE cont4_db    !INTENT(IN OUT)
IMPLICIT NONE
!     Dummy arguments
INTEGER (kind=4), INTENT(IN) :: iwric                      !INTENT(IN)
REAL (kind=8), INTENT(IN) :: ttime,emass(:,:)              !INTENT(IN)
REAL (kind=8), POINTER :: coord(:,:),coorb(:,:),coort(:,:) !INTENT(IN)
REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)                !INTENT(OUT)
END SUBROUTINE celmn4
