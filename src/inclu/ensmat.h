SUBROUTINE ensmat(nvarl,lm,locst,glost,elem)
!*************************************************************************
!
!     assembles a local symmetric matrix into the global symmetric matrix
!     both matrices are stored as arrays of the upper triangle only
!
!*************************************************************************
!USE kinc_db, ONLY : nn,maxa,npsdf,nesdf,ftsdf,maxav
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl)
REAL (kind=8),INTENT(IN) :: locst(1:*)
REAL (kind=8),INTENT(IN OUT) :: glost(1:*)
INTEGER (kind=4), OPTIONAL :: elem
END SUBROUTINE ensmat
