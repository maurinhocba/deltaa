 SUBROUTINE tr4to8(gausv,newgv,iv,nv,ng,mv)
 IMPLICIT NONE
 INTEGER , INTENT(IN) :: iv,ng,nv
 REAL(kind=8) :: gausv(:,:)
 REAL(kind=8) :: newgv(6,8)
 REAL(kind=8), OPTIONAL :: mv
 END SUBROUTINE tr4to8
