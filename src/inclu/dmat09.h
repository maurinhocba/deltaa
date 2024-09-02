 SUBROUTINE dmat09( nucom,rr,volfr,newmt, cm,prop, gausv,d)

 IMPLICIT NONE

 INTEGER(kind=4), INTENT(IN) :: nucom
 REAL (kind=8), POINTER ::  rr(:,:),volfr(:)
 LOGICAL, INTENT(IN OUT) :: newmt
 REAL (kind=8), POINTER, OPTIONAL:: cm(:,:),prop(:,:)
 REAL (kind=8), OPTIONAL, INTENT(IN) :: gausv(:)
 REAL (kind=8), INTENT(OUT), OPTIONAL :: d(:,:)
 END SUBROUTINE dmat09
