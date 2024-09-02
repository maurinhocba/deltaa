 SUBROUTINE stre02_1 (mat,e,stres,hyper,newm,j)

 ! computes the internal nodal forces  1D (truss elements)

 USE mat_dba, ONLY : mater
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: e
 REAL (kind=8), INTENT(IN OUT) :: stres(:)
 REAL (kind=8), INTENT(OUT) :: j
 TYPE (mater), POINTER :: mat
 LOGICAL, INTENT(IN) :: hyper
 LOGICAL, INTENT(IN OUT) :: newm
 END SUBROUTINE stre02_1
