 SUBROUTINE corrb2(st11,st22,efpst,c,prop,b,d,ierr,dstpl)
 !-------------------------------------------------------------------
 !
 !     Planar and transversal Anisotropy
 !
 !-------------------------------------------------------------------
 USE lispa0
 IMPLICIT NONE

 REAL (kind=8),INTENT(IN) :: prop(:),   & !material properties
                             c(:),      & !elasticity matrix (orthotropic)
                             b(:),      & !flow rule matrix
                             d(:)         !yield function derivative
 REAL (kind=8),INTENT(IN OUT) :: st11,st22   !trial and corrected stresses
 REAL (kind=8),INTENT(IN OUT) :: efpst    !effective plastic strain
                                          !(IN) present   (OUT) increment
 REAL (kind=8),INTENT(OUT) :: dstpl(2)    !increment in plastic strains
 INTEGER (kind=4), INTENT(OUT) :: ierr    !error flag (0: O.K., 1:error)
 END SUBROUTINE corrb2
