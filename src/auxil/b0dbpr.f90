SUBROUTINE b0dbpr(b0,b,d,e,ncol,nrow)
!**********************************************************************
!                              T
!****THIS ROUTINE PERFORMS E=B0 D B WITH B A GENERAL MATRIX
!    AND DMATX SYMMETRIC IS DEFINED AS AN ARRAY AND THE RESULTING
!    EMATX IS GIVEN AS AN ARRAY TOO.
!               T                                                T
!             B0 = NCOL * NROW     NOTICE THAT THE ARGUMENT IS B0
!             B  = NCOL * NROW     NOTICE THAT THE ARGUMENT IS B-T
!             D  = NROW * NROW     the full matrix is pased
!             E  = NCOL * NCOL     the full matriz IS CONSTRUCTED
!             A  = NROW            AUXILIAR VECTOR
!
!**********************************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nrow, ncol
REAL (kind=8),INTENT(IN) :: b0(ncol,nrow),b(ncol,nrow),d(nrow,nrow)
REAL (kind=8),INTENT(IN OUT) :: e(ncol,ncol)

INTEGER (kind=4) :: i, j
REAL (kind=8) :: a(nrow)

!                       MATRIX E IS NOT MADE ZERO INITIALLY

DO i=1,ncol                                           !row i of e
  !                         A =  (row I of B0(T) * D
  a = MATMUL(b0(i,:),d)
  !             E = (B0*D) * B
  DO j=1,ncol                                         !column J of E
    e(i,j) = e(i,j)+DOT_PRODUCT(a,b(j,:))
  END DO                                              ! J=1,NCOL

END DO                                                ! I=1,NCOL
RETURN

END SUBROUTINE b0dbpr
