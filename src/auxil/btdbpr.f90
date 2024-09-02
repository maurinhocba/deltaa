SUBROUTINE btdbpr(b,d,e,ncol,nrow)
!**********************************************************************
!                             T
!****THIS ROUTINE PERFORMS E=B D B WITH B A GENERAL MATRIX
!    AND DMATX SYMMETRIC IS DEFINED AS AN ARRAY AND THE RESULTING
!    EMATX IS GIVEN AS AN ARRAY TOO.
!              T                                                T
!             B  = NCOL * NROW     NOTICE THAT THE ARGUMENT IS B
!             D  = NROW * NROW     ONLY THE UPPER TRIANGLE IS PASSED
!             E  = NCOL * NCOL     ONLY THE UPPER TRIANGLE IS CONSTRUCTED
!             A  = NROW            AUXILIAR VECTOR
!
!**********************************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nrow, ncol
REAL (kind=8),INTENT(IN) :: b(ncol,nrow), d(*)
REAL (kind=8),INTENT(IN OUT) :: e(*)

INTEGER (kind=4) :: i, j, k, l, i_j, k_l
REAL (kind=8) :: a(nrow)

!                       MATRIX E IS NOT MADE ZERO INITIALLY
i_j = 0
DO i=1,ncol                                           !row i of e

  !                         A =  (column I of B(T) * D

  !                         LOWER TRIANGLE AND DIAGONAL OF D MATRIX
  k_l = 0
  DO l=1,nrow                                         !column L of D
    a(l) = 0
    DO k=l,nrow                                       !row K of D
      k_l = k_l + 1
      a(l) = a(l) +  b(i,k)*d(k_l)
    END DO                                            ! k = l,nrow
  END DO                                              ! l = 1,nrow
  !                         UPPER TRIANGLE OF D MATRIX
  k_l = 0
  DO k=1,nrow-1                                       !row K of D
    k_l = k_l + 1                                     !skip diagonal term
    DO l=k+1,nrow                                     !column L of D
      k_l = k_l + 1
      a(l) = a(l) +  b(i,k) * d(k_l)
    END DO                                            !K= L+1,NROW
  END DO                                              ! L= 1,NROW-1
  !             E = tp(B) * A * B
  DO j=i,ncol                                         !column J of E
    i_j = i_j+1
    DO l=1,nrow
      e(i_j) = e(i_j)+a(l)*b(j,l)
    END DO
  END DO                                              ! J=1,NCOL

END DO                                                ! I=1,NCOL
RETURN

END SUBROUTINE btdbpr
