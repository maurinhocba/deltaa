SUBROUTINE iteinv(a,b,m,x,lcr,rtol,n,iter,iw,ji)
!
!     Computes the least (absolute value) eigenvalue of a generalized
!     problem           A x = lcr * B * x
!
!     input variables:
!        a (:)  : matrix (n*n) expresed in compacted form (factorized)
!        b (:)  : matrix (n*n) expresed in compacted form
!        m (n+1)   : array indicating the position of diagonal elements
!                    in "a" and "b"
!        x (n)     : iegenvector prediction (if ji = 1)
!        n         : system order
!        iter      : maximum number of iterations allowed
!        rtol      : tolerance in convergence for eigenvalues
!        iw        : standard output unit
!        ji        : flag
!                 ji = 0 initializes eigenvector with ones
!                 ji = 1 x includes prediction for eigenvector
!
!     output
!        lcr       : least (absolute value) eigenvalue
!        x         : associated eigenvector
!
!     Reference  :
!             Bathe ,K.J. and Wilson , E.L. ," Numerical methods in
!             finite element analysis ", Prentice Hall , 1976, pag 420
!
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: n,iter,iw,m(n+1),ji
REAL (kind=8),INTENT(IN) :: b(:),rtol
REAL (kind=8),INTENT(IN OUT) :: a(:),lcr,x(n)

END SUBROUTINE iteinv
