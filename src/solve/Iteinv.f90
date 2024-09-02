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

 INTERFACE
   INCLUDE 'colsol.h'
   INCLUDE 'mulmat.h'
 END INTERFACE

 INTEGER (kind=4) :: l
 REAL (kind=8) :: c,d,e,aut,aut1
 REAL (kind=8),ALLOCATABLE :: y(:)

 ALLOCATE( y(n) )
 IF (ji == 0) x = 1d0   ! initializes eigenvector x
 aut = 1d0              ! initializes eingevalue
 CALL mulmat (b,m,y,x,n)                  !y = B*x
 l = 0                                    !initializes counter
 e = 2d0*rtol                             !initializes error
 DO WHILE(l < iter .AND. rtol-e < 0)
    l = l+1
    x = y
    CALL colsol(a,m,n,2,iw,1,0,v=x)
    c = DOT_PRODUCT(y,x)
    CALL mulmat(b,m,y,x,n)
    d = DOT_PRODUCT(y,x)
    aut1 = c/d
    d = 1d0/SQRT(ABS(d))
    y = y*d
    e = ABS((aut1-aut)/aut)
    aut = aut1
    WRITE(6,8500) e,aut,aut+lcr,l
 END DO
 IF(l == iter)WRITE(iw,8000)iter,rtol,e,aut+lcr
 d = MAXVAL(x)            !positive value
 c = MINVAL(x)            !negative value
 IF( ABS(c) > d )  d = c
 x = x/d                  !normalizes vector x
 lcr = lcr + aut

 DEALLOCATE( y )
 RETURN

 8000 FORMAT(' no convergence for the maximum number of iterations ',  &
      &'allowed'/'  number of iterations ',i5/                         &
      &'  tolerance in convergence.',e15.5/'  last convergence factor',&
      &   e15.5,/,'  last computed eigenvalue ',e15.5/)
 8500 FORMAT('+',t12,'  E  = ',e12.5,'  Eig= ',e12.5,' Eig+L0 = ',e12.5,' it=',i4)

 END SUBROUTINE iteinv
