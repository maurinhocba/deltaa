SUBROUTINE mulmat(b,m,y,x,n)

!     subrutina que realiza la multiplicacion matricial de una matriz
!     cuadrada con una matriz columna   y = b * x   donde b es una
!     matriz simetrica escrita como vector en los elementos bajo el
!     skyline
!         variables de entrada
!             b (:)          : matriz n*n expresada como vector
!             y (:)          : vector columna
!             x (:)          : vector columna
!             m (:)          : arreglo que guarda las posiciones de
!                              los elementos diagonales en b
!             n              : orden de las matriz

IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: n,m(:)
REAL (kind=8),INTENT(IN) :: b(:),x(:)
REAL (kind=8),INTENT(OUT):: y(:)

INTEGER (kind=4) :: i,k,k1,j,j1

DO i = 1,n
  k  = m(i)         !  k  = posicion del elemento diagonal
  k1 = m(i+1) - 1   !  k1 = posicion del 1er elem. no-nulo de la fila
  y(i) = b(k)*x(i)
  j1 = i
  DO j =k+1, k1
    j1    = j1-1
    y(i)  = y(i)  + b(j)*x(j1)
    y(j1) = y(j1) + b(j)*x(i)
  END DO
END DO
RETURN

END SUBROUTINE mulmat
