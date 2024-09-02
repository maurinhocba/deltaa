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

END SUBROUTINE mulmat
