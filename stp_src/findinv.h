    SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
            IMPLICIT NONE
            !Declarations
            INTEGER, INTENT(IN) :: n
            INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
            REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
            REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
    END SUBROUTINE FINDInv

