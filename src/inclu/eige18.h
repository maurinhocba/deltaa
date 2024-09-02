 SUBROUTINE eige18(str,r,lb,flag)

 ! compute eigenvalues and eigenvector of symmetric tensor STR

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT (IN) :: str(6)     !Symmetric tensor
 REAL(kind=8), INTENT (OUT) :: r(3,3), & !components of eigevectors
                               lb(3)     !eigenvalues of STR
 INTEGER (kind=4), INTENT(OUT) :: flag   !flag

 END SUBROUTINE eige18
