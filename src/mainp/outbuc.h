SUBROUTINE outbuc(npoin,ndofn,istep,ifpre,x,lambd,buckl,flag)
!     normalizes eigevector and print into output file

IMPLICIT NONE
INTEGER(kind=4),INTENT(IN)::npoin,ndofn,istep,ifpre(:,:),flag
REAL (kind=8),INTENT(IN) :: lambd,buckl
REAL (kind=8),INTENT(IN OUT) :: x(:)
END SUBROUTINE outbuc
