 SUBROUTINE coract(ndime,npoin,neulr,ndofn,lambd,coord,ddisp,coora,euler)

 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: ndime,npoin,ndofn,neulr
 REAL (kind=8),INTENT(IN) :: ddisp(:),lambd,coord(:,:)
 REAL (kind=8),INTENT(IN OUT) :: coora(:,:)
 REAL (kind=8), POINTER :: euler(:,:)

 END SUBROUTINE coract
