SUBROUTINE resvpl (istop,ttime,gvect,resid,flag)
!********************************************************************
!
!***   evaluation of integral (b)**t*(sigma)
!
!********************************************************************
!USE npoi_db, ONLY : ndime,ndofn,npoin,coora,coor1,locsy,locs1
!USE kin2_db, ONLY : ifpre,id
!USE cont_db, ONLY : ncont,fcont
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN OUT) :: istop
REAL (kind=8),INTENT(IN) :: ttime
REAL (kind=8),INTENT(IN OUT) :: gvect(:,:),resid(:)
LOGICAL, INTENT(IN), OPTIONAL :: flag

END SUBROUTINE resvpl
