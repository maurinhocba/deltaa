SUBROUTINE dynast(istop,actio)
!********************************************************************
!
! *** iterative dynamic process routine
!
!********************************************************************
!USE lispa0
!USE curv_db
!USE npoi_db, ONLY : ndime,neulr,ndofn,npoin, &
!                    coord,coora,coorc,euler,locsy
!USE solv_db
!USE dyna_db
!USE kin0_db, ONLY : ndepd,naris
!USE kin1_db, ONLY : neq,maxa,maxav
!USE kin2_db, ONLY : nvelr,ifpre,velor,nvfix
!USE load_db, ONLY : nload,force,loadv,loass
IMPLICIT NONE
INTEGER (kind=4), INTENT(IN OUT) :: istop
CHARACTER(len=*),INTENT(IN):: actio
END SUBROUTINE dynast
