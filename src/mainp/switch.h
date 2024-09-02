SUBROUTINE switch(istop,actio)
!********************************************************************
!
! *** Branch switching iterative routine
!
!********************************************************************
USE c_input
USE curv_db
USE ctrl_db, ONLY : ndime,neulr,ndofn,npoin,neq,maxa,nload
USE npo_db, ONLY :coord,coora,coorc,coor1,euler,locsy,locs1,ifpre,force,loadv,loass
USE solv_db
USE kinc_db, ONLY : ndepd,naris,maxav,nvelr,velor,nvfix
IMPLICIT NONE
INTEGER (kind=4), INTENT(IN OUT) :: istop
CHARACTER(len=*),INTENT(IN):: actio
END SUBROUTINE switch
