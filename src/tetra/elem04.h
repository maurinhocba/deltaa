 SUBROUTINE elem04(task ,ttime,iload,iforce,                       &
                   coora,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 !USE load_db, ONLY : igrav,gravy,gv,loadv
 !USE npoi_db, ONLY : npoin,label,oldlb,coord
 !USE ele04_db
 IMPLICIT NONE
 CHARACTER(len=6), INTENT(IN) :: task

 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 END SUBROUTINE elem04
