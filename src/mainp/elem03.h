 SUBROUTINE elem03(task ,ttime,iload,iforce,                       &
                   coora,locsy,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 !USE ctrl_db, ONLY: ndime, npoin, top, bottom
 !USE ele03_db
 !USE npo_db, ONLY : ifpre, coord, euler, coort, coorb, ifact, loadv
 !USE loa_db, ONLY : gravy, gv, igrav
 IMPLICIT NONE
 CHARACTER(len=6), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime,coora(:,:),locsy(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL  ::  flag1,flag2

 END SUBROUTINE elem03
