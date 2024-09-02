 SUBROUTINE elemt2(task ,ttime,iload,iforce,                       &
                   coord,sumat,mass,emass,resid,gstif,             &
                   istop,flag1,flag2,iwrit)

 IMPLICIT NONE
 ! Dummy arguments
 CHARACTER(len=*) task
 INTEGER (kind=4), INTENT(IN), OPTIONAL  :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL  :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL  :: ttime
 REAL (kind=8), INTENT(OUT), OPTIONAL  :: sumat
 REAL (kind=8), INTENT(IN), OPTIONAL  :: coord(:,:)
 REAL (kind=8), INTENT(IN OUT), OPTIONAL  :: mass(:),emass(:,:),resid(:,:),gstif(:)
 LOGICAL, OPTIONAL  :: flag1,flag2

 END SUBROUTINE elemt2
