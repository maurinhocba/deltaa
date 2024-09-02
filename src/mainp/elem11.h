 SUBROUTINE elem11(task ,ttime,iload,iforce,                &
                   coora,sumat,mass,emass,resid,gstif,       &
                   istop,flag1,flag2,iwrit)

 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: task
 INTEGER (kind=4), INTENT(IN), OPTIONAL :: iwrit,iload,iforce
 INTEGER (kind=4), INTENT(IN OUT), OPTIONAL :: istop
 REAL (kind=8), INTENT(IN), OPTIONAL :: ttime,coora(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL :: sumat
 REAL (kind=8), INTENT(IN OUT), OPTIONAL :: mass(:),resid(:,:),gstif(:),     &
                                  emass(:,:)
 LOGICAL, OPTIONAL ::  flag1,flag2
 END SUBROUTINE elem11
