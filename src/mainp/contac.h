 SUBROUTINE contac(itask,dtime,ttime,toutd,iwrit,neq,  &
                   veloc,force,stiff,ustif,sname,maxve)

 !     main contac routine

 !USE solv_db, ONLY : nsymm
 IMPLICIT NONE

 !        Dummy arguments
 CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
 INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,maxve,neq
 REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,toutd,veloc(:)
 REAL (kind=8),INTENT(OUT), OPTIONAL :: force(:),stiff(:),ustif(:)
 CHARACTER(len=6), OPTIONAL, INTENT(IN) :: sname
 END SUBROUTINE contac
