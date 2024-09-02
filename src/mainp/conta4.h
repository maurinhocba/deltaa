SUBROUTINE conta4(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif, &
                  sname,maxve)

!     main contac routine (ALGORITHM 4)

!USE lispa0
!USE npoi_db, ONLY : ndofn, npoin, label(:), oldlb(:), emass(:,:), &
!                    coora(:,:), coord(:,:)
!USE kin1_db
!USE kin2_db
!USE cont_db ! iwric, bottom, top, coorb(:,:), coort(:,:) , fcont(:,:)
!USE cont4_db

IMPLICIT NONE
   !        Dummy arguments
CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
CHARACTER(len=6),INTENT(IN), OPTIONAL :: sname  !surface name
INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,nsymm,maxve
REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,veloc(:)
REAL (kind=8), INTENT(IN OUT), OPTIONAL :: force(:),stiff(:),ustif(:)

END SUBROUTINE conta4
