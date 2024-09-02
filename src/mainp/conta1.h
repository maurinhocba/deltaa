SUBROUTINE conta1(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif, &
               sname,maxve)
!     main contac routine (ALGORITHM 1) 2D problems

!USE npoi_db, ONLY : ndofn, npoin, label, oldlb, emass, coora, coord,
!                    coor1
!USE kin1_db
!USE kin2_db
!USE cont_db   ! bottom, top, fcont(:,:), coorb(:,:), coort(:,:)
!USE cont1_db

IMPLICIT NONE
   !        Dummy arguments
CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
CHARACTER(len=6),INTENT(IN), OPTIONAL :: sname ! surface name
INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,nsymm,maxve
REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,veloc(:)
REAL (kind=8), INTENT(IN OUT), OPTIONAL :: force(:),stiff(:),ustif(:)
END SUBROUTINE conta1
