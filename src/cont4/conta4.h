SUBROUTINE conta4(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif)

!     main contac routine (ALGORITHM 4)

!USE npo_db ! botton, top, ndofn, npoin,
!           ! label(:,:), oldlb(:,:), emass(:,:), fcont(:,:),
!           ! coora(:,:), coord(:,:), coorb(:,:), coort(:,:)
!USE kin1_db
!USE kin2_db
!USE cont_db
!USE cont4_db

IMPLICIT NONE
   !        Dummy arguments
CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,nsymm
REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,veloc(:)
REAL (kind=8), INTENT(IN OUT), OPTIONAL :: force(:),stiff(:),ustif(:)
END SUBROUTINE conta4
