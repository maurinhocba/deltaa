SUBROUTINE stiffm (istop,ttime,nullm,iforce,stiff,force,nsymm, &
                   ustif,flag)
!********************************************************************
!
!***   evaluation of tangent stiffness matrix
!
!********************************************************************
!USE npoi_db, ONLY : coora,coor1,locsy,locs1
!USE cont_db, ONLY : ncont
IMPLICIT NONE
LOGICAL, INTENT(IN) :: nullm
INTEGER (kind=4),INTENT(IN) :: iforce
INTEGER (kind=4),INTENT(IN OUT) :: istop,nsymm
REAL (kind=8),INTENT(IN) :: ttime
REAL (kind=8),INTENT(OUT) :: stiff(:),force(:,:),ustif(:)
LOGICAL, INTENT(IN), OPTIONAL :: flag

END SUBROUTINE stiffm
