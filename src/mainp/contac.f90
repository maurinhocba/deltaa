 SUBROUTINE contac(itask,dtime,ttime,toutd,iwrit,neq,  &
                   veloc,force,stiff,ustif,sname,maxve)

 !     main contac routine

 USE solv_db, ONLY : nsymm
 IMPLICIT NONE

 !        Dummy arguments
 CHARACTER(len=6),INTENT(IN) :: itask  !task to perform
 INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,maxve,neq
 REAL (kind=8),INTENT(IN), OPTIONAL :: dtime,ttime,toutd,veloc(:)
 REAL (kind=8),INTENT(OUT), OPTIONAL :: force(:),stiff(:),ustif(:)
 CHARACTER(len=6), OPTIONAL, INTENT(IN) :: sname

! INTERFACE
!   INCLUDE 'conta1.h'
!   INCLUDE 'conta4.h'
! END INTERFACE
!
! SELECT CASE (ncont)
!
! CASE (1)
!   CALL conta1(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif, &
!               sname,maxve)
!
! CASE (4)
!   CALL conta4(itask,ttime,iwrit,veloc,dtime,nsymm,force,stiff,ustif, &
!               sname,maxve)
!
! END SELECT

 RETURN
 END SUBROUTINE contac
