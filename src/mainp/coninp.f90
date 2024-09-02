      SUBROUTINE coninp(itask,ncont,dtcal,ttime,toutd,iwrit,  &
     &                  neq,veloc,nvfix,maxve,velor,sname)

!     main contac input routine

      IMPLICIT NONE

!        Dummy arguments
      CHARACTER(len=*),INTENT(IN):: itask  !task to perform
      INTEGER(kind=4),INTENT(IN):: ncont
      INTEGER(kind=4),INTENT(IN),OPTIONAL:: iwrit,neq,nvfix
      INTEGER(kind=4),INTENT(IN),OPTIONAL:: maxve
      REAL(kind=8),INTENT(IN),OPTIONAL:: dtcal,ttime,toutd,veloc(:),velor(:)
      CHARACTER(len=*),OPTIONAL,INTENT(IN):: sname
!        Local variables

      INTERFACE
        INCLUDE 'conta1.h'
        !INCLUDE 'conta2.h'
        !INCLUDE 'conta3.h'
        INCLUDE 'conta4.h'
        !INCLUDE 'conta5.h'
        !INCLUDE 'conta6.h'
        !INCLUDE 'dbea3d.h'
      END INTERFACE

      SELECT CASE (ncont)

      !CASE (1)
      !  CALL conta1(itask,dtcal,ttime,iwrit,neq,veloc,nvfix,maxve,velor,sname)
      !CASE (2)
      !  CALL conta2(itask,dtcal,ttime,toutd,sname)
      !CASE (3)
      !  CALL conta3(itask,dtcal,ttime,toutd,sname)
      CASE (4)
        CALL conta4(itask,dtcal,ttime,iwrit,neq,veloc,nvfix,maxve,velor,sname)
      !CASE (5)
      !  CALL conta5(itask,dtcal,sname,ttime)
      !CASE (6)
      !  CALL conta6(itask,dtcal,sname,ttime)
      !CASE (7) ! special contact algorithm - drawbeads
      !  CALL dbea3d(itask,dtcal,ttime,toutd)
      CASE DEFAULT
        CALL runend (' Wrong contact algorithm code')

      END SELECT

      RETURN
      END SUBROUTINE coninp
