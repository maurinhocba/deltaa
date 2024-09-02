      SUBROUTINE coninp(itask,ncont,dtcal,ttime,toutd,iwrit,  &
     &                  neq,veloc,nvfix,maxve,velor,sname)

!     main contac input routine

      IMPLICIT NONE

!        Dummy arguments
      CHARACTER(len=*),INTENT(IN) :: itask  !task to perform
      INTEGER (kind=4),INTENT(IN) :: ncont
      INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,neq,nvfix
      INTEGER (kind=4),INTENT(IN), OPTIONAL :: maxve
      REAL (kind=8),INTENT(IN), OPTIONAL :: dtcal,ttime,toutd,veloc(:), &
     &                                      velor(:)
      CHARACTER(len=*), OPTIONAL, INTENT(IN) :: sname
      END SUBROUTINE coninp
