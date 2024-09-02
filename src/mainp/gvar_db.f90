MODULE gvar_db
IMPLICIT NONE

  LOGICAL :: actchk=.FALSE.       !Flag for check data input file
  !File identification to import element sets
  INTEGER(kind=4) :: fimpo=70, lab1=0
  LOGICAL :: renum=.FALSE., seque=.FALSE., inter=.TRUE., overw=.FALSE.

  INTEGER (kind=4), POINTER :: nodset(:)       ! (numnp), nodes of the eset
  INTEGER (kind=4) :: numpo                    ! number of points in old mesh set
  LOGICAL          :: ksnxt=.FALSE.            ! Compute Stiffnes in next iteration

END MODULE gvar_db
