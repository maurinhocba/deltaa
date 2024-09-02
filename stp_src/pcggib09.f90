SUBROUTINE pcggib09( )
!  Write mesh information for GiD for SHREV elements (FGF)
USE data_db,ONLY: shrev, shrev_nvarg, shrev_head, shrev_sets
IMPLICIT NONE

  TYPE(shrev),POINTER:: e
  INTEGER(kind=4):: etype,iset,ngaus
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss =  shrev_nvarg > 0
  IF (.NOT.gauss) RETURN

  e => shrev_head                 !point to first element set
  DO iset=1,shrev_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    etype = 2
    ngaus = e%ngaus
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),ngaus,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib09
