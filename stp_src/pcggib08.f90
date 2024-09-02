SUBROUTINE pcggib08( )
!  Write mesh information for GiD for BEAME elements (Simo Theory)
USE data_db,ONLY: beame, beame_nvarg, beame_head, beame_sets
IMPLICIT NONE

  TYPE(beame),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss = beame_nvarg > 0

  e => beame_head                 !point to first element set
  DO iset=1,beame_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    etype = 2
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),1,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
 END DO

RETURN
END SUBROUTINE pcggib08
