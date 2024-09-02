SUBROUTINE pcggib01( )
!  Write mesh information for GiD for SPOT elements
USE data_db,ONLY: spot, spot_nvarg, spot_head, spot_sets
IMPLICIT NONE

  TYPE(spot),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss = spot_nvarg > 0         !Gauss variables to be print
  IF (.NOT.gauss) RETURN

  e => spot_head                 !point to first element set
  DO iset=1,spot_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    etype = 2
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),1,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib01
