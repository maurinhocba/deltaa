SUBROUTINE pcggib06( )
!  Write mesh information for GiD for 3-D-Shells (SIMO theory)
USE data_db,ONLY: shl3d, shl3d_nvarg, shl3d_head, shl3d_sets
IMPLICIT NONE

  TYPE(shl3d),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss =  shl3d_nvarg > 0
  IF (.NOT.gauss) RETURN

  e => shl3d_head                 !point to first element set
  DO iset=1,shl3d_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    IF( e%nnode == 6 .OR. e%nnode == 3)THEN
      etype = 3
    ELSE
      etype = 4
    END IF
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),1,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib06
