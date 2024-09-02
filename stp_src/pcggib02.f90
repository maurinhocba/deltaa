SUBROUTINE pcggib02( )
!  Write mesh information for GiD for TRUSS elements
USE data_db,ONLY: truss, truss_nvarg, truss_head, truss_sets
IMPLICIT NONE

  TYPE(truss),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss = truss_nvarg > 0         !Gauss variables to be print
  IF (.NOT.gauss) RETURN

  e => truss_head                 !point to first element set
  DO iset=1,truss_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    etype = 2
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),1,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib02
