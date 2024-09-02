SUBROUTINE pcggib12( )
!  Write mesh information for GiD for 3-D-Shell (BST)
USE data_db,ONLY: bst, bst_nvarg, bst_head, bst_sets
IMPLICIT NONE

  TYPE(bst),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname

  gauss =  bst_nvarg > 0
  IF (.NOT.gauss) RETURN

  e => bst_head                 !point to first element set
  DO iset=1,bst_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)
    etype = e%nnode
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),1,0,1)
    CALL GID_ENDGAUSSPOINT()
    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib12
