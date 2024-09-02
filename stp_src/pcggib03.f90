SUBROUTINE pcggib03( )
!  Write mesh information for GiD for 2-D-Solids
USE data_db,ONLY: sol2d, sol2d_nvarg, sol2d_head, sol2d_sets
IMPLICIT NONE

  TYPE(sol2d),POINTER:: e
  INTEGER(kind=4):: etype,iset
  LOGICAL :: gauss
  CHARACTER(len=32):: sname
  CHARACTER (len=34) :: gpname

  gauss =  sol2d_nvarg > 0
  IF (.NOT.gauss) RETURN

  e => sol2d_head                 !point to first element set
  DO iset=1,sol2d_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)

    IF( e%nnode == 3 )THEN
      etype = 3
    ELSE
      etype = 4
    END IF
    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),e%ngaus,0,1)
    CALL GID_ENDGAUSSPOINT()

    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib03
