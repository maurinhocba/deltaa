SUBROUTINE pcmgib01( )
!  Write mesh information for GiD for SPOT elements
USE data_db,ONLY: spot, spot_head, spot_sets, spot_nodes, label
IMPLICIT NONE

  TYPE(spot),POINTER:: e
  INTEGER(kind=4):: iset,iel
  INTEGER(kind=4),ALLOCATABLE:: idmat(:)
  CHARACTER(len=32):: sname

 INTERFACE
   SUBROUTINE prtcob_spot(nelem,nnode,iset,sname,x,fl,lnods,matno)
     INTEGER (kind=4), INTENT(IN) :: nelem,nnode,fl,iset,matno
     CHARACTER (len=32), INTENT(IN) :: sname
     INTEGER(kind=4), POINTER :: lnods(:,:)
     REAL (kind=8), POINTER :: x(:,:,:)
   END SUBROUTINE prtcob_spot
 END INTERFACE

  e => spot_head                 !point to first element set
  DO iset=1,spot_sets            !for each element set
    sname = TRIM(e%sname)
    CALL prtcob_spot(e%nelem,e%nnode,iset,sname,e%a_x,e%first_l,e%lnods,e%matno)        !write coordinates and headers

    e => e%next
  END DO

RETURN
END SUBROUTINE pcmgib01
