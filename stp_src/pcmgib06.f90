SUBROUTINE pcmgib06( )
!  Write mesh information for GiD for 3-D-Shells (SIMO theory)
USE data_db,ONLY: shl3d, shl3d_head, shl3d_sets, shl3d_nodes, label
IMPLICIT NONE

  TYPE(shl3d),POINTER:: e
  INTEGER(kind=4):: iset,iel
  INTEGER(kind=4),ALLOCATABLE:: idmat(:)
  CHARACTER (len=32) :: sname

  e => shl3d_head                 !point to first element set
  DO iset=1,shl3d_sets            !for each element set
    sname = TRIM(e%sname)
    CALL prtcob(e%nnode,2,sname,shl3d_nodes(1,2),iset)        !write coordinates and headers

    CALL GID_BEGINELEMENTS()
    ALLOCATE(idmat(e%nnode+1))
    DO iel = 1,e%nelem            !for each element in the set
      idmat(1:e%nnode) = label(e%lnods(1:e%nnode,iel))
      idmat(e%nnode+1) = e%matno(iel)
      CALL GID_WRITEELEMENTMAT(iel,idmat)
    END DO
    DEALLOCATE(idmat)
    CALL GID_ENDELEMENTS()

    e => e%next
  END DO

RETURN
END SUBROUTINE pcmgib06
