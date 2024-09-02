SUBROUTINE pcmgib12( )
!  Write mesh information for GiD for 3-D-Shell (BST)
USE data_db,ONLY: bst, bst_head, bst_sets, bst_nodes, label
IMPLICIT NONE

  TYPE(bst),POINTER:: e
  INTEGER(kind=4):: iset,iel,nn
  INTEGER(kind=4),ALLOCATABLE:: idmat(:)
  CHARACTER(len=32):: sname

  e => bst_head                 !point to first element set
  DO iset=1,bst_sets            !for each element set
    nn = e%nnode
    sname = TRIM(e%sname)
    CALL prtcob(nn,2,sname,bst_nodes(1,2),iset)        !write coordinates and headers

    CALL GID_BEGINELEMENTS()
    ALLOCATE(idmat(nn+1))
    DO iel = 1,e%nelem            !for each element in the set
      idmat(1:nn) = label(e%lnods(1:nn,iel))
      idmat(nn+1) = e%matno(iel)
      CALL GID_WRITEELEMENTMAT(iel,idmat)
    END DO
    DEALLOCATE(idmat)
    CALL GID_ENDELEMENTS()

    e => e%next
  END DO

RETURN
END SUBROUTINE pcmgib12
