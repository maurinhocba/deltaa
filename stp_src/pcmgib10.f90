SUBROUTINE pcmgib10( )
!  Write mesh information for GiD for Rigid Solids
USE data_db,ONLY: ndime, rigid, rigid_head, rigid_sets, rigid_nodes, label
IMPLICIT NONE

  TYPE(rigid),POINTER:: e
  INTEGER(kind=4):: iset,iel,ndim,n
  INTEGER(kind=4),ALLOCATABLE:: idmat(:)
  CHARACTER(len=32):: sname

  e => rigid_head                 !point to first element set
  DO iset=1,rigid_sets            !for each element set
    ndim = ndime                  !element dimension
    IF( e%ntype == 2 ) ndim = ndim - 1 !for boundaries (contact surfaces)
    IF( e%ntype == 1 ) ndim = 0 !for points
    sname = TRIM(e%sname)
    CALL prtcob(e%nnode,ndim,sname,rigid_nodes(1,2),iset)        !write coordinates and headers

    CALL GID_BEGINELEMENTS()
    ALLOCATE(idmat(e%nnode+1))
    DO iel = 1,e%nelem            !for each element in the set
      idmat(1:e%nnode) = label(e%lnods(1:e%nnode,iel))
      idmat(e%nnode+1) = e%matno
      CALL GID_WRITEELEMENTMAT(iel,idmat)
    END DO
    DEALLOCATE(idmat)
    CALL GID_ENDELEMENTS()
    IF( e%nmast > 0 )THEN
      n = LEN_TRIM(e%sname)
      sname = '"MN_'//e%sname(1:n)//'"'
      CALL prtcob(1,0,sname,rigid_nodes(1,2),10000+iset)        !write coordinates and headers

      CALL GID_BEGINELEMENTS()
      ALLOCATE(idmat(2))
      idmat(1) = label(e%nmast)
      idmat(2) = e%matno
      CALL GID_WRITEELEMENTMAT(1,idmat)
      DEALLOCATE(idmat)
      CALL GID_ENDELEMENTS()
    END IF


    e => e%next
  END DO

RETURN
END SUBROUTINE pcmgib10
