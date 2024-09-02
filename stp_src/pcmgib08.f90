SUBROUTINE pcmgib08( )
!  Write mesh information for GiD for BEAME elements (Simo Theory)
USE data_db,ONLY: beame, beame_head, beame_sets, beame_nodes, label
IMPLICIT NONE

  TYPE(beame),POINTER:: e
  INTEGER(kind=4):: iset,iel
  INTEGER(kind=4),ALLOCATABLE:: idmat(:)
  CHARACTER(len=32):: sname

 e => beame_head                 !point to first element set
 DO iset=1,beame_sets            !for each element set
   sname = TRIM(e%sname)
   CALL prtcob(e%nnode,1,sname,beame_nodes(1,2),iset)        !write coordinates and headers


    CALL GID_BEGINELEMENTS()
    ALLOCATE(idmat(e%nnode+1))
    IF( e%nnode == 2)THEN
      DO iel = 1,e%nelem            !for each element in the set
        idmat(1:e%nnode) = label(e%lnods(1:e%nnode,iel))
        idmat(e%nnode+1) = e%matno
        CALL GID_WRITEELEMENTMAT(iel,idmat)
      END DO
    ELSE !e%nnode == 3
      DO iel = 1,e%nelem            !for each element in the set
        idmat(1) = label(e%lnods(1,iel))
        idmat(2) = label(e%lnods(3,iel))
        idmat(3) = label(e%lnods(2,iel))
        idmat(e%nnode+1) = e%matno
        CALL GID_WRITEELEMENTMAT(iel,idmat)
      END DO
    END IF
    DEALLOCATE(idmat)
    CALL GID_ENDELEMENTS()

   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcmgib08
