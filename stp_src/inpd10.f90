 SUBROUTINE inpd10( nel, sname, flag)
 !
 !  read data for a RIGID BODIES sets
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN OUT) :: nel
 LOGICAL, INTENT(IN) :: flag
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,iaux,n
 TYPE (rigid), POINTER :: eset

 ALLOCATE (eset)             !get memory for this set
 NULLIFY (eset%next)         !nullify pointer to next set

 IF( rigid_sets == 0 )THEN   !for the first 2-D SOLID set
   rigid_head => eset
 ELSE                        !for subsequent sets
   rigid_tail%next => eset
 END IF
 rigid_tail => eset          !last set position

 rigid_sets = rigid_sets + 1 !increase number of sets

 eset%sname = sname            !set name

 IF( flag )THEN  !boundary (contact surface)
   READ(18) eset%nnode         !number of nodes per element and
   eset%nmast = 0              !no master node
   eset%ntype = 2              !a surface
   eset%set = rigid_sets+100   !set position (possibly unnecessary)
 ELSE            !rigid body
   eset%set = nsets            !set position (possibly unnecessary)
   READ(17) eset%nnode,iaux,eset%nmast         !number of nodes per element and
   SELECT CASE ( iaux )
   CASE (0)
     eset%ntype = 1
     nel = eset%nnode
     eset%nnode = 1
   CASE (1,2,3)
     eset%ntype = 3
   CASE (4,5,6)
     eset%ntype = 2
   END SELECT
 END IF

 eset%nelem = nel              !number of elements in the set
 ALLOCATE ( eset%lnods(eset%nnode,nel) )   !Connectivities

 ! read connectivities

 IF( flag )THEN  !boundary (contact surface)
   READ(18) ((eset%lnods(n,ielem),n=1,eset%nnode),ielem=1,nel)
   eset%matno = 100+eset%set
 ELSE            !whole body (RIGID)
   READ(17) eset%matno,((eset%lnods(n,ielem),n=1,eset%nnode),ielem=1,nel)
 END IF

 RETURN
 END SUBROUTINE inpd10
