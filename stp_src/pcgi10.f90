 SUBROUTINE pcgi10 ( )
 !
 !  Write mesh information for GiD for Rigid Solids
 !  No material is written for rigid solids
 !  An option may be to write only if MATNO /= 0
 !
 USE data_db
 IMPLICIT NONE
 TYPE( rigid ), POINTER :: e
 INTEGER :: iset,iel,n,ndim
 CHARACTER (len=32) :: sname

 e => rigid_head                 !point to first element set

 DO iset=1,rigid_sets            !for each element set
   ndim = ndime                  !element dimension
   IF( e%ntype == 1 ) ndim = 0        !for points
   IF( e%ntype == 2 ) ndim = ndim - 1 !for boundaries (contact surfaces)
   n = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:n)//'"'
   CALL prtcoo(e%nnode,ndim,sname,rigid_nodes(1,2),iset)        !write coordinates and headers

   DO iel = 1,e%nelem            !for each element in the set
     WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode) !,e%matno  or e%set
   END DO
   IF( e%nmast > 0 )THEN
     n = LEN_TRIM(e%sname)
     sname = '"MN_'//e%sname(1:n)//'"'
     CALL prtcoo(1,0,sname,rigid_nodes(1,2),10000+iset)  !write coordinates and headers
     WRITE(11,"(10i8)") 1,label(e%nmast) !,e%matno  or     e%set
   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi10
