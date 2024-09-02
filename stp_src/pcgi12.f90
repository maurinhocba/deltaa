 SUBROUTINE pcgi12 ( )
 !
 !  Write mesh information for GiD for 3-D-Shell (BST)
 !
 USE data_db
 IMPLICIT NONE
 TYPE( bst ), POINTER :: e
 INTEGER :: iset,iel,n,k
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname

 gauss =  bst_nvarg > 0

 e => bst_head                 !point to first element set

 DO iset=1,bst_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo(e%nnode,2,sname,bst_nodes(1,2),iset)        !write coordinates and headers

   DO iel = 1,e%nelem            !for each element in the set
     WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno(iel) !e%set
   END DO
   IF( gauss )THEN
     IF( e%nnode == 3 ) etype = 'Triangle     '
     IF( e%nnode == 4 ) etype = 'Quadrilateral'
     WRITE(13,"('GaussPoints ',a,' Elemtype ',a13,x,a)")  &
                gpname(1:k),etype,sname
     WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
             &  '  Natural Coordinates: internal',/,     &
             &  'End GaussPoints')")e%ngaus
   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi12
