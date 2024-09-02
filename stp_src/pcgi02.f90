 SUBROUTINE pcgi02 ( )
 !
 !  Write mesh information for GiD for TRUSS elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE( truss ), POINTER :: e
 INTEGER :: iset,iel,n,k
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname

 gauss = truss_nvarg > 0         !Gauss variables to be print

 e => truss_head                 !point to first element set

 DO iset=1,truss_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo(e%nnode,1,sname,truss_nodes(1,2),iset)        !write coordinates and headers
   DO iel = 1,e%nelem            !for each element in the set
     IF( e%nnode == 2)THEN
       WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno !e%set
     ELSE !e%nnode == 3
       WRITE(11,"(10i8)") iel,label(e%lnods(1,iel)),label(e%lnods(3,iel)), &
                              label(e%lnods(2,iel)),e%set !e%matno
     END IF
   END DO
   IF( gauss )THEN
     etype = 'Linear       '
     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname
     WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
             &  '  Natural Coordinates: internal',/,     &
             &  'End GaussPoints')")e%ngaus
   END IF
   e => e%next
 END DO


 RETURN

 END SUBROUTINE pcgi02
