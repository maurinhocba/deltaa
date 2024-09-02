 SUBROUTINE pcgi09 ( )
 !
 !  Write mesh information for GiD for SHREV elements (FGF)
 !
 USE data_db
 IMPLICIT NONE
 TYPE( shrev ), POINTER :: e
 INTEGER :: iset,iel,n,k,n1,l,i
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname

 gauss =  shrev_nvarg > 0

 e => shrev_head                 !point to first element set

 DO iset=1,shrev_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   IF( ntype /= 4 )THEN
     CALL prtcoo(e%nnode,1,sname,shrev_nodes(1,2),iset)        !write coordinates and headers

     DO iel = 1,e%nelem            !for each element in the set
       IF( e%nnode == 2)THEN
         WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno !e%set
       ELSE !e%nnode == 3
         WRITE(11,"(10i8)") iel,label(e%lnods(1,iel)),label(e%lnods(3,iel)), &
                                label(e%lnods(2,iel)),e%matno  !e%set
       END IF
     END DO
   ELSE
     WRITE(11,"('End Elements')")
     WRITE(11,"('MESH ',a,' dimension =   2',' ElemType     Linear   ', &
     &          ' Nnode =  2',/,'Coordinates')")sname(1:k)
     DO i=1,ngp
       WRITE(11,"(i8,3e18.10)")labea(i),coora(1:ndime,i)
     END DO
     WRITE(11,"('End coordinates')")
     WRITE(11,"('Elements')")
     n = 0
     DO iel = 1,e%nelem            !for each element in the set
       n1 = labea(e%lnods(3,iel))
       DO l=1,e%ngaus-1
         n = n+1
         WRITE(11,"(10i8)") n,n1,n1+1,e%matno+10 !e%set
         n1 = n1+1
       END DO
     END DO
   END IF
   IF( gauss )THEN
     etype = 'Linear       '
     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname
     IF( ntype /= 4 )THEN
       WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
               &  '  Nodes not included',/,     &
               &  '  Natural Coordinates: internal',/,     &
               &  'End GaussPoints')")e%ngaus
     ELSE
       WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
               &  '  Nodes included               ',/,     &
               &  '  Natural Coordinates: internal',/,     &
               &  'End GaussPoints')")2 !e%ngaus
     END IF
   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi09
