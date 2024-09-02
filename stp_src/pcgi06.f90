 SUBROUTINE pcgi06 ( )
 !
 !  Write mesh information for GiD for 3-D-Shells (SIMO theory)
 !
 USE data_db
 IMPLICIT NONE
 TYPE( shl3d ), POINTER :: e
 INTEGER :: iset,iel,n,k
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname

 gauss =  shl3d_nvarg > 0

 e => shl3d_head                 !point to first element set

 DO iset=1,shl3d_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo(e%nnode,2,sname,shl3d_nodes(1,2),iset)        !write coordinates and headers

   DO iel = 1,e%nelem            !for each element in the set
     WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno(iel) !e%set
   END DO
   IF( gauss )THEN
     IF( e%nnode == 6 .OR.  e%nnode == 3 )THEN
       etype = 'Triangle     '
     ELSE
       etype = 'Quadrilateral'
     END IF
     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname
     IF( e%ngaus > 0 )THEN
       WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
               &  '  Natural Coordinates: internal',/,     &
               &  'End GaussPoints')")e%ngaus
     ELSE
       WRITE(13,"('  Number Of Gauss Points:',i3,/,        &
               &  '  Natural Coordinates: Given')")ABS(e%ngaus)
       WRITE(13,"(3e20.12)") 1d0/6d0,1d0/6d0
       WRITE(13,"(3e20.12)") 2d0/3d0,1d0/6d0
       WRITE(13,"(3e20.12)") 1d0/6d0,2d0/3d0
       WRITE(13,"('End GaussPoints')")
     END IF
   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi06
