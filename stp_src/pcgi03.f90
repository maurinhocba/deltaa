 SUBROUTINE pcgi03 ( )
 !
 !  Write mesh information for GiD for 2-D-Solids
 !
 USE data_db
 IMPLICIT NONE
 TYPE( sol2d ), POINTER :: e
 INTEGER :: iset,iel,n,k
 LOGICAL :: gauss
 CHARACTER (len=13) :: etype
 CHARACTER (len=32) :: sname
 CHARACTER (len=34) :: gpname


 gauss =  sol2d_nvarg > 0

 e => sol2d_head                 !point to first element set

 DO iset=1,sol2d_sets            !for each element set
   k = LEN_TRIM(e%sname)
   sname = '"'//e%sname(1:k)//'"'
   gpname = '"GP'//sname(2:32)
   k = k + 4
   CALL prtcoo(e%nnode,2,sname,sol2d_nodes(1,2),iset)        !write coordinates and headers

   DO iel = 1,e%nelem            !for each element in the set
     WRITE(11,"(10i8)") iel,(label(e%lnods(n,iel)),n=1,e%nnode),e%matno(iel) !e%set
   END DO
   IF( gauss )THEN
     IF( e%nnode == 3 .OR. e%nnode == 6 )THEN
       etype = 'Triangle     '
     ELSE
       etype = 'Quadrilateral'
     END IF
     WRITE(13,"('GaussPoints  ',a,' Elemtype  ',a13,3x,a)") &
                gpname(1:k),etype,sname
     WRITE(13,"('  Number Of Gauss Points:',i3)") e%ngaus
     IF( e%nnode  == 6 .AND. e%ngaus == 3 .AND. e%ver == 1)THEN
        WRITE(13,"('  Natural Coordinates: Given')")
        WRITE(13,"(2e20.10)") 1d0/6D0, 1D0/6D0, 2d0/3D0, 1d0/6D0, 1d0/6D0, 2d0/3D0
     ELSE
        !WRITE(13,"('  Natural Coordinates: Given')")
        !WRITE(13,"(2e20.10)") 0.5, 0.0, 0.5 , 0.5, 0.0, 0.5
        !WRITE(13,"(2e20.10)")  0.5 , 0.5, 0.0, 0.5, 0.5, 0.0 
        !WRITE(13,"(2e20.10)") 0.0, 0.5, 0.5, 0.0, 0.5 , 0.5 
        WRITE(13,"('  Natural Coordinates: internal')")
     END IF
     WRITE(13,"('End GaussPoints')")
   END IF
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcgi03
