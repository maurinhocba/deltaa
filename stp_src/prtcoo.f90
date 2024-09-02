 SUBROUTINE prtcoo(nnode,dim,snam,nodes,iset)

 ! print coordinates for GiD

 USE data_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nnode,dim,nodes(npoin),iset
 CHARACTER (len=32), INTENT(IN) :: snam

 CHARACTER (len=13) elmt
 INTEGER (kind=4) i,l
 CHARACTER (len=32) :: sname

 sname = snam
 SELECT CASE (dim)           !element dimension
 CASE (0)
   IF( nnode == 1 ) elmt = 'Point        '
 CASE (1)
   IF( nnode == 2 .OR. nnode == 3 ) elmt = 'Linear       '
 CASE (2)
   IF( nnode == 2 ) elmt = 'Linear       '
   IF( nnode == 3 .OR. nnode == 6) elmt = 'Triangle     '
   IF( nnode == 4 .OR. nnode == 8) elmt = 'Quadrilateral'
 CASE (3)
   IF( nnode == 4) elmt = 'Tetrahedra   '
   IF( nnode == 6 .OR. nnode == 15 .OR. nnode == 18) elmt = 'Prism        '
   IF( nnode == 8 .OR. nnode == 20) elmt = 'Hexahedra    '
 END SELECT

 IF(first)THEN
   first = .FALSE.
 ELSE
   WRITE(11,"('End Elements')")
 END IF
 l = LEN_TRIM(sname)
 WRITE(11,"('MESH ',a,' dimension =',i2,' ElemType ',a13, &
&          ' Nnode = ',i2,/,'Coordinates')")sname(1:l),ndime,elmt,nnode
 IF( .NOT.wsdisp .AND. wtdisp )THEN  !write original coordinates
   DO i=1,npoin
     IF ( nodes(i) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coord(1:ndime,i)*tdisp_f
     IF ( nodes(i) > 0 ) meshn(i) = .TRUE.
   END DO
 ELSE !Write stage coordinates
   DO i=1,npoin
     IF ( nodes(i) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coors(1:ndime,i)*sdisp_f
     IF ( nodes(i) > 0 ) meshn(i) = .TRUE.
   END DO
 END IF
 WRITE(11,"('End coordinates')")
 WRITE(11,"('Elements')")

 RETURN
 END SUBROUTINE prtcoo
