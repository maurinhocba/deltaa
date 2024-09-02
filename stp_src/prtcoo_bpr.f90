 SUBROUTINE prtcoo_bpr (nnode,nelem,snam,nodes,iset,lnods,ilab)

 ! print coordinates for GiD

 USE data_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem,nnode,nodes(npoin),iset,ilab
 CHARACTER (len=32), INTENT(IN) :: snam
 INTEGER(kind=4) :: lnods(nnode,nelem)

 CHARACTER (len=13) elmt
 INTEGER (kind=4) i,j,l
 CHARACTER (len=32) :: sname
 REAL(kind=8) :: x(3)

 sname = snam
 elmt = 'Prism        '

 IF(first)THEN
   first = .FALSE.
 ELSE
   WRITE(11,"('End Elements')")
 END IF
 l = LEN_TRIM(sname)
 WRITE(11,"('MESH ',a,' dimension =',i2,' ElemType ',a13, &
&          ' Nnode = ',i2,/,'Coordinates')")sname(1:l),3,elmt,nnode

 IF( .NOT.wsdisp .AND. wtdisp )THEN  !write original coordinates
   DO i=1,npoin
     IF ( nodes(i) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coord(1:ndime,i)*tdisp_f
     IF ( nodes(i) > 0 ) meshn(i) = .TRUE.
   END DO
   l = ilab
   DO i=1,nelem
     DO j=1,3
       x = (coord(:,lnods(j,i)) + coord(:,lnods(j+3,i)))/2d0
       l = l+1
       WRITE(11,"(i8,3e18.10)")l,x*tdisp_f
     END DO
   END DO
 ELSE !Write stage coordinates
   DO i=1,npoin
     IF ( nodes(i) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coors(1:ndime,i)*sdisp_f
     IF ( nodes(i) > 0 ) meshn(i) = .TRUE.
   END DO
   l = ilab
   DO i=1,nelem
     DO j=1,3
       x = (coors(:,lnods(j,i)) + coors(:,lnods(j+3,i)))/2d0
       l = l+1
       WRITE(11,"(i8,3e18.10)")l,x*tdisp_f
     END DO
   END DO
 END IF
 WRITE(11,"('End coordinates')")
 WRITE(11,"('Elements')")

 RETURN
 END SUBROUTINE prtcoo_bpr
