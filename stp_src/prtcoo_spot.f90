 SUBROUTINE prtcoo_spot(nelem,nnode,iset,sname,x,fl,lnods,matno)

 ! print coordinates for GiD

 USE data_db, ONLY : ndime,npoin,wsdisp,wtdisp,first,spot_nodes,label,coord,meshn,coors,tdisp_f,sdisp_f
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem,nnode,fl,iset,matno
 CHARACTER (len=32), INTENT(IN) :: sname
 INTEGER(kind=4), POINTER :: lnods(:,:)
 REAL (kind=8), POINTER :: x(:,:,:)

 CHARACTER (len=13) :: elmt = 'Linear       '
 INTEGER (kind=4) i,l,n,naux,iel
 LOGICAL :: fauxi

 IF(first)THEN
   first = .FALSE.
 ELSE
   WRITE(11,"('End Elements')")
 END IF
 l = LEN_TRIM(sname)
 WRITE(11,"('MESH ',a,' dimension =',i2,' ElemType ',a13, &
 &          ' Nnode = ',i2,/,'Coordinates')")sname(1:l),ndime,elmt,2
 fauxi = ASSOCIATED( x )
 IF( fauxi )THEN
   naux = SIZE(x,2)
 ELSE
   naux = 0
 END IF

 IF( .NOT.wsdisp .AND. wtdisp )THEN  !write original coordinates
   IF( nnode == 2 )THEN
     DO i=1,npoin
       IF ( spot_nodes(i,2) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coord(1:ndime,i)*tdisp_f
       IF ( spot_nodes(i,2) > 0 ) meshn(i) = .TRUE.
     END DO
   END IF
   DO i=1,naux                 !
     WRITE(11,"(i8,3e18.10)")fl+i,x(1:ndime,i,1)*tdisp_f
   END DO
 ELSE                                !Write stage coordinates
   IF( nnode == 2 )THEN
     DO i=1,npoin
       IF ( spot_nodes(i,2) == iset ) WRITE(11,"(i8,3e18.10)")label(i),coors(1:ndime,i)*sdisp_f
       IF ( spot_nodes(i,2) > 0 ) meshn(i) = .TRUE.
     END DO
     DO i=1,naux
       WRITE(11,"(i8,3e18.10)")fl+i,x(1:ndime,i,1)*sdisp_f   !always original coordinates
     END DO
   ELSE
     DO i=1,naux
       WRITE(11,"(i8,3e18.10)")fl+i,x(1:ndime,i,2)*sdisp_f
     END DO
   END IF
 END IF
 WRITE(11,"('End coordinates')")
 WRITE(11,"('Elements')")
 IF( nnode == 2 )THEN
   n = fl
   DO iel = 1,nelem            !for each element in the set
     IF( lnods(2,iel) > 0)THEN
       WRITE(11,"(10i8)") iel,label(lnods(:,iel)),matno !iset
     ELSE
       n = n+1
       WRITE(11,"(10i8)") iel,label(lnods(1,iel)),n,matno !iset
     END IF
   END DO
 ELSE
   n = fl + 1
   DO iel = 1,nelem            !for each element in the set
     WRITE(11,"(10i8)") iel,n,n+1,matno !iset
     n = n+2
   END DO
 END IF

 RETURN
 END SUBROUTINE prtcoo_spot
