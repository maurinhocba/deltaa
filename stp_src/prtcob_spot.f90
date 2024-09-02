 SUBROUTINE prtcob_spot(nelem,nnode,iset,sname,x,fl,lnods,matno)

 ! print coordinates for GiD

 USE data_db, ONLY : ndime,npoin,wsdisp,wtdisp,spot_nodes,label,coord,coors,meshn,tdisp_f,sdisp_f
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem,nnode,fl,iset,matno
 CHARACTER (len=32), INTENT(IN) :: sname
 INTEGER(kind=4), POINTER :: lnods(:,:)
 REAL (kind=8), POINTER :: x(:,:,:)

 INTEGER (kind=4), PARAMETER :: elmt = 2, nn = 2
 INTEGER (kind=4) i,l,n,naux,iel,idmat(3)
 REAL (kind=8) :: y(3)
 LOGICAL :: fauxi

 fauxi = ASSOCIATED( x )
 IF( fauxi )THEN
   naux = SIZE(x,2)
 ELSE
   naux = 0
 END IF

 IF(ndime == 2)y(3) = 0D0
 CALL GID_BEGINMESH(TRIM(sname),ndime,elmt,nn)
 CALL GID_BEGINCOORDINATES()
 IF( .NOT.wsdisp .AND. wtdisp )THEN  !write original coordinates
   IF( nnode == 2 )THEN
     DO i=1,npoin
       IF ( spot_nodes(i,2) == iset ) THEN
         y(1:ndime) = coord(:,i)*tdisp_f
         CALL GID_WRITECOORDINATES(label(i),y(1),y(2),y(3))
         meshn(i) = .TRUE.
       END IF
     END DO
   END IF
   DO i=1,naux                 !
     y(1:ndime) = x(:,i,1)*tdisp_f
     CALL GID_WRITECOORDINATES(fl+i,y(1),y(2),y(3))
   END DO
 ELSE                                !Write stage coordinates
   IF( nnode == 2 )THEN
     DO i=1,npoin
       IF ( spot_nodes(i,2) == iset ) THEN
         y(1:ndime) = coors(:,i)*tdisp_f
         CALL GID_WRITECOORDINATES(label(i),y(1),y(2),y(3))
         meshn(i) = .TRUE.
       END IF
     END DO
     DO i=1,naux                 !
       y(1:ndime) = x(:,i,1)*sdisp_f
       CALL GID_WRITECOORDINATES(fl+i,y(1),y(2),y(3))      !always original coordinates
     END DO
   ELSE
     DO i=1,naux                 !
       y(1:ndime) = x(:,i,1)*sdisp_f
       CALL GID_WRITECOORDINATES(fl+i,y(1),y(2),y(3))      !always original coordinates
     END DO
   END IF
 END IF
 CALL GID_ENDCOORDINATES()
 CALL GID_BEGINELEMENTS()
 IF( nnode == 2 )THEN
   n = fl
   DO iel = 1,nelem            !for each element in the set
     IF( lnods(2,iel) > 0 )THEN
       idmat = (/label(lnods(1:2,iel)),matno/)
     ELSE
       n = n+1
       idmat = (/ label(lnods(1,iel)),n,matno/)
     END IF
     CALL GID_WRITEELEMENTMAT(iel,idmat)
   END DO
 ELSE
   n = fl + 1
   DO iel = 1,nelem            !for each element in the set
     idmat = (/n,n+1,matno/)
     CALL GID_WRITEELEMENTMAT(iel,idmat)
     n = n+2
   END DO
 END IF
 CALL GID_ENDELEMENTS()

 RETURN
 END SUBROUTINE prtcob_spot
