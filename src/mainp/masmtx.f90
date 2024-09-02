 SUBROUTINE masmtx(istop)
 !******************************************************************
 !
 ! *** calculates lumped mass for each element
 !
 !******************************************************************
 USE c_input
 USE cms_db
 USE ctrl_db, ONLY : ndime,ndofn,nrotd,npoin,numct,ndyna
 USE npo_db, ONLY : emass,label,ifpre,mass
 USE outp_db, ONLY : iwrit
 USE kinc_db, ONLY : maxav,ndumn
 USE ndp_db, ONLY : nndp

 IMPLICIT NONE
 INTEGER (kind=4), INTENT(OUT):: istop

 INTEGER (kind=4) :: i,ipoin,chnode
 REAL    (kind=8) :: xcmas(6),sumat
 LOGICAL          :: lumpd

 INTERFACE
   INCLUDE 'elemnt.h'
   INCLUDE 'contac.h'
 END INTERFACE

 emass = 0d0
 mass  = 0d0
 sumat = 0d0
 lumpd = (MOD(ndyna,2) == 0)

 CALL elemnt ('MASMTX', istop=istop, sumat=sumat, flag1=lumpd, &
                        mass=mass )
 IF(istop == 1) RETURN

 WRITE(lures,"(//'Total Mass of the System (Without Point Masses):'  &
               ,e15.7)",ERR=9999)sumat

 !     concentrated masses

 IF (nconm > 0) THEN
   IF(iwrit == 1) WRITE(lures,"(5x,'Concentrated Masses')",ERR=9999)

   DO i = 1,nconm
     ipoin = nodcms(i)
     ipoin = chnode(ipoin)
     xcmas(1:nrotd) = cmass(1:nrotd,i)
     WRITE(lures,"(i8,6f10.3)",ERR=9999) label(ipoin),xcmas
     IF(lumpd)THEN
       emass(1:nrotd,ipoin) = emass(1:nrotd,ipoin) + xcmas
     ELSE
       CALL ensmat(nrotd,ifpre(1,ipoin),xcmas(1),mass(1))
     END IF
   END DO
 END IF

 IF( numct > 0)CALL contac ('MASMTX')
 CALL elemnt ('RIGBDY')
 IF( ndumn > 0) CALL rigbdc (nndp)   !nodal dependencies of slave DOFs


 IF(iwrit == 1) THEN
   IF(ndime == 2)WRITE(lures,"(//,'  Lumped Mass Matrix' /,          &
                              ' Node',6X,'X',11X,'Y',11X,'a')",ERR=9999)
   IF(ndime == 3)WRITE(lures,"(//,'  Lumped Mass Matrix' /,' Node',  &
               6X,'X',11X,'Y',11X,'z',10x,'RX',10X,'RY',10X,'RZ')",ERR=9999)
   DO ipoin=1,npoin
     WRITE(lures,"(i5,8e12.4)",ERR=9999)label(ipoin),emass(1:ndofn,ipoin)
   END DO
 END IF
 RETURN
 9999 CALL runen2(' ')
 END SUBROUTINE masmtx
