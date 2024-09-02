 SUBROUTINE coract(ndime,npoin,neulr,ndofn,lambd,coord,ddisp, &
                   coora,euler)

 !     UPDATES CONFIGURATION AND NODAL SYSTEMS
 USE ctrl_db, ONLY : nrotd
 USE kinc_db, ONLY : nn,nn1,nm,nescv,ndepd,naris,&
                     nndpd,distd,npsdf,nesdf,ftsdf,nvelr,velor
 USE npo_db, ONLY : naeul,ifpre,psia
 USE nsld_db, ONLY : nnsld
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: ndime,npoin,ndofn,neulr
 REAL (kind=8),INTENT(IN) :: ddisp(:),lambd,coord(:,:)
 REAL (kind=8),INTENT(IN OUT) :: coora(:,:)
 REAL (kind=8), POINTER :: euler(:,:)

 INTEGER (kind=4) :: ipoin,idofn,ieq,i,ie,ib,iq,iv,na,idofa
 REAL    (kind=8) :: theta(3),auxil
 !      REAL (kind=8),PARAMETER :: pi = 3.1415926535898_8

 INTERFACE
   INCLUDE 'actrot.h'
   INCLUDE 'arisfi.h'
   INCLUDE 'fixdpd.h'
   INCLUDE 'fixsld.h'
 END INTERFACE

 iv = nvelr+1
 !     translational degrees of freedom
 DO ipoin=1,npoin
   DO idofn=1,ndime
     ieq = ifpre(idofn,ipoin)
     SELECT CASE (ieq)
     CASE (:-nn1)           !fixed DOF
       coora(idofn,ipoin) = coord(idofn,ipoin)+velor(-ieq-nn,iv)*lambd
     CASE (-nn)
       coora(idofn,ipoin) = coord(idofn,ipoin)	
     CASE (-nm:-1)          !slave DOF
       ib = npsdf(-ieq)
       ie = npsdf(-ieq+1)-1
       auxil = 0d0
       DO i = ib,ie
         iq = nesdf(i)
         IF(iq > 0)THEN
           auxil = auxil + ddisp(iq)*ftsdf(i)
         ELSE IF(iq < -nn)THEN
           auxil = auxil + velor(-iq-nn,iv)*ftsdf(i)*lambd
         END IF
       END DO
       coora(idofn,ipoin) = coord(idofn,ipoin) + auxil
     CASE (0)                 !null DOF
       coora(idofn,ipoin) = coord(idofn,ipoin)	
     CASE (1:)                !active DOF
       coora(idofn,ipoin) = coord(idofn,ipoin) + ddisp(ieq)
     END SELECT
   END DO
 END DO
 !     rotational degrees of freedom
 IF (neulr == 1) THEN
   idofa = 3
   na = ndofn - idofa
   idofn=ndime+1
   DO ipoin=1,npoin
     IF( .NOT. naeul(ipoin) )CYCLE
     ieq = ifpre(idofn,ipoin)
     SELECT CASE (ieq)
     CASE (:-nn1)                !fixed DOF
       euler(1,ipoin) = euler(1,ipoin)+lambd*velor(-ieq-nn,iv)
     CASE (-nn:-1)               !slave DOF
       ib = npsdf(-ieq)
       ie = npsdf(-ieq+1)-1
       auxil = 0d0
       DO i = ib,ie
         iq = nesdf(i)
         IF(iq > 0)THEN
           auxil = auxil + ddisp(iq)*ftsdf(i)
         ELSE IF(iq < -nn)THEN
           auxil = auxil + velor(-iq-nn,iv)*ftsdf(i)*lambd
         END IF
       END DO
       euler(1,ipoin) = euler(1,ipoin)+auxil
 !         CASE (0)                     !null DOF
 !           euler(1,ipoin) = euler(1,ipoin)
     CASE (1:)                     !active  DOF
       euler(1,ipoin) = euler(1,ipoin)+ddisp(ieq)
     END SELECT
   END DO
 ELSE IF (neulr == 9) THEN
   theta(3) = 0d0
   idofa = 6
   na = ndofn - idofa
   DO ipoin=1,npoin
     IF( .NOT. naeul(ipoin) )CYCLE
     DO idofn=ndime+1,nrotd
       ieq = ifpre(idofn,ipoin)
       SELECT CASE (ieq)
       CASE (:-nn1)        !fixed DOF
         theta(idofn-ndime) = lambd*velor(-ieq-nn,iv)
       CASE (-nn:-1)        !slave DOF
         ib = npsdf(-ieq)
         ie = npsdf(-ieq+1)-1
         auxil = 0d0
         DO i = ib,ie
           iq = nesdf(i)
           IF(iq > 0)THEN
             auxil = auxil + ddisp(iq)*ftsdf(i)
           ELSE IF(iq < -nn)THEN
            auxil = auxil + velor(-iq-nn,iv)*ftsdf(i)*lambd
           END IF
         END DO
         theta(idofn-ndime) = auxil
       CASE (0)               !null DOF
         theta(idofn-ndime) = 0d0
       CASE (1:)              !active or DOF
         theta(idofn-ndime) = ddisp(ieq)
       END SELECT
     END DO
     !write(58,"(2e15.6)")theta(1:2)
     CALL actrot(euler(1:9,ipoin),theta)
     !write(58,"(i7,9e15.6)")ipoin,euler(7:9,ipoin)
   END DO
 END IF
 IF(na > 0) THEN
   DO ipoin=1,npoin
     DO idofn=1,na
       ieq = ifpre(idofn+idofa,ipoin)
       SELECT CASE (ieq)
       CASE (:-nn1)        !fixed DOF
         psia(idofn,ipoin) = psia(idofn,ipoin) + lambd*velor(-ieq-nn,iv)
       CASE (-nn:-1)        !slave DOF
         ib = npsdf(-ieq)
         ie = npsdf(-ieq+1)-1
         auxil = 0d0
         DO i = ib,ie
           iq = nesdf(i)
           IF(iq > 0)THEN
             auxil = auxil + ddisp(iq)*ftsdf(i)
           ELSE IF(iq < -nn)THEN
            auxil = auxil + velor(-iq-nn,iv)*ftsdf(i)*lambd
           END IF
         END DO
         psia(idofn,ipoin) = psia(idofn,ipoin) + auxil
       CASE (0)               !null DOF
         !nothing
       CASE (1:)              !active or DOF
         psia(idofn,ipoin) = psia(idofn,ipoin) + ddisp(ieq)
       END SELECT
     END DO
   END DO
 END IF

 IF(ndepd > 0) CALL fixdpd(coora,euler)
 IF(naris > 0) CALL arisfi(naris,coora,nndpd(1:3,ndepd+1:ndepd+naris),euler)
 IF(nnsld > 0) CALL fixsld(coora,euler)

 RETURN

 END SUBROUTINE coract
