 SUBROUTINE rest_kc (ndofn,npoin)

 !   dumps kinematic condition databases

 USE kin0_db
 USE kin1_db
 USE kin2_db
 USE ifix_db
 USE nes_db
 USE nar_db
 USE ndp_db
 USE rve_db
      USE nsld_db, ONLY : nnsld,nsldat,nsldfg
 IMPLICIT NONE
 INTEGER (kind=4) :: ndofn,npoin

 INTEGER (kind=4) :: i,j,nrve
 TYPE (rve_set), POINTER :: rves
 TYPE (rve_nod), POINTER :: rven
 !--------------------------------------------------------------

 CALL rest_kin0( )                     !kinematic constraints
 CALL rest_kin1( )                     !assemble data for kinematic constraints
 CALL rest_kin2(ndofn,npoin)           !equation information and restrains
 CALL rest_ifx                         !boundary restrains data base

 READ (51) nnes
 IF (nnes > 0) THEN
   ALLOCATE (nesdat(nnes))
   DO i=1,nnes
     READ (51) nesdat(i)%node, nesdat(i)%ndof, nesdat(i)%factor
   END DO
 END IF

 IF (naris > 0) THEN
   ALLOCATE (nardat(3,naris))
   DO i=1,naris
     READ (51) nardat(1:3,i)
   END DO
 END IF

 READ (51) nndp
 IF (nndp > 0) THEN
   ALLOCATE (ndpdat(2,nndp))
   DO i=1,nndp
     READ (51) ndpdat(1:2,i)
   END DO
 END IF

 IF (nvelr > 0) THEN
   DO i=1,nvelr
     ALLOCATE ( rves)
     NULLIFY(rves%head)
     READ (51) rves%lc, rves%factor, rves%nrv, rves%dspflg
     nrve = rves%nrv
     DO j=1,nrve
       ALLOCATE (rven)
       READ (51) rven%node, rven%v(1:ndofn)
       CALL add_rven (rven,rves%head,rves%tail)
     END DO
     CALL add_rves (rves, headv, tailv)
   END DO
 END IF

 READ (51) nnsld
 IF (nnsld > 0) THEN
   ALLOCATE( nsldat(2,nnsld), nsldfg(3,nnsld) )
   DO i=1,nnsld
     READ (51) nsldat(1:2,i),nsldfg(1:3,i)
   END DO
 END IF

 RETURN

 END SUBROUTINE rest_kc
