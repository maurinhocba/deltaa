 SUBROUTINE wridva( )
 USE data_db
 IMPLICIT NONE

 INTEGER (kind=4) idime,ipoin,j,i
 REAL (kind=8), PARAMETER :: valor1=0.999D+99, valor2 =0.999D-99
 REAL(kind=8), POINTER :: dsp(:,:)
 CHARACTER, PARAMETER :: tmp(3,3)= (/'N',' ',' ','B','T',' ','N','B','T' /)
 CHARACTER (len=15) :: var1,var2
 REAL (kind=8) :: x(3)

 INTERFACE
   SUBROUTINE wridva_spot(d,f)
     REAL(kind=8), POINTER :: d(:,:)
     REAL(kind=8) :: f
   END SUBROUTINE wridva_spot
   SUBROUTINE wridva_bpr(d,f,type)
     REAL(kind=8), POINTER :: d(:,:)
     REAL(kind=8) :: f
     INTEGER(kind=4) :: type
   END SUBROUTINE wridva_bpr
   INCLUDE 'cab_gid.h'
   INCLUDE 'cab_gid_bin.h'
 END INTERFACE

 x = 0d0
 !                    Print displacements
 IF( wtdisp )THEN                 !total displacements
   IF(ip == 2)THEN
     CALL cab_gid(2,1,tdisp_l(1),tdisp_l(2),ndime,loadstep,ttime,units=tdisp_u)
     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),displ(1:ndime,ipoin)*tdisp_f
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(displ,tdisp_f )
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(displ,tdisp_f,1)
     IF( ntype == 4 .AND. ngp > 0 )THEN
       DO ipoin=1,ngp
         WRITE(13,2003)labea(ipoin),dispa(1:ndime,ipoin)*tdisp_f
       END DO
     END IF
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5 ) THEN
     CALL cab_gid_bin(2,1,tdisp_l(1),tdisp_l(2),ndime,loadstep,ttime,units=tdisp_u)
     IF (ndime == 2) x(3) = 0d0
     DO ipoin=1,npoin
        IF (meshn(ipoin)) THEN
          x(1:ndime) = displ(:,ipoin)*tdisp_f
          CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
        END IF
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(displ,tdisp_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(displ,tdisp_f,1)
     CALL GID_ENDRESULT()

   END IF

 END IF

 IF( wsdisp )THEN               !stage displacements
   ALLOCATE(dsp(ndime,npoin))
   DO ipoin=1,npoin
     DO idime=1,ndime
       dsp(idime,ipoin) = displ(idime,ipoin) &
                         + coord(idime,ipoin) - coors(idime,ipoin)
     END DO
   END DO
   IF(ip == 2)THEN
     CALL cab_gid(2,1,sdisp_l(1),sdisp_l(2),ndime,loadstep,ttime,units=sdisp_u)
     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),dsp(1:ndime,ipoin)*sdisp_f
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(dsp,sdisp_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(dsp,sdisp_f,1)
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     CALL cab_gid_bin(2,1,sdisp_l(1),sdisp_l(2),ndime,loadstep,ttime,units=sdisp_u)
     DO ipoin=1,npoin
       IF (meshn(ipoin)) THEN
         x(1:ndime)=dsp(1:ndime,ipoin)*sdisp_f
         CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
       END IF
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(dsp,sdisp_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(dsp,sdisp_f,1)
     CALL GID_ENDRESULT()

   END IF
 END IF

 IF( widisp )THEN
   IF(ip == 2)THEN
     CALL cab_gid(2,1,idisp_l(1),idisp_l(2),ndime,loadstep,ttime,units=idisp_u)
     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),dispi(1:ndime,ipoin)*idisp_f
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(dispi,idisp_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(dispi,idisp_f,1)
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5 ) THEN
     CALL cab_gid_bin(2,1,idisp_l(1),idisp_l(2),ndime,loadstep,ttime,units=idisp_u)
     DO ipoin=1,npoin
       IF (meshn(ipoin))THEN
         x(1:ndime)=dispi(1:ndime,ipoin)*idisp_f
         CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
       END IF
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(dispi,idisp_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(dispi,idisp_f,1)
     CALL GID_ENDRESULT()

   END IF

 END IF
 !                    Print Velocities
 IF( wveloc )THEN
   DO ipoin=1,npoin
     DO idime=1,ndime
       IF(ABS(veloc(idime,ipoin)) > valor1)    &
              veloc(idime,ipoin) = SIGN(valor1,veloc(idime,ipoin))
       IF(ABS(veloc(idime,ipoin)) < valor2 .AND. veloc(idime,ipoin) /= 0.0) &
              veloc(idime,ipoin) = SIGN(valor2,veloc(idime,ipoin))
     END DO
   END DO
   IF(ip == 2)THEN
     CALL cab_gid(2,1,veloc_l(1),veloc_l(2),ndime,loadstep,ttime,units=veloc_u)

     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),veloc(1:ndime,ipoin)*veloc_f
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(veloc,veloc_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(veloc,veloc_f,1)
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     CALL cab_gid_bin(2,1,veloc_l(1),veloc_l(2),ndime,loadstep,ttime,units=veloc_u)
     DO ipoin=1,npoin
       IF (meshn(ipoin))THEN
         x(1:ndime)=veloc(1:ndime,ipoin)*idisp_f
         CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
       END IF
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(veloc,veloc_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(veloc,veloc_f,1)
     CALL GID_ENDRESULT()
   END IF

 END IF

 !                    Print Accelerations
 IF( waccel )THEN
   DO ipoin=1,npoin
     DO idime=1,ndime
       IF(ABS(accel(idime,ipoin)) > valor1)    &
              accel(idime,ipoin) = SIGN(valor1,accel(idime,ipoin))
        IF(ABS(accel(idime,ipoin)) < valor2 .AND. accel(idime,ipoin) /= 0.0) &
               accel(idime,ipoin) = SIGN(valor2,accel(idime,ipoin))
     END DO
   END DO
   IF(ip == 2)THEN
     CALL cab_gid(2,1,accel_l(1),accel_l(2),ndime,loadstep,ttime,units=accel_u)

     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),accel(1:ndime,ipoin)*accel_f
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(accel,accel_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(accel,accel_f,1)
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     CALL cab_gid_bin(2,1,accel_l(1),accel_l(2),ndime,loadstep,ttime,units=accel_u)
     DO ipoin=1,npoin
       IF (meshn(ipoin))THEN
         x(1:ndime)=accel(1:ndime,ipoin)*idisp_f
         CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
       END IF
     END DO
     IF( spot_auxiliar_nodes > 0 )CALL wridva_spot(accel,accel_f)
     IF( sol3d_auxiliar_nodes > 0 )CALL wridva_bpr(accel,accel_f,1)
     CALL GID_ENDRESULT()

   END IF

 END IF

 IF ( weuler ) THEN   ! prints local axes

   x = 0d0
   IF(ip == 2)THEN
     IF( ndime == 2 )THEN
       CALL cab_gid(1,1,euler_l(1),euler_l(2),1,loadstep,ttime,units=euler_u)
     ELSE ! (ndime == 3 )THEN
       CALL cab_gid(2,1,euler_l(1),euler_l(2),3,loadstep,ttime,units=euler_u)
     END IF

     DO ipoin=1,npoin
       IF(meshn(ipoin))THEN
         x(1:neulr) = euler(1:neulr,ipoin)*euler_f
         WRITE(13,2003)label(ipoin),x
       END IF
     END DO
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     IF (neulr == 3) THEN
       CALL cab_gid_bin(2,1,euler_l(1),euler_l(2),3,loadstep,ttime,units=euler_u)
       DO ipoin=1,npoin
         IF(meshn(ipoin))THEN
           x(1:neulr) = euler(1:neulr,ipoin)*euler_f
           CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
         END IF
       END DO
     ELSE ! IF (neulr == 1) THEN
       CALL cab_gid_bin(1,1,euler_l(1),euler_l(2),1,loadstep,ttime,units=euler_u)
       DO ipoin=1,npoin
          IF (meshn(ipoin)) CALL GID_WRITESCALAR(label(ipoin),euler(1,ipoin)*euler_f)
       END DO
     END IF
     CALL GID_ENDRESULT()

   END IF
 END IF

 IF ( wangve ) THEN   ! prints local axes velocities

   x = 0d0
   IF(ip == 2)THEN
     IF( ndime == 2 )THEN
       CALL cab_gid(1,1,angve_l(1),angve_l(2),1,loadstep,ttime,units=angve_u)
     ELSE ! IF( ndime == 3 )THEN
       CALL cab_gid(2,1,angve_l(1),angve_l(2),3,loadstep,ttime,units=angve_u)
     END IF

     DO ipoin=1,npoin
       SELECT CASE (nrotd)
       CASE (1)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anvel(1,ipoin)*angve_f
       CASE (2)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anvel(1:2,ipoin)*angve_f,0d0
       CASE (3)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anvel(1:3,ipoin)*angve_f
       END SELECT
     END DO
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     SELECT CASE (nrotd)
     CASE (2:3)
       CALL cab_gid_bin(2,1,angve_l(1),angve_l(2),3,loadstep,ttime,units=angve_u)
       DO ipoin=1,npoin
         IF (meshn(ipoin)) THEN
           x(1:nrotd) = anvel(1:nrotd,ipoin)*angve_f
           CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
         END IF
       END DO
     CASE (1)
       CALL cab_gid_bin(1,1,euler_l(1),euler_l(2),1,loadstep,ttime,units=angve_u)
       DO ipoin=1,npoin
         IF (meshn(ipoin)) CALL GID_WRITESCALAR(label(ipoin),anvel(1,ipoin)*angve_f)
       END DO
     END SELECT
     CALL GID_ENDRESULT()


   END IF
 END IF

 IF ( wangac ) THEN   ! prints local axes accelerations

   x = 0d0
   IF(ip == 2)THEN

     IF( ndime == 2 )THEN
       CALL cab_gid(1,1,angac_l(1),angac_l(2),1,loadstep,ttime,units=angac_u)
     ELSE ! IF( ndime == 3 )THEN
       CALL cab_gid(2,1,angac_l(1),angac_l(2),3,loadstep,ttime,units=angac_u)
     END IF

     DO ipoin=1,npoin
       SELECT CASE (nrotd)
       CASE (1)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anacc(1,ipoin)*angac_f
       CASE (2)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anacc(1:2,ipoin)*angac_f,0d0
       CASE (3)
         IF(meshn(ipoin))WRITE(13,2003)label(ipoin),anacc(1:3,ipoin)*angac_f
       END SELECT
     END DO
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5) THEN
     SELECT CASE (nrotd)
     CASE (2:3)
       CALL cab_gid_bin(2,1,angac_l(1),angac_l(2),3,loadstep,ttime,units=angac_u)
       DO ipoin=1,npoin
         IF (meshn(ipoin))THEN
           x(1:nrotd) = anacc(1:nrotd,ipoin)*angac_f
           CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
         END IF
       END DO
     CASE (1)
       CALL cab_gid_bin(1,1,angac_l(1),angac_l(2),1,loadstep,ttime,units=angac_u)
       DO ipoin=1,npoin
         IF (meshn(ipoin)) CALL GID_WRITESCALAR(label(ipoin),anacc(1,ipoin)*angac_f)
       END DO
     END SELECT
     CALL GID_ENDRESULT()


   END IF
 END IF

 !                    Print nodal temperatures
 IF( wtempe )THEN
   IF(ip == 2)THEN
     DO j=1,ndoft
       CALL cab_gid(1,1,TRIM(tempe_l(1))//tmp(j,ndoft),TRIM(tempe_l(2))//tmp(j,ndoft),1,loadstep,ttime,units=tempe_u)
       DO ipoin=1,npoin
         IF(meshn(ipoin))WRITE(13,2003) label(ipoin),tempe(j,ipoin)*tempe_f
       END DO
       WRITE(13,"('End Values')")
     END DO
   ELSE IF (ip == 4 .OR. ip == 5) THEN
     DO j=1,ndoft
       i = MIN(13,LEN_TRIM(tempe_l(1)))
       var1 = tempe_l(1)(1:i)//'_'//tmp(j,ndoft)
       var2 = TRIM(tempe_l(2))//'_'//tmp(j,ndoft)
       CALL cab_gid_bin(1,1,var1,var2,1,loadstep,ttime,units=tempe_u)
       DO ipoin=1,npoin
         IF (meshn(ipoin)) CALL GID_WRITESCALAR(label(ipoin),tempe(j,ipoin)*tempe_f)
       END DO
       CALL GID_ENDRESULT()
     END DO
   END IF

 END IF

 !                    Print additional DOFs
 IF( addof .AND. waddof )THEN
   IF(ip == 2)THEN
     CALL cab_gid(2,1,addof_l(1),addof_l(2),ndime,loadstep,ttime,units=addof_u)
     DO ipoin=1,npoin
       IF(meshn(ipoin))WRITE(13,2003)label(ipoin),psi(1:nadd,ipoin)*addof_f,0d0
     END DO
     WRITE(13,"('End Values')")

   ELSE IF (ip == 4 .OR. ip == 5 ) THEN
     CALL cab_gid_bin(2,1,addof_l(1),addof_l(2),ndime,loadstep,ttime,units=addof_u)
     x(2:3) = 0d0
     DO ipoin=1,npoin
        IF (meshn(ipoin)) THEN
          x(1:nadd) = psi(:,ipoin)*addof_f
          CALL GID_WRITEVECTOR(label(ipoin),x(1),x(2),x(3))
        END IF
     END DO
     CALL GID_ENDRESULT()

   END IF

 END IF

 RETURN

 2000 FORMAT(1x,i2,2x,a8,5i5,a3)
 2001 FORMAT(1x,i2,i5,3e12.5)
 2002 FORMAT(a15,i5,e12.4,3i5)
 2003 FORMAT(i8,3e15.7)

 END SUBROUTINE wridva
