 SUBROUTINE prgi08( )
 !
 !  print results for GiD for beame elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE (beame), POINTER :: eset
 INTEGER :: ielem,g,iv,is,i,ng
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname
 INTEGER, PARAMETER :: o(3,3) = RESHAPE( (/ 1, 0, 0,  1, 2, 0, 1, 3, 2 /), &
                                         (/ 3, 3 /) )

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => beame_head          !point to first set

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values
   ng = eset%ngaus

   IF(beame_force > 1)THEN
     CALL cab_gid(2,2,beame_force_l(5),beame_force_l(6),3,loadstep,ttime,gpname,units=beame_force_u)
     iv = iv+3
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*beame_force_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(beame_momen > 1)THEN
     CALL cab_gid(2,2,beame_momen_l(5),beame_momen_l(6),3,loadstep,ttime,gpname,units=beame_momen_u)
     iv = iv+3
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*beame_momen_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(beame_eqpst > 1)THEN
!     CALL cab_gid(1,2,beame_eqpst_l(3),beame_eqpst_l(4),3,loadstep,ttime,gpname,units=beame_eqpst_u)
     CALL cab_gid(1,2,beame_eqpst_l(3),beame_eqpst_l(4),1,loadstep,ttime,gpname,units=beame_eqpst_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (MAX(eset%elvar(iv,o(g,ng),ielem),0.)*beame_eqpst_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   !IF(beame_vmise > 1)THEN
   ! CALL cab_gid(1,2,beame_vmise_l(3),beame_vmise_l(4),3,loadstep,ttime,gpname,units=beame_vmise_u))
   !  iv = iv+1
   !  DO ielem=1,eset%nelem
   !    WRITE(13,"(i8,(e15.5))")ielem, &
   !              (eset%elvar(iv,o(g,ng),ielem)*beame_vmise_f, g=1,ng)
   !  END DO
   ! WRITE(13,"('End Values')")
   !END IF

   IF(beame_arrat > 1)THEN
!     CALL cab_gid(1,2,beame_arrat_l(3),beame_arrat_l(4),3,loadstep,ttime,gpname,units=beame_arrat_u))
     CALL cab_gid(1,2,beame_arrat_l(3),beame_arrat_l(4),1,loadstep,ttime,gpname,units=beame_arrat_u)  !SCALAR ???
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*beame_arrat_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

 !        print Nodal variables

 IF( beame_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => beame_nodes
 vargs => beame_vargs
 iv = 0

 IF(MOD(beame_force,2) == 1)THEN
   CALL cab_gid(2,1,beame_force_l(1),beame_force_l(2),3,loadstep,ttime,gpname,units=beame_force_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-2:iv,is)*beame_force_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(beame_momen,2) == 1)THEN
   CALL cab_gid(2,1,beame_momen_l(1),beame_momen_l(2),3,loadstep,ttime,gpname,units=beame_momen_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-2:iv,is)*beame_momen_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(beame_eqpst,2) ==  1)THEN
   CALL cab_gid(1,1,beame_eqpst_l(1),beame_eqpst_l(2),1,loadstep,ttime,gpname,units=beame_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*beame_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 !IF(MOD(beame_vmise,2) ==  1)THEN
!   CALL cab_gid(1,1,beame_vmise_l(1),beame_vmise_l(2),1,loadstep,ttime,gpname,units=beame_vmise_u)
 !  iv = iv+1
 !  is = 0
 !  DO i=1,npoin
 !   IF(nodes(i,1) == 0)CYCLE
 !    is = is + 1
 !    WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*beame_vmise_f
 !  END DO
 !  WRITE(13,"('End Values')")
 !END IF

 IF(MOD(beame_arrat,2) ==  1)THEN
   CALL cab_gid(1,1,beame_arrat_l(1),beame_arrat_l(2),1,loadstep,ttime,gpname,units=beame_arrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*beame_arrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi08
