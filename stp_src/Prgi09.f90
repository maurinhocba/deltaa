 SUBROUTINE prgi09(etype)
 !
 !  print results for GiD for shrev elements
 !
 USE data_db
 IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: etype

 TYPE (shrev), POINTER :: eset
 INTEGER :: ielem,g,iv,is,i,ng,nst,iv0,nstre
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname,auxvl
 INTEGER, PARAMETER :: o(3,3) = RESHAPE( (/ 1, 0, 0,  1, 2, 0, 1, 3, 2 /), &
                                         (/ 3, 3 /) )

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 nst = 2
 IF( ntype == 1 )nst = 1

 eset => shrev_head          !point to first set
 nstre = eset%nstre

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values
   ng = eset%ngaus

   IF(shrev_force > 1)THEN
     CALL cab_gid(1,2,shrev_force_l(5),shrev_force_l(5),nst,loadstep,ttime,gpname,units=shrev_force_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shrev_force_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
     IF(ntype > 1 .AND. ntype < 5)THEN
       CALL cab_gid(1,2,shrev_force_l(6),shrev_force_l(6),1,loadstep,ttime,gpname,units=shrev_force_u)
       iv = iv+1
       DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                   (eset%elvar(iv,o(g,ng),ielem)*shrev_force_f, g=1,ng)
       END DO
       WRITE(13,"('End Values')")
     END IF
   END IF

   IF(shrev_momen > 1)THEN
     CALL cab_gid(1,2,shrev_momen_l(5),shrev_momen_l(5),1,loadstep,ttime,gpname,units=shrev_momen_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shrev_momen_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
     IF(ntype > 1 .AND. ntype < 5)THEN
       CALL cab_gid(1,2,shrev_momen_l(6),shrev_momen_l(6),1,loadstep,ttime,gpname,units=shrev_momen_u)
       iv = iv+1
       DO ielem=1,eset%nelem
         WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                   (eset%elvar(iv,o(g,ng),ielem)*shrev_momen_f, g=1,ng)
       END DO
       WRITE(13,"('End Values')")
     END IF
   END IF

   IF(shrev_shear > 1)THEN
     CALL cab_gid(1,2,shrev_shear_l(3),shrev_shear_l(4),1,loadstep,ttime,gpname,units=shrev_shear_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shrev_shear_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   ! additional stresses at Gauss points Only
   IF( (ntype == 1 .AND. nstre == 5).OR. nstre == 7 )THEN
     IF( shrev_momen > 1 )THEN   !Gauss points moments desired?
       CALL cab_gid(1,2,'Add_Moment','M_phi',1,loadstep,ttime,gpname,units=shrev_momen_u)
       iv = iv+1
       DO ielem=1,eset%nelem
         WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                   (eset%elvar(iv,o(g,ng),ielem)*shrev_momen_f, g=1,ng)
       END DO
       WRITE(13,"('End Values')")
     END IF
     IF(shrev_shear > 1)THEN
       CALL cab_gid(1,2,'Add_shear','Q_phi',1,loadstep,ttime,gpname,units=shrev_shear_u)
       iv = iv+1
       DO ielem=1,eset%nelem
         WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem, &
                   (eset%elvar(iv,o(g,ng),ielem)*shrev_shear_f, g=1,ng)
       END DO
       WRITE(13,"('End Values')")
     END IF
   END IF

   IF(shrev_eqpst > 1)THEN
     CALL cab_gid(1,2,shrev_eqpst_l(3),shrev_eqpst_l(4),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
     iv = iv+1
     DO ielem=1,eset%nelem
         WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem,  &
                 (MAX(eset%elvar(iv,o(g,ng),ielem),0.)*shrev_eqpst_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shrev_vmise > 1)THEN
     CALL cab_gid(1,2,shrev_vmise_l(3),shrev_vmise_l(4),1,loadstep,ttime,gpname,units=shrev_vmise_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem,  &
                 (eset%elvar(iv,o(g,ng),ielem)*shrev_vmise_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shrev_thrat > 1)THEN
     CALL cab_gid(1,2,shrev_thrat_l(3),shrev_thrat_l(4),1,loadstep,ttime,gpname,units=shrev_thrat_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")ielem,  &
                 (eset%elvar(iv,o(g,ng),ielem)*shrev_thrat_f, g=1,ng)
     END DO
     WRITE(13,"('End Values')")
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

 !        print Nodal variables

 IF( shrev_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => shrev_nodes
 vargs => shrev_vargs
 iv = 0
 iv0 = 1

 IF(MOD(shrev_force,2) == 1)THEN
   CALL cab_gid(nst,1,shrev_force_l(1),shrev_force_l(2),nst,loadstep,ttime,gpname,units=shrev_force_u)
   iv = iv+nst
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,2e15.5,(/,8x,2e15.5))")label(i),vargs(iv0:iv,is)*shrev_force_f
   END DO
   WRITE(13,"('End Values')")
   iv0 = iv0 + nst
 END IF

 IF(MOD(shrev_momen,2) == 1)THEN
   CALL cab_gid(nst,1,shrev_momen_l(1),shrev_momen_l(2),nst,loadstep,ttime,gpname,units=shrev_momen_u)
   iv = iv+nst
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,2e15.5,(/,8x,2e15.5)))")label(i),vargs(iv0:iv,is)*shrev_momen_f
   END DO
   WRITE(13,"('End Values')")
   iv0 = iv0 + nst
 END IF

 IF(MOD(shrev_shear,2) == 1 )THEN
   CALL cab_gid(1,1,shrev_shear_l(1),shrev_shear_l(2),1,loadstep,ttime,gpname,units=shrev_shear_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,e15.5,(/,8x,e15.5))")label(i),vargs(iv,is)*shrev_shear_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shrev_eqpst,2) ==  1)THEN
   iv = iv+1
   auxvl = shrev_eqpst_l(1)
   IF( etype == 11 ) auxvl = TRIM(auxvl)//'_1'
   CALL cab_gid(1,1,auxvl,shrev_eqpst_l(2),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,e15.5,(/,8x,e15.5))")label(i),vargs(iv,is)*shrev_eqpst_f
   END DO
   WRITE(13,"('End Values')")
   iv0 = iv0+1
   IF( etype == 11 )THEN
     iv = iv+1
     auxvl = shrev_eqpst_l(1)
     auxvl = TRIM(auxvl)//'_n'
     CALL cab_gid(1,1,auxvl,shrev_eqpst_l(2),1,loadstep,ttime,gpname,units=shrev_eqpst_u)
     is = 0
     DO i=1,npoin
       IF(nodes(i,1) == 0)CYCLE
       is = is + 1
       WRITE(13,"(i8,e15.5,(/,8x,e15.5))")label(i),vargs(iv,is)*shrev_eqpst_f
     END DO
     WRITE(13,"('End Values')")
     iv0 = iv0+1
   END IF
 END IF

 IF(MOD(shrev_vmise,2) ==  1 .AND. etype == 9)THEN
   CALL cab_gid(1,1,shrev_vmise_l(1),shrev_vmise_l(2),1,loadstep,ttime,gpname,units=shrev_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,e15.5,(/,8x,e15.5))")label(i),vargs(iv,is)*shrev_vmise_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shrev_thrat,2) ==  1)THEN
   CALL cab_gid(1,1,shrev_thrat_l(1),shrev_thrat_l(2),1,loadstep,ttime,gpname,units=shrev_thrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,e15.5,(/,8x,e15.5))")label(i),vargs(iv,is)*shrev_thrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi09
