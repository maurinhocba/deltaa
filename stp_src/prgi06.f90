 SUBROUTINE prgi06( )
 !
 !  print results for GiD for shl3d elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE (shl3d), POINTER :: eset
 TYPE (udmat), POINTER :: postd
 INTEGER :: ielem,g,iv,is,i,ng,ngaus,nelem,naux,matno,nvar,j,dim
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname
 !         non-standard Gauss positions for GiD
 INTEGER, PARAMETER :: o(4,3) = RESHAPE( (/ 1, 0, 0, 0, 1, 2, 3, 0,  1, 3, 4, 2 /), &
                                          (/ 4, 3 /) )

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => shl3d_head          !point to first set

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values
   ngaus = ABS(eset%ngaus)    !number of Gauss points per element
   nelem = eset%nelem    !number of elements
   ng = MAX(ngaus-1,1)   !1-2 for triangles 3 for quads

   IF(shl3d_force > 1)THEN
     CALL cab_gid(4,2,shl3d_force_l(5),shl3d_force_l(6),3,loadstep,ttime,gpname,units=shl3d_force_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_force_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_momen > 1)THEN
     CALL cab_gid(4,2,shl3d_momen_l(5),shl3d_momen_l(6),3,loadstep,ttime,gpname,units=shl3d_momen_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_momen_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_shear > 1)THEN
     CALL cab_gid(1,2,shl3d_shear_l(5),shl3d_shear_l(5),1,loadstep,ttime,gpname,units=shl3d_shear_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shl3d_shear_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")

     CALL cab_gid(1,2,shl3d_shear_l(6),shl3d_shear_l(6),1,loadstep,ttime,gpname,units=shl3d_shear_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shl3d_shear_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_logst > 1)THEN
     CALL cab_gid(2,2,shl3d_logst_l(5),shl3d_logst_l(6),3,loadstep,ttime,gpname,units=shl3d_logst_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*shl3d_logst_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_eqpst > 1)THEN
     CALL cab_gid(1,2,shl3d_eqpst_l(3),shl3d_eqpst_l(4),1,loadstep,ttime,gpname,units=shl3d_eqpst_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (MAX(eset%elvar(iv,o(g,ng),ielem),0.)*shl3d_eqpst_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_vmise > 1)THEN
     CALL cab_gid(1,2,shl3d_vmise_l(3),shl3d_vmise_l(4),1,loadstep,ttime,gpname,units=shl3d_vmise_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shl3d_vmise_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shl3d_thrat > 1)THEN
     CALL cab_gid(1,2,shl3d_thrat_l(3),shl3d_thrat_l(4),1,loadstep,ttime,gpname,units=shl3d_thrat_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*shl3d_thrat_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   !  user defined internal variables
   naux = eset%nstre - 11
   IF( naux > 0 )THEN
     matno = eset%matno(1)
     postd => udmat_head
     DO
       IF( postd%matno == matno )EXIT
       postd => postd%next
     END DO
     nvar = postd%nvar
     DO i=1,nvar
       SELECT CASE (postd%type(i))
       CASE (0) !scalar
         CALL cab_gid(1,2,postd%name(1,i),postd%name(2,i),1,            &
                      loadstep,ttime,gpname)
         iv = iv+1
         DO ielem=1,nelem
           WRITE(13,"(i8,(e15.5))")ielem,(eset%elvar(iv,o(g,ng),ielem),g=1,ngaus)
         END DO
       CASE (1) !vector
         CALL cab_gid(2,2,postd%name(1,i),postd%name(2,i),postd%dim(i), &
                      loadstep,ttime,gpname)
         j  = iv+1
         iv = iv+postd%dim(i)
         DO ielem=1,nelem
           WRITE(13,"(i8,(2e15.5))")ielem, (eset%elvar(j:iv,o(g,ng),ielem),g=1,ngaus)
         END DO
       CASE (2) !tensor
         dim = postd%dim(i)
         CALL cab_gid(4,2,postd%name(1,i),postd%name(2,i),dim, &
                      loadstep,ttime,gpname)
         j  = iv+1
         iv = iv+dim
         DO ielem=1,nelem
           WRITE(13,"(i8,(4e15.5))")ielem, (eset%elvar(j:iv,o(g,ng),ielem),g=1,ngaus)
         END DO
       END SELECT
       WRITE(13,"('End Values')")
     END DO
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

 !        print Nodal variables

 IF( shl3d_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => shl3d_nodes
 vargs => shl3d_vargs
 iv = 0

 IF(MOD(shl3d_force,2) == 1)THEN
   CALL cab_gid(4,1,shl3d_force_l(1),shl3d_force_l(2),3,loadstep,ttime,gpname,units=shl3d_force_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-2:iv,is)*shl3d_force_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_momen,2) == 1)THEN
   CALL cab_gid(4,1,shl3d_momen_l(1),shl3d_momen_l(2),3,loadstep,ttime,gpname,units=shl3d_momen_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-2:iv,is)*shl3d_momen_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_shear,2) == 1)THEN
   CALL cab_gid(1,1,shl3d_shear_l(2),shl3d_shear_l(2),1,loadstep,ttime,gpname,units=shl3d_shear_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*shl3d_shear_f
   END DO
   WRITE(13,"('End Values')")

   CALL cab_gid(1,1,shl3d_shear_l(3),shl3d_shear_l(3),1,loadstep,ttime,gpname,units=shl3d_shear_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*shl3d_shear_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_logst,2) == 1)THEN
   CALL cab_gid(2,1,shl3d_logst_l(1),shl3d_logst_l(2),3,loadstep,ttime,gpname,units=shl3d_logst_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-2:iv,is)*shl3d_logst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_eqpst,2) ==  1)THEN
   CALL cab_gid(1,1,shl3d_eqpst_l(1),shl3d_eqpst_l(2),1,loadstep,ttime,gpname,units=shl3d_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*shl3d_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_vmise,2) ==  1)THEN
   CALL cab_gid(1,1,shl3d_vmise_l(1),shl3d_vmise_l(2),1,loadstep,ttime,gpname,units=shl3d_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*shl3d_vmise_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(shl3d_thrat,2) ==  1)THEN
   CALL cab_gid(1,1,shl3d_thrat_l(1),shl3d_thrat_l(2),1,loadstep,ttime,gpname,units=shl3d_thrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*shl3d_thrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi06

