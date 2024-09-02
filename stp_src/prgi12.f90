 SUBROUTINE prgi12( )
 !
 !  print results for GiD for
 !  3-D BST Shell elements
 !
 USE data_db
 IMPLICIT NONE

 ! local variables
 TYPE (bst), POINTER :: eset
 TYPE (udmat), POINTER :: postd
 INTEGER :: ielem,i,iv,nelem,nnode,nstre,ngaus,is,ic, &
            nvar,naux,matno,j,dim
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname
 CHARACTER (len=17) :: saux,caux

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => bst_head

 !        print Gaussian variables

 DO
   gpname = 'GP'//eset%sname
   nelem = eset%nelem
   nnode = eset%nnode
   nstre = eset%nstre
   ngaus = eset%ngaus
   iv = 0       !initializes pointer

   IF(bst_force > 1)THEN
     CALL cab_gid(4,2,bst_force_l(5),bst_force_l(6),3,loadstep,ttime,gpname,units=bst_force_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem,eset%elvar(iv-2:iv,1,ielem)*bst_force_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_momen > 1)THEN
     CALL cab_gid(4,2,bst_momen_l(5),bst_momen_l(6),3,loadstep,ttime,gpname,units=bst_momen_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem,eset%elvar(iv-2:iv,1,ielem)*bst_momen_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_shear > 1)THEN
     CALL cab_gid(2,2,bst_shear_l(4),bst_shear_l(5),2,loadstep,ttime,gpname,units=bst_shear_u)
     iv = iv+2
     DO ielem=1,nelem
       WRITE(13,"(i8,(2e15.5))")ielem, eset%elvar(iv-1:iv,1,ielem)*bst_shear_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_logst > 1)THEN
     CALL cab_gid(3,2,bst_logst_l(5),bst_logst_l(6),3,loadstep,ttime,gpname,units=bst_logst_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, eset%elvar(iv-2:iv,1,ielem)*bst_logst_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_curva > 1)THEN
     CALL cab_gid(4,2,bst_curva_l(5),bst_curva_l(6),3,loadstep,ttime,gpname,units=bst_curva_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, eset%elvar(iv-2:iv,1,ielem)*bst_curva_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_thrat > 1)THEN
     CALL cab_gid(1,2,bst_thrat_l(3),bst_thrat_l(4),1,loadstep,ttime,gpname,units=bst_thrat_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,eset%elvar(iv,1,ielem)*bst_thrat_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_eqpst > 1)THEN
     is = LEN_TRIM(bst_eqpst_l(3))
     saux = bst_eqpst_l(3)(1:is)//'_1'
     ic = LEN_TRIM(bst_eqpst_l(4))
     caux = bst_eqpst_l(4)(1:ic)//'_1'
     CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,MAX(eset%elvar(iv,1,ielem),0.)*bst_eqpst_f
     END DO
     WRITE(13,"('End Values')")

     saux = bst_eqpst_l(3)(1:is)//'_N'
     caux = bst_eqpst_l(4)(1:ic)//'_N'
     CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,MAX(eset%elvar(iv,1,ielem),0.)*bst_eqpst_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_vmise > 1)THEN
     is = LEN_TRIM(bst_vmise_l(3))
     saux = bst_vmise_l(3)(1:is)//'_1'
     ic = LEN_TRIM(bst_vmise_l(4))
     caux = bst_vmise_l(4)(1:ic)//'_1'
     CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,eset%elvar(iv,1,ielem)*bst_vmise_f
     END DO
     WRITE(13,"('End Values')")

     saux = bst_vmise_l(3)(1:is)//'_N'
     caux = bst_vmise_l(4)(1:ic)//'_N'
     CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,eset%elvar(iv,1,ielem)*bst_vmise_f
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(bst_fldma > 1)THEN
     CALL cab_gid(1,2,bst_fldma_l(3),bst_fldma_l(4),1,loadstep,ttime,gpname)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem,eset%elvar(iv,1,ielem)
     END DO
     WRITE(13,"('End Values')")
   END IF
   IF (bst_wfldFZ) THEN
     CALL cab_gid(1,2,bst_fldFZ_l(1),bst_fldFZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_fr%name))
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))") ielem, eset%elvar(iv,1,ielem)
     END DO
     WRITE(13,"('End Values')")
   END IF
   IF (bst_wfldSZ) THEN
     CALL cab_gid(1,2,bst_fldSZ_l(1),bst_fldSZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_sf%name))
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))") ielem, eset%elvar(iv,1,ielem)
     END DO
     WRITE(13,"('End Values')")
   END IF

   !  user defined internal variables
   naux = eset%nstre - 13
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
           WRITE(13,"(i8,(e15.5))")ielem,eset%elvar(iv,1,ielem)
         END DO
       CASE (1) !vector
         CALL cab_gid(2,2,postd%name(1,i),postd%name(2,i),postd%dim(i), &
                      loadstep,ttime,gpname)
         j  = iv+1
         iv = iv+postd%dim(i)
         DO ielem=1,nelem
           WRITE(13,"(i8,(3e15.5))")ielem, eset%elvar(j:iv,1,ielem)
         END DO
       CASE (2) !tensor
         dim = postd%dim(i)
         CALL cab_gid(4,2,postd%name(1,i),postd%name(2,i),dim, &
                      loadstep,ttime,gpname)
         j  = iv+1
         iv = iv+dim
         DO ielem=1,nelem
           WRITE(13,"(i8,(6e15.5))")ielem, eset%elvar(j:iv,1,ielem)
         END DO
       END SELECT
       WRITE(13,"('End Values')")
     END DO
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED(eset) )EXIT

 END DO

 !        print Nodal variables

 IF( bst_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => bst_nodes
 vargs => bst_vargs
 iv = 0       !initializes pointer to vargs

 IF(MOD(bst_force,2) == 1)THEN
   CALL cab_gid(4,1,bst_force_l(1),bst_force_l(2),3,loadstep,ttime,gpname,units=bst_force_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*bst_force_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_momen,2) == 1)THEN
   CALL cab_gid(4,1,bst_momen_l(1),bst_momen_l(2),3,loadstep,ttime,gpname,units=bst_momen_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*bst_momen_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_shear,2) == 1)THEN
   CALL cab_gid(2,1,bst_shear_l(1),bst_shear_l(2),2,loadstep,ttime,gpname,units=bst_shear_u)
   iv = iv+2
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(3e15.5))")label(i),vargs(iv-1:iv,is)*bst_shear_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_logst,2) == 1)THEN
   CALL cab_gid(2,1,bst_logst_l(1),bst_logst_l(2),3,loadstep,ttime,gpname,units=bst_logst_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*bst_logst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_curva,2) == 1)THEN
   CALL cab_gid(4,1,bst_curva_l(1),bst_curva_l(2),3,loadstep,ttime,gpname,units=bst_curva_u)
   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*bst_curva_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_thrat,2) == 1)THEN
   CALL cab_gid(1,1,bst_thrat_l(1),bst_thrat_l(2),1,loadstep,ttime,gpname,units=bst_thrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*bst_thrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_eqpst,2) == 1)THEN
   is = LEN_TRIM(bst_eqpst_l(1))
   saux = bst_eqpst_l(1)(1:is)//'_1'
   ic = LEN_TRIM(bst_eqpst_l(2))
   caux = bst_eqpst_l(2)(1:ic)//'_1'
   CALL cab_gid(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*bst_eqpst_f
   END DO
   WRITE(13,"('End Values')")

   is = LEN_TRIM(bst_eqpst_l(1))
   saux = bst_eqpst_l(1)(1:is)//'_N'
   ic = LEN_TRIM(bst_eqpst_l(2))
   caux = bst_eqpst_l(2)(1:ic)//'_N'
   CALL cab_gid(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*bst_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_vmise,2) == 1)THEN
   is = LEN_TRIM(bst_vmise_l(1))
   saux = bst_vmise_l(1)(1:is)//'_1'
   ic = LEN_TRIM(bst_vmise_l(2))
   caux = bst_vmise_l(2)(1:ic)//'_1'
   CALL cab_gid(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*bst_vmise_f
   END DO
   WRITE(13,"('End Values')")

   is = LEN_TRIM(bst_vmise_l(1))
   saux = bst_vmise_l(1)(1:is)//'_N'
   ic = LEN_TRIM(bst_vmise_l(2))
   caux = bst_vmise_l(2)(1:ic)//'_N'
   CALL cab_gid(1,1,saux,caux,1,loadstep,ttime,gpname,units=bst_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*bst_vmise_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(bst_fldma,2) == 1)THEN
   CALL cab_gid(1,1,bst_fldma_l(1),bst_fldma_l(2),1,loadstep,ttime,gpname)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi12
