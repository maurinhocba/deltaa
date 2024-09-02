 SUBROUTINE prgi03( )
 !
 !  print results for GiD for
 !  2-D solid elements
 !
 USE data_db
 IMPLICIT NONE

 ! local variables
 TYPE (sol2d), POINTER :: eset
 INTEGER :: ielem,g,i,iv,nelem,nnode,ngaus,is,ng
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8) , POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname
 !     non standart GiD Gauss point positions (INTERNAL)
 INTEGER, PARAMETER :: o(9,4) = RESHAPE( (/ 1, 0, 0, 0, 0, 0, 0, 0, 0,   &
                                            1, 3, 4, 2, 0, 0, 0, 0, 0,   & !quad4 
                                            1, 7, 9, 3, 4, 8, 6, 2, 5,   & !quad8 & quad9
                                            1, 2, 3, 0, 0, 0, 0, 0, 0     /),& !trian-6
                               (/ 9,4 /) )

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => sol2d_head

 !        print Gaussian variables

 DO
   gpname = 'GP'//eset%sname
   nelem = eset%nelem
   nnode = eset%nnode
   ngaus = eset%ngaus
   ng = NINT(SQRT(REAL(ngaus)))
   IF( ngaus == 3 ) ng = 4
   iv = 0       !initializes pointer

   IF(sol2d_stres > 1)THEN
     CALL cab_gid(4,2,sol2d_stres_l(6),sol2d_stres_l(7),4,loadstep,ttime,gpname,units=sol2d_stres_u)
     iv = iv+4
     DO ielem=1,nelem
       WRITE(13,"(i8,(4e15.5))")ielem, &
                 (eset%elvar(iv-3:iv,o(g,ng),ielem)*sol2d_stres_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_logst > 1)THEN
     CALL cab_gid(2,2,sol2d_logst_l(5),sol2d_logst_l(6),3,loadstep,ttime,gpname,units=sol2d_logst_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*sol2d_logst_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_shtst > 1)THEN
     CALL cab_gid(2,2,sol2d_shtst_l(5),sol2d_shtst_l(6),3,loadstep,ttime,gpname,units=sol2d_shtst_u)
     iv = iv+3
     DO ielem=1,nelem
       WRITE(13,"(i8,(3e15.5))")ielem, &
                 (eset%elvar(iv-2:iv,o(g,ng),ielem)*sol2d_shtst_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_thrat > 1)THEN
     CALL cab_gid(1,2,sol2d_thrat_l(3),sol2d_thrat_l(4),1,loadstep,ttime,gpname,units=sol2d_thrat_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*sol2d_thrat_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_eqpst > 1)THEN
     CALL cab_gid(1,2,sol2d_eqpst_l(3),sol2d_eqpst_l(4),1,loadstep,ttime,gpname,units=sol2d_eqpst_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (MAX(eset%elvar(iv,o(g,ng),ielem),0.)*sol2d_eqpst_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_vmise > 1)THEN
     CALL cab_gid(1,2,sol2d_vmise_l(3),sol2d_vmise_l(4),1,loadstep,ttime,gpname,units=sol2d_vmise_u)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem)*sol2d_vmise_f, g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(sol2d_fldma > 1)THEN
     CALL cab_gid(1,2,sol2d_fldma_l(3),sol2d_fldma_l(4),1,loadstep,ttime,gpname)
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,o(g,ng),ielem), g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF
   IF (sol2d_wfldFZ) THEN
     CALL cab_gid(1,2,sol2d_fldFZ_l(1),sol2d_fldFZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_fr%name))
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))") ielem, &
                 (eset%elvar(iv,o(g,ng),ielem), g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF
   IF (sol2d_wfldSZ) THEN
     CALL cab_gid(1,2,sol2d_fldSZ_l(1),sol2d_fldSZ_l(1),1,loadstep,ttime,gpname,TRIM(rrt_sf%name))
     iv = iv+1
     DO ielem=1,nelem
       WRITE(13,"(i8,(e15.5))") ielem, &
                 (eset%elvar(iv,o(g,ng),ielem), g=1,ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF


   eset => eset%next
   IF( .NOT.ASSOCIATED(eset) )EXIT
 END DO

 !        print Nodal variables

 IF( sol2d_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => sol2d_nodes
 vargs => sol2d_vargs
 iv = 0       !initializes pointer to vargs

 IF(MOD(sol2d_stres,2) == 1)THEN
   CALL cab_gid(4,1,sol2d_stres_l(1),sol2d_stres_l(2),4,loadstep,ttime,gpname,units=sol2d_stres_u)

   iv = iv+4
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-3:iv,is)*sol2d_stres_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_logst,2) == 1)THEN
   CALL cab_gid(2,1,sol2d_logst_l(1),sol2d_logst_l(2),3,loadstep,ttime,gpname,units=sol2d_logst_u)

   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*sol2d_logst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_shtst,2) == 1)THEN
   CALL cab_gid(2,1,sol2d_shtst_l(1),sol2d_shtst_l(2),3,loadstep,ttime,gpname,units=sol2d_shtst_u)

   iv = iv+3
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv-2:iv,is)*sol2d_shtst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_thrat,2) == 1)THEN
   CALL cab_gid(1,1,sol2d_thrat_l(1),sol2d_thrat_l(2),1,loadstep,ttime,gpname,units=sol2d_thrat_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol2d_thrat_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_eqpst,2) == 1)THEN
   CALL cab_gid(1,1,sol2d_eqpst_l(1),sol2d_eqpst_l(2),1,loadstep,ttime,gpname,units=sol2d_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol2d_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_vmise,2) == 1)THEN
   CALL cab_gid(1,1,sol2d_vmise_l(1),sol2d_vmise_l(2),1,loadstep,ttime,gpname,units=sol2d_vmise_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(4e15.5))")label(i),vargs(iv,is)*sol2d_vmise_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(sol2d_fldma,2) == 1)THEN
   CALL cab_gid(1,1,sol2d_fldma_l(1),sol2d_fldma_l(2),1,loadstep,ttime,gpname)
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
 END SUBROUTINE prgi03
