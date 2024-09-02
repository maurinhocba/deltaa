 SUBROUTINE prgi02( )
 !
 !  print results for GiD for truss elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE (truss), POINTER :: eset
 INTEGER :: ielem,g,iv,is,i
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => truss_head          !point to first set

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values

   IF(truss_force > 1)THEN
     CALL cab_gid(1,2,truss_force_l(3),truss_force_l(4),1,loadstep,ttime,gpname,units=truss_force_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*truss_force_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(truss_stres > 1)THEN
     CALL cab_gid(1,2,truss_stres_l(3),truss_stres_l(4),1,loadstep,ttime,gpname,units=truss_stres_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*truss_stres_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(truss_eqpst > 1)THEN
     CALL cab_gid(1,2,truss_eqpst_l(3),truss_eqpst_l(4),1,loadstep,ttime,gpname,units=truss_eqpst_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*truss_eqpst_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

 !        print Nodal variables

 IF( truss_nvarn == 0 )RETURN    !no variables to print, => exit

 nodes => truss_nodes
 vargs => truss_vargs
 iv = 0

 IF(MOD(truss_force,2) == 1)THEN
   CALL cab_gid(1,1,truss_force_l(1),truss_force_l(2),1,loadstep,ttime,gpname,units=truss_force_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*truss_force_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(truss_stres,2) == 1)THEN
   CALL cab_gid(1,1,truss_stres_l(1),truss_stres_l(2),1,loadstep,ttime,gpname,units=truss_stres_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*truss_stres_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 IF(MOD(truss_eqpst,2) ==  1)THEN
   CALL cab_gid(1,1,truss_eqpst_l(1),truss_eqpst_l(2),1,loadstep,ttime,gpname,units=truss_eqpst_u)
   iv = iv+1
   is = 0
   DO i=1,npoin
     IF(nodes(i,1) == 0)CYCLE
     is = is + 1
     WRITE(13,"(i8,(e15.5))")label(i),vargs(iv,is)*truss_eqpst_f
   END DO
   WRITE(13,"('End Values')")
 END IF

 RETURN
 END SUBROUTINE prgi02
