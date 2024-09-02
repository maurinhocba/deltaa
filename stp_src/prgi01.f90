 SUBROUTINE prgi01( )
 !
 !  print results for GiD for spot elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE (spot), POINTER :: eset
 INTEGER :: ielem,g,iv,is,i
 INTEGER, POINTER :: nodes(:,:)
 REAL(kind=8), POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname

 INTERFACE
   INCLUDE 'cab_gid.h'
 END INTERFACE

 eset => spot_head          !point to first set

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values

   IF(spot_force > 1)THEN
     CALL cab_gid(1,2,spot_force_l(1),spot_force_l(2),1,loadstep,ttime,gpname,units=spot_force_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*spot_force_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(spot_shear > 1)THEN
     CALL cab_gid(1,2,spot_shear_l(1),spot_shear_l(2),1,loadstep,ttime,gpname,units=spot_shear_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*spot_shear_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(spot_eqpst > 1)THEN
     CALL cab_gid(1,2,spot_eqpst_l(1),spot_eqpst_l(2),1,loadstep,ttime,gpname,units=spot_eqpst_u)
     iv = iv+1
     DO ielem=1,eset%nelem
       WRITE(13,"(i8,(e15.5))")ielem, &
                 (eset%elvar(iv,g,ielem)*spot_eqpst_f, g=1,eset%ngaus)
     END DO
     WRITE(13,"('End Values')")
   END IF

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

 RETURN
 END SUBROUTINE prgi01
