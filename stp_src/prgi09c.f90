 SUBROUTINE prgi09c( )
 !
 !  print results for GiD for beam elements
 !
 USE data_db
 IMPLICIT NONE
 TYPE (shrev), POINTER :: eset
 INTEGER :: ielem,g,iv,is,i,ng
 REAL, POINTER :: vargs(:,:)
 CHARACTER (len=34) :: gpname

 eset => shrev_head          !point to first set

 DO                          !loop for all the sets
   gpname = 'GP'//eset%sname   !gauss point name string "GPsname"
   iv = 0  !initializes pointer to Gauss values
   ng = eset%ngaus

   IF(shrev_force > 1)THEN

     WRITE(13,"('Result ""Beam_N"" ""',a10,'"" ',e12.4,  &
      &             ' Scalar OnGaussPoints  ""',a,'""')") &
      &             loadstep,ttime,TRIM(gpname)
     WRITE(13,"('Values')")
     iv = iv+1
     is = 0
     DO ielem=1,eset%nelem
       DO g=1,ng-1
         is = is + 1
         WRITE(13,"(i8,(e15.5))")is, &
                 (eset%elvar(iv,i,ielem), i=g,g+1)
       END DO
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shrev_shear > 1)THEN
     WRITE(13,"('Result ""Beam_Q"" ""',a10,'"" ',e12.4,  &
      &             ' Scalar OnGaussPoints  ""',a,'""')") &
      &             loadstep,ttime,TRIM(gpname)
     WRITE(13,"('Values')")
     iv = iv+1
     is = 0
     DO ielem=1,eset%nelem
       DO g=1,ng-1
         is = is + 1
         WRITE(13,"(i8,(e15.5))")is, &
                 (eset%elvar(iv,i,ielem), i=g,g+1)
       END DO
     END DO
     WRITE(13,"('End Values')")
   END IF

   IF(shrev_momen > 1)THEN
     WRITE(13,"('Result ""Beam_M"" ""',a10,'"" ',e12.4,  &
      &             ' Scalar OnGaussPoints  ""',a,'""')") &
      &             loadstep,ttime,TRIM(gpname)
     WRITE(13,"('Values')")
     iv = iv+1
     is = 0
     DO ielem=1,eset%nelem
       DO g=1,ng-1
         is = is + 1
         WRITE(13,"(i8,(e15.5))")is, &
                 (eset%elvar(iv,i,ielem), i=g,g+1)
       END DO
     END DO
     WRITE(13,"('End Values')")
   END IF

   WRITE(13,"('Result ""Beam_u"" ""',a10,'"" ',e12.4,  &
    &             ' Scalar OnGaussPoints  ""',a,'""')") &
    &             loadstep,ttime,TRIM(gpname)
   WRITE(13,"('Values')")
   iv = iv+1
   is = 0
   DO ielem=1,eset%nelem
     DO g=1,ng-1
       is = is + 1
       WRITE(13,"(i8,(e15.5))")is, &
               (eset%elvar(iv,i,ielem), i=g,g+1)
     END DO
   END DO
   WRITE(13,"('End Values')")

   WRITE(13,"('Result ""Beam_v"" ""',a10,'"" ',e12.4,  &
    &             ' Scalar OnGaussPoints  ""',a,'""')") &
    &             loadstep,ttime,TRIM(gpname)
   WRITE(13,"('Values')")
   iv = iv+1
   is = 0
   DO ielem=1,eset%nelem
     DO g=1,ng-1
       is = is + 1
       WRITE(13,"(i8,(e15.5))")is, &
               (eset%elvar(iv,i,ielem), i=g,g+1)
     END DO
   END DO
   WRITE(13,"('End Values')")

   eset => eset%next
   IF( .NOT.ASSOCIATED (eset) )EXIT
 END DO

! !        print Nodal variables
!
! IF( shrev_nvarn == 0 )RETURN    !no variables to print, => exit
!
! vargs => shrev_vargs
!
! IF(MOD(shrev_force,2) == 1)THEN
!
!   WRITE(13,"('Result ""N"" ""',a10,'"" ',e12.4,  &
!    &             ' Scalar OnNodes')") &
!    &             loadstep,ttime
!   WRITE(13,"('Values')")
!   DO i=1,ngp
!     WRITE(13,"(i8,(e15.5))")labea(i),vargs(1,i)
!   END DO
!   WRITE(13,"('End Values')")
! END IF
!
! IF(MOD(shrev_shear,2) == 1)THEN
!   WRITE(13,"('Result ""Q"" ""',a10,'"" ',e12.4,  &
!    &             ' Scalar OnNodes')") &
!    &             loadstep,ttime
!   WRITE(13,"('Values')")
!   DO i=1,ngp
!     WRITE(13,"(i8,(e15.5))")labea(i),vargs(2,i)
!   END DO
!   WRITE(13,"('End Values')")
! END IF
!
! IF(MOD(shrev_momen,2) == 1)THEN
!   WRITE(13,"('Result ""M"" ""',a10,'"" ',e12.4,  &
!    &             ' Scalar OnNodes')") &
!    &             loadstep,ttime
!   WRITE(13,"('Values')")
!   DO i=1,ngp
!     WRITE(13,"(i8,(e15.5))")labea(i),vargs(3,i)
!   END DO
!   WRITE(13,"('End Values')")
! END IF
!
!
!   WRITE(13,"('Result ""u"" ""',a10,'"" ',e12.4,  &
!    &             ' Scalar OnNodes')") &
!    &             loadstep,ttime
!   WRITE(13,"('Values')")
!   DO i=1,ngp
!     WRITE(13,"(i8,(e15.5))")labea(i),vargs(4,i)
!   END DO
!   WRITE(13,"('End Values')")
!
!
!
!   WRITE(13,"('Result ""v"" ""',a10,'"" ',e12.4,  &
!    &             ' Scalar OnNodes')") &
!    &             loadstep,ttime
!   WRITE(13,"('Values')")
!   DO i=1,ngp
!     WRITE(13,"(i8,(e15.5))")labea(i),vargs(5,i)
!   END DO
!   WRITE(13,"('End Values')")
!

 RETURN
 END SUBROUTINE prgi09c
