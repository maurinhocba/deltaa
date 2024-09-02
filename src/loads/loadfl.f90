 SUBROUTINE loadfl (istop) ! (ttime,dtime)

 ! calculates and applies follower load

 USE ctrl_db, ONLY: ndime,ndofn,nload,ntype,npoin,istep,nstep,neq,ttime,dtime
 USE loa_db
 USE npo_db
 USE hydf_db
 IMPLICIT NONE
 INTEGER(kind=4), INTENT(OUT) :: istop
 !REAL (kind=8), INTENT(IN) :: ttime,dtime
 INTERFACE
    INCLUDE 'aplfol.h'
 END INTERFACE

 ! Local
 INTEGER (kind=4) :: il,i !,j,k
 TYPE (loa_set), POINTER :: loas
 !TYPE(foll_seg),POINTER:: seg
 !---------------------------------------------------------------


 loas => headls  ! points to first set of loads
 DO il = 1,nload ! loop over load sets
   ! if set contains follower load data
   IF (loas%numfl > 0) THEN ! if the set has follower loads
     force(1:neq,il) = 0d0  !initializes assembled load vector
     ! assembles conservative load if present
     IF(loas%iplod + loas%nedge + loas%nsurf + loas%igrav > 0)THEN
       DO i=1,npoin         !for each node in the mesh
         ! assemble into force vectors
         IF(ANY(loadv(1:ndofn,i,il) /= 0d0)) & !only if load exists
         CALL ensvec(ndofn,ifpre(1,i),loadv(1,i,il),force(1,il))
       END DO
     END IF
     ! compute and assemble follower loads
     CALL aplfol(ndime,coora,force(1:neq,il),ntype,loas%headf)

   END IF
   !IF (ASSOCIATED (loas%hydfl) ) & ! if set contains hydroforming data
   !  CALL hydflo (loas%hydfl, loass(il), loas%factor, ttime)
   loas => loas%next
 END DO

 END SUBROUTINE loadfl
