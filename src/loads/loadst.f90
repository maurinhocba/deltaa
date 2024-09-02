 SUBROUTINE loadst (ttime,nsymm,force,stiff,ustif,coora)

 ! calculates and applies follower load

 USE npo_db, ONLY : loass
 USE loa_db, ONLY : loa_set,headls
 USE ctrl_db, ONLY : ndime,ntype,nload

 !USE hydf_db
 IMPLICIT NONE
 ! dummy arguments
   REAL (kind=8), INTENT(IN) :: ttime,     & !to compute load factor
                                force(:),  & !assembled forces
                                stiff(:),  & !stiffness matrix (symmetric)
                                ustif(:),  & !stiffness matrix (antisymmetric)
                                coora(:,:)   !present coordinates

 ! Local
 INTEGER (kind=4) :: il
 TYPE (loa_set), POINTER :: loas
   INTEGER(kind=4):: nsymm    !0: symmetric matrix only, 1:include unsymmetric part

 INTERFACE
   INCLUDE 'folost.h'
 END INTERFACE
 !---------------------------------------------------------------
 !for 2-D problems in plane stress/strain, symmetric matrix is null
 IF( nsymm == 0 .AND. ndime == 2 .AND. ntype /= 3 )RETURN

 loas => headls  ! points to first set of loads
 DO il = 1,nload ! loop over load sets
   ! if set contains follower load data
   IF (loas%numfl > 0) THEN ! if the set has follower loads
     CALL folost(ttime,nsymm,force,stiff,ustif,coora,loas%headf,loas%factor,loass(il))

   END IF
   loas => loas%next
 END DO

 END SUBROUTINE loadst
