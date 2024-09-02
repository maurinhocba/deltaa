 SUBROUTINE actcd( )

 !keeps the coordinates at the begining of the stage ==> COORS
 !applies initial conditions on displacements and local systems

 USE npo_db, ONLY : coora,coors,coorc,psia,psic,euler,eulei,naeul,cpx,locsy
 USE ctrl_db, ONLY: ndime, npoin,ndofn, initial_displacements, initial_rotations
 IMPLICIT NONE

 INTERFACE
   INCLUDE 'inrotm.h'
 END INTERFACE

 !Local variable
 INTEGER (kind=4) :: i
 REAL(kind=8) :: dlb(3,3),lb(3,3)

 IF( initial_displacements )THEN
   IF( ASSOCIATED(cpx))THEN   !incremental displacement at control points
     DO i=1,npoin
       IF(cpx(1,i) == 0 )CYCLE
       coorc(:,i) = 2d0*coorc(:,i) -(coorc(:,cpx(2,i)) + coorc(:,cpx(3,i)))/2d0
     END DO
   END IF
   DO i=1,npoin
     coora(:,i) = coorc(:,i) + coora(:,i)  !incremental displacement at the begining of the strategy
   END DO
   IF( ndofn == 8 .OR. ndofn == 4 )THEN
     DO i=1,npoin
       psic(:,i) = psic(:,i) + psia(:,i)   !incremental displacement at the begining of the strategy
     END DO
   END IF
   initial_displacements = .FALSE.
 END IF
 IF( initial_rotations )THEN
   DO i=1,npoin
     IF(.NOT.naeul(i))CYCLE
     IF(ANY(eulei(:,i) /= 0d0))THEN
       IF( ndime == 2 )THEN
         !euler(1,i) = euler(1,i) +eulei(1,i)
         locsy(1,i) = euler(1,i) +eulei(1,i)
       ELSE
         CALL inrotm(eulei(:,i),dlb)        !compute incremental rotation matrix
         lb = RESHAPE(euler(:,i),(/3,3/))   !recover previous rotation matrix
         lb = MATMUL(lb,dlb)                !new rotation matrix
         !euler(:,i) = RESHAPE( lb, (/ 9 /))
         locsy(:,i) = RESHAPE( lb, (/ 9 /))  !changed 18/9/20 
       END IF
     END IF
   END DO
   initial_rotations = .FALSE.
   DEALLOCATE( eulei )
 END IF
 IF (ASSOCIATED(coors)) DEALLOCATE(coors)
 ALLOCATE(coors(ndime,npoin))
 DO i=1,npoin
   coors(:,i) = coora(:,i)
   coorc(:,i) = coora(:,i)  !no incremental displacement at the begining of the strategy
   IF( ndofn == 8 .OR. ndofn == 4 ) psia(:,i) = psic(:,i)      !add displacement at the begining of the strategy
 END DO
 RETURN
 END SUBROUTINE actcd
