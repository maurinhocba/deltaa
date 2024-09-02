 SUBROUTINE stiffm (istop,ttime,nullm,iforce,stiff,force,nsymm, &
                    ustif,flag)
 !********************************************************************
 !
 !***   evaluation of tangent stiffness matrix
 !
 !********************************************************************
 USE npo_db, ONLY : coora,coor1,locsy,locs1
 USE ctrl_db, ONLY : numct,npoin,maxa,neulr
 USE loa_db, ONLY : ifloa
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: nullm
 INTEGER (kind=4),INTENT(IN) :: iforce
 INTEGER (kind=4),INTENT(IN OUT) :: istop,nsymm
 REAL (kind=8),INTENT(IN) :: ttime
 REAL (kind=8),INTENT(OUT) :: stiff(:),force(:,:),ustif(:)
 LOGICAL, INTENT(IN), OPTIONAL :: flag

 !local variables
 INTEGER(kind=4) :: i
 LOGICAL :: fc

 INTERFACE
   INCLUDE 'timing.h'
   INCLUDE 'elemnt.h'
   INCLUDE 'contac.h'
   INCLUDE 'loadst.h'
 END INTERFACE

 !     zero symmetric and asymmetric components of stiffness matrices
 fc = .NOT.ASSOCIATED(coor1)
 IF( fc )THEN   !default values for coordinates
   coor1 => coora
   locs1 => locsy
 ELSE IF( PRESENT(flag) )THEN
   IF(flag) THEN
     DO i=1,npoin
       coor1(:,i) = coora(:,i)
     END DO
     IF( neulr > 0 )THEN
       DO i=1,npoin
         locs1(:,i) = locsy(:,i)
       END DO
     END IF
   END IF
 END IF

 CALL timing(5,.TRUE.)
 IF(nullm)THEN
   stiff = 0d0
   IF(ABS(nsymm) == 1)THEN
     ustif = 0d0
     nsymm = 1
   END IF
 END IF

 !***   Stiffness matrices due to structural elements

 CALL elemnt('STIFFM', ttime=ttime, stiff=stiff, ustif=ustif, iforce=iforce, istop=istop, &
              resid=force )

 CALL timing(5,.FALSE.)


 IF(numct > 0) THEN
   CALL timing(7,.TRUE.)
 !  CALL contac('STIFFM',ttime=ttime,nsymm=nsymm,force=force(:,iforce), &
 !               stiff=stiff,ustif=ustif)
   CALL timing(7,.FALSE.)
 END IF

 ! follower forces
 CALL timing(6,.TRUE.)
 IF (ifloa > 0) CALL loadst(ttime,nsymm,force(:,iforce),stiff,ustif,coor1)
 CALL timing(6,.FALSE.)

 !     adds symmetric + asymmetric components
 !IF(nsymm == 1) THEN
 !  CALL timing(5,.TRUE.)
 !  IF(ANY(ustif /= 0d0))THEN
 !    DO i=1,maxa
 !      stiff(i) = stiff(i) + ustif(i)      !lower triangle
 !      ustif(i) = stiff(i) - 2d0*ustif(i)  !upper triangle
 !    END DO
 !  ELSE
 !    nsymm = -1
 !  END IF
 !  CALL timing(5,.FALSE.)
 !END IF

 IF( fc ) NULLIFY( coor1, locs1 )

 RETURN

 END SUBROUTINE stiffm
