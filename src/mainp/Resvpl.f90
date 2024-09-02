 SUBROUTINE resvpl (istop,ttime,resid,gvect,flag)
 !********************************************************************
 !
 !***   evaluation of integral (b)**t*(sigma)
 !
 !********************************************************************
 USE ctrl_db, ONLY : ndime,ndofn,npoin,numct,neulr
 USE npo_db, ONLY : coora,coor1,locsy,locs1,ifpre,fcont,id
 USE loa_db, ONLY : ifloa
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN OUT) :: istop
 REAL (kind=8),INTENT(IN) :: ttime
 REAL (kind=8),INTENT(IN OUT) :: gvect(:),resid(:,:)
 LOGICAL, INTENT(IN), OPTIONAL :: flag
 !local variables
 INTEGER (kind=4) :: i
 LOGICAL :: fc

 INTERFACE
   INCLUDE 'timing.h'
   INCLUDE 'elemnt.h'
   INCLUDE 'contac.h'
 END INTERFACE

 CALL timing(6,.TRUE.)

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

 gvect = 0d0     !initializes internal forces
 resid = 0d0

 !***   loop over all the elemnt sets
 CALL elemnt('RESVPL', resid=resid, istop=istop, ttime=ttime)
 CALL ensvec(ndofn*npoin,ifpre(1,1),resid(1,1),gvect(1))

 CALL timing(6,.FALSE.)

 ! compute contact forces
 IF(numct > 0) THEN
   CALL timing(7,.TRUE.)
   fcont = 0d0
   CALL contac('FORCES',ttime=ttime)
   CALL ensvec(ndime*npoin,id(1,1),fcont(1,1),gvect(1))
   CALL timing(7,.FALSE.)
 END IF

 CALL timing(6,.TRUE.)
 gvect = - gvect

 !update follower forces
 IF (ifloa > 0) CALL loadfl (istop) !(ttime,0d0)

 IF( fc ) NULLIFY( coor1, locs1 )
 CALL timing(6,.FALSE.)

 RETURN
 END SUBROUTINE resvpl
