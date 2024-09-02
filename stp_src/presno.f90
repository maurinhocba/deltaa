 SUBROUTINE presno( )
 !
 !  Compute pressure due to Blank_holder and distance to Tools
 !
 USE cont_db
 USE data_db
 IMPLICIT NONE

 INTEGER(kind=4) :: isurf,ncnod,nsegm,i,j,k,ll
 REAL (kind=8) :: a3,t1(3),t2(3),t3(3)
 TYPE (surf_db), POINTER :: surf                      ! CONTACT surface
 INTEGER (kind=4), POINTER :: lcseg(:,:),lcnod(:)     ! connectivities and list of nodes
 REAL (kind=8), ALLOCATABLE :: vector(:), & !(ncnod)
                               ar(:,:),   & !(2,npoin)
                               al(:)        !(npoin)
 REAL (kind=8), PARAMETER  :: pi=3.141592653589793

 INTERFACE
   INCLUDE 'cab_gid.h'
   INCLUDE 'cab_gid_bin.h'
 END INTERFACE

 surf => shead                                     !point to first surface
 IF(wpress) ALLOCATE ( ar(2,npoin))                !reserve memory, for every node 2 values 1)?  2)?
 ar = 0d0                                          !initializes
 IF(wpress) presn = 0d0                            !initializes  nodal presss
 IF(wwrink) wrinkn= 0d0                            !initializes
 DO isurf=1,nsurf                                  !for each surface
   ncnod = surf%ncnod                              !number of nodes in the surf
   nsegm = surf%nsegm                              !number of segments in surfa
   lcseg => surf%lcseg                             !connectivities
   lcnod => surf%lcnod                             !list of nodes
   ll = 1                                          !binder pressure
   IF( surf%cpress ) ll = 2                        !contact pressure
   IF( surf%press .AND. wpress )THEN               !if press desired
     ALLOCATE(al(npoin))
     al = 0d0                                      !initializes array
     DO i=1,nsegm
       t1(1:ndime) = coorf(:,lcseg(2,i)) - coorf(:,lcseg(1,i))   !side 3
       IF( ndime == 3 )THEN
         t2 = coorf(:,lcseg(3,i)) - coorf(:,lcseg(1,i))   !-side 2
         t3(1) = t1(2)*t2(3) - t1(3)*t2(2)                !triangle normal
         t3(2) = t1(3)*t2(1) - t1(1)*t2(3)
         t3(3) = t1(1)*t2(2) - t1(2)*t2(1)
         a3 = SQRT(DOT_PRODUCT(t3,t3))/6d0                !area/3
       ELSE
         a3 = SQRT(DOT_PRODUCT(t1(1:ndime),t1(1:ndime)))/2d0                !length/2
         IF( ntype == 3 ) a3 = a3*pi*(coorf(1,lcseg(2,i)) + coorf(1,lcseg(1,i)))
       END IF
       al(lcseg(:,i)) = al(lcseg(:,i)) + a3             !add no nodal area
     END DO
   END IF
   ALLOCATE( vector(ncnod) )                       !reserve memory
   IF( surf%press )THEN                            !if press written
     IF( wpress )THEN                              !if press desired
       READ(44) (vector(i),i=1,ncnod)              !read data
       DO i=1,ncnod                                !for each node
         j = lcnod(i)                              !associated node
         presn(ll,j) = presn(ll,j) + vector(i)/al(j)     !pressure
       END DO
       ar(ll,:) = ar(ll,:) + al
       DEALLOCATE(al)
     ELSE
       READ(44)                                    !skip data
     END IF
   END IF
   IF( surf%wrink )THEN                            !if gaps written
     DO k=1,2                                      !for each value
       IF( wwrink )THEN                            !if gaps desired
         READ(44) (vector(i),i=1,ncnod)            !read data
         DO i=1,ncnod                              !for each node
           j = lcnod(i)                            !associated node
           IF( wrinkn(k,j) == 0d0) THEN            !no previous value
             wrinkn(k,j) = vector(i)
           ELSE                                    !previous value exist
             wrinkn(k,j) = MIN(wrinkn(k,j),vector(i))  !take minimum value
           END IF
         END DO
       ELSE
         READ(44)                                  !skip data
       END IF
     END DO
   END IF
   surf => surf%next                               !point to next value
   DEALLOCATE( vector )                            !release memory
 END DO

 IF( wpress )THEN          !if press desired
   IF( ip == 2 )THEN       !for GiD
     IF( ANY(ar(1,:) > 0d0) )THEN  !binder press
       CALL cab_gid(1,1,press_l(2),press_l(2),1,loadstep,ttime,units=press_u)
       DO i=1,npoin          !for each point
         IF(ar(1,i) > 0) WRITE(13,"(i8,e13.5)")label(i),presn(1,i)*press_f
       END DO
       WRITE(13,"('End Values')")
     END IF
     IF( ANY(ar(2,:) > 0d0 ))THEN  !contact press
       CALL cab_gid(1,1,press_l(3),press_l(3),1,loadstep,ttime,units=press_u)
       DO i=1,npoin          !for each point
         IF(ar(2,i) > 0) WRITE(13,"(i8,e13.5)")label(i),presn(2,i)*press_f
       END DO
       WRITE(13,"('End Values')")
     END IF
   ELSE IF (ip == 4 .OR. ip == 5) THEN      !for GiD binary
     IF( ANY(ar(1,:) > 0d0) )THEN  !binder press
       CALL cab_gid_bin(1,1,press_l(2),press_l(2),1,loadstep,ttime,units=press_u)
       DO i=1,npoin          !for each point
         IF (ar(1,i) > 0) CALL GID_WRITESCALAR(label(i),presn(1,i)*press_f )
       END DO
       CALL GID_ENDRESULT()
     END IF
     IF( ANY(ar(2,:) > 0d0 ))THEN  !contact press
       CALL cab_gid_bin(1,1,press_l(3),press_l(3),1,loadstep,ttime,units=press_u)
       DO i=1,npoin          !for each point
         IF (ar(2,i) > 0) CALL GID_WRITESCALAR(label(i),presn(2,i)*press_f )
       END DO
       CALL GID_ENDRESULT()
     END IF
   END IF
   DEALLOCATE( ar)   !release memory
 END IF

 IF( wwrink )THEN          !if gaps information desired
   IF( ip == 2 )THEN       !for GiD
 !    CALL cab_gid(1,2,wrink_l(1),wrink_l(2),3,loadstep,ttime,units=wrink_u)
     CALL cab_gid(2,1,wrink_l(1),wrink_l(2),3,loadstep,ttime,units=wrink_u)   ! scalar ???
     DO i=1,npoin         !for each Point
       WRITE(13,"(i8,3e13.5)")label(i),wrinkn(1:2,i)*wrink_f ,(wrinkn(1,i)-wrinkn(2,i))*wrink_f
     END DO
     WRITE(13,"('End Values')")
   ELSE IF( ip == 4 .OR. ip == 5 )THEN       !for GiD binary
     CALL cab_gid_bin(2,1,wrink_l(1),wrink_l(2),3,loadstep,ttime,units=wrink_u)
     DO i=1,npoin          !for each point
       !IF (ar(i) > 0) CALL GID_WRITEVECTOR(label(i),wrinkn(1,i)*wrink_f ,wrinkn(2,i)*wrink_f ,   &
       CALL GID_WRITEVECTOR(label(i),wrinkn(1,i)*wrink_f ,wrinkn(2,i)*wrink_f ,   &
         (wrinkn(1,i)-wrinkn(2,i))*wrink_f )
     END DO
     CALL GID_ENDRESULT()
   END IF
 END IF

 RETURN
 END SUBROUTINE presno
