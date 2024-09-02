 SUBROUTINE elemnt (TASK, name, dtime, ttime, istop, flag1, flag2, &
                    sumat, iload, iforce, mass, stiff, ustif, resid)
 !********************************************************************
 !
 !*** standard ELEMENT routine
 !
 !********************************************************************
 USE esets_db   !element data base
 USE ctrl_db, ONLY : top,bottom,npoin,neulr
 USE npo_db, ONLY : coor1,coord,euler,locs1,coort,coorb,ifact, &    !point information
                    emass
 USE surf_db    !to get a surface from element data
 USE outp_db, ONLY : iwrit

 IMPLICIT NONE
 !                        Dummy (compulsory) arguments
 CHARACTER(len=*),INTENT(IN):: task
 !                        Dummy (optional) arguments
 CHARACTER(len=*),OPTIONAL:: name
 INTEGER (kind=4), OPTIONAL :: istop,iload,iforce
 REAL (kind=8), OPTIONAL :: dtime,ttime,sumat,resid(:,:),stiff(:),ustif(:),mass(:)
 LOGICAL, OPTIONAL :: flag1,flag2

 !      Local variables
 LOGICAL :: fc
 INTEGER (kind=4) :: iset,i

 INTERFACE
   INCLUDE 'elemt1.h'
   INCLUDE 'elemt2.h'
   INCLUDE 'elem03.h'
   INCLUDE 'elem04.h'
   INCLUDE 'elem05.h'
   INCLUDE 'elemt6.h'
   INCLUDE 'elemt7.h'
   INCLUDE 'elemt8.h'
   INCLUDE 'elemt9.h'
   INCLUDE 'elem10.h'
   INCLUDE 'elem11.h'
   INCLUDE 'elem12.h'
   INCLUDE 'elem13.h'
   INCLUDE 'elem14.h'
   INCLUDE 'elem15.h'
   INCLUDE 'elem16.h'
   INCLUDE 'elem17.h'
   INCLUDE 'elem18.h'
   INCLUDE 'elem19.h'
   INCLUDE 'elem20.h'
   INCLUDE 'elem24.h'
   INCLUDE 'elem25.h'
   INCLUDE 'elem26.h'
   INCLUDE 'elem27.h'
   INCLUDE 'elem29.h'
   INCLUDE 'elem30.h'
 END INTERFACE

 fc = ASSOCIATED(coor1)
 IF( task == 'GAUSSC' .OR. task == 'MASMTX' .OR. task == 'RIGBDY' )THEN
   IF(fc)THEN
     DO i=1,npoin
       coor1(:,i) = coord(:,i)
     END DO
     IF( neulr > 0 )THEN
       DO i=1,npoin
         locs1(:,i) = euler(:,i)
       END DO
     END IF
   ELSE
     coor1 => coord
     locs1 => euler
   END IF
 END IF

 IF( task == 'RESVPL' )THEN
   IF( bottom ) coorb = 0d0
   IF( top    ) coort = 0d0
   IF( top .OR. bottom   )ifact = 0
 END IF

 IF (PRESENT(flag2) .AND. TASK /= 'OUTPUT') flag2=.FALSE.

 DO iset=1,esets

   SELECT CASE ( eset(iset) )
   CASE (1)
     CALL elemt1(TASK, ttime,iload,iforce,                   &
                 coor1,locs1,sumat,mass,emass,resid,stiff,   &
                 istop,flag1,flag2,iwrit)

   CASE (2)
     CALL elemt2(TASK, ttime,iload,iforce,           &
                 coor1,sumat,mass,emass,resid,stiff, &
                 istop,flag1,flag2,iwrit)

   CASE (3)
     CALL elem03(TASK, ttime,iload,iforce,           &
                 coor1,locs1,sumat,mass,emass,resid, &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (4)
     CALL elem04(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,ustif,istop,flag1,flag2,iwrit)

   CASE (5)
     CALL elem05(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (6)
     CALL elemt6(TASK, ttime,iload,iforce,           &
                 coor1,locs1,sumat,mass,emass,resid, &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (7)
     CALL elemt7(TASK, ttime,iload,iforce,           &
                 coor1,locs1,sumat,mass,emass,resid, &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (8)
     CALL elemt8(TASK, ttime,iload,iforce,           &
                 coor1,locs1,sumat,mass,emass,resid, &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (9)
     CALL elemt9(TASK, ttime,iload,iforce,           &
                 coor1,locs1,sumat,mass,emass,resid, &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (10)
     CALL elem10(TASK, name, ttime,iload,iforce,     &
                 coor1,locs1,sumat,mass,emass,       &
                 istop,flag1,flag2,iwrit)

   CASE (11)
     CALL elem11(TASK, ttime,iload,iforce,           &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (12)
     CALL elem12(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,ustif,istop,flag1,flag2,iwrit)

   CASE (13)
     CALL elem13(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (14)
     CALL elem14(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (15)
     CALL elem15(TASK, ttime,iload,iforce,           &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (16)
     CALL elem16(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (17)
     CALL elem17(TASK, ttime,iload,iforce,           &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,ustif,istop,flag1,flag2,iwrit)

   CASE (18)
     CALL elem18(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,ustif,istop,flag1,flag2,iwrit)

   CASE (19)
     CALL elem19(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (20)
     CALL elem20(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (24)
     CALL elem24(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (25)
     CALL elem25(TASK, ttime,iload,iforce,           &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (26)
     CALL elem26(TASK, name,ttime,iload,iforce,      &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,ustif,istop,flag1,flag2,iwrit)

   CASE (27)
     CALL elem27(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (29)
     CALL elem29(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   CASE (30)
     CALL elem30(TASK, name, ttime,iload,iforce,     &
                 coor1,sumat,mass,emass,resid,       &
                 stiff,istop,flag1,flag2,iwrit)

   END SELECT
 END DO

 IF( task == 'RESVPL' )THEN
   IF( bottom ) THEN
     DO i=1,npoin
       IF(ifact(i) > 0)THEN
        coorb(:,i) = coorb(:,i)/ifact(i) + coor1(:,i)
       ELSE
        coorb(:,i) = coor1(:,i)
       END IF
     END DO
   END IF
   IF( top ) THEN
     DO i=1,npoin
       IF(ifact(i) > 0)THEN
        coort(:,i) = coort(:,i)/ifact(i) + coor1(:,i)
       ELSE
        coort(:,i) = coor1(:,i)
       END IF
     END DO
   END IF
 END IF

 IF( .NOT.fc) NULLIFY(coor1,locs1)

 RETURN

 END SUBROUTINE elemnt
