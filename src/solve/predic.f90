 SUBROUTINE predic(istep,lauto,ecdis,neq,arcln,dlamb,disax,displ,   &
                   ddisp,delta,karcl,piter,diter,newtv,ncdis)
 !***********************************************************************
 !
 !*** this routine predicts displacement and load step increments
 !    according to selected path and previous increments
 !
 !***********************************************************************
 USE ctrl_db, ONLY : ndime,npoin
 USE npo_db, ONLY : ifpre
 IMPLICIT NONE
 !       routine arguments
 INTEGER (kind=4),INTENT(IN) :: istep,neq,lauto,karcl,newtv
 INTEGER (kind=4),INTENT(IN OUT) :: ecdis,ncdis
 REAL (kind=8),INTENT(IN) :: displ(:),piter,diter
 REAL (kind=8),INTENT(IN OUT) :: dlamb,arcln,disax(:,:),ddisp(:),delta
 !       local variables
 INTEGER (kind=4) :: ecphi,neq1,i,j,n
 REAL    (kind=8) :: a,b,c,signr,d(16),d1(4,4),s1,s2,s3,s,st

 INTERFACE
   INCLUDE 'invmtx.h'
 END INTERFACE

 neq1 = neq+1
 !     New equation for control displacement
 IF(karcl == 4 .AND. MOD(lauto,2) == 1 .AND. istep > 1)THEN
   ecphi = ecdis
   a = 0d0   !initializes the largest displacement increment
   DO n=1,npoin        !search for the node and DOF
     DO j=1,ndime
       i = ifpre(j,n)            !associated DOF
       IF(i > 0)THEN             !If an active DOF
         IF( ABS(displ(i)) >= a ) THEN        !if displacement is greater
           a = ABS(displ(i))                  !updates maximum Translational DOF
           ncdis = 10*n+j                     !updates node and DOF
           ecdis = i                          !updates ACTIVE DOF
         END IF
       END IF
     END DO
   END DO
   arcln = displ(ecdis)    !new ARC-LENGTH
   IF(ecdis /= ecphi)WRITE(55,"(' change TO',i8,i3)")ncdis/10,MOD(ncdis,10)
 END IF

 IF( lauto <= 1 .OR. istep == 1)THEN
   !     first step or standard prediction
   IF( karcl == 0 .OR. karcl == 5) THEN   ! load or tools control
     delta = dlamb
     IF(lauto == 1) delta = delta*SQRT(diter/piter) ! variable load_inc
   ELSE IF(ABS(arcln) > 0d0) THEN                   ! exists length
     IF(lauto == 1) arcln = arcln*SQRT(diter/piter) ! variable arc_length
     SELECT CASE (karcl)
     CASE (1:3)           ! arc_lenght methods
       IF(istep == 1) THEN                          ! displ is unknown
         signr = 1d0
       ELSE                                         ! previous displ
         signr = DOT_PRODUCT(disax(1:neq,1),displ)
         signr = signr/ABS(signr)
       END IF
       delta = signr*arcln/ SQRT(DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1)))
     CASE (4)             ! displ. control
       delta = (arcln-ddisp(ecdis))/disax(ecdis,1)
     END SELECT  ! karcl == ?
   ELSE ! length = 0.0d0
     delta = dlamb
     SELECT CASE (karcl)
     CASE (1:3)        ! arc_length methods
       arcln=SQRT(DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1)))*dlamb
     CASE (4)            ! displ. control
       arcln = disax(ecdis,1)*dlamb
     END SELECT
   END IF ! ABS(arcln) > 0d0
   !     keep a vector for comparison & compute ddisp
   IF(istep == 1) THEN
     disax(1:neq,2) = disax(1:neq,1)*delta          ! keep tangent vector
     ddisp = ddisp + delta * disax(1:neq,1)
   ELSE
     disax(1:neq,2) = displ                         ! keep last increment
     IF(karcl == 5) THEN                            ! fixed tools
       a = 0.90d0
       b = (1d0-a)*delta
       IF(lauto == 1) a = a*SQRT(diter/piter)       ! variable arc_length
       ddisp = ddisp + displ * a + disax(1:neq,1) * b
     ELSE
       ddisp = ddisp + delta * disax(1:neq,1)
     END IF
   END IF

 ELSE
   !     prediction based on previous increments
   disax(1:neq1,7) = disax(1:neq1,6)
   disax(1:neq1,6) = disax(1:neq1,5)
   disax(1:neq,5) = displ
   disax(neq+1,5) = dlamb
   IF( lauto == 3) arcln = arcln*SQRT(diter/piter)
   IF( lauto == 3) dlamb = dlamb*SQRT(diter/piter)
   SELECT CASE (karcl)
   CASE (0,5)
     s  = dlamb
     s1 = disax(neq+1,5)
     s2 = disax(neq+1,6)
     s3 = disax(neq+1,7)
     st = 1d0
   CASE (1:3)
     s  = arcln
     s1 = SQRT(DOT_PRODUCT(disax(1:neq,5),disax(1:neq,5)))
     s2 = SQRT(DOT_PRODUCT(disax(1:neq,6),disax(1:neq,6)))
     s3 = SQRT(DOT_PRODUCT(disax(1:neq,7),disax(1:neq,7)))
     st = SQRT(DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1)))
     disax(1:neq,2) = displ                         ! keep last increment
   CASE (4)
     s  = arcln
     s1 = disax(ecdis,5)
     s2 = disax(ecdis,6)
     s3 = disax(ecdis,7)
     st = disax(ecdis,1)
   END SELECT
   c = s+s1
   IF( s2 == 0d0) THEN                             !second step
     IF( newtv == 0 ) THEN
       a = s/s1
       ddisp(1:neq) = ddisp(1:neq) + disax(1:neq,5) * a
       delta        = disax(neq+1,5)*a
     ELSE
       a = s*(c/s1)/st
       b = - (s/s1)**2
       ddisp = ddisp + disax(1:neq,1)*a + disax(1:neq,5)*b
       delta = a + disax(neq+1,5)*b
     END IF
   ELSE IF (s3 == 0d0 .AND. newtv == 0) THEN  !third step
     a = -c*s/(s2*(s1+s2))
     b =  c*(c+s2)/(s1*(s1+s2)) - 1d0
     ddisp = ddisp + disax(1:neq,6)*a + disax(1:neq,5)*b
     delta = disax(neq+1,6)*a + disax(neq+1,5)*b
   ELSE                                   !default step
     IF(newtv == 0) THEN
       a = -(s3+s2)
       b = -s2
       d = (/ 1d0,   a, a**2, a**3, 1d0,  b,  b**2, b**3,   &
              1d0, 0d0,  0d0,  0d0, 1d0, s1, s1**2, s1**3 /)
       CALL invmtx ( d, d1, signr  ,4)
       d(5:8) = (/ 1d0, c, c**2, c**3 /)
       d(1:4) = MATMUL( d1 , d(5:8))
       a = - d(1)
       b = - d(1) - d(2)
       c =   d(4) - 1d0
       ddisp = ddisp + disax(1:neq,7)*a+disax(1:neq,6)*b +disax(1:neq,5)*c
       delta = disax(neq+1,7)*a+disax(neq+1,6)*b+disax(neq+1,5)*c
     ELSE
       b = -s2
       d = (/ 0d0, 1d0, 2d0*s1, 3d0*s1**2, 1d0,  b,  b**2, b**3,    &
              1d0, 0d0,    0d0,       0d0, 1d0, s1, s1**2, s1**3 /)
       CALL invmtx ( d, d1, signr  ,4)
       d(5:8) = (/ 1d0, c, c**2, c**3 /)
       d(1:4) = MATMUL( d1 , d(5:8))
       a =  d(1)/st
       b = -d(2)
       c =  d(4) - 1d0
       ddisp = ddisp + disax(1:neq,1)*a+disax(1:neq,6)*b + disax(1:neq,5)*c
       delta  = a + disax(neq+1,6)*b + disax(neq+1,5)*c
     END IF
   END IF
 END IF
 IF(karcl == 3)arcln = SQRT(DOT_PRODUCT(ddisp,ddisp))
 RETURN
 END SUBROUTINE predic
