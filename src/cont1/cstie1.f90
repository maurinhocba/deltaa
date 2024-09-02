 SUBROUTINE cstie1(rssdb,stiff,cprop,issdb,facto,nsymm,astif)

 !.... COMPUTE STIFFNESS MATRIX FOR A THREE-NODE 2-D CONTACT ELEMENT

 IMPLICIT NONE
 !       arguments
 INTEGER (kind=4), INTENT(IN) :: nsymm, & !1: asymmetric matrices allowed  0: symmetric only
                                 issdb(:) !node slave data base
 REAL (kind=8), INTENT(IN) :: rssdb(:), & !node slave data base
                              cprop(:), & !contact pair parameters
                              facto       !contact factor  2d0 maux/dtime**2
 REAL (kind=8), INTENT(OUT) :: stiff(:),& !element stiffness (symmetric part)
                               astif(:)   !element stiffness (asymmetric part)
 !       local variables
 LOGICAL frict,slip
 INTEGER (kind=4), PARAMETER :: nee = 6
 INTEGER (kind=4) i,j,k
 REAL (kind=8) gap,r,shp1,shp2,b(6),pnlty,vn1,vn2,vt1,vt2,     &
               bs(6),vnm(6),cofri,pnltc,iter         !,dxm,g_dxm

 !.... READ surface properties and normal direction

 pnlty = cprop(1) * facto
 iter = DBLE(issdb(4)+1)/4d0
 IF( iter < 1d0 ) pnlty = pnlty*iter

 !slip = (ABS(rssdb(5) - rssdb(4))  >  .1e-11)
 slip = issdb(3) == 2
 IF( slip )THEN  !sliding
   cofri = cprop(4)         !kinet
 ELSE                     !stuck
   cofri = cprop(3)         !static
 END IF

 IF(cofri > 0 ) THEN
   frict = .TRUE.
   pnltc = cprop(2) * facto
 ELSE
   frict = .FALSE.
 END IF

 gap = rssdb(1)
 vn1 = rssdb(2)
 vn2 = rssdb(3)
 r   = rssdb(4)
 vt1 =  vn2
 vt2 = -vn1

 !     FORM SHAPE FUNCTIONS
 shp1 = 1d0-r
 shp2 = r
 !dxm = rssdb(5)
 !g_dxm = gap/dxm

 !     FORM B OPERATOR
 b(1) = vn1
 b(2) = vn2
 b(3) = -shp1*vn1
 b(4) = -shp1*vn2
 b(5) = -shp2*vn1
 b(6) = -shp2*vn2
 !     FORM TS OPERATOR
 bs(1) = vt1
 bs(2) = vt2
 bs(3) = -shp1*vt1
 bs(4) = -shp1*vt2
 bs(5) = -shp2*vt1
 bs(6) = -shp2*vt2
 !     FORM NM OPERATOR
 vnm(1) = 0.0
 vnm(2) = 0.0
 vnm(3) = -vn1
 vnm(4) = -vn2
 vnm(5) =  vn1
 vnm(6) =  vn2
 !     compute stiffness matrix
 IF(nsymm == 1)astif(1:15) = 0d0
 k = 0
 DO j = 1,nee      !for each column
   DO i = j,nee      !for each row
     k = k+1
     ! Compute tangent operator for linearized kinematics
     stiff(k) = b(i)*b(j)*pnlty
     ! add extra terms to compute consistent tangent operator for
     ! fully non-linear kinematics
     !stiff(k) = stiff(k) - pnlty*g_dxm* &
     !  (bs(i)*vnm(j) + vnm(i)*bs(j) + g_dxm*vnm(i)*vnm(j))
     ! ADD TERMS DUE TO FRICTION
     if(frict) then
       IF(slip) THEN !is this OK, MUST BE REVISED
         IF(nsymm == 1) THEN
           stiff(k) = stiff(k) - (bs(i)*b(j)+b(i)*bs(j))/2d0*pnlty*cofri 
           astif(k-j) = -(bs(i)*b(j)-b(i)*bs(j))/2d0*pnlty*cofri
         ELSE
           stiff(k) = stiff(k) + bs(i)*bs(j)*pnltc*cofri
         END IF
       ELSE
         stiff(k) = stiff(k) + bs(i)*bs(j)*pnltc
       END IF
     END IF
   END DO
 END DO

 RETURN
 END SUBROUTINE cstie1
