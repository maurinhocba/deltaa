SUBROUTINE projt4(xs,nearn,x,lcseg,   &
                  nhseg,issdb,prdat,cutof,gapin,curms,cu)

!.... project the slave node onto a 3-d master triangular surface

IMPLICIT NONE
!     arguments
      LOGICAL, INTENT(IN) :: curms
INTEGER (kind=4), INTENT(IN) :: lcseg(:,:),nhseg(:,:)
INTEGER (kind=4), INTENT(IN OUT) :: nearn,issdb(:)
REAL (kind=8), INTENT (IN) :: x(:,:),cutof,gapin,xs(:),cu(:,:)
REAL (kind=8), INTENT (OUT) :: prdat(:)

!     local variables
INTEGER (kind=4) :: imn,jmn,kmn,nchck,loop,nnear,i,o(0:40)
REAL (kind=8) :: vd(3),vj(3),vk(3),vt(3),auxi(3),area2,  &
                 gap,dist,y(3),mgap,pgap,                &
                 xita,eta,xita1,eta1,zeta1,              &
                 vn(3),vnn(3),vns(3),gaps,bx,be,gapin0,cutof0
LOGICAL ::  chck1,chck2,sig
REAL (kind=8), PARAMETER :: toler = 0.01d0, tole1 = 1.01d0

!REAL (kind=8), SAVE :: maxg = -1d-6  !variable for debug
INTEGER (kind=4) :: nn(41,3)     !variable for statistic records
COMMON /proj3/ nn

INTERFACE
  INCLUDE 'projtr.h'
END INTERFACE

!.... initialize to zero some control parameters
pgap = MIN( prdat(9), 0d0)                  !previos penetration
gapin0 = MIN(gapin,pgap/2d0)                !maximum incremental gap
cutof0 = cutof + pgap/2d0                   !maximum penetration
cutof0 = MAX(cutof0, 3d0*cutof )            !maximum penetration
nchck = 0            !initializes number of candidates
o(0)  = nearn        !keep nearn
sig = .FALSE.        !initializes to NO full projection
mgap = 1d10          !initializes max gap found
vnn = 0d0            !initializes weighted normal
loop = 0             !initializes loop search
!.... loop until projection
search : DO
  !....   identify global master nodes
  imn = lcseg(1,nearn)        !1st master node
  jmn = lcseg(2,nearn)        !2nd master node
  kmn = lcseg(3,nearn)        !3rd master node
  vd = xs(1:3) - x(1:3,imn)          !distance vector
  dist = SQRT(DOT_PRODUCT(vd,vd))    !distance
  mgap = MIN(mgap,dist)              !for boundaries
  IF(dist == 0d0) THEN               !IF nodes coincide
    nchck = 0                        !no penetration
    sig = .TRUE.                     !keep the node
    mgap = 0d0                       !update max gap
    EXIT search                      !exit search
  END IF
  !....   form tangent vectors to master segment
  vj = x(1:3,jmn) - x(1:3,imn)       !side i-->j
  vk = x(1:3,kmn) - x(1:3,imn)       !side i-->k
  ! compute normal vector  Vn = Vj x Vk
  vn(1) = vj(2)*vk(3) - vj(3)*vk(2)
  vn(2) = vj(3)*vk(1) - vj(1)*vk(3)
  vn(3) = vj(1)*vk(2) - vj(2)*vk(1)
  ! compute twice the area of the triangle and normalizes Vn
  area2 = SQRT(DOT_PRODUCT(vn,vn))
  vn = vn/area2                      !unit normal
  !....   project distance vector over the normal ==> gap
  gap = DOT_PRODUCT(vd,vn)           !penetration
  ! compute the projection of the distance vector to the tangent plane
  vt = vd - gap*vn                   !projection point relative to master
  !  check IF the projection is between the two tangent vectors
  !        auxi = Vj x Vt
  auxi(1) = vj(2)*vt(3) - vj(3)*vt(2)
  auxi(2) = vj(3)*vt(1) - vj(1)*vt(3)
  auxi(3) = vj(1)*vt(2) - vj(2)*vt(1)
  eta1 = DOT_PRODUCT(auxi,vn)/area2          !2nd local coordinate
  !        auxi = Vt x Vk
  auxi(1) = vt(2)*vk(3) - vt(3)*vk(2)
  auxi(2) = vt(3)*vk(1) - vt(1)*vk(3)
  auxi(3) = vt(1)*vk(2) - vt(2)*vk(1)
  xita1 = DOT_PRODUCT(vn,auxi)/area2         !1st local coordinate
  zeta1 = 1d0 - eta1 - xita1                 !3rd local coordinate

  chck1 = eta1 >= 0d0 .AND. xita1 >= 0d0 .AND. zeta1 >= 0d0  !projects?

  IF( curms )THEN !modify gap and normal according to surface curvature
    gap = gap - xita1*eta1 *cu(1,nearn)/2d0  &
              - eta1 *zeta1*cu(2,nearn)/2d0  &
              - zeta1*xita1*cu(3,nearn)/2d0

    bx = (eta1 *(cu(3,nearn)-cu(2,nearn))+(zeta1-xita1)*cu(1,nearn))/2d0
    be = (xita1*(cu(3,nearn)-cu(1,nearn))+(zeta1-eta1 )*cu(2,nearn))/2d0
    vj = vj + bx * vn
    vk = vk + be * vn
    vn(1) = vj(2)*vk(3) - vj(3)*vk(2)
    vn(2) = vj(3)*vk(1) - vj(1)*vk(3)
    vn(3) = vj(1)*vk(2) - vj(2)*vk(1)
    vn = vn/SQRT(DOT_PRODUCT(vn,vn))
  END IF

  IF( gap < cutof0 ) THEN          ! discard segment if GAP is too large
    IF( loop == 40 )EXIT search   !to many iterations
    ! Find the corresponding side element if
    IF( xita1 > zeta1 .AND. eta1 >= zeta1 )THEN
      nnear = nhseg(1,nearn)
    ELSE IF( eta1 > xita1 .AND. zeta1 >= xita1 )THEN
      nnear = nhseg(2,nearn)
    ELSE
      nnear = nhseg(3,nearn)
    END IF

    IF( nnear ==  0 )EXIT search   !No side element, EXIT
    DO i=0,loop          !see if element not checked yet
      IF( nnear == o(i) )EXIT search !element already checked
    END DO
    loop = loop + 1      !increase counter
    nearn = nnear        !new element to check
    o(loop) = nnear      !store element as checked
    CYCLE search         !new search
  END IF

  mgap = MIN(gap,mgap)             ! Compare with previous gap
  chck2 = (gap < 0d0 )             ! gap < 0  == penetration

  IF( chck1 )THEN                            !IF Full projection
    sig = .TRUE.                             !keep full projection
    xita  = xita1
    eta   = eta1
    IF ( chck2 )THEN                         !IF Penetration
      gaps  = gap                        !keep gap
      nchck = 1                          !number of efective penetrations
    END IF

  ELSE ! projection is out of element boundary

    IF( xita1 > zeta1 .AND. eta1 >= zeta1 )THEN
      nnear = nhseg(1,nearn)
    ELSE IF( eta1 > xita1 .AND. zeta1 >= xita1 )THEN
      nnear = nhseg(2,nearn)
    ELSE
      nnear = nhseg(3,nearn)
    END IF

    !check if projects over the expanded segments
    IF ( chck2 ) THEN         ! if penetrated
      !correct local coord. and check if corrected coord. success
      IF( eta1  > -toler .AND. eta1  < tole1 .AND. &
          xita1 > -toler .AND. xita1 < tole1 .AND. &
          zeta1 > -toler .AND. zeta1 < tole1       )THEN

        !compute corrected coordinates (next six lines could be commented)
        IF(eta1  < 0d0 ) eta1  = 0d0
        IF(eta1  > 1d0 ) eta1  = 1d0
        IF(xita1 < 0d0 ) xita1 = 0d0
        IF(xita1 > 1d0 ) xita1 = 1d0
        IF(zeta1 < 0d0 ) zeta1 = 0d0
        IF(zeta1 > 1d0 ) zeta1 = 1d0

        nchck = nchck + 1   !increase number of projected elements
        IF(nchck == 1 .OR. gap > gaps) THEN
          ! if the first or nearer save as segment (ims) for projection
          vns   = vn         ! segment normal
          gaps  = gap        ! gap candidate
          xita = xita1
          eta  = eta1
        END IF  ! nchck = 1 or gap > gaps
        vnn = vnn - vn/gap   !weigthed normal
      END IF  ! projection over expanded triangle

    END IF ! if gap is negative (effective penetration)

    IF( nnear > 0 )THEN       !if a new segment to search exists
      DO i=0,loop             !see if it has not been checked yet
        IF( nnear == o(i) ) EXIT search  !if checked EXIT search
      END DO
      loop = loop + 1         ! not checked, new loop
      IF( loop == 41 )EXIT search  !too many check, EXIT search
      nearn = nnear           ! new segment to check
      o(loop) = nnear         ! remember segment number
      CYCLE search            ! go for a new check
    END IF

  END IF ! if full projection
  EXIT  search !segment OK, or no side element EXIT search
END DO search      !projection loop

! next lines for statistics only
IF( loop > 0 )THEN
  IF( chck1 )THEN                  !full projection
    IF( chck2 )THEN                   !penetrated
      nn(loop,3) = nn(loop,3) + 1    !Full projection and penetration
    ELSE
      nn(loop,2) = nn(loop,2) + 1    !Full projection but no penetration
    END IF
  ELSE
    nn(loop,1) = nn(loop,1) + 1      !No projection
  END IF
END IF
! end of statistics lines

!.... SELECT projection type

IF(nchck >= 1) THEN
  issdb(1) = nearn          !actual segment
  issdb(4) = issdb(4) + 1   !increase iterative penetration
  IF(nchck == 1) THEN       !only one effective penetration
    !....   CASE 1: project onto one master segment
    gap = gaps              !gap
    prdat(2:4) = vn         !normal
  ELSE                      !two or more effective penetration
    !....   CASE 2: weigthed projection onto segment
    gap = SQRT(DOT_PRODUCT(vnn,vnn))   !auxiliar
    vnn = vnn/gap           !normalizes weigthed normal
    gap = gaps/ DOT_PRODUCT(vnn,vns)  !definite gap
    y = xs(1:3) - x(1:3,imn) - vnn * gap  !equivalent position
    !   projection in local system
    CALL projtr(lcseg(1,nearn),x,lcseg(:,nearn),xita,eta,y)
    prdat(2:4) = vnn                  !store normal
    ! WRITE(55,"(e12.4,8f7.4)")gap,xs,vnn,xita,eta
    IF(gap - pgap < gapin0)THEN
      WRITE(55,"('mg',3i6,2e13.4)")lcseg(:,nearn),gap,prdat(9)
      gap = pgap + gapin0
    END IF
  END IF
  prdat(1) = gap       !penetration
  prdat(5) = xita      !first area coordinate
  prdat(6) = eta       !second area coordinate
  ! debug print
  !IF( gap < maxg )THEN
  !  WRITE(55,"(' mgap ',3i6,2e15.3)")lcseg(:,nearn),gap,maxg
  !  maxg = gap
  !END IF
  ! end debug print
ELSE
  prdat(1) = MAX(mgap,0d0)       ! gap >= 0
  issdb(4) = MAX( issdb(4)-1, 0) ! iterative penetration set to 0
  IF(sig)THEN
    issdb(1) =  nearn         !keep segment with projection
  ELSE
    issdb(1) =  o(0)          !original nearest segment
  END IF
END IF

RETURN
END SUBROUTINE projt4
