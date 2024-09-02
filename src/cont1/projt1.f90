      SUBROUTINE projt1(xs,nearn,x,isn,lcseg,   &
                        nhseg,issdb,prdat,cutof,gapin,curms,cu)

!.... project the slave node onto a 2-d master segment

      IMPLICIT NONE
!     arguments
      LOGICAL, INTENT(IN) :: curms
      INTEGER (kind=4), INTENT(IN) :: isn,lcseg(:,:),nhseg(:,:)
      INTEGER (kind=4), INTENT(IN OUT) :: nearn,issdb(:)
      REAL (kind=8), INTENT (IN) :: x(:,:),cutof,gapin,xs(:),cu(:)
      REAL (kind=8), INTENT (IN OUT) :: prdat(:)

!     local variables
      INTEGER (kind=4) :: imn,jmn,nchck,masts,loop,nnear,i,o(0:2)
      REAL (kind=8) :: vd(2),vj(2),vt(2),lengt,gap,dist,y(2),mgap, &
                       xita,xita1,vn(2),vnn(2),vns(2),gaps,bx,     &
                       pgap,gapin0,cutof0
      LOGICAL ::  chck1,chck2,sig
      REAL (kind=8), PARAMETER :: toler = 0.01d0, tole1 = 1.01d0

      REAL (kind=8), SAVE :: maxg = -.25d-3  !variable for debug

!.... initialize to zero some control parameters
      pgap = MIN( prdat(6), 0d0)                  !previos converged penetration
      gapin0 = MIN(gapin,pgap/2d0)                !maximum incremental gap
      cutof0 = cutof + pgap/2d0                   !maximum penetration
      cutof0 = MAX(cutof0, 3d0*cutof )            !maximum penetration
      nchck = 0            !initializes number of candidates
      sig = .FALSE.        !initializes to NO full projection
      mgap = 1d10          !initializes max gap found
      vnn = 0d0            !initializes weighted normal
      masts = 0            !initializes to no segment
      loop = 0             !initializes loop search
      o(loop) = nearn
!.... loop until projection
      search : DO
        !....   identify global master nodes
        imn = lcseg(1,nearn)        !1st master node
        jmn = lcseg(2,nearn)        !2nd master node
        vd = xs(1:2) - x(1:2,imn)          !distance vector
        dist = SQRT(DOT_PRODUCT(vd,vd))    !distance
        mgap = MIN(mgap,dist)              !for boundaries
        IF(dist == 0d0) THEN               !IF nodes coincide
          nchck = 0                        !no penetration
          sig = .TRUE.                     !keep the node
          mgap = 0d0                       !update max gap
          masts = nearn
          EXIT search                      !exit search
        END IF
        !....   form tangent vectors to master segment
        vj = x(1:2,jmn) - x(1:2,imn)       !side i-->j
        ! compute segment length and normalizes Vj
        lengt = SQRT(DOT_PRODUCT(vj,vj))
        vj = vj/lengt                      !unit normal
        ! compute normal vector
        vn(1) = -vj(2)
        vn(2) =  vj(1)
        !....   project distance vector over the normal ==> gap
        gap = DOT_PRODUCT(vd,vn)           !penetration
        ! compute the projection of the distance vector to the tangent plane
        vt = vd - gap*vn                   !projection point relative to master
        !  check IF the projection is over the segment
        xita1 = DOT_PRODUCT(vt,vj)/lengt         ! local coordinate

        chck1 = xita1 >= 0d0 .AND. xita1 <= 1d0    !projects?

        IF( curms )THEN !modify gap and normal according to surface curvature
          gap = gap - xita1*(1d0-xita1) *cu(nearn)/2d0
          bx = (0.5d0-xita1)*cu(nearn)/2d0
          vn(1) = -(vj(2)*lengt+bx*vn(2))
          vn(2) =  (vj(1)*lengt+bx*vn(1))
          vn = vn/SQRT(DOT_PRODUCT(vn,vn))
        END IF

        IF( gap < cutof0 ) THEN          ! discard segment if GAP is too large
          IF(chck1 .AND. gap/cutof0 < 3d0) &
            WRITE(55,"('ct',3i6,2e13.4)")isn,lcseg(:,nearn),gap,prdat(2)
          IF( loop == 2 )EXIT search    ! too many segments considered
          ! Find the corresponding side element if
          IF( xita1 <= 0.5d0 )THEN      !previous segment
            nnear = nhseg(1,nearn)
          ELSE                          !next segment
            nnear = nhseg(2,nearn)
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
          masts = nearn                            !segment candidate
          xita  = xita1
          IF ( chck2 )THEN                         !IF Penetration
            gaps  = gap                        !keep gap
            nchck = 1                          !number of efective penetrations
          END IF

        ELSE ! projection is out of element boundary

          IF( xita1 < 0d0 )THEN       ! previous segment
            nnear = nhseg(1,nearn)
          ELSE                        ! next segment
            nnear = nhseg(2,nearn)
          END IF

          !check if projects over the expanded segments
          IF ( chck2 ) THEN         ! if penetrated
            !correct local coord. and check if corrected coord. success
            IF( xita1 > -toler .AND. xita1 < tole1 )THEN

              !compute corrected coordinates (next two lines could be commented)
              IF(xita1 < 0d0 ) xita1 = 0d0
              IF(xita1 > 1d0 ) xita1 = 1d0

              nchck = nchck + 1   !increase number of projected elements
              IF(nchck == 1 .OR. gap > gaps) THEN
                ! if the first or nearer save as segment (ims) for projection
                masts = nearn      ! new segment candidate
                vns   = vn         ! segment normal
                gaps  = gap        ! gap candidate
                xita = xita1
              END IF  ! nchck = 1 or gap > gaps
              vnn = vnn - vn/gap   !weigthed normal
            END IF  ! projection over expanded segment

          END IF ! if gap is negative (effective penetration)

          IF( nnear > 0 )THEN       !if a new segment to search exists
            DO i=0,loop             !see if it has not been checked yet
              IF( nnear == o(i) ) EXIT search  !if checked EXIT search
            END DO
            loop = loop + 1         ! not checked, new loop
            IF( loop == 3 )EXIT search  !too many check, EXIT search
            nearn = nnear           ! new segment to check
            o(loop) = nnear         ! remember segment number
            CYCLE search            ! go for a new check
          END IF

        END IF ! if full projection
        EXIT  search !segment OK, or no side element EXIT search
      END DO search      !projection loop

!.... SELECT projection type

      IF(nchck >= 1) THEN
        issdb(1) = nearn          !actual segment
        issdb(4) = issdb(4) + 1   !increase iterative penetration
        IF(nchck == 1) THEN       !only one effective penetration
          !....   CASE 1: project onto one master segment
          gap = gaps              !gap
          prdat(2:3) = vn         !normal
        ELSE                      !two or more effective penetration
          !....   CASE 2: weigthed projection onto segment
          gap = SQRT(DOT_PRODUCT(vnn,vnn))   !auxiliar
          vnn = vnn/gap           !normalizes weighted normal
          gap = gaps/ DOT_PRODUCT(vnn,vns)  !definite gap
          y = xs(1:2) - x(1:2,imn) - vnn * gap  !equivalent position
          !  recompute local coordinate
          xita = DOT_PRODUCT(y,vj)/lengt
          prdat(2:3) = vnn                  !store normal
          WRITE(55,"(e12.4,3f9.6)")gap,vnn,xita  !debug print
        END IF
        IF(gap - pgap < gapin0)THEN
          WRITE(55,"('mg',3i6,2e13.4)")isn,lcseg(:,nearn),gap,prdat(2)
          gap = pgap + gapin0
        END IF
        prdat(1) = gap       !penetration
        prdat(4) = xita      !local coordinate
        ! debug print
        IF( gap < maxg )THEN
          WRITE(55,"(' mgap ',2e15.3)")gap,maxg
          maxg = gap
        END IF
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
      END SUBROUTINE projt1
