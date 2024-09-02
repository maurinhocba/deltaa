      SUBROUTINE force1(isn,issdb,rssdb,lcseg,indco,cprop,x,dtime,  &
                        emass,fcont,surtf,cmptf,cursl,tn,icnod,press, &
                        wear,wwear,npres)

!.... form contact residual force for a node-segment 2-d contact element

      USE c_input, ONLY : lures,runend
      IMPLICIT NONE
!     dummy arguments
      LOGICAL, INTENT(IN) :: cmptf,cursl,press,wear
      INTEGER (kind=4), INTENT (IN) :: indco,isn,lcseg(:,:),icnod
      INTEGER (kind=4), INTENT (IN OUT) :: issdb(:)
      REAL (kind=8), INTENT (IN) :: cprop(:),x(:,:),emass(:,:), &
                     dtime,tn(:,:)
      REAL (kind=8), INTENT (IN OUT) :: rssdb(:),surtf(:),fcont(:,:),wwear(:)
      REAL (kind=8), INTENT (OUT) :: npres
!     local variables
      LOGICAL :: last,same
      INTEGER (kind=4) :: n1,n2,n1o,n2o,nearn,onear,n,node
      REAL (kind=8) :: r,s,ro,so,ds,pgap,pnltc,cofri,m,maux,iter,         &
                       tslip,eslip,fricm,frict,lengt,auxil,xfact,         &
                       vns(2),t1(2),ts(2),eresf(2,3),pslip

      INTERFACE
        INCLUDE 'vecuni.h'
      END INTERFACE

!.... read gap, normal direction and isoparametric coordinates of the
!.... projection point, from real-type slave surface database
      IF( rssdb(1) < cprop(4) )THEN
        WRITE(55,"(' gap > cutof, node:',i6,2e15.4)")isn,rssdb(1),rssdb(6)
        WRITE(* ,"(' gap > cutof, node:',i6,2e15.4)")isn,rssdb(1),rssdb(6)
      END IF
      iter = DBLE(issdb(4))/4d0
      IF(iter > 0.8d0)THEN
        !pgap = -cprop(1)*rssdb(1)                         !normal penalty * gap
        pgap = -cprop(1)*rssdb(1)*(1d0+rssdb(1)/cprop(4)) !normal penalty * gap
      ELSE
        !pgap = -cprop(1)*(rssdb(6) + (rssdb(1) - rssdb(6))*iter)
        xfact= rssdb(6) + (rssdb(1) - rssdb(6))*iter
        pgap = -cprop(1)*xfact*(1d0+xfact/cprop(4)) !normal penalty * gap
      END IF

      vns = rssdb(2:3)          !normal direction
      IF( cursl .AND. indco /= 1)THEN
        IF(indco == 0)THEN
          t1 = (-tn(:,icnod)+vns)/2d0       !average normal
          vns = t1/SQRT(DOT_PRODUCT(t1,t1))
        ELSE         ! IF(indco == 2)THEN
          vns = -tn(:,icnod)              !node normal direction
        END IF
        rssdb(2:3) = vns
      END IF
      IF( pgap <= 0d0 )RETURN
      vns = vns*pgap               !normal direction(scaled)
                                   !segment linear shape functions
      s     = rssdb(4)             !node 2 (xita)
      r     = 1d0 - s              !node 1 (1-xita)
      nearn = issdb(1)             !segment

!.... form residual force vector
      IF( indco == 0 .OR. indco == 1) THEN   !slave surface
        eresf(1:2,1)  = vns     !normal residual forces
      END IF
      IF( indco == 0 .OR. indco == 2) THEN   !master surface
        !....   form ns operator
        eresf(1:2,2) = -r*vns   !normal residual forces node1
        eresf(1:2,3) = -s*vns   !normal residual forces node2
      END IF

!     friction treatment
      pslip = 0d0
      last  = issdb(3) > 0            !contact at previous step
      pnltc = cprop(2)                !tangential penalty parameter
      IF( issdb(3) == 1)THEN          !no-slip in previous step
        cofri = cprop(3)              !Static friction coefficient
      ELSE
        cofri = cprop(6)              !Kinetic friction coefficient
      END IF
      IF(cofri /= 0D0) THEN           !if no friction, that's all
        !  computes friction forces  tslip: total slip;  eslip : elastic slip
        IF( .NOT.last ) THEN          ! IF no penetration in previous step
          frict = 0d0
        ELSE                          ! Penetration in previous step

          !  computes friction forces  tslip: total slip;  eslip : elastic slip
          fricm = cofri*pgap              !maximum friction force,
          n1 = lcseg(1,nearn)             !nodes defining master segment
          n2 = lcseg(2,nearn)
          t1 = x(:,n2) - x(:,n1)        !side vectors n1->n2
          so = rssdb(5)                 !old onset coordinates
          onear = issdb(2)              !previous segment
          same = onear == nearn         !same of different segments
          IF( same )THEN  ! IF penetration in previous step was in the same segment
            ds = s - so                 !differences in natural coordinates
            ts = t1*ds                  !slip vector
            CALL vecuni(2,ts,tslip)     !ts = unit vector, tslip = total slip
            frict = pnltc*tslip         !tangential force (+)
            IF(frict > fricm) THEN      !compare with maximum force
              ! IF friction force > maximum friction force ==> slip
              frict = fricm             !assign maximum force
              eslip = frict/pnltc       !maximum slip
              rssdb(4) = s - ds*eslip/tslip   !updates onset natural coordinates
            END IF
          ELSE  ! segment in previous step was in other segment
            n1o= lcseg(1,onear)     !nodes defining previous master segment
            n2o= lcseg(2,onear)
            ro = 1d0-so             !previous r coord.
            !vector between actual and onset point
            ts= r *x(:,n1) + s *x(:,n2) -ro*x(:,n1o) -so*x(:,n2o)
            CALL vecuni(2,t1,lengt)       ! side length
            auxil = DOT_PRODUCT(ts,vns)   ! projects vector over normal
            ts = ts - auxil*vns           ! orthogonal projection
            CALL vecuni(2,ts,tslip)       ! ts = unit vector, tslip = total slip
            frict = pnltc*tslip           ! friction force
            IF(frict > fricm) THEN        ! compare with maximum friction force
              ! IF friction force > maximum friction force ==> slip
              eslip = fricm/pnltc         ! maximum relative slip
              frict = fricm               ! assign maximum friction force
            ELSE
              eslip = tslip               ! all tslip is elastic
            END IF  !maximum tangential force exceded
            auxil = DOT_PRODUCT(t1,ts)*eslip/lengt  ! slip Ds coordinate
            rssdb(4) = s - auxil        ! store onset natural coordinate
          END IF  !same or different segment
          IF( cursl .AND. indco /= 1 .AND. frict /= 0d0)THEN
            auxil = DOT_PRODUCT(ts,vns)
            ts = ts - auxil*vns
            ts = ts/SQRT(DOT_PRODUCT(ts,ts))
          END IF

          ts = -ts*frict
          IF( indco  == 0 .OR. indco == 1) THEN     ! slave surface
            eresf(1:2,1) = eresf(1:2,1) + ts   ! add to normal forces
          END IF
          IF( indco == 0 .OR. indco == 2 )THEN      ! master surface
            eresf(1:2,2) = eresf(1:2,2) - r*ts      ! add to normal forces node 1
            eresf(1:2,3) = eresf(1:2,3) - s*ts      ! add to normal forces node 2
          END IF
        END IF !if projection in previous step
      END IF !if friction

      !  computes equivalent mass
      m = 0d0  !initializes
      maux = MAXVAL(emass(1:2,isn))   ! mass of slave node
      IF( maux /= 0d0) m = 1d0/maux
      IF( indco  == 0) THEN           ! IF forces in both surfaces
        !....   consider master nodes mass
        DO n =1,2                         !For each master node
          node = lcseg(n,nearn)           !global node number
          maux = MAXVAL(emass(1:2,node))  !max mass in any direction
          IF(maux > 0)THEN                ! add to equivalent mass
            SELECT CASE (n)
            CASE (1)                      !first node
               m = m + r**2/maux            !use r
            CASE (2)                      !second node
               m = m + s**2/maux            !use s
            END SELECT
          END IF
        END DO
      END IF
      IF(m == 0d0) THEN
        WRITE(*,    *)'  Error in nodal mass ',isn,lcseg(1:2,nearn)
        WRITE(lures,*)'  Error in nodal mass, Internal numbers are: Slave ',isn, &
                    & '  Master segment ',lcseg(1:2,nearn)
        WRITE(lures,*)'  Check data, at least one surface Must have mass'
        CALL runend('CONTA1: zero masses in contact pair')
      END IF
      xfact = 2d0/(m*dtime**2)     !computes xfact

!     Add forces into global contact forces

      IF(indco == 0 .OR. indco == 1) THEN           !forces on SLAVE
        fcont(1:2,isn) = fcont(1:2,isn) - eresf(1:2,1)*xfact
        IF( cmptf ) surtf = surtf - eresf(1:2,1)*xfact
        IF( press ) npres = pgap*xfact*dtime
      END IF
      IF(indco == 0 .OR. indco == 2)THEN           !forces on MASTER
        DO n=1,2                 !for each node
          node = lcseg(n,nearn)  !global number
          fcont(1:2,node) = fcont(1:2,node) - eresf(1:2,1+n)*xfact
        END DO
        IF( cmptf .AND. indco == 2 ) &
          surtf = surtf + (eresf(1:2,2)+eresf(1:2,3))*xfact
      END IF
      IF( wear .AND. pslip /= 0d0)THEN
        pslip = ABS( frict*pslip*xfact )
        wwear(isn) = wwear(isn) + pslip
        node = lcseg(1,nearn)  !global number
        wwear(node) = wwear(node) + pslip*r
        node = lcseg(2,nearn)  !global number
        wwear(node) = wwear(node) + pslip*s
      END IF

      RETURN
      END SUBROUTINE force1
