SUBROUTINE force4(isn,issdb,rssdb,lcseg,indco,cprop,x,  &
                  dtime,emass,fcont,surtf,cmptf,cursl,tn,icnod,press, &
                  wear,wwear,npres )

!.... form contact residual force for a node-triangle 3-d contact element

USE c_input, ONLY : runend,lures
IMPLICIT NONE
!     dummy arguments
LOGICAL, INTENT(IN) :: cmptf,cursl,press,wear  !compute total contact forces?
INTEGER (kind=4), INTENT (IN) :: indco,isn,lcseg(:,:),icnod
INTEGER (kind=4), INTENT (IN OUT) :: issdb(:)
REAL (kind=8), INTENT (IN) :: cprop(:),x(:,:),emass(:,:),dtime,tn(:,:)
REAL (kind=8), INTENT (IN OUT) :: rssdb(:),surtf(:),fcont(:,:),wwear(:)
REAL (kind=8), INTENT (OUT) :: npres
!     local variables
LOGICAL :: same,last
INTEGER (kind=4) :: n1,n2,n3,n1o,n2o,n3o,nearn,onear,n,node
REAL (kind=8) :: r,s,t,ro,so,to,ds,dt,pgap,pnltc,cofri,m,maux,    &
                 tslip,eslip,fricm,frict,a2,auxil,xfact,pslip,    &
                 vns(3),t1(3),t2(3),t3(3),ts(3),ta(3),eresf(3,4),iter

INTERFACE
  INCLUDE 'vecuni.h'
  INCLUDE 'vecpro.h'
END INTERFACE

!.... read gap, normal direction and isoparametric coordinates of the
!.... projection point, from real-type slave surface database
IF( rssdb(1) < cprop(4) )THEN
  WRITE(55,"(' gap > cutof, node:',i6,2e15.4)")isn,rssdb(1),rssdb(9)
  WRITE(* ,"(' gap > cutof, node:',i6,2e15.4)")isn,rssdb(1),rssdb(9)
END IF
iter = DBLE(issdb(4))/4d0
IF(iter > 0.8d0)THEN
  !pgap = -cprop(1)*rssdb(1)   !normal penalty * gap
  pgap = -cprop(1)*rssdb(1)*(1d0+rssdb(1)/cprop(4)) !normal penalty * gap
ELSE
  !pgap = -cprop(1)*(rssdb(9) + (rssdb(1) - rssdb(9))*iter)
  xfact= rssdb(9) + (rssdb(1) - rssdb(9))*iter
  pgap = -cprop(1)*xfact*(1d0+xfact/cprop(4)) !normal penalty * gap
END IF

vns = rssdb(2:4)          !normal direction
IF( cursl .AND. indco /= 1)THEN
  IF(indco == 0)THEN
    t1  = (-tn(:,icnod)+vns)/2d0       !average normal
    vns = t1/SQRT(DOT_PRODUCT(t1,t1))
  ELSE         ! IF(indco == 2)THEN
    vns = -tn(1:3,icnod)              !node normal direction
  END IF
  rssdb(2:4) = vns
END IF
IF( pgap <= 0d0 )RETURN
vns = vns*pgap               !normal direction(scaled)
                             !triangular area coordinates
s     = rssdb(5)             !node 2 (xita)
t     = rssdb(6)             !node 3 (eta)
r     = 1d0 - s - t          !node 1 (zeta)
nearn = issdb(1)             !segment

!.... form residual force vector
IF( indco == 0 .OR. indco == 1) THEN   !slave surface
  eresf(1:3,1)  = vns     !normal residual forces
END IF
IF( indco == 0 .OR. indco == 2) THEN   !master surface
  !....   form ns operator
  eresf(1:3,2) = -r*vns   !normal residual forces node1
  eresf(1:3,3) = -s*vns   !normal residual forces node2
  eresf(1:3,4) = -t*vns   !normal residual forces node3
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
IF(cofri /= 0D0 ) THEN  ! if friction
  !  computes friction forces  tslip: total slip;  eslip : elastic slip
  IF( .NOT.last ) THEN          ! IF no penetration in previous step
    frict = 0d0
  ELSE                          ! Penetration in previous step

    fricm = cofri*pgap              !maximum friction force,
    n1 = lcseg(1,nearn)             !nodes defining master segment
    n2 = lcseg(2,nearn)
    n3 = lcseg(3,nearn)
    t1 = x(1:3,n2) - x(1:3,n1)    !side vectors n1->n2
    t2 = x(1:3,n3) - x(1:3,n1)    !             n1->n3
    so = rssdb(7)                 !old onset coordinates
    to = rssdb(8)
    onear = issdb(2)              !previous segment
    same = onear == nearn         !same of different segments
    IF( same )THEN  ! IF penetration in previous step was in the same segment
      ds = s - so                 !differences in natural coordinates
      dt = t - to
      ts = t1*ds + t2*dt          !slip vector
      CALL vecuni(3,ts,tslip)     !ts = unit vector, tslip = total slip
      frict = pnltc*tslip         !tangential force
      IF(frict > fricm) THEN      !compare with maximum force
        ! IF friction force > maximum friction force ==> slip
        frict = fricm             !assign maximum force
        eslip = frict/pnltc       !maximum slip
        rssdb(5) = s - ds*eslip/tslip   !updates onset natural coordinates
        rssdb(6) = t - dt*eslip/tslip
      END IF
    ELSE  ! segment in previous step was in other segment
      n1o= lcseg(1,onear)     !nodes defining previous master segment
      n2o= lcseg(2,onear)
      n3o= lcseg(3,onear)
      ro = 1d0-so-to          !previous r coord.
      !vector between actual and onset point
      ts= r *x(:,n1) + s *x(:,n2) + t *x(:,n3)   &
         -ro*x(:,n1o)- so*x(:,n2o)- to*x(:,n3o)
      CALL vecpro(t1,t2,t3)         ! t3 = t1 x t2
      CALL vecuni(3,t3,a2)          ! a2 = twice the area
      auxil = DOT_PRODUCT(ts,t3)    ! projects vector over normal
      ts = ts - auxil*t3            ! orthogonal projection
      CALL vecuni(3,ts,tslip)       ! ts = unit vector, tslip = total slip
      frict = pnltc*tslip           ! friction force
      IF(frict > fricm) THEN        ! compare with maximum friction force
        ! IF friction force > maximum friction force ==> slip
        eslip = fricm/pnltc         ! maximum relative slip
        frict = fricm               ! assign maximum friction force
      ELSE
        eslip = tslip               ! all tslip is elastic
      END IF  !maximum tangential force exceded
      CALL vecpro(ts,t2,ta)       ! ta = ts x t2
      auxil = DOT_PRODUCT(ta,t3)/a2 *eslip  ! slip Ds coordinate
      rssdb(5) = s - auxil        ! store onset natural coordinate
      CALL vecpro(t1,ts,ta)       ! ta = t1 x ts
      auxil = DOT_PRODUCT(ta,t3)/a2 *eslip  ! slip Dt coordinate
      rssdb(6) = t - auxil        ! store onset natural coordinate
    END IF  !same or different segment
    IF( cursl .AND. indco /= 1 .AND. frict /= 0d0)THEN
      a2 = DOT_PRODUCT(ts,vns)
      ts = ts - a2*vns
      ts = ts/SQRT(DOT_PRODUCT(ts,ts))
    END IF
    ts = -ts*frict
    IF( indco  == 0 .OR. indco == 1) THEN     ! slave surface
      eresf(1:3,1) = eresf(1:3,1) + ts   ! add to normal forces
    END IF
    IF( indco == 0 .OR. indco == 2 )THEN      ! master surface
      eresf(1:3,2) = eresf(1:3,2) - r*ts      ! add to normal forces node 1
      eresf(1:3,3) = eresf(1:3,3) - s*ts      ! add to normal forces node 2
      eresf(1:3,4) = eresf(1:3,4) - t*ts      ! add to normal forces node 3
    END IF
  END IF ! projection in previous step
END IF !if friction

!  computes equivalent mass
m = 0d0  !initializes
maux = MAXVAL(emass(1:3,isn))   ! mass of slave node
IF( maux /= 0d0) m = 1d0/maux
IF( indco == 0 ) THEN           ! IF forces in both surfaces
  !....   consider master nodes mass
  DO n =1,3                         !For each master node
    node = lcseg(n,nearn)           !global node number
    maux = MAXVAL(emass(1:3,node))  !max mass in any direction
    IF(maux > 0)THEN                ! add to equivalent mass
      SELECT CASE (n)
      CASE (1)                      !first node
         m = m + r**2/maux            !use r
      CASE (2)                      !second node
         m = m + s**2/maux            !use s
      CASE (3)                      !third node
         m = m + t**2/maux            !use t
      END SELECT
    END IF
  END DO
END IF
IF(m == 0d0) THEN
  WRITE(*    ,*)'  Error in nodal mass ',isn,lcseg(1:3,nearn)
  WRITE(lures,*)'  Error in nodal mass ',isn,lcseg(1:3,nearn)
  CALL runend('FORCE4: zero mass in a contact pair')
END IF
xfact = 2d0/(m*dtime**2)     !computes xfact

!     Add forces into global contact forces

IF( press ) npres = pgap*xfact*dtime
IF(indco == 0 .OR. indco == 1) THEN           !forces on SLAVE
  fcont(1:3,isn) = fcont(1:3,isn) - eresf(1:3,1)*xfact
  IF( cmptf ) surtf = surtf - eresf(1:3,1)*xfact
END IF
IF(indco == 0 .OR. indco == 2)THEN           !forces on MASTER
  DO n=1,3                 !for each node
    node = lcseg(n,nearn)  !global number
    fcont(1:3,node) = fcont(1:3,node) - eresf(1:3,1+n)*xfact
  END DO
  IF( cmptf .AND. indco == 2 ) &
    surtf = surtf + (eresf(1:3,2)+eresf(1:3,3)+eresf(1:3,4))*xfact
END IF

IF( wear .AND. pslip /= 0d0)THEN
  pslip = ABS( frict*pslip*xfact )
  wwear(isn) = wwear(isn) + pslip
  node = lcseg(1,nearn)  !global number
  wwear(node) = wwear(node) + pslip*r
  node = lcseg(2,nearn)  !global number
  wwear(node) = wwear(node) + pslip*s
  node = lcseg(3,nearn)  !global number
  wwear(node) = wwear(node) + pslip*t
END IF
RETURN
END SUBROUTINE force4
