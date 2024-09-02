SUBROUTINE cscf1a(lcseg,nhseg,ncnod,lcnod,xc,            &
                  issdb,rssdb,coors,coorm,               &
                  fcont,surtf,emass,indco,cprop,         &
                  dtime,cmptf,nsegm,freq,                &
                  cursl,curms,tn,cu,press,presn,         &
                  wear,wwear,wrink,mingp)

! perform search for contact interface
! and compute contact forces for a given master/slave pair

IMPLICIT NONE
!   Dummy arguments
LOGICAL, INTENT(IN) :: cmptf,cursl,curms,press,wear,wrink
INTEGER (kind=4), INTENT(IN) :: ncnod,indco,nsegm,freq,               &
                                lcseg(:,:),nhseg(:,:),lcnod(:)
INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)
REAL (kind=8), INTENT(IN) :: coors(:,:),coorm(:,:),emass(:,:),        &
                             dtime,cprop(4),xc(:,:),cu(:),tn(:,:)
REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:),fcont(:,:),surtf(:),      &
                                 presn(:),wwear(:),mingp(:)

!  local variables
INTEGER (kind=4) :: icnod,nearn,isn,nnp
REAL (kind=8) :: npres
INTERFACE
  INCLUDE 'nearst.h'
  INCLUDE 'projt1.h'
  INCLUDE 'force1.h'
END INTERFACE


IF( cmptf ) surtf = 0d0
!.... loop over all slave nodes
DO icnod = 1, ncnod

  !... step 1 ... look for the nearest segment

  isn = lcnod(icnod) ! global number of the current slave node
  nearn = issdb(1,icnod) !previous nearest segment
  nnp   =-issdb(3,icnod) !previous time of effective penetration

  IF(nearn == 0 .OR. nnp > freq ) THEN   !if no previous or too many steps
    !global search in ALL master segment the nearest segment to ISN
    CALL nearst(xc,nsegm,nearn,coors(:,isn))
    issdb(3,icnod) = 0   ! re-initializes counter

  ELSE IF( nearn < 0)THEN           !If a search is advisable
    nearn = -nearn
    IF(nearn > 1000000 )CYCLE      !proj. seems not necessary yet
  END IF

  !.... step 2 .... project slave node onto the master surface

  !....   two-dimensional 2-node master segments
  CALL projt1(coors(:,isn),nearn,coorm,isn,lcseg, &
              nhseg,issdb(:,icnod),rssdb(:,icnod),cprop(4),cprop(5),curms,cu)
  IF( issdb(4,icnod) == 0 ) CYCLE   ! no penetration

  !.... step 3 .... Compute contact force

  CALL force1(isn,issdb(:,icnod),rssdb(:,icnod),lcseg,indco,cprop, &
              coorm,dtime,emass,fcont,surtf,cmptf,cursl,tn,icnod,  &
              press,wear,wwear,npres )
  IF( press ) presn(icnod) = presn(icnod) + npres

END DO
RETURN
END SUBROUTINE cscf1a
