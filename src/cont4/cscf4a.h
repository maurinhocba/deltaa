SUBROUTINE cscf4a(lcseg,nhseg,ncnod,lcnod,xc,            &
                  issdb,rssdb,coors,coorm,               &
                  fcont,surtf,emass,indco,cprop,         &
                  dtime,cmptf,nsegm,freq,                &
                  cursl,curms,tn,cu,press,presn,         &
                  wear,wwear)

! perform search for contact interface
! and compute contact forces for a given master/slave pair

IMPLICIT NONE
!   Dummy arguments
LOGICAL, INTENT(IN) :: cmptf,cursl,curms,press,wear
INTEGER (kind=4), INTENT(IN) :: ncnod,indco,nsegm,freq,   &
                  lcseg(:,:),nhseg(:,:),lcnod(:)
INTEGER (kind=4), INTENT(IN OUT) :: issdb(:,:)
REAL (kind=8), INTENT(IN) :: coors(:,:),coorm(:,:),emass(:,:),        &
                             dtime,cprop(:),xc(:,:),tn(:,:),cu(:,:)
REAL (kind=8), INTENT(IN OUT) :: rssdb(:,:),fcont(:,:),surtf(:),      &
                                 presn(:),wwear(:)

END SUBROUTINE cscf4a
