SUBROUTINE force4(isn,issdb,rssdb,lcseg,indco,cprop,x,  &
                  dtime,emass,fcont,surtf,cmptf,cursl,tn,icnod,press, &
                  wear,wwear,npres )

!.... form contact residual force for a node-triangle 3-d contact element

!USE c_input, ONLY : runend,lures
IMPLICIT NONE
!     dummy arguments
LOGICAL, INTENT(IN) :: cmptf,cursl,press,wear  !compute total contact forces?
INTEGER (kind=4), INTENT (IN) :: indco,isn,lcseg(:,:),icnod
INTEGER (kind=4), INTENT (IN OUT) :: issdb(:)
REAL (kind=8), INTENT (IN) :: cprop(:),x(:,:),emass(:,:),dtime,tn(:,:)
REAL (kind=8), INTENT (IN OUT) :: rssdb(:),surtf(:),fcont(:,:),wwear(:)
REAL (kind=8), INTENT (OUT) :: npres
END SUBROUTINE force4
