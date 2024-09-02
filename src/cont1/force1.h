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
      END SUBROUTINE force1
