SUBROUTINE update(istep,itera,lauto,ecdis,neq,arcln,lambd,dlamb,  &
                  disax,displ,ddisp,delta,karcl,piter,diter,newtv,ncdis)
!***********************************************************************
!
!***  this routine calculates the change in load step according to
!     selected path and updates load level and displacement vectors
!
!***********************************************************************
!USE lispa0
IMPLICIT NONE
!       routine arguments
INTEGER (kind=4),INTENT(IN) :: istep,itera,neq,lauto,karcl,newtv
INTEGER (kind=4),INTENT(IN OUT) :: ecdis,ncdis
REAL (kind=8),INTENT(IN) :: piter,diter
REAL (kind=8),INTENT(IN OUT) :: dlamb,lambd,arcln,disax(:,:),displ(:),ddisp(:)
REAL (kind=8),INTENT(OUT) :: delta

END SUBROUTINE update
