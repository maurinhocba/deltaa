SUBROUTINE predic(istep,lauto,ecdis,neq,arcln,dlamb,disax,displ,   &
                  ddisp,delta,karcl,piter,diter,newtv,ncdis)
!***********************************************************************
!
!*** this routine predicts displacement and load step increments
!    according to selected path and previous increments
!
!***********************************************************************
!USE npoi_db, ONLY : ndime,npoin
!USE kin2_db, ONLY : ifpre
IMPLICIT NONE
!       routine arguments
INTEGER (kind=4),INTENT(IN) :: istep,neq,lauto,karcl,newtv
INTEGER (kind=4),INTENT(IN OUT) :: ecdis,ncdis
REAL (kind=8),INTENT(IN) :: displ(:),piter,diter
REAL (kind=8),INTENT(IN OUT) :: dlamb,arcln,disax(:,:),ddisp(:),delta
END SUBROUTINE predic
