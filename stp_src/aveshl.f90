SUBROUTINE aveshl( )
!
! interpolate smoothed values for quadratic elements
! 3-D Shell (Simo) elements
!
USE data_db
IMPLICIT NONE

! local variables
TYPE (shl3d), POINTER :: eset
INTEGER :: ielem,i,n,np,nn,nelem,nnode,nvert,lnods(6)
INTEGER, POINTER :: nodes(:,:)
REAL(kind=8), POINTER :: vargs(:,:),accpn(:)

vargs => shl3d_vargs         !pointers to use short names
accpn => shl3d_accpn
nodes => shl3d_nodes

eset => shl3d_head           !point to first set

DO
  IF( .NOT.ASSOCIATED(eset) )EXIT
  nelem = eset%nelem      !number of elements in the set
  nnode = eset%nnode      !number of nodes per element
  IF( nnode > 4 )THEN     !only for quadratic elements
    nvert = nnode/2       !number of lower order nodes

    DO ielem=1,nelem                             !for each element
      lnods(1:nnode) = eset%lnods(:,ielem)       !element connectivities
      lnods(1:nnode) = nodes(lnods(1:nnode),1)     !local nodes
      DO i=nvert+1,2*nvert                   !for each mid side node
        n = lnods(i)                         !associated node
        IF( accpn(n) /= 0) CYCLE
        np = lnods(i-nvert)               !previous node
        nn = lnods(MOD(i-nvert,nnode)+1)  !next node
        vargs(:,n) = ( vargs(:,np) + vargs(:,nn) )/2d0  !average value
      END DO
    END DO
  END IF
  eset => eset%next
END DO

RETURN
END SUBROUTINE aveshl
