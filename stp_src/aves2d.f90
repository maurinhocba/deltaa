SUBROUTINE aves2d( )
!
! interpolate smoothed values for quadratic elements
! 2-D solid elements
!
USE data_db
IMPLICIT NONE

! local variables
TYPE (sol2d), POINTER :: eset
INTEGER :: ielem,i,n,np,nn,nelem,nnode,nvert,lnods(9)
INTEGER, POINTER :: nodes(:,:)
REAL(kind=8), POINTER :: vargs(:,:),accpn(:)

vargs => sol2d_vargs         !pointers to use short names
accpn => sol2d_accpn
nodes => sol2d_nodes

eset => sol2d_head           !point to first set

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
      IF( nnode == 9 )THEN
        n = lnods(9)                  !associated node
        IF( accpn(n) /= 0) CYCLE
        vargs(:,n) = ( vargs(:,lnods(1)) + vargs(:,lnods(2))        &
                     + vargs(:,lnods(3)) + vargs(:,lnods(4)) )/4d0
      END IF
    END DO
  END IF
  eset => eset%next
END DO

RETURN
END SUBROUTINE aves2d
