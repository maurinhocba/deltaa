SUBROUTINE mastd4(nsegm,nnseg,lcseg,nhseg,confor,x)
!
!   Generates data base for surfaces that will act as master in some pairs.
!   This data includes the neighbours segments
!
!.... input
!....   nsegm = number of segments of the master surface considered
!....   nnseg = number of nodes per segment of the surface considered
!....   lcseg(inseg,icseg) = global node number for the local element node
!....                        [inseg] of the segment [icseg]
!....   confor : .TRUE.   conforming surface is considered
!....            .FALSE.  non-conforming surface is considered
!....   x     = nodal coordinates
!.... output
!....   nhseg(nnseg,nsegm) = neighbours segments to each side

IMPLICIT NONE
!     Dummy arguments
INTEGER (kind=4), INTENT(IN) :: nsegm,nnseg,lcseg(:,:)
INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: confor
REAL (kind=8), INTENT(IN), OPTIONAL :: x(:,:)
END SUBROUTINE mastd4
