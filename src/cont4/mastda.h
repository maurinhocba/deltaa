SUBROUTINE mastda (ns,ni,x,nhseg)
!
!     Compares sides to detect internal interfaces. If matched sides
!     found assign neighbour segments for contact search
!
!     INPUT
!       ns    : number of boundary sides
!       ni(6,ns): boundary sides data
!         (1) : associated segment
!         (2) : associated side
!         (3) : first node of the side
!         (4) : second node of the side
!       x(3,:): nodal coordinates
!
!    OUTPUT
!       nhseg(3,:) neighbours segments (conected segments to each side)

IMPLICIT NONE

!     Dummy arguments
INTEGER (kind=4), INTENT(IN) :: ns,ni(:,:)
INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
REAL (kind=8), INTENT (IN) :: x(:,:)
END SUBROUTINE mastda
