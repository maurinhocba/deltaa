      SUBROUTINE mastd1 (ns,ni,x,nhseg)
!
!     Compares sides to detect internal interfaces. If matched sides
!     found assign neighbour segments for contact search
!
!     INPUT
!       ns    : number of boundary sides
!       ni(4,ns): boundary sides data
!         (1) : associated segment
!         (2) : associated side
!         (3) : node of the side
!         (4) : opposite node of the side
!       x(2,:): nodal coordinates
!
!    OUTPUT
!       nhseg(2,:) neighbours segments (connected segments to each side)

      IMPLICIT NONE

!     Dummy arguments
      INTEGER (kind=4), INTENT(IN) :: ns,ni(:,:)
      INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
      REAL (kind=8), INTENT (IN) :: x(:,:)

      END SUBROUTINE mastd1
