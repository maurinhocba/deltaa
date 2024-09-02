SUBROUTINE smot08(ngaus,nstre,nnode,dvolu,strsg,strea,shape,facta)
!***********************************************************************
!
!****this routine computes the elemental contribution to nodal stresses
!       -----> strea
!***********************************************************************
IMPLICIT NONE

!     routine parameters

INTEGER, INTENT(IN) :: ngaus, & !number of Gauss points per element
                       nstre, & !number of variables to smooth
                       nnode    !number of nodes per element

REAL (kind=8), INTENT(IN) :: strsg(nstre,ngaus), & !Gaussian Variables
                             dvolu(ngaus),       & !Gauss point lengths
                             shape(nnode,ngaus)    !shape functions
REAL (kind=8), INTENT(OUT):: facta,              & !nodal factor
                             strea(nstre,nnode)    !nodal contributions


!     local variables

INTEGER (kind=4) i,g,n

REAL (kind=8)  auxva

!     evaluates the length

facta = 1d0/SUM(dvolu)  !inverse of length

strea = 0d0             !initializes

DO n=1,nnode            !for each node
  DO g=1,ngaus          !from each Gauss point
    auxva = shape(n,g)  !*dvolu(g)
    DO i=1,nstre
      strea(i,n) = strea(i,n) + auxva*strsg(i,g)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE smot08
