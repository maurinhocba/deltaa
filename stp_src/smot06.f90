SUBROUTINE smot06(ngaus,nstre,dvolu,strsg,strea, &
                  shape,facta,nvert)
!***********************************************************************
!
!****this routine computes the elemental contribution to nodal stresses
!    for two dimensional elements
!       -----> strea
!***********************************************************************
IMPLICIT NONE

!     routine parameters

INTEGER, INTENT(IN) ::  ngaus, & !Number of Gauss points
                        nstre, & !Number of variables to smooth
                        nvert    !Number of Vertex per element

REAL(kind=8), INTENT(IN) :: strsg(nstre,ngaus), & !Gauss variables
                            dvolu(ngaus),       & !Gauss point volumes
                            shape(nvert,ngaus)    !Nodal shape functions
REAL(kind=8), INTENT(OUT) :: facta,              & !nodal factor
                             strea(nstre,nvert)    !nodal contributions

!     local variables

INTEGER (kind=4) i,g,n

REAL (kind=8)  auxva,emass(16),streb(nstre,nvert)

!     evaluates the area

facta = 1d0/SUM(dvolu)

CALL smom06(ngaus,emass,dvolu,shape,nvert)
!     integrate nodal weighted stresses and compute nodal contributions

streb = 0d0

DO n=1,nvert
  DO g=1,ngaus
    auxva = shape(n,g)*dvolu(g)
    DO i=1,nstre
      streb(i,n) = streb(i,n) + auxva*strsg(i,g)
    END DO
  END DO
END DO

CALL proma2(strea,streb,emass,nstre,nvert,nvert)

RETURN
END SUBROUTINE smot06
