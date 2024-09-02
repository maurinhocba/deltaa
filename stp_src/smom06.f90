SUBROUTINE smom06(ngaus,emass,dvolu,shape,nvert)
!**********************************************************************
!
 !****this routine evaluates the element smoothing matrix and inverts it
!    for element 06
!**********************************************************************
IMPLICIT NONE

!     routine variables

INTEGER, INTENT(IN) :: ngaus, & !number of integration points
                       nvert    !number of vertex nodes
REAL(kind=8), INTENT(IN) :: dvolu(ngaus),      & !Gauss point volumes
                            shape(nvert,ngaus)   !Shape functions
REAL(kind=8), INTENT(OUT) :: emass(nvert,nvert)  !inverse of mass matrix

!     local variables

INTEGER (kind=4) i,j,g
REAL    (kind=8) a(4,4),auxva,t1,t2,t3,t4,deter,denom

!***consistent smoothing matrix

IF(nvert == 3) THEN     !for triangles
  auxva = SUM(dvolu)    !total volume
  a(1,1) = auxva/6d0
  a(1,2) = auxva/12d0
  a(1,3) = auxva/12d0
  a(2,2) = auxva/6d0
  a(2,3) = auxva/12d0
  a(3,3) = auxva/6d0

ELSE IF(nvert == 4)THEN  !for quadrilaterals

  !        compute standard mass matrix
  a = 0d0                !initializes

  DO g=1,ngaus           !for each gauss point
    DO i=1,nvert         !for each vertex
      auxva = shape(i,g)*dvolu(g)
      DO j=i,nvert
        a(i,j) = a(i,j) + auxva*shape(j,g)
      END DO
    END DO
  END DO

END IF

!              assign the symmetric part
DO i=1,nvert
  DO j=1,i-1
    a(i,j) = a(j,i)
  END DO
END DO

IF(nvert == 4) THEN  ! For quadrilaterals

  !*** inverse of a 4*4 matrix

  t1 =   a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)    &
       + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)    &
       - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
  t2 = - a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)    &
       - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)    &
       + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
  t3 = + a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)    &
       + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)    &
       - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
  t4 = - a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)    &
       - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)    &
       + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
  deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4

  denom=1./deter
  emass(1,1) = t1*denom
  emass(2,1) = t2*denom
  emass(3,1) = t3*denom
  emass(4,1) = t4*denom
  emass(1,2)=(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)         &
              - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)         &
              + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
  emass(2,2)=(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)         &
              + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)         &
              - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
  emass(3,2)=(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)         &
              - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)         &
              + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
  emass(4,2)=(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)         &
              + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)         &
              - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
  emass(1,3)=(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)         &
              + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)         &
              - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
  emass(2,3)=(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)         &
              - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)         &
              + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
  emass(3,3)=(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)         &
              + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)         &
              - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
  emass(4,3)=(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)         &
              - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)         &
              + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
  emass(1,4)=(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)         &
              - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)         &
              + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
  emass(2,4)=(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)         &
              + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)         &
              - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
  emass(3,4)=(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)         &
              - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)         &
              + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
  emass(4,4)=(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)         &
              + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)         &
              - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

ELSE   !for triangles

  deter = a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)                 &
         +a(3,1)*a(2,1)*a(2,3)-a(1,3)*a(2,2)*a(3,1)                 &
         -a(2,3)*a(3,2)*a(1,1)-a(3,3)*a(1,2)*a(2,1)
  denom = 1./deter

  emass(1,1) = denom*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
  emass(2,1) = denom*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
  emass(3,1) = denom*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
  emass(1,2) = denom*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
  emass(2,2) = denom*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
  emass(3,2) = denom*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
  emass(1,3) = denom*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
  emass(2,3) = denom*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
  emass(3,3) = denom*(a(1,1)*a(2,2)-a(2,1)*a(1,2))
END IF

END SUBROUTINE smom06
