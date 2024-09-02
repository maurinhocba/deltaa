SUBROUTINE update(istep,itera,lauto,ecdis,neq,arcln,lambd,dlamb,  &
                  disax,displ,ddisp,delta,karcl,piter,diter,newtv,ncdis)
!***********************************************************************
!
!***  this routine calculates the change in load step according to
!     selected path and updates load level and displacement vectors
!
!***********************************************************************
USE lispa0
IMPLICIT NONE
!       routine arguments
INTEGER (kind=4),INTENT(IN) :: istep,itera,neq,lauto,karcl,newtv
INTEGER (kind=4),INTENT(IN OUT) :: ecdis,ncdis
REAL (kind=8),INTENT(IN) :: piter,diter
REAL (kind=8),INTENT(IN OUT) :: dlamb,lambd,arcln,disax(:,:),displ(:),ddisp(:)
REAL (kind=8),INTENT(OUT) :: delta
!      local variables
INTEGER (kind=4) :: i
REAL (kind=8) :: a,b,c,v1,v2,v3,v4,v5,v6,angl1,angl2,root1,root2

INTERFACE
  INCLUDE 'predic.h'
END INTERFACE

IF(itera == 1) THEN                                ! prediction

  ! *** for the first iteration compute delta y ddisp
  CALL predic(istep,lauto,ecdis,neq,arcln,dlamb,                  &
              disax,displ,ddisp,delta,karcl,piter,diter,newtv,ncdis)

  displ = 0d0                        ! initialize increm. displacement

ELSE                                               ! correction
  !  *** for any other iteration compute delta according to path
  SELECT CASE (karcl)

  CASE (0,5)                               ! fixed load or tools

    delta = 0d0

  CASE (1)                                 ! normal plane path

    delta = - DOT_PRODUCT(ddisp,disax(1:neq,1))/                  &
              DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1))

  CASE (2)                                 ! updated normal plane path

    delta = - DOT_PRODUCT(ddisp,displ)/                           &
              DOT_PRODUCT(disax(1:neq,1),displ)

  CASE (3)                                 ! spherical path

    a = 0d0
    b = 0d0
    c = 0d0
    DO i=1,neq
      v1 = disax(i,1)
      v3 = displ(i) + ddisp(i)
      a = a + v1*v1
      b = b + v3*v1
      c = c + v3*v3
    END DO
    c = c - arcln * arcln
    v1 = b * b - a * c
    IF(v1 >= 0d0) THEN               ! REAL roots
      v1 = SQRT(v1)
      root1 = (-b+v1)/a
      root2 = (-b-v1)/a
      angl1 = 0d0
      angl2 = 0d0
      DO i=1,neq                     ! compute angles with previous step
        v1 = displ(i)
        v2 = ddisp(i)
        v3 = disax(i,1)
        v4 = v1 + v2 + root1*v3
        v5 = v1 + v2 + root2*v3
        v6 = disax(i,2)                ! displ of previous step
        angl1 = angl1 + v4*v6
        angl2 = angl2 + v5*v6
      END DO ! i=1,neq
      IF(angl1 >= angl2) THEN          ! choose smaller angle (greater cosine)
        delta = root1
      ELSE
        delta = root2
      END IF
      IF(angl1 < 0d0 .AND. angl2 < 0d0) WRITE(lures,930) ! two negative angles
    ELSE                             ! no solution, take minimum distance
      WRITE(lures,920) SQRT(-v1)/b
      delta = -b/a
      !arcln = arcln/2d0
      !delta = -dlamb + arcln/SQRT(DOT_PRODUCT(disax(1:neq,1),disax(1:neq,1)))
      !ddisp = -displ
    END IF

  CASE (4)                                 ! displacement control

    delta = -ddisp(ecdis)/disax(ecdis,1)

  END SELECT ! (karcl == ?)

  ddisp = ddisp + delta*disax(1:neq,1)
  dlamb = dlamb + delta
  ! correct if excesive correction in second iteration
  ! IF(itera == 2)THEN
  !   a = DOT_PRODUCT(ddisp,ddisp)/DOT_PRODUCT(displ,displ)
  !   IF(a > 1d0)THEN
  !     WRITE(lures,"(' correction in second iteration bigger than ', &
  !             & 'prediction',/)")
  !     ddisp = -displ/2d0
  !     delta = -d/2d0
  !   END IF
  ! END IF
END IF ! (itera == ?)

lambd = lambd + delta                       ! correct load factor
RETURN

!  WRITE(lures,"(' step',i4,' itera',i3,' L',e12.4,' dL',e12.4)")
!        &            istep,      itera,     lambd,      delta
920 FORMAT(1x,'**warning** UPDATE: two complex roots into spherical',  &
    &      ' path.',/,21x,'imaginary part / REAL part =',e15.5,/,      &
    &      21x,'only REAL part taken in root as minimum distance')
930 FORMAT(1x,'**warning** UPDATE: two negatives angles in ',         &
    &  'spherical path' )

END SUBROUTINE update
