SUBROUTINE chckar(nsegm,nnseg,lnods,x,label)
!
!     Checks aspect ratio of triangles defining a master surface
!     writes warnings when a specified value is exceeded
!
IMPLICIT NONE
INTEGER (KIND=4), INTENT(IN) :: nsegm,nnseg,lnods(:,:),label(:)
REAL (KIND=8), INTENT(IN) :: x(:,:)

INTEGER i,j,k,jm(1),jj
REAL (KIND=8) :: y(3,nnseg),s(3,nnseg),l(nnseg),lmax,h(3),ar,proy

IF(nnseg /= 3) RETURN    !only for triangles

DO i=1,nsegm                        !for each segment
  y=x(1:3,lnods(1:nnseg,i))         !nodal coordinates
  DO j=1,nnseg                      !for each node of the triangle
    k = MOD(j,nnseg) + 1               !next node
    s(1:3,j) = y(1:3,k) - y(1:3,j)     !vector j->k
    l(j) = SQRT(DOT_PRODUCT(s(1:3,j),s(1:3,j)))  !length of side j->k
  END DO
  jm   = MAXLOC(l) !position of longest side of the triangle
  jj = jm(1)       !
  k = MOD(jj,nnseg) + 1  !node k
  lmax = l(jj)           !length of longest side
  proy = DOT_PRODUCT(s(1:3,jj),s(1:3,k))/lmax**2 !Proyection of side k
  h = s(1:3,k) - proy*s(1:3,jj)      !Orthogonal projection
  ar = lmax/SQRT(DOT_PRODUCT(h,h))   !aspect ratio
  IF( ar > 20d0 )THEN
     WRITE(55,"(' WARNING: Aspect ratio excessive in segment',          &
              &  4i7,' Value: ',f7.2)")i,label(lnods(1:3,i)),ar
     WRITE(55,"(20X,'Internal ',3i7)")lnods(1:3,i)  !internal numeration
     WRITE(*,"(' WARNING: Aspect ratio excessive in segment',           &
              &  4i7,' Value: ',f7.2)")i,label(lnods(1:3,i)),ar
  END IF
END DO
RETURN
END SUBROUTINE chckar
