SUBROUTINE readcw (areas,coord,ndime,ntype)
!
!   read surface data to compute surface contact wearing
!

IMPLICIT NONE
INTEGER (kind=4), INTENT(IN) :: ndime,ntype
REAL (kind=8), INTENT(IN) :: coord(:,:)
REAL (kind=8), INTENT(OUT) :: areas(:)


INTEGER (kind=4) :: i,j,k,nsegm
INTEGER (kind=4), ALLOCATABLE :: lcseg(:,:)
REAL (kind=8) :: x(ndime,ndime),t(ndime,ndime),a3
REAL (kind=8), PARAMETER :: pi = 3.141592654
LOGICAL :: iscod

areas = 0d0
DO
  READ(45) nsegm,iscod
  IF( nsegm == 0 ) EXIT
  ALLOCATE( lcseg(ndime,nsegm) )
  READ(45) (lcseg(:,i),i=1,nsegm)
  DO i=1,nsegm
    x = coord(:,lcseg(:,i))
    t(:,1) = x(:,2) - x(:,1)
    IF( ndime == 3 )THEN
      t(:,2) = x(:,3) - x(:,1)
      t(1,3) = t(2,1)*t(3,2) - t(3,1)*t(2,2)
      t(2,3) = t(3,1)*t(1,2) - t(1,1)*t(3,2)
      t(3,3) = t(1,1)*t(2,2) - t(2,1)*t(1,2)
      a3 = SQRT(DOT_PRODUCT(t(:,3),t(:,3)))/6d0
    ELSE
      a3 = SQRT(DOT_PRODUCT(t(:,1),t(:,1)))/2d0
    END IF
    IF(ndime == 2 .AND. ntype == 3) a3 = a3*pi*(x(1,2) + x(1,1))
    DO j=1,ndime
      k = lcseg(j,i)
      IF( iscod )THEN
        areas(k) = areas(k)-a3
      ELSE
        areas(k) = areas(k)+a3
      END IF
    END DO
  END DO
  DEALLOCATE( lcseg )
END DO

RETURN
END SUBROUTINE readcw
