 SUBROUTINE stran6(d0,d,t,j0,stran)

 ! compute the metric tensor in the deformed configuration

 IMPLICIT NONE
 REAL(kind=8), INTENT(IN) :: d0(3,2), & !covariant vectors (original)
                             d(3,2),  & !covariant vectors (deformed)
                             t(3),    & !original normal vector
                             j0         !original jacobian
 REAL(kind=8), INTENT(OUT) :: stran(3)  !surface metric tensor
 !local variables
 REAL(kind=8) :: t1(3),t2(3),l,ji(2,2)

 l = SQRT(DOT_PRODUCT(d0(:,1),d0(:,1)))  !length of first original vector
 t1 = d0(:,1)/l                          !tangent unit vector along 1st local dir
 CALL vecpro(t(1),t1(1),t2(1))           !tangent unit vector along 2nd local dir 
 ji(1,1) =  DOT_PRODUCT(d0(:,2),t2)/j0   !d(xita)/dx
 ji(1,2) = -DOT_PRODUCT(d0(:,2),t1)/j0   !d(xita)/dy
 ji(2,1) = -DOT_PRODUCT(d0(:,1),t2)/j0   !d(eta)/dx
 ji(2,2) =  DOT_PRODUCT(d0(:,1),t1)/j0   !d(eta)/dy
 t1 = d(:,1)*ji(1,1) + d(:,2)*ji(2,1)    !tangent vector (x) in the deformed conf.
 t2 = d(:,1)*ji(1,2) + d(:,2)*ji(2,2)    !tangent vector (y) in the deformed conf.
 stran(1) = DOT_PRODUCT(t1,t1)           !g11
 stran(2) = DOT_PRODUCT(t2,t2)           !g22
 stran(3) = DOT_PRODUCT(t1,t2)           !g12
 RETURN
 END SUBROUTINE stran6
