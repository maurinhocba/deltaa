 SUBROUTINE corr20 (ntype,t,lb3,gausv,newm,is,ierr,props,gm,km,stres,sigma,np,curve)
   !
   ! updates stresses for von Mises yield function
   ! plane strain and axilsymmetric problems
   !
   USE lispa0
   USE mat_dba, ONLY : inte_cr
   IMPLICIT NONE
   REAL (kind=8), INTENT(IN) :: t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                props(5), & !plastic properties (hardening)
                                gm,km       !elastic properties (km = 3K)
   REAL (kind=8), INTENT(IN OUT) :: gausv(6)  !Internal variables (1:5) =(Fp)^(-1)  6:=ep
   REAL (kind=8), INTENT(OUT) :: stres(4),  & !Kirchhoff stresses (for post-process)
                                 sigma(4)     !2nd Piola-Kirchhoff (for internal forces)
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np,       & !number of points defining curve (is = 5)
                                  ntype       !1:plane stress 2:plane strain 3:axil
   INTEGER(kind=4), INTENT(OUT) :: ierr     !error flag (1)
   LOGICAL, INTENT(IN OUT) :: newm          !new material flag
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve

   ! local variables
   INTEGER(kind=4) :: i
   REAL (kind=8) :: te(2,2),bet(3),delta,r1,r2,stran(4),lc(3), &
                    yield,aprim,efpst,f(5),a,b,c,d,al,de,se,dvs(3),j3

   REAL (kind=8), PARAMETER :: r32 = 1.224744871d0, toler = 1d-4
   REAL (kind=8), SAVE :: c0,c1,c2,c3,m2,m3  !material parameters

   IF( newm )THEN
     newm = .FALSE.
     c0 = props(1)     !c0 constant or Initial yield
     c1 = props(2)     !linear hardening or efref
     c2 = props(3)     !exponent for non-linear hardening
     c3 = props(4)     !residual flow stress
     m2 = gm*2.0d0     !2*Mu to compute stresses
     m3 = gm*3.0d0     !3*Mu to use as auxiliar
   END IF
   !     setup initial yield FUNCTION radius
   efpst = gausv(6)          !initial (old) Equivalent Plastic Strain
   IF( is == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     yield = inte_cr (curve,np,efpst,i)    !s_y
     aprim = curve(3,i)                    !A'
   ELSE
     CALL isoha14(is,yield,aprim,efpst,c0,c1,c2,c3)  !compute s_y and A'
   END IF
   ! compute elastic trial gradient  te + lc(3) from t + lb3
   f = gausv(1:5)        !previous plastic gradient (inverse)
   te(1,1) = t(1,1)*f(1) + t(1,2)*f(2)
   te(2,1) = t(2,1)*f(1) + t(2,2)*f(2)
   te(1,2) = t(1,1)*f(3) + t(1,2)*f(4)
   te(2,2) = t(2,1)*f(3) + t(2,2)*f(4)
   lc(3)   = lb3   *f(5)
   !computes eigen-decomposition and log strains
   stran(1) = DOT_PRODUCT(te(:,1),te(:,1))  !Ce
   stran(2) = DOT_PRODUCT(te(:,2),te(:,2))
   stran(3) = DOT_PRODUCT(te(:,1),te(:,2))
   CALL eige20(stran(1),r1,r2,lc(1),j3,ierr) !lc = deviatoric   j3=J^1/3
   IF( ierr == 1 )RETURN              !error, stop program
   delta = stran(4)*km                !p
   bet   = m2*stran(1:3)              !deviatoric log stresses
   de = 0d0                           !initializes consistency parameter
   se = SQRT(bet(1)*bet(1)+bet(2)*bet(2)+bet(3)*bet(3))*r32  !vMstress
   yield = se - yield       !f
   IF( yield/se > toler )THEN !plastic
     ! radial return (no iteration) in terms of A' at eps old
     de =  yield /(aprim + m3)    !increment in eps
     al = 1d0 -  m3*de/se         !alpha multiplier
     bet = al*bet                 !corrected stresses
     lc  = (lc**al)*j3            !new elastic eigenvalues of U
   ELSE
     lc  = lc*j3               !new elastic eigenvalues of U
   END IF
   ! diagonal stress (2nd P-k) tensor including mean pressure (at NEW intermediate configuration)
   DO i=1,3
     bet(i) = (bet(i)+delta)/lc(i)**2
   END DO
   ! 2nd Piola-kirchoff stress tensor at NEW intermediate configuration
   stres(1) = r1*bet(1)*r1 + r2*bet(2)*r2    !(1,1)
   stres(2) = r2*bet(1)*r2 + r1*bet(2)*r1    !(2,2)
   stres(3) = r1*(bet(1)-bet(2))*r2          !(1,2)
   stres(4) = bet(3)                         !(3,3)
   IF( de > 0d0 ) THEN !updates internal variables
     ! lb'_i = exp[-(1-al)*Ln(lb_i'(trial)) ]
     al = al - 1d0          ! - (1-alpha)
     DO i=1,3
       lc(i) = EXP(al*stran(i))    !exp( -alpha*Ln(lb'_i))
     END DO
     ! exp[(1-al)*li' * Ni x Ni
     a = r1*r1*lc(1) + r2*r2*lc(2)  !11
     b = r1*r2*(lc(1)-lc(2))        !12 = 21
     c = r2*r2*lc(1) + r1*r1*lc(2)  !22
     d = lc(3)                      !33
     ! new plastic gradient (inverse)
     gausv(1) = f(1)*a + f(3)*b     !11
     gausv(2) = f(2)*a + f(4)*b     !21
     gausv(3) = f(1)*b + f(3)*c     !12
     gausv(4) = f(2)*b + f(4)*c     !22
     gausv(5) = f(5)*d              !33
     gausv(6) = gausv(6) + de       !update equivalent plastic strain
     f = gausv(1:5)                 !new plastic gradient
   END IF
   ! 2nd Piola Kirchhoff stress tensor at original configuration
   sigma(1) = stres(1)*f(1)*f(1)+stres(2)*f(3)*f(3)+2d0*stres(3)*f(3)*f(1)
   sigma(2) = stres(1)*f(2)*f(2)+stres(2)*f(4)*f(4)+2d0*stres(3)*f(2)*f(4)
   sigma(3) = stres(1)*f(1)*f(2)+stres(2)*f(3)*f(4)+stres(3)*(f(1)*f(4)+f(2)*f(3))
   sigma(4) = stres(4)*f(5)*f(5)
   ! Kirchhoff stresses at spatial configuration
   stres(1) = sigma(1)*t(1,1)*t(1,1)+sigma(2)*t(1,2)*t(1,2)+2d0*sigma(3)*t(1,2)*t(1,1)
   stres(2) = sigma(1)*t(2,1)*t(2,1)+sigma(2)*t(2,2)*t(2,2)+2d0*sigma(3)*t(2,1)*t(2,2)
   stres(3) = sigma(1)*t(1,1)*t(2,1)+sigma(2)*t(1,2)*t(2,2)+sigma(3)*(t(1,1)*t(2,2)+t(1,2)*t(2,1))
   stres(4) = sigma(4)*lb3*lb3
 RETURN
 END SUBROUTINE corr20

 SUBROUTINE eige20 (stran,r1,r2,lb,j,ierr)

 ! compute eigen-decomposition
 ! and log strains

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT (IN OUT) :: lb(3),  & !IN: eigenvalues of U   OUT: deviatoric eigenvalues of U
                                  stran(4)  !(IN) C  (OUT) Log deviatoric strains & jac
 REAL(kind=8), INTENT (OUT) :: r1,r2,     & !components of first eigenvector
                               j            !j^(1/3)
 INTEGER (kind=4), INTENT(OUT) :: ierr      !error flag

 ! local variables
 REAL (kind=8) :: c,r,l1,l2

 !compute eigenvalues
 c = (stran(1)+stran(2))/2d0                        !center of circle
 r = SQRT((stran(1)-stran(2))**2/4d0+stran(3)**2)   !circle radius
 l1 = c+r                                 !first (maximum) eigenvalue
 l2 = c-r                                 !second (minimum) eigenvalue
 !compute eigenvectors
 c = stran(1) - l1                        !first diagonal element
 IF( c < -1d-15)THEN                      !check
   r = -stran(3)/c                        ! v = (r,1)
   c = SQRT(r**2+1d0)                     ! eigenvector length
   r1 = r/c                               !first component of eigevector
   r2 = 1d0/c                             !first component of eigevector
 ELSE
   r1 = 1d0                               !local direction is the vector
   r2 = 0d0
 END IF
 !Compute eigenvalues and principal strains
 IF( l2 < 0d0 .OR. lb(3) < 0d0 )THEN
   WRITE(55,*)" Too distorted mesh, negative eigen-value detected"
   WRITE(55,*)' U^2 = ',stran
   WRITE(55,*)' l1 = ',l1,' l2 = ',l2
   WRITE(55,*)' CALLED FROM CORR20'
   ierr = 1
   RETURN
 END IF
 lb(1) = SQRT(l1)     !previous values where eigenvalues of U^2
 lb(2) = SQRT(l2)
 j = (lb(1)*lb(2)*lb(3))**(1d0/3d0)   !elastic jacobian^(1/3)
 lb = lb/j                            !deviatoric eigenvalues
 stran(1) = LOG(lb(1))                !dev log strains
 stran(2) = LOG(lb(2))
 stran(3) = LOG(lb(3))
 stran(4) = LOG(j)                    !Ln(jac^(1/3))
 RETURN
 END SUBROUTINE eige20
