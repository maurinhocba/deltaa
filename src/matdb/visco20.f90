 SUBROUTINE visco20(ms,m2,ap,yi,vpar,stra,dtime,g,ier)
 !
 ! compute viscoplastic consistency parameter defined by a power law
 ! using a Newton-Raphson schemme combined with bisection algorithm
 !
 !USE ctrl_db, ONLY: vefac !factor for strain rate smoothing
 IMPLICIT NONE
 !dummy variables
 REAL (kind=8), INTENT(IN)    :: ms,      & !Mises stress
                                 m2,      & !2 * mu
                                 ap,      & !hardening curve tangent
                                 yi,      & !yield stress
                                 vpar(3), & !viscoplastic model parameters
                                 stra(2), & !efective plastic strain & evp strain rate
                                 dtime       !time increment
 REAL (kind=8), INTENT(OUT)   :: g    !viscoplas. consistency parameter
 INTEGER(kind=4), INTENT(OUT) :: ier  !if no convergence then is TRUE
 !local variables
 REAL (kind=8), PARAMETER :: vefac = 0d0
 REAL (kind=8) :: gmax,  & !radial return plastic consist. parameter (maximum)
                  fbar,  & !mises regularized function
                  fprm,  & !1st derivative of mises regularized function
                  dg,    & !increment in viscoplastic consist parameter
                  e_vp,  & !viscoplastic effective strain
                  r_vp     !viscoplastic effective strain rate
 REAL (kind=8) :: ftrl,m2a,r_vpi,go !aux
 INTEGER(kind=4) :: cont   !step counter

 REAL (kind=8)            :: a=0d0,b=1d0,c=1d0,d=0d0 ! deriv parameters
 REAL (kind=8), PARAMETER :: r23 = 0.816496580927726d0, & ! sqrt(2/3)
                             tol = 1d-6
 INTEGER(kind=4), PARAMETER :: ncon = 30

 !elastic-plastic consistency parameter (maximum) from radial return
 gmax = 3d0*(ms - r23*yi)/(2d0*ap+3d0*m2) ! 3*f_trial / (2*A'+6*mu)

 ier  = 0        ! initialize
 e_vp = stra(2)  ! old equivalent strain
 r_vp = stra(1)  ! old equivalent strain rate 
 ftrl = ms - r23*yi        ! f* = |s*| - sqrt(2/3)*sy0
 m2a  = m2 + (2d0/3d0)*ap  ! 2*mu + (2/3)*A'
 g  = 0d0
 go = HUGE(1d0)
 cont = 1

 DO
   ! regularized mises criteria
   fbar = ftrl - m2a*e_vp - vpar(1)*(e_vp**vpar(2))*(r_vp**vpar(3))
   ! first derivative of regularized mises function 
   ! IF( ABS( fbar / ftrl ) <= tol )EXIT
   IF( ABS(g-go/gmax) <= tol .OR. cont > ncon )EXIT
   ! derivative of hardening term
   IF( vpar(2) > 0d0 )THEN
     a = 2d0*r23*vpar(2)*(e_vp**(vpar(2)-1d0))
     c = e_vp**vpar(2)
   END IF
   ! derivative of rate sensivity term
   IF( vpar(3) > 0d0 )THEN
     b = r_vp**vpar(3)
     d = (r23/dtime)*vpar(3)*(r_vp**(vpar(3)-1d0))
   END IF
   ! evaluates first derivative
   fprm = - m2 - (2d0/3d0)*ap - vpar(1)*(a*b+c*d) !tangent to curve fbar(vp)
   dg = fbar/fprm !increment in viscoplastic parameter
   go = g         !save previous parameter
   g  = g - dg     !corrected viscoplastic parameter
   
   DO     ! bisection if g is greather gmax
     IF( g < gmax )EXIT 
     g = g / 2d0  
   END DO
      
   !update variables
   e_vp  = stra(1) + r23 * g  !new equivalent strain
   r_vpi = (r23 * g) / dtime  !new equivalent strain rate
   r_vp  = (r_vpi + vefac*r_vp) / (1d0 + vefac) !new smoothed eq. strain rate

   cont = cont + 1 
 END DO

 IF(g > gmax) g = gmax

 !IF (ier /= 0) WRITE(55,999) ncon,gmax,g
 !999 FORMAT('Viscoplastic parameter correction no convergence',/,  &
 !          'iter = ',     i3,/,                                   &
 !          'gmax = ', f16.12,/,                                   &
 !          'g_vp = ', f16.12,/,                                   &
 !          'CALLED FROM VISCO20',/)

 RETURN
 END SUBROUTINE visco20
