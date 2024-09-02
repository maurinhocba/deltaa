 SUBROUTINE corr20 (eulrf,t,lb3,gausv,is,props,propv,gm,km, &
                    stres,sigma,np,curve,ierr,tcont,dg,fac,alpha,dther)
   !
   ! updates stresses for von Mises yield function
   ! plane strain and axilsymmetric problems
   !
   USE lispa0
   USE ctrl_db, ONLY : dtime,therm
   USE mat_dba, ONLY : inte_cr
   IMPLICIT NONE
   !dummy variables for mechanical analysis
   LOGICAL, INTENT(IN) :: eulrf  !TRUE use spatial configuration else use intermediate configuration
   REAL (kind=8), INTENT(IN) :: t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                props(5), & !plastic properties (hardening)
                                propv(3), & !viscoplastic properties
                                gm,km       !elastic properties (km = 3K)
   REAL (kind=8), INTENT(IN OUT) :: gausv(:), & !Internal variables (1:5) =(Fp)^(-1)  6:=ep  7: ep_dot
                                    stres(4)    !Kirchhoff stresses (for post-process)
   REAL (kind=8), INTENT(OUT)    :: sigma(4)    !for internal forces evaluation
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np          !number of points defining curve (is = 5)
   INTEGER(kind=4), INTENT(OUT) :: ierr     !error flag (1)
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve
   !dummy variables for thermo-mechanical analysis
   REAL (kind=8), OPTIONAL, INTENT(OUT) :: tcont(3)    !temperature dependant contributions for coupled thermo-mechanical
   REAL (kind=8), OPTIONAL, INTENT(IN)  :: fac,      & !(1+nu) for plane strain
                                           alpha,    & !thermal expansion coefficient
                                           dther,    & !element temperature increment
                                           dg          !derivative of shear modulus with respect to temperature  

 END SUBROUTINE corr20
