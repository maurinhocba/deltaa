 SUBROUTINE corr20 (eulrf,t,lb3,gausv,is,props,propv,gm,km, &
                    stres,sigma,np,curve,ierr,tcont,dg,fac,alpha,dther)
   !
   ! plastic correction algorithm for linear triangle in plane problems
   ! with or without coupled thermo-mechanical analysis
   !
   USE lispa0
   !USE ctrl_db, ONLY : dtime,therm
   USE mat_dba, ONLY : inte_cr
   IMPLICIT NONE
   !dummy variables for mechanical analysis
   LOGICAL, INTENT(IN) :: eulrf  !TRUE use spatial configuration else use intermediate configuration
   REAL (kind=8), INTENT(IN) :: t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                props(5), & !plastic properties (hardening)
                                propv(3), & !viscoplastic properties
                                gm,km       !elastic properties (km = 3K)
   REAL (kind=8), INTENT(IN OUT) :: gausv(:), & !Internal variables TLF: (1:5)=Fp^-1 (6)=ep
                                                !                   ULF: (1:4)=be^-1 (5)=ep_dot (6)=ep
                                    stres(4)    !Kirchhoff stresses (for post-process)
   REAL (kind=8), INTENT(OUT)    :: sigma(4)    !for internal forces evaluation
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np          !number of points defining curve (is = 5)
   INTEGER(kind=4), INTENT(OUT) :: ierr     !error flag (1)
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve
   !dummy variables for thermo-mechanical analysis only
   REAL (kind=8), OPTIONAL, INTENT(OUT) :: tcont(3)    !temperature dependant contributions for coupled thermo-mechanical
   REAL (kind=8), OPTIONAL, INTENT(IN)  :: fac,      & !(1+nu) for plane strain and axilsymmetry
                                           alpha,    & !thermal expansion coefficient
                                           dther,    & !element temperature increment
                                           dg          !derivative of shear modulus with respect to temperature

   ! local variables
   INTEGER(kind=4) :: i
   REAL (kind=8) :: yield,aprim,efps0,ethm,estr(4),pstr(4),treo,tren,mises
   LOGICAL :: visc !.TRUE. for elastic-visco-plastic material

   !     setup initial yield FUNCTION radius
   efps0 = gausv(6)          !initial (old) Equivalent Plastic Strain
   IF( is == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     yield = inte_cr (curve,np,efps0,i)    !s_y
     aprim = curve(3,i)                    !A'
   ELSE
     CALL isoha14(is,yield,aprim,efps0,props(1),props(2),props(3),props(4))  !compute s_y and A'
   END IF

   IF(eulrf)THEN !update lagrangian formulation (actual configuration) from Garcia Garino tesis
     ethm = 0d0
     !IF( therm )THEN
     !  treo = stres(1)+stres(2)+stres(4) !previous kirchhoff stress trace
     !  ethm = 2d0*fac*alpha*dther        !thermal strain trace * 2
     !END IF

     visc = ( propv(1) > 0d0 .AND. ANY(propv(2:) > 0d0) ) !is TRUE when viscoplast material
     !plastic corrector
     CALL corr20_a (visc,t,lb3,ethm,gausv,aprim,yield, &
                    propv,gm,km,stres,pstr,0d0,ierr)
     sigma = stres !for internal forces evaluation

     !IF( therm )THEN !evaluates thermal dissipation contributions
     !  tren = stres(1)+stres(2)+stres(4) !new kirchhoff stress trace
     !  estr(1) = stres(1)-tren/3d0 ! new deviatoric kirchhoff stress
     !  estr(2) = stres(2)-tren/3d0 !
     !  estr(3) = stres(3)          !
     !  estr(4) = stres(4)-tren/3d0 !
     !  mises = SQRT(estr(1)*estr(1)+estr(2)*estr(2)+estr(4)*estr(4)+ &
     !               2d0*estr(3)*estr(3)) !von mises stress
     !  estr =  estr/gm/2d0  !new deviatoric elastic strain tensor
     !  !store coupled thermo-mechanical factors
     !  tcont    = 0d0 !
     !  tcont(1) = mises*gausv(6)       ! plastic work (s_ij : e_ij)
     !  tcont(2) = (tren-treo)/(km/3d0) ! volume elastic Gough-Joule effect
     !  tcont(3) = dther*dg*(estr(1)*pstr(1)+estr(2)*pstr(2)+ &  !change in plastic work due temperature
     !             estr(4)*pstr(4)+2d0*estr(3)*pstr(3))          !DT * dgm/dT * e_ij : ep_ij
     !END IF

   ELSE !total lagrangian formulation (intermediate configuration) from Crisfield - not thermal yet
     CALL corr20_i (t,lb3,gausv,aprim,yield,gm,km,stres,sigma,ierr)

   END IF

   RETURN
 END SUBROUTINE corr20


