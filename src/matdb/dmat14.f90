 SUBROUTINE dmat14(c,prop,chi,stres,defps,efpst,d,newmt)

 IMPLICIT NONE

 LOGICAL, INTENT(IN OUT) :: newmt
 REAL (kind=8), INTENT(IN) :: c(4),prop(5),chi(12),stres(3),defps,efpst
 REAL (kind=8), INTENT(OUT) :: d(3,3)

 INTEGER (kind=4), SAVE :: is
 REAL (kind=8), SAVE :: de(3,3)=0d0,di11,di12,di22,di33,c0,expn, &
                        aprim,e0,rfs,al(4)
 REAL (kind=8) :: a(3),s(3),yield,g2,beta,th11,th12,th22,th33,  &
                  ti11,ti12,ti22,ti33,aux,dlamb


 !First TASK : Set Constants for a new material

 IF(newmt)THEN
   de(1,1) = c(1)
   de(1,2) = c(2)
   de(2,2) = c(3)
   de(3,3) = c(4)
   aux = 1d0/(c(1)*c(3) - c(2)*c(2))
   di11 = c(3)*aux
   di12 =-c(2)*aux
   di22 = c(1)*aux
   di33 = 1d0/c(4)

   is   = INT(prop(5))   !Isotropic hardening model
   IF( is > 0 )THEN
     c0   = prop(1)        !Initial yield or C0 constant
     expn = prop(3)        !exponent for non-linear hardening
     rfs  = prop(4)        !residual flow stress

     SELECT CASE (is)
     CASE (1)
       aprim = 0d0
     CASE (2)
       aprim = prop(2)
     CASE DEFAULT
       e0    = prop(2)        !non-linear hardening
     END SELECT

     al(1) = chi(2)+chi(3)  !g+h
     al(2) = -chi(3)        !-h
     al(3) = chi(1)+chi(3)  !f+h
     al(4) = 2d0*chi(6)     !2n
   END IF
   newmt = .FALSE.
   RETURN
 END IF

 ! Second TASK : compute tangent constitutive matrix

 IF( defps == 0d0)THEN
   d = de
 ELSE
   ! equivalente stress
   SELECT CASE (is)
   CASE (1)
     yield = c0
     !aprim = 0d0
   CASE (2)
     yield = c0 + aprim*efpst     !linear hardening
     !aprim = constant
   CASE (3)
     yield = c0*(e0+efpst)**expn            !non-linear (exponential) hardening
     aprim = expn*c0/(e0+efpst)**(1d0-expn) !derivative
   CASE (4)
     yield = c0+e0*efpst+(rfs-c0)*(1d0-1d0/EXP(expn*efpst)) !linear + saturation law hardening
     aprim = e0 + (rfs-c0)*expn/EXP(expn*efpst)             !derivative
   END SELECT

   ! Derivative or yield function
   s(1) = (al(1)*stres(1) + al(2)*stres(2))/yield
   s(2) = (al(2)*stres(1) + al(3)*stres(2))/yield
   s(3) = al(4)*stres(3)/yield
   dlamb = defps/yield
   g2 = 1d0 - aprim*dlamb/1.5d0                 !gamma2
   beta = 4d0/9d0*aprim/g2                      !beta/yield^2
   ti11 = di11 + dlamb/1.5d0                    !theta^(-1)
   ti12 = di12 - dlamb/3.0d0
   ti22 = di22 + dlamb/1.5d0
   ti33 = di33 + dlamb*2d0
   aux = ti11*ti22 - ti12*ti12                  !theta
   th11 = ti22/aux
   th12 =-ti12/aux
   th22 = ti22/aux
   th33 = ti33/aux
   a(1) = th11*s(1) + th12*s(2)                 !theta*s
   a(2) = th22*s(2) + th12*s(1)
   a(3) = th33*s(3)
   aux = a(1)*s(1)+a(2)*s(2)+a(3)*s(3) + beta   !s*theta*s + b
   d(1,1) = th11 - a(1)*a(1)/aux
   d(1,2) = th12 - a(1)*a(2)/aux
   d(1,3) =      - a(1)*a(3)/aux
   d(2,2) = th22 - a(2)*a(2)/aux
   d(2,3) =      - a(2)*a(3)/aux
   d(3,3) = th33 - a(3)*a(3)/aux
 END IF
 RETURN
 END SUBROUTINE dmat14
