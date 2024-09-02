 SUBROUTINE dmatx6(c,prop,stres,defps,efpst,d,newmt)

 IMPLICIT NONE

 LOGICAL, INTENT(IN OUT) :: newmt
 REAL (kind=8), INTENT(IN) :: c(5),prop(17),stres(5),defps,efpst
 REAL (kind=8), INTENT(OUT) :: d(5,5)

 REAL (kind=8), SAVE :: de(5,5)=0d0,di11,di12,di22,di33,di44,c0,expn, &
                        aprim,e0,rfs,al(6)
 REAL (kind=8) :: a(5),s(5),syield,g2,beta,th11,th12,th22,th33,th44,  &
                  ti11,ti12,ti22,ti33,ti44,aux,dlamb


 !First TASK : Set Constants for a new material

 IF(newmt)THEN
   de(1,1) = c(1)
   de(1,2) = c(2)
   de(2,2) = c(3)
   de(3,3) = c(4)
   de(4,4) = c(5)
   de(5,5) = c(5)
   aux = 1d0/(c(1)*c(3) - c(2)*c(2))
   di11 = c(3)*aux
   di12 =-c(2)*aux
   di22 = c(1)*aux
   di33 = 1d0/c(4)
   di44 = 1d0/c(5)
   c0    = prop(1)
   expn  = prop(3)
   rfs   = prop(4)
   IF( expn == 0d0)THEN
     aprim = prop(2)
   ELSE
     e0 = prop(2)
   END IF
   al(1:6) = prop(12:17)
   !newmt = .FALSE.
   RETURN
 END IF

 ! Second TASK : compute tangent constitutive matrix

 IF( defps == 0d0)THEN
   d = de
 ELSE
   ! equivalent stress
   IF( expn == 0d0)THEN
     syield = c0 + aprim*efpst     !linear hardening
   ELSE IF( rfs == 0d0 )THEN
     syield = c0*(e0+efpst)**expn            !non-linear (exponential) hardening
     aprim = expn*c0/(e0+efpst)**(1d0-expn)  !derivative
   ELSE
     syield = c0+e0*efpst+(rfs-c0)*(1d0-1d0/EXP(expn*efpst)) !linear + saturation law hardening
     aprim = e0 + (rfs-c0)*expn/EXP(expn*efpst)            !derivative
   END IF
   ! Derivative or yield function
   s(1) = (al(1)*stres(1) + al(2)*stres(2))/syield
   s(2) = (al(2)*stres(1) + al(3)*stres(2))/syield
   s(3) = al(4)*stres(3)/syield
   s(4) = al(5)*stres(4)/syield
   s(5) = al(6)*stres(5)/syield
   dlamb = defps/syield
   g2 = 1d0 - aprim*dlamb/1.5d0                 !gamma2
   beta = 4d0/9d0*aprim/g2                      !beta/syield^2
   ti11 = di11 + dlamb/1.5d0                    !theta^(-1)
   ti12 = di12 - dlamb/3.0d0
   ti22 = di22 + dlamb/1.5d0
   ti33 = di33 + dlamb*2d0
   ti44 = di44 + dlamb*2d0
   aux = ti11*ti22 - ti12*ti12                  !theta
   th11 = ti22/aux
   th12 =-ti12/aux
   th22 = ti22/aux
   th33 = ti33/aux
   th44 = ti44/aux
   a(1) = th11*s(1) + th12*s(2)                 !theta*s
   a(2) = th22*s(2) + th12*s(1)
   a(3) = th33*s(3)
   a(4) = th44*s(4)
   a(5) = th44*s(5)
   aux = a(1)*s(1)+a(2)*s(2)+a(3)*s(3)+a(4)*s(4)+a(5)*s(5) + beta   !s*theta*s + b
   d(1,1) = th11 - a(1)*a(1)/aux
   d(1,2) = th12 - a(1)*a(2)/aux
   d(1,3) =      - a(1)*a(3)/aux
   d(2,2) = th22 - a(2)*a(2)/aux
   d(2,3) =      - a(2)*a(3)/aux
   d(3,3) = th33 - a(3)*a(3)/aux
   d(4,4) = th44 - a(4)*a(4)/aux
   d(5,5) = th44 - a(5)*a(5)/aux
 END IF
 RETURN
 END SUBROUTINE dmatx6
