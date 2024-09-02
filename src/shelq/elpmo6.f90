 SUBROUTINE elpmo6(matty,amb,prop,dmatx,str,stp,ehist,a1,a2,tgp,d)
 !***********************************************************************
 !
 !**** this routine compute the elasto plastic moduli for shells (elem.6)
 !
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4) matty
 REAL    (kind=8) amb,prop(*),dmatx(*),str(8),ehist(5),stp(8),    &
                  a1(3),a2(3),tgp(3),d(8,8)
 !***  local variables
 INTEGER (kind=4) i
 REAL   (kind=8) pymat(3),khard,epstr,expo,yield,kp,f1,coefm,coenm,&
                 st(8),dsp(8),spc,ab1(3),ab2(3),ap1(3),ap2(3),rfs, &
                 t11,t12,t21,t22,deter,l11,l12,l21,l22,def
 SAVE pymat
 DATA pymat /0.3333333333333333d0,1d0,2d0/
 !***
 def = amb
 spc = ehist(4)+ehist(5)
 IF(matty == 0) THEN  !large strain
   !       curvilinear system with unit metric (schmidt)
   ab1 = a1
   CALL vecuni(3,ab1,t22)
   t11 = DOT_PRODUCT(ab1,a2)
   DO i = 1,3
     ab2(i) = a2(i) - t11*ab1(i)
   END DO
   CALL vecuni(3,ab2,t22)
   !       evaluates contravariant base
   t11  =  ab2(2)*tgp(3) - ab2(3)*tgp(2)
   t22  = -ab1(2)*tgp(3) + ab1(3)*tgp(2)
   t12  =  ab1(2)*ab2(3) - ab1(3)*ab2(2)
   deter =  ab1(1)*t11 + ab2(1)*t22 + tgp(1)*t12
   ap1(1) = t11/deter
   ap1(2) = (-ab2(1)*tgp(3)  + ab2(3)*tgp(1)) /deter
   ap1(3) = ( ab2(1)*tgp(2)  - ab2(2)*tgp(1)) /deter
   ap2(1) = t22/deter
   ap2(2) = ( ab1(1)*tgp(3)  - ab1(3)*tgp(1)) /deter
   ap2(3) = (-ab1(1)*tgp(2)  + ab1(2)*tgp(1)) /deter

   l11 = DOT_PRODUCT(ap1,a1)
   l12 = DOT_PRODUCT(ap1,a2)
   l21 = DOT_PRODUCT(ap2,a1)
   l22 = DOT_PRODUCT(ap2,a2)

   deter = l11*l22-l12*l21

   t11 =  l22/deter
   t12 = -l12/deter
   t21 = -l21/deter
   t22 =  l11/deter
 END IF

 IF(spc > 0) THEN
   f1    = prop(6)*def**2
   coefm = prop(8)/def**2
   coenm = prop(9)/def
   IF(matty == 0) THEN
     !         plastic Almansi strain (step n)
     dsp(1)=stp(1)*t11**2 +stp(2)*t21**2 +stp(3)*2d0*t11*t21
     dsp(2)=stp(1)*t12**2 +stp(2)*t22**2 +stp(3)*2d0*t12*t22
     dsp(3)=stp(1)*t11*t12+stp(2)*t21*t22+stp(3)*(t11*t22 +t21*t12)
     dsp(4)=stp(4)*t11**2 +stp(5)*t21**2 +stp(6)*2d0*t11*t21
     dsp(5)=stp(4)*t12**2 +stp(5)*t22**2 +stp(6)*2d0*t12*t22
     dsp(6)=stp(4)*t11*t12+stp(5)*t21*t22+stp(6)*(t11*t22 +t21*t12)
     dsp(7)=stp(7)*t11/def+stp(8)*t21/def
     dsp(8)=stp(7)*t12/def+stp(8)*t22/def
     !                       Kirchhoff stresses
     st(1)=str(1)*l11**2 +str(2)*l21**2 +2*str(3)*l11*l21
     st(2)=str(1)*l12**2 +str(2)*l22**2 +2*str(3)*l12*l22
     st(3)=str(1)*l11*l12+str(2)*l21*l22+str(3)*(l11*l22 +l21*l12)
     st(4)=str(4)*l11**2 +str(5)*l21**2 +2*str(6)*l11*l21
     st(5)=str(4)*l12**2 +str(5)*l22**2 +2*str(6)*l12*l22
     st(6)=str(4)*l11*l12+str(5)*l21*l22+str(6)*(l11*l22 +l21*l12)
     st(7)=str(7)*l11*def+str(8)*l21*def
     st(8)=str(7)*l12*def+str(8)*l22*def
   ELSE
     dsp = stp
     st  = str
     !
     IF(matty == 1 .OR. matty == 2)THEN
       st(1) = st(1)* 6.45d0    !3.4d0
       st(5) = st(5)* 0.13d0    !0.525D0
     END IF
   END IF
   !       IF kinematic hardening
   khard = prop(5)
   IF(khard > 0) THEN
     st(1) = st(1) - khard*dsp(1)
     st(2) = st(2) - khard*dsp(2)
     st(3) = st(3) - khard*dsp(3)
     st(4) = st(4) - khard*dsp(4)*f1
     st(5) = st(5) - khard*dsp(5)*f1
     st(6) = st(6) - khard*dsp(6)*f1
     st(7) = st(7) - khard*dsp(7)*def
     st(8) = st(8) - khard*dsp(8)*def
   END IF
   epstr = ehist(1)
   expo  = prop(3)
   rfs   = prop(4)
   IF(expo == 0d0) THEN
     yield = prop(1) + prop(2)*epstr
     kp = prop(2)
   ELSE IF (rfs == 0d0 )THEN
     IF(epstr > 0.02d0) THEN
       yield = prop(1)*(prop(2)+epstr)**expo
       kp = prop(1)*expo/(prop(2)+epstr)**(1d0-expo)
     ELSE
       kp = prop(1)*expo/(prop(2)+0.02d0)**(1d0-expo)
       yield = prop(1)*(prop(2)+0.02d0)**expo-kp*(0.02d0-epstr)
     END IF
   ELSE
     yield  = prop(1)+prop(2)*epstr+(rfs-prop(1))*(1d0-1d0/EXP(expo*epstr)) !linear + saturation law hardening
     kp = prop(2) + (rfs-prop(1))*expo/EXP(expo*epstr)            !derivative
   END IF
   CALL j2smt9( 8, 3,yield,kp,khard,f1,prop(7),coefm,             &
               coenm,prop(10),prop(13),d,st,ehist(4),pymat)
 ELSE
   CALL modps6(matty,d,dmatx,prop,prop(1))
 END IF
 IF(matty == 0) CALL modul6(d,t11,t12,t21,t22,def,spc)
 RETURN
 END SUBROUTINE elpmo6
