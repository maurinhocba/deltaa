 SUBROUTINE elpms6(nlayr,amb,dm,stp,stp2,str,a1,a2,dt1,dt2, &
                   dmatx,thick)
 !***********************************************************************
 !
 !**** this routine compute the elasto plastic moduli for shells (elem.6)
 !     Layered model
 !
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4) nlayr
 REAL    (kind=8) amb,dm(*),str(5,nlayr),stp(6,nlayr), &
                  stp2(6,nlayr),a1(3),a2(3),dt1(3),dt2(3),dmatx(8,8),thick
 !***  local variables
 INTEGER (kind=4) i,j,l
 REAL   (kind=8) f2,lb,zl,zl2,dt,tn,efpst,                  &
                 stra1(6),stran(3),u2(3),aux,aux1,d(5,5)
 REAL (kind=8), SAVE :: cm(5,5) = 0d0
 LOGICAL :: plast,dum

 !***
 tn = thick*amb     !present thickness
 dt = tn/nlayr      !layer thickness
 zl = -(tn+dt)/2d0  !initializes
 plast = ANY(stp(6,:) > 0d0 )      !if step is plastic in some layer
 ! mid-surface metric tensor
 stra1(1) = DOT_PRODUCT(a1,a1)
 stra1(2) = DOT_PRODUCT(a2,a2)
 stra1(3) = DOT_PRODUCT(a1,a2)
 stra1(4) = DOT_PRODUCT(dt1,a1)
 stra1(5) = DOT_PRODUCT(dt2,a2)
 stra1(6) = DOT_PRODUCT(dt1,a2) + DOT_PRODUCT(dt2,a1)

 IF( plast )THEN
   DO l=1,nlayr                        ! for each layer
     zl = zl + dt                      ! Z coordinate of the layer
     zl2 = zl*zl                       ! Z coord squared
     stran = stra1(1:3)+2d0*stra1(4:6)*zl  ! layer U^2
     ! inverse of U^2
     aux1  = (stran(1)*stran(2)-stran(3)**2)
     u2(1) = stran(2)/aux1
     u2(2) = stran(1)/aux1
     u2(3) = -stran(3)/aux1

     ! compute layer constitutive matrix (tangent algorithmic)
     lb = stp(6,l)              !
     efpst = stp2(6,l) + lb       !Total Eff. plastic strain
     CALL dmatx6(dum,dum,str(1,l),lb,efpst,d(1,1),.FALSE.)

     cm(1,1) = d(1,1)*u2(1)*u2(1) + 2d0*d(1,3)*u2(1)*u2(3) &
             + d(3,3)*u2(3)*u2(3)

     cm(1,2) = d(1,3)*u2(1)*u2(3) + d(1,2)*u2(1)*u2(2)     &
             + d(3,3)*u2(3)*u2(3) + d(2,3)*u2(3)*u2(2)

     cm(1,3) = ( d(1,1)*u2(1)*u2(3) + d(1,2)*u2(1)*u2(3)   &
             + d(1,3)*( u2(1)*u2(1) + u2(1)*u2(2) + u2(3)*u2(3) )   &
             + d(2,3)*u2(3)*u2(3)                                   &
             + d(3,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/2d0

     cm(2,2) = d(2,2)*u2(2)*u2(2) + 2d0*d(2,3)*u2(2)*u2(3) &
             + d(3,3)*u2(3)*u2(3)

     cm(2,3) = ( d(1,2)*u2(2)*u2(3) + d(1,3)*u2(3)*u2(3)   &
             + d(2,2)*u2(2)*u2(3)                                   &
             + d(2,3)*( u2(1)*u2(2) + u2(2)*u2(2) +u2(3)*u2(3) )    &
             + d(3,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/2d0

     cm(3,3) = (d(3,3)*(u2(1)*u2(1)                        &
             + 2d0*u2(1)*u2(2)+u2(2)*u2(2))                         &
             + d(1,1)*u2(3)*u2(3) + 2d0*d(1,2)*u2(3)*u2(3)          &
             + 2d0*d(1,3)*( u2(1)*u2(3) + u2(2)*u2(3) )             &
             + d(2,2)*u2(3)*u2(3)                                   &
             + 2d0*d(2,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/4d0

     cm(4,4) = d(4,4)
     cm(5,5) = d(5,5)

     DO j=1,3
       DO i=j,3
         dmatx(i,j)     = dmatx(i,j)     + cm(i,j)        !membrane
         dmatx(i,j+3)   = dmatx(i,j+3)   + cm(i,j)*zl     !coupling
         dmatx(i+3,j+3) = dmatx(i+3,j+3) + cm(i,j)*zl2    !bending
       END DO
     END DO
     dmatx(7,7) = dmatx(7,7) + cm(4,4)
     dmatx(8,8) = dmatx(8,8) + cm(5,5)
   END DO        !l=1,nlayr

   dmatx(2,4) = dmatx(1,5)
   dmatx(3,4) = dmatx(1,6)
   dmatx(3,5) = dmatx(2,6)

   ! scale D matrix
   aux = thick/nlayr     ! Initial thickness / number of layer
   dmatx = aux*dmatx
 ELSE  !the point is elastic
   f2 = tn*tn/12d0
   dmatx(1,1) = dm(1)*thick
   dmatx(2,1) = dm(2)*thick
   dmatx(2,2) = dm(3)*thick
   dmatx(3,3) = dm(4)*thick
   dmatx(4,4) = dmatx(1,1)*f2
   dmatx(5,4) = dmatx(2,1)*f2
   dmatx(5,5) = dmatx(2,2)*f2
   dmatx(6,6) = dmatx(3,3)*f2
   dmatx(7,7) = dm(5)*thick
   dmatx(8,8) = dmatx(7,7)
   !CALL modul6(d,t11,t12,t21,t22,def,spc)
 END IF

 RETURN
 END SUBROUTINE elpms6
