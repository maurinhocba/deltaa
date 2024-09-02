 SUBROUTINE shapes(nnode,ngaus,shape,deriv,posgp,weigh,ndim,etype)
 !
 !  compute Shape functions
 !
 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT(IN) :: nnode,ngaus,ndim,etype
 REAL (kind=8), INTENT(OUT) :: shape(nnode,ngaus),posgp(ndim,ngaus), &
                               deriv(nnode,ndim,ngaus),weigh(ngaus)
 !     local variables
 INTEGER (kind=4) g

 REAL (kind=8), PARAMETER :: a = 0.577350269189626, b = 0.774596669241483
 REAL (kind=8) :: c,d,e,cd,c2,d2,cc,dd,ccd,cdd,cd2
 REAL (kind=8) :: s,t,q,sm,tm,zm,sp,tp,zp,bf, &
                  omg, omh, omr, opg, oph, opr, tpgphpr, tmgphpr, tmgmhpr, tpgmhpr, &
                  tpgphmr, tmgphmr, tmgmhmr, tpgmhmr
 !REAL (kind=8) :: xita,eta,zeta,zzet,l1,l2

 INTERFACE
      SUBROUTINE gaussq    (ngaus ,posgp ,weigp )
      IMPLICIT NONE
      INTEGER (kind=4) ngaus
      REAL    (kind=8) posgp(ngaus),weigp(ngaus)
      END SUBROUTINE gaussq
 END INTERFACE

 SELECT CASE (ndim)  !element dimension
 CASE (1)            !Linear
   SELECT CASE (ngaus)
   CASE (1)
     posgp(1,1) = 0d0
     weigh(1) = 2d0
   CASE (2)
     posgp(1,1:2) = (/ -a, a /)
     weigh(1:2) = 1d0
   CASE (3)
     posgp(1,1:3) = (/ -b, 0d0, b /)
     weigh(1:3) = (/ 5d0/9d0, 8d0/9d0, 5d0/9d0 /)
   END SELECT
   SELECT CASE (nnode)
   CASE (2)
     DO g=1,ngaus
       c = (1d0 - posgp(1,g))/2d0
       d =  1d0 - c
       shape(1:2,g) = (/ c, d /)
       deriv(1,1,g) = -0.5d0
       deriv(2,1,g) =  0.5d0
     END DO
   CASE (3)
     DO g=1,ngaus
       shape(1,g)   =  posgp(1,g)*(-1d0+posgp(1,g))/2d0
       shape(2,g)   =  1d0-posgp(1,g)**2
       shape(3,g)   =  posgp(1,g)*( 1d0+posgp(1,g))/2d0
       deriv(1,1,g) =  posgp(1,g)-0.5d0
       deriv(2,1,g) = -posgp(1,g)*2d0
       deriv(3,1,g) =  posgp(1,g)+0.5d0
     END DO
   END SELECT

 CASE (2)
   SELECT CASE (nnode)
   CASE (4,8,9)  !Quadrilateral
     SELECT CASE (ngaus)
     CASE (1)
       posgp = 0d0
       weigh(1) = 4d0
     CASE (4)
       posgp = RESHAPE((/ -a,-a, -a, a, a,-a , a,a  /),(/2,4/))
       weigh(1:4) = 1d0
     CASE (9)
       posgp = RESHAPE((/  -b, -b,    -b, 0d0,    -b,  b, &
                          0d0, -b,   0d0, 0d0,   0d0,  b, &
                            b, -b,     b, 0d0,     b,  b  /),(/2,9/))
       !weigh(1:9) = ??
     END SELECT

     SELECT CASE (nnode)
     CASE (4)
       DO g=1,ngaus
         c = posgp(1,g)
         d = posgp(2,g)
         shape(1,g) = (1d0-c)*(1d0-d)/4d0
         shape(2,g) = (1d0+c)*(1d0-d)/4d0
         shape(3,g) = (1d0+c)*(1d0+d)/4d0
         shape(4,g) = (1d0-c)*(1d0+d)/4d0
         deriv(1,1,g) = (-1d0+d)/4d0
         deriv(2,1,g) = -deriv(1,1,g)
         deriv(3,1,g) = (1d0+d)/4d0
         deriv(4,1,g) = -deriv(3,1,g)
         deriv(1,2,g) = (-1d0+c)/4d0
         deriv(2,2,g) = (-1d0-c)/4d0
         deriv(3,2,g) = -deriv(2,2,g)
         deriv(4,2,g) = -deriv(1,2,g)
       END DO

     CASE (8)
       DO g=1,ngaus
         c = posgp(1,g)
         d = posgp(2,g)
         cd = c*d
         c2 = c*2d0
         d2 = d*2d0
         cc = c*c
         dd = d*d
         ccd = cc*d
         cdd = cd*d
         cd2 = cd*2d0

         !     ***  shape functions

         shape(1,g) = (-1d0+ cd+ cc+ dd- ccd- cdd)/4d0
         shape(2,g) = (-1d0- cd+ cc+ dd- ccd+ cdd)/4d0
         shape(3,g) = (-1d0+ cd+ cc+ dd+ ccd+ cdd)/4d0
         shape(4,g) = (-1d0- cd+ cc+ dd+ ccd- cdd)/4d0
         shape(5,g) = ( 1d0- d- cc+ ccd)/2d0
         shape(6,g) = ( 1d0+ c- dd- cdd)/2d0
         shape(7,g) = ( 1d0+ d- cc- ccd)/2d0
         shape(8,g) = ( 1d0- c- dd+ cdd)/2d0

         !     *** and derivatives

         deriv(1,1,g) = ( d+ c2- cd2- dd)/4d0
         deriv(2,1,g) = (-d+ c2- cd2+ dd)/4d0
         deriv(3,1,g) = ( d+ c2+ cd2+ dd)/4d0
         deriv(4,1,g) = (-d+ c2+ cd2- dd)/4d0
         deriv(5,1,g) = -c+ cd
         deriv(6,1,g) = (1d0- dd)/2d0
         deriv(7,1,g) = -c- cd
         deriv(8,1,g) = (-1d0+ dd)/2d0
         deriv(1,2,g) = ( c+ d2- cc- cd2)/4d0
         deriv(2,2,g) = (-c+ d2- cc+ cd2)/4d0
         deriv(3,2,g) = ( c+ d2+ cc+ cd2)/4d0
         deriv(4,2,g) = (-c+ d2+ cc- cd2)/4d0
         deriv(5,2,g) = (-1d0+ cc)/2d0
         deriv(6,2,g) =  -d- cd
         deriv(7,2,g) = (1d0- cc)/2d0
         deriv(8,2,g) =  -d+ cd
       END DO

     END SELECT

   CASE (3,6)  !triangle
     SELECT CASE (ngaus)
     CASE (1)
       posgp(1:2,1) = 1d0/3d0
       weigh(1) = 0.5d0
     CASE (3)
       c = 1d0/6d0
       weigh(1:3) = c
       IF( etype == 3 )THEN
         d = 2d0/3d0
         posgp(1:2,1:3) = RESHAPE((/ c,d, c,c, d,c  /),(/2,3/))
       ELSE
         d = 1d0/2d0
         posgp(1:2,1:3) = RESHAPE((/ d,0d0, d,d, 0d0,d  /),(/2,3/))
       END IF
     END SELECT

     SELECT CASE (nnode)
     CASE (3)
       DO g=1,ngaus
         shape(1,g) = 1d0 - posgp(1,g) - posgp(2,g)
         shape(2,g) = posgp(1,g)
         shape(3,g) = posgp(2,g)
         deriv(1,1,g) = -1.d0
         deriv(1,2,g) = -1.d0
         deriv(2,1,g) =  1.d0
         deriv(2,2,g) =  0.d0
         deriv(3,1,g) =  0.d0
         deriv(3,2,g) =  1.d0
       END DO
     CASE (6)
       DO g=1,ngaus
         c = posgp(1,g)
         d = posgp(2,g)
         e = 1d0 -c -d
         shape(1,g) = (2d0*e-1d0)*e
         shape(2,g) = (2d0*c-1d0)*c
         shape(3,g) = (2d0*d-1d0)*d
         shape(4,g) = 4d0*e*c
         shape(5,g) = 4d0*c*d
         shape(6,g) = 4d0*e*d

         deriv(1,1,g) =  1d0-4d0*e
         deriv(2,1,g) =  4d0*c-1d0
         deriv(3,1,g) =  0d0
         deriv(4,1,g) =  4d0*(e-c)
         deriv(5,1,g) =  4d0*d
         deriv(6,1,g) = -4d0*d
         deriv(1,2,g) =  1d0-4d0*e
         deriv(2,2,g) =  0d0
         deriv(3,2,g) =  4d0*d-1d0
         deriv(4,2,g) = -4d0*c
         deriv(5,2,g) =  4d0*c
         deriv(6,2,g) =  4d0*(e-d)
       END DO
     END SELECT

   END SELECT

 CASE (3)         !3 dimensional-solids

   SELECT CASE (nnode)
   CASE (4)  !Tetrahedra
     SELECT CASE (ngaus)
     CASE (1)
       posgp(1:3,1) = 1d0/4d0
       weigh(1) = 1d0/6d0
     END SELECT
     DO g=1,ngaus

       s = posgp(1,g)
       t = posgp(2,g)
       q = posgp(3,g)

       shape(   1,g) = 1.-s-t-q
       deriv(1, 1,g) =-1.
       deriv(1, 2,g) =-1.
       deriv(1, 3,g) =-1.

       shape(   2,g) = s
       deriv(2, 1,g) = 1.
       deriv(2, 2,g) = 0.
       deriv(2, 3,g) = 0.

       shape(   3,g) = t
       deriv(3, 1,g) = 0.
       deriv(3, 2,g) = 1.
       deriv(3, 3,g) = 0.

       shape(   4,g) = q
       deriv(4, 1,g) = 0.
       deriv(4, 2,g) = 0.
       deriv(4, 3,g) = 1.
     END DO

   CASE (6)  !Linear Prism

     IF(etype == 27 )THEN  !Solid-Shell with six integration points
        posgp(:,1) = (/ 0.5d0, 0.0d0, -1.0d0 /) ; posgp(:,3) = (/ 0.5d0, 0.5d0, -1.0d0 /) ; posgp(:,5) = (/ 0.0d0, 0.5d0, -1.0d0 /)
        posgp(:,2) = (/ 0.5d0, 0.0d0,  1.0d0 /) ; posgp(:,4) = (/ 0.5d0, 0.5d0,  1.0d0 /) ; posgp(:,6) = (/ 0.0d0, 0.5d0,  1.0d0 /)
        weigh = 1d0/6d0
     ELSE
       posgp(1:2,:) = 1d0/3d0  ! one integration point in the central axis
       IF(etype == 5 .OR. etype == 12 )THEN  !Solid-Shell
         posgp(3,:) = (/ -1d0, 1d0 /)  ! data is on external surfaces only
         weigh      = 1d0
       ELSE
         CALL gaussq( ngaus,posgp(3,:),weigh)
       END IF
     END IF


     DO g=1,ngaus
       s = posgp(1,g)
       t = posgp(2,g)
       q = posgp(3,g)
       sp = (1d0-q)/2d0   !N1
       tp = (1d0+q)/2d0   !N2
       zp = 1d0-s-t

       !*** shape functions

       shape(1,g) = zp * sp
       shape(2,g) = s  * sp
       shape(3,g) = t  * sp
       shape(4,g) = zp * tp
       shape(5,g) = s  * tp
       shape(6,g) = t  * tp

       !*** and derivatives

       deriv(1,1,g) = -sp
       deriv(2,1,g) =  sp
       deriv(3,1,g) =  0d0
       deriv(4,1,g) = -tp
       deriv(5,1,g) =  tp
       deriv(6,1,g) =  0d0

       deriv(1,2,g) = -sp
       deriv(2,2,g) =  0d0
       deriv(3,2,g) =  sp
       deriv(4,2,g) = -tp
       deriv(5,2,g) =  0d0
       deriv(6,2,g) =  tp

       deriv(1,3,g) = -zp/2d0
       deriv(2,3,g) = -s /2d0
       deriv(3,3,g) = -t /2d0
       deriv(4,3,g) =  zp/2d0
       deriv(5,3,g) =  s /2d0
       deriv(6,3,g) =  t /2d0

     END DO

   !CASE (12)  !Bezier Prism
   !  posgp(:,1) = (/ 0.5d0, 0.0d0, -1.0d0 /) ; posgp(:,2) = (/ 0.5d0, 0.5d0, -1.0d0 /) ; posgp(:,3) = (/ 0.0d0, 0.5d0, -1.0d0 /)
   !  posgp(:,4) = (/ 0.5d0, 0.0d0,  1.0d0 /) ; posgp(:,5) = (/ 0.5d0, 0.5d0,  1.0d0 /) ; posgp(:,6) = (/ 0.0d0, 0.5d0,  1.0d0 /)
   !  weigh = 1d0/6d0
   !  DO g=1,ngaus
   !    xita = posgp(1,g)
   !    eta  = posgp(1,g)
   !    zeta = posgp(1,g)
   !    zzet = 1d0-xita-eta  !third area coordinate
   !    l1 = (1d0-zeta)/2d0  !linear functions L1
   !    l2 = (1d0+zeta)/2d0  !                 L2
   !    !*** shape functions
   !
   !    shape(1) = zzet * zzet
   !    shape(2) = xita * xita
   !    shape(3) = eta  * eta
   !    shape(4) = 2d0* xita * zzet
   !    shape(5) = 2d0* eta  * xita
   !    shape(6) = 2d0* zzet * eta
   !
   !    !*** and derivatives
   !
   !    deriv(1,1) = -2d0*zzet
   !    deriv(2,1) =  2d0*xita
   !    deriv(3,1) =  0d0
   !    deriv(4,1) =  2d0*(zzet-xita)
   !    deriv(5,1) =  2d0*eta
   !    deriv(6,1) = -2d0*eta
   !
   !    deriv(1,2) = -2d0*zzet
   !    deriv(2,2) =  0d0
   !    deriv(3,2) =  2d0*eta
   !    deriv(4,2) = -2d0*xita
   !    deriv(5,2) =  2d0*xita
   !    deriv(6,2) =  2d0*(zzet-eta)
   !
   !    deriv(1: 6,3) = -shape(1:6)/2d0
   !    deriv(7:12,3) =  shape(1:6)/2d0
   !
   !    shape(7:12) = shape(1:6) * l2
   !    shape(1: 6) = shape(1:6) * l1
   !    deriv(7:12,1:2) = deriv(1:6,1:2)*l2
   !    deriv(1: 6,1:2) = deriv(1:6,1:2)*l1
   !  END DO

   CASE (15)  !Quadratic Prism
     c=0.333333333333333d0
     d=0.666666666666667d0
     e=0.166666666666667d0

     SELECT CASE (ngaus)
     CASE ( 1 )
       posgp(1:3,1:1) = RESHAPE((/  c, c, 0d0 /), (/3,1/))
       weigh = 1.0d0
     CASE ( 2 )
       posgp(1:3,1:2) = RESHAPE((/  c, c, -a, c, c, a /), (/3,2/))
       weigh = 0.5d0
     CASE ( 6 )
       posgp(1:3,1:6) = RESHAPE((/  e,e,-a, d,e,-a, e,d,-a, &
                                    e,e, a, d,e, a, e,d, a  /), (/3,6/))
       weigh = e
     END SELECT

     DO g=1,ngaus
       s = posgp(1,g)
       t = posgp(2,g)
       q = posgp(3,g)
       sp = (1d0-q)/2d0   !N1
       tp = (1d0+q)/2d0   !N2
       zp = 1d0-s-t


       bf = 1d0 - q**2                      !bubble function
       shape(1,g) = zp*(2d0*zp-1d0)*sp - zp*bf/2d0
       shape(2,g) = s*(2d0*s-1d0)*sp - s*bf/2d0
       shape(3,g) = t *(2d0*t -1d0)*sp - t *bf/2d0
       shape(4,g) = zp*(2d0*zp-1d0)*tp - zp*bf/2d0
       shape(5,g) = s*(2d0*s-1d0)*tp - s*bf/2d0
       shape(6,g) = t *(2d0*t -1d0)*tp - t *bf/2d0
       shape(7,g) = 4d0*zp*s*sp
       shape(8,g) = 4d0*s*t *sp
       shape(9,g) = 4d0*t *zp*sp
       shape(10,g) = zp*bf
       shape(11,g) = s*bf
       shape(12,g) = t *bf
       shape(13,g) = 4d0*zp*s*tp
       shape(14,g) = 4d0*s*t *tp
       shape(15,g) = 4d0*t *zp*tp

       deriv(1,1,g) = (-4d0*zp+1d0)*sp + bf/2d0
       deriv(2,1,g) =  (4d0*s-1d0)*sp - bf/2d0
       deriv(3,1,g) =  0d0
       deriv(4,1,g) = (-4d0*zp+1d0)*tp + bf/2d0
       deriv(5,1,g) =  (4d0*s-1d0)*tp - bf/2d0
       deriv(6,1,g) =  0d0
       deriv(7,1,g) =  4d0*(zp-s)*sp
       deriv(8,1,g) =  4d0*t *       sp
       deriv(9,1,g) = -4d0*t *       sp
       deriv(10,1,g) = -bf
       deriv(11,1,g) =  bf
       deriv(12,1,g) =  0d0
       deriv(13,1,g) =  4d0*(zp-s)*tp
       deriv(14,1,g) =  4d0*t *       tp
       deriv(15,1,g) = -4d0*t *       tp

       deriv(1,2,g) = (-4d0*zp+1d0)*sp + bf/2d0
       deriv(2,2,g) =  0d0
       deriv(3,2,g) =  (4d0*t -1d0)*sp - bf/2d0
       deriv(4,2,g) = (-4d0*zp+1d0)*tp + bf/2d0
       deriv(5,2,g) =  0d0
       deriv(6,2,g) =  (4d0*t -1d0)*tp - bf/2d0
       deriv(7,2,g) = -4d0*s*       sp
       deriv(8,2,g) =  4d0*s*       sp
       deriv(9,2,g) =  4d0*(-t+zp)*sp
       deriv(10,2,g) = -bf
       deriv(11,2,g) =  0d0
       deriv(12,2,g) =  bf
       deriv(13,2,g) = -4d0*s*      tp
       deriv(14,2,g) =  4d0*s*      tp
       deriv(15,2,g) = 4d0*(-t+zp)*tp

       deriv(1,3,g) = -zp*(2d0*zp-1d0)/2d0 + zp*q
       deriv(2,3,g) = -s*(2d0*s-1d0)/2d0 + s*q
       deriv(3,3,g) = -t *(2d0*t -1d0)/2d0 + t *q
       deriv(4,3,g) =  zp*(2d0*zp-1d0)/2d0 + zp*q
       deriv(5,3,g) =  s*(2d0*s-1d0)/2d0 + s*q
       deriv(6,3,g) =  t *(2d0*t -1d0)/2d0 + t *q
       deriv(7,3,g) = -2d0*zp*s
       deriv(8,3,g) = -2d0*s*t
       deriv(9,3,g) = -2d0*t *zp
       deriv(10,3,g) = -zp*2d0*q
       deriv(11,3,g) = -s*2d0*q
       deriv(12,3,g) = -t *2d0*q
       deriv(13,3,g) = 2d0*zp*s
       deriv(14,3,g) = 2d0*s*t
       deriv(15,3,g) = 2d0*t *zp
     END DO

   CASE (18)  !Quadratic Prism
     c=0.333333333333333d0
     d=0.666666666666667d0
     e=0.166666666666667d0

     SELECT CASE (ngaus)
     CASE ( 1 )
       posgp(1:3,1:1) = RESHAPE((/  c, c, 0d0 /), (/3,1/))
       weigh = 1.0d0
     CASE ( 2 )
       posgp(1:3,1:2) = RESHAPE((/  c, c, -a, c, c, a /), (/3,2/))
       weigh = 0.5d0
     CASE ( 6 )
       posgp(1:3,1:6) = RESHAPE((/  e,e,-a, d,e,-a, e,d,-a, &
                                    e,e, a, d,e, a, e,d, a  /), (/3,6/))
       weigh = e
     CASE ( 8 )
       posgp(1:3,1:6) = RESHAPE((/  e,e,-a, d,e,-a, e,d,-a, &
                                    e,e, a, d,e, a, e,d, a  /), (/3,6/))
       weigh = e
     CASE ( 9 )
       posgp(1:3,1:9) = RESHAPE((/  e,e,-b, d,e,-b, e,d,-b, &
                                    e,e, 0, d,e, 0, e,d, 0, &
                                    e,e, b, d,e, b, e,d, b  /), (/3,9/))
       weigh = (/ 5d0/9d0, 5d0/9d0, 5d0/9d0, 8d0/9d0, 8d0/9d0, 8d0/9d0, 5d0/9d0, 5d0/9d0,  5d0/9d0  /)/2d0
     END SELECT

     DO g=1,ngaus
       s = posgp(1,g)
       t = posgp(2,g)
       q = posgp(3,g)
       sp = -q*(1d0-q)/2d0  !quadratic functions N1
       tp =  q*(1d0+q)/2d0  !quadratic functions N3
       zp = 1d0-s-t

       bf = 1d0 - q**2                      !bubble function

       shape( 1,g) = zp*(2d0*zp-1d0) *sp
       shape( 2,g) =  s*(2d0*s-1d0)  *sp
       shape( 3,g) =  t*(2d0*t -1d0) *sp
       shape( 7,g) = 4d0*zp*s        *sp
       shape( 8,g) = 4d0*s *t        *sp
       shape( 9,g) = 4d0*t *zp       *sp
       shape(10,g) = zp*(2d0*zp-1d0) *bf
       shape(11,g) = s*(2d0*s-1d0)   *bf
       shape(12,g) = t *(2d0*t -1d0) *bf
       shape(16,g) = 4d0*zp*s        *bf
       shape(17,g) = 4d0*s *t        *bf
       shape(18,g) = 4d0*t *zp       *bf
       shape( 4,g) = zp*(2d0*zp-1d0) *tp
       shape( 5,g) = s *(2d0*s -1d0) *tp
       shape( 6,g) = t *(2d0*t -1d0) *tp
       shape(13,g) = 4d0*zp*s        *tp
       shape(14,g) = 4d0*s *t        *tp
       shape(15,g) = 4d0*t *zp       *tp

       deriv( 1,1,g) =-(4d0*zp-1d0)  *sp
       deriv( 2,1,g) = (4d0*s-1d0)   *sp
       deriv( 3,1,g) =  0d0
       deriv( 7,1,g) =  4d0*(zp-s)   *sp
       deriv( 8,1,g) =  4d0*t        *sp
       deriv( 9,1,g) = -4d0*t        *sp
       deriv(10,1,g) =-(4d0*zp-1d0)  *bf
       deriv(11,1,g) = (4d0*s -1d0)  *bf
       deriv(12,1,g) =  0d0
       deriv(16,1,g) =  4d0*(zp-s )  *bf
       deriv(17,1,g) =  4d0*t        *bf
       deriv(18,1,g) = -4d0*t        *bf
       deriv( 4,1,g) =-(4d0*zp-1d0)  *tp
       deriv( 5,1,g) = (4d0*s -1d0)  *tp
       deriv( 6,1,g) =  0d0
       deriv(13,1,g) =  4d0*(zp-s )  *tp
       deriv(14,1,g) =  4d0*t        *tp
       deriv(15,1,g) = -4d0*t        *tp

       deriv( 1,2,g) =-(4d0*zp-1d0)  *sp
       deriv( 2,2,g) =  0d0
       deriv( 3,2,g) = (4d0*t -1d0)  *sp
       deriv( 7,2,g) = -4d0*s        *sp
       deriv( 8,2,g) =  4d0*s        *sp
       deriv( 9,2,g) =  4d0*(zp-t )  *sp
       deriv(10,2,g) =-(4d0*zp-1d0)  *bf
       deriv(11,2,g) =  0d0
       deriv(12,2,g) = (4d0*t -1d0)  *bf
       deriv(16,2,g) = -4d0*s        *bf
       deriv(17,2,g) =  4d0*s        *bf
       deriv(18,2,g) =  4d0*(zp-t )  *bf
       deriv( 4,2,g) =-(4d0*zp-1d0)  *tp
       deriv( 5,2,g) =  0d0
       deriv( 6,2,g) = (4d0*t -1d0)  *tp
       deriv(13,2,g) = -4d0*s        *tp
       deriv(14,2,g) =  4d0*s        *tp
       deriv(15,2,g) =  4d0*(zp-t )  *tp

       deriv( 1,3,g) =  zp*(2d0*zp-1d0) *(q-0.5d0)
       deriv( 2,3,g) =  s *(2d0*s -1d0) *(q-0.5d0)
       deriv( 3,3,g) =  t *(2d0*t -1d0) *(q-0.5d0)
       deriv( 7,3,g) =  4d0*zp*s        *(q-0.5d0)
       deriv( 8,3,g) =  4d0*s *t        *(q-0.5d0)
       deriv( 9,3,g) =  4d0*t *zp       *(q-0.5d0)
       deriv(10,3,g) =  zp*(2d0*zp-1d0) *(-2d0*q)
       deriv(11,3,g) =  s *(2d0*s -1d0) *(-2d0*q)
       deriv(12,3,g) =  t *(2d0*t -1d0) *(-2d0*q)
       deriv(16,3,g) =  4d0*zp*s        *(-2d0*q)
       deriv(17,3,g) =  4d0*s *t        *(-2d0*q)
       deriv(18,3,g) =  4d0*t *zp       *(-2d0*q)
       deriv( 4,3,g) =  zp*(2d0*zp-1d0) *(q+0.5d0)
       deriv( 5,3,g) =  s *(2d0*s -1d0) *(q+0.5d0)
       deriv( 6,3,g) =  t *(2d0*t -1d0) *(q+0.5d0)
       deriv(13,3,g) =  4d0*zp*s        *(q+0.5d0)
       deriv(14,3,g) =  4d0*s *t        *(q+0.5d0)
       deriv(15,3,g) =  4d0*t *zp       *(q+0.5d0)
     END DO

   CASE (8,20)  !Hexahedra
     SELECT CASE (ngaus)
     CASE (1)         !reduced integration
       posgp(1:3,1) = 0d0
       weigh(1) = 8d0
     CASE (2)         !solid-shell
       posgp(1:3,1) = (/ 0d0,0d0,-1d0/)
       posgp(1:3,2) = (/ 0d0,0d0,+1d0/)
       weigh(1:2) = 4d0
     CASE (4)         !spetial rule
       posgp(1:3,1:4) = RESHAPE((/  a,-a,-a, -a, a,-a, -a,-a, a ,  a, a, a /), &
       !posgp(1:3,1:4) = RESHAPE((/  a, a, a, -a,-a, a, -a, a,-a ,  a,-a,-a /), &
                                   (/3,4/))
       weigh(1:4) = 2d0
     CASE (8)         !standard rule
       posgp(1:3,1:8) = RESHAPE((/ -a,-a,-a,  a,-a,-a, -a, a,-a ,  a, a,-a,    &
                                   -a,-a, a,  a,-a, a, -a, a, a ,  a, a, a /), &
                                   (/3,8/))
       !posgp(1:3,1:8) = RESHAPE((/ -a,-a,-a, -a,-a, a, -a, a,-a , -a, a, a,    &
       !                             a,-a,-a,  a,-a, a,  a, a,-a ,  a, a, a /), &
       !                            (/3,8/))
       weigh(1:8) = 1d0
     END SELECT

     DO g=1,ngaus
       s = posgp(1,g)
       t = posgp(2,g)
       q = posgp(3,g)

       IF( nnode == 8 )THEN
         !       auxiliar values
         sm = 0.5d0*(1d0-s)
         tm = 0.5d0*(1d0-t)
         zm = 0.5d0*(1d0-q)
         sp = 0.5d0*(1d0+s)
         tp = 0.5d0*(1d0+t)
         zp = 0.5d0*(1d0+q)
         !   shape functions

         shape(1,g) = sm*tm*zm
         shape(2,g) = sp*tm*zm
         shape(3,g) = sp*tp*zm
         shape(4,g) = sm*tp*zm
         shape(5,g) = sm*tm*zp
         shape(6,g) = sp*tm*zp
         shape(7,g) = sp*tp*zp
         shape(8,g) = sm*tp*zp

         !     derivatives

         deriv(1,1,g) = -0.5d0*tm*zm
         deriv(2,1,g) = +0.5d0*tm*zm
         deriv(3,1,g) = +0.5d0*tp*zm
         deriv(4,1,g) = -0.5d0*tp*zm
         deriv(5,1,g) = -0.5d0*tm*zp
         deriv(6,1,g) = +0.5d0*tm*zp
         deriv(7,1,g) = +0.5d0*tp*zp
         deriv(8,1,g) = -0.5d0*tp*zp
         deriv(1,2,g) = -0.5d0*sm*zm
         deriv(2,2,g) = -0.5d0*sp*zm
         deriv(3,2,g) = +0.5d0*sp*zm
         deriv(4,2,g) = +0.5d0*sm*zm
         deriv(5,2,g) = -0.5d0*sm*zp
         deriv(6,2,g) = -0.5d0*sp*zp
         deriv(7,2,g) = +0.5d0*sp*zp
         deriv(8,2,g) = +0.5d0*sm*zp
         deriv(1,3,g) = -0.5d0*sm*tm
         deriv(2,3,g) = -0.5d0*sp*tm
         deriv(3,3,g) = -0.5d0*sp*tp
         deriv(4,3,g) = -0.5d0*sm*tp
         deriv(5,3,g) = +0.5d0*sm*tm
         deriv(6,3,g) = +0.5d0*sp*tm
         deriv(7,3,g) = +0.5d0*sp*tp
         deriv(8,3,g) = +0.5d0*sm*tp

       ELSE IF (nnode == 20 )THEN
         !       auxiliar values

         omg=1.d0-s
         omh=1.d0-t
         omr=1.d0-q
         opg=1.d0+s
         oph=1.d0+t
         opr=1.d0+q
         tpgphpr=2.d0+s+t+q
         tmgphpr=2.d0-s+t+q
         tmgmhpr=2.d0-s-t+q
         tpgmhpr=2.d0+s-t+q
         tpgphmr=2.d0+s+t-q
         tmgphmr=2.d0-s+t-q
         tmgmhmr=2.d0-s-t-q
         tpgmhmr=2.d0+s-t-q

         !       shape functions

         shape( 1,g)=-omg*omh*omr*tpgphpr/8.d0
         shape( 2,g)=-opg*omh*omr*tmgphpr/8.d0
         shape( 3,g)=-opg*oph*omr*tmgmhpr/8.d0
         shape( 4,g)=-omg*oph*omr*tpgmhpr/8.d0
         shape( 5,g)=-omg*omh*opr*tpgphmr/8.d0
         shape( 6,g)=-opg*omh*opr*tmgphmr/8.d0
         shape( 7,g)=-opg*oph*opr*tmgmhmr/8.d0
         shape( 8,g)=-omg*oph*opr*tpgmhmr/8.d0
         shape( 9,g)=omg*opg*omh*omr/4.d0
         shape(10,g)=omh*oph*opg*omr/4.d0
         shape(11,g)=omg*opg*oph*omr/4.d0
         shape(12,g)=omh*oph*omg*omr/4.d0
         shape(13,g)=omr*opr*omg*omh/4.d0
         shape(14,g)=omr*opr*opg*omh/4.d0
         shape(15,g)=omr*opr*opg*oph/4.d0
         shape(16,g)=omr*opr*omg*oph/4.d0
         shape(17,g)=omg*opg*omh*opr/4.d0
         shape(18,g)=omh*oph*opg*opr/4.d0
         shape(19,g)=omg*opg*oph*opr/4.d0
         shape(20,g)=omh*oph*omg*opr/4.d0

         !       local derivatives of the shape functions: xi-derivative

         deriv( 1,1,g)=omh*omr*(tpgphpr-omg)/8.d0
         deriv( 2,1,g)=(opg-tmgphpr)*omh*omr/8.d0
         deriv( 3,1,g)=(opg-tmgmhpr)*oph*omr/8.d0
         deriv( 4,1,g)=oph*omr*(tpgmhpr-omg)/8.d0
         deriv( 5,1,g)=omh*opr*(tpgphmr-omg)/8.d0
         deriv( 6,1,g)=(opg-tmgphmr)*omh*opr/8.d0
         deriv( 7,1,g)=(opg-tmgmhmr)*oph*opr/8.d0
         deriv( 8,1,g)=oph*opr*(tpgmhmr-omg)/8.d0
         deriv( 9,1,g)=(omg-opg)*omh*omr/4.d0
         deriv(10,1,g)=omh*oph*omr/4.d0
         deriv(11,1,g)=(omg-opg)*oph*omr/4.d0
         deriv(12,1,g)=-omh*oph*omr/4.d0
         deriv(13,1,g)=-omr*opr*omh/4.d0
         deriv(14,1,g)=omr*opr*omh/4.d0
         deriv(15,1,g)=omr*opr*oph/4.d0
         deriv(16,1,g)=-omr*opr*oph/4.d0
         deriv(17,1,g)=(omg-opg)*omh*opr/4.d0
         deriv(18,1,g)=omh*oph*opr/4.d0
         deriv(19,1,g)=(omg-opg)*oph*opr/4.d0
         deriv(20,1,g)=-omh*oph*opr/4.d0

         !       local derivatives of the shape functions: eta-derivative

         deriv( 1,2,g)=omg*omr*(tpgphpr-omh)/8.d0
         deriv( 2,2,g)=opg*omr*(tmgphpr-omh)/8.d0
         deriv( 3,2,g)=opg*(oph-tmgmhpr)*omr/8.d0
         deriv( 4,2,g)=omg*(oph-tpgmhpr)*omr/8.d0
         deriv( 5,2,g)=omg*opr*(tpgphmr-omh)/8.d0
         deriv( 6,2,g)=opg*opr*(tmgphmr-omh)/8.d0
         deriv( 7,2,g)=opg*(oph-tmgmhmr)*opr/8.d0
         deriv( 8,2,g)=omg*(oph-tpgmhmr)*opr/8.d0
         deriv( 9,2,g)=-omg*opg*omr/4.d0
         deriv(10,2,g)=(omh-oph)*opg*omr/4.d0
         deriv(11,2,g)=omg*opg*omr/4.d0
         deriv(12,2,g)=(omh-oph)*omg*omr/4.d0
         deriv(13,2,g)=-omr*opr*omg/4.d0
         deriv(14,2,g)=-omr*opr*opg/4.d0
         deriv(15,2,g)=omr*opr*opg/4.d0
         deriv(16,2,g)=omr*opr*omg/4.d0
         deriv(17,2,g)=-omg*opg*opr/4.d0
         deriv(18,2,g)=(omh-oph)*opg*opr/4.d0
         deriv(19,2,g)=omg*opg*opr/4.d0
         deriv(20,2,g)=(omh-oph)*omg*opr/4.d0

         !       local derivatives of the shape functions: zeta-derivative

         deriv( 1,3,g)=omg*omh*(tpgphpr-omr)/8.d0
         deriv( 2,3,g)=opg*omh*(tmgphpr-omr)/8.d0
         deriv( 3,3,g)=opg*oph*(tmgmhpr-omr)/8.d0
         deriv( 4,3,g)=omg*oph*(tpgmhpr-omr)/8.d0
         deriv( 5,3,g)=omg*omh*(opr-tpgphmr)/8.d0
         deriv( 6,3,g)=opg*omh*(opr-tmgphmr)/8.d0
         deriv( 7,3,g)=opg*oph*(opr-tmgmhmr)/8.d0
         deriv( 8,3,g)=omg*oph*(opr-tpgmhmr)/8.d0
         deriv( 9,3,g)=-omg*opg*omh/4.d0
         deriv(10,3,g)=-omh*oph*opg/4.d0
         deriv(11,3,g)=-omg*opg*oph/4.d0
         deriv(12,3,g)=-omh*oph*omg/4.d0
         deriv(13,3,g)=(omr-opr)*omg*omh/4.d0
         deriv(14,3,g)=(omr-opr)*opg*omh/4.d0
         deriv(15,3,g)=(omr-opr)*opg*oph/4.d0
         deriv(16,3,g)=(omr-opr)*omg*oph/4.d0
         deriv(17,3,g)=omg*opg*omh/4.d0
         deriv(18,3,g)=omh*oph*opg/4.d0
         deriv(19,3,g)=omg*opg*oph/4.d0
         deriv(20,3,g)=omh*oph*omg/4.d0

       END IF
     END DO

   END SELECT

 END SELECT

 END SUBROUTINE shapes
