 SUBROUTINE shape4 (nnode, shape, deriv ,xita, eta, zeta, bezier, order )
 !********************************************************************
 !
 !***  calculates shape functions and their derivatives for 3d
 !     prismatic elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER(Kind=4) :: nnode
 REAL (kind=8)    deriv(nnode,3),shape(nnode),xita,eta,zeta
 LOGICAL, INTENT(IN), OPTIONAL :: bezier, order

 REAL(kind=8), PARAMETER :: fq =0.25d0, fl = 0.75d0
 REAL (kind=8) l1,l2,zzet,bf,l12,l22,f1,f2,fb
 LOGICAL :: bez,ord

   bez =  PRESENT(bezier)
   IF( bez )  bez = bezier
   ord =  PRESENT(order)
   IF( ord ) ord = order

   l1 = (1d0-zeta)/2d0  !linear functions L1
   l2 = (1d0+zeta)/2d0  !                 L2
   zzet = 1d0-xita-eta  !third area coordinate

   SELECT CASE (nnode)
   CASE (6)
     !*** shape functions

     shape(1) = zzet * l1
     shape(2) = xita * l1
     shape(3) = eta  * l1
     shape(4) = zzet * l2
     shape(5) = xita * l2
     shape(6) = eta  * l2

     !*** and derivatives

     deriv(1,1) = -l1
     deriv(2,1) =  l1
     deriv(3,1) =  0d0
     deriv(4,1) = -l2
     deriv(5,1) =  l2
     deriv(6,1) =  0d0

     deriv(1,2) = -l1
     deriv(2,2) =  0d0
     deriv(3,2) =  l1
     deriv(4,2) = -l2
     deriv(5,2) =  0d0
     deriv(6,2) =  l2

     deriv(1,3) = -zzet/2d0
     deriv(2,3) = -xita/2d0
     deriv(3,3) = -eta /2d0
     deriv(4,3) =  zzet/2d0
     deriv(5,3) =  xita/2d0
     deriv(6,3) =  eta /2d0

   CASE (12)

     IF(bez) THEN  !Bezier Prism
       !*** shape functions   (middle surface)
       shape(1) = zzet * zzet
       shape(2) = xita * xita
       shape(3) = eta  * eta
       shape(4) = 2d0* xita * zzet
       shape(5) = 2d0* eta  * xita
       shape(6) = 2d0* zzet * eta
       !*** and derivatives   (in plane)
       deriv(1,1) = -2d0*zzet
       deriv(2,1) =  2d0*xita
       deriv(3,1) =  0d0
       deriv(4,1) =  2d0*(zzet-xita)
       deriv(5,1) =  2d0*eta
       deriv(6,1) = -2d0*eta

       deriv(1,2) = -2d0*zzet
       deriv(2,2) =  0d0
       deriv(3,2) =  2d0*eta
       deriv(4,2) = -2d0*xita
       deriv(5,2) =  2d0*xita
       deriv(6,2) =  2d0*(zzet-eta)
     ELSE
       !*** shape functions   (middle surface)
       shape(1) = (2d0*xita-1d0)*xita
       shape(2) = (2d0*eta-1d0)*eta
       shape(3) = (2d0*zzet-1d0)*zzet
       shape(4) = 4d0*xita*eta
       shape(5) = 4d0*eta*zzet
       shape(6) = 4d0*xita*zzet
       !*** and derivatives  (in plane)
       deriv(1,1) =  1d0-4d0*xita
       deriv(2,1) =  4d0*eta-1d0
       deriv(3,1) =  0d0
       deriv(4,1) =  4d0*(xita-eta)
       deriv(5,1) =  4d0*zzet
       deriv(6,1) = -4d0*zzet

       deriv(1,2) =  1d0-4d0*xita
       deriv(2,2) =  0d0
       deriv(3,2) =  4d0*zzet-1d0
       deriv(4,2) = -4d0*eta
       deriv(5,2) =  4d0*eta
       deriv(6,2) =  4d0*(xita-zzet)
     END IF
     ! derivatives (out of plane)
     deriv(1: 6,3) = -shape(1:6)/2d0
     deriv(7:12,3) =  shape(1:6)/2d0
     ! Shape functions
     shape(7:12) = shape(1:6) * l2  !top surface
     shape(1: 6) = shape(1:6) * l1  !bottom surface
     ! derivatives (in plane)
     deriv(7:12,1:2) = deriv(1:6,1:2)*l2 !top surface
     deriv(1: 6,1:2) = deriv(1:6,1:2)*l1 !bottom surface

     IF( ord ) THEN
       shape(:)   = shape((/1:3,7:9,4:6,10:12/))
       deriv(:,1) = deriv((/1:3,7:9,4:6,10:12/),1)
       deriv(:,2) = deriv((/1:3,7:9,4:6,10:12/),2)
       deriv(:,3) = deriv((/1:3,7:9,4:6,10:12/),3)
     END IF

   CASE (15)
     IF( bez )THEN
       ! bezier quadratic polynomials
       l12 = l1*l1                                            !l12 = l1*l1
       l22 = l2*l2                                            !l22 = l2*l2
       bf  = 1d0 - zeta*zeta    !twice the polynomial         !bf  = (1d0 - zeta*zeta)/2d0
                                                              !f1  = fq*l12+fl*l1
                                                              !f2  = fq*l22+fl*l2
                                                              !fb  = fq*bf
                                                              !
       shape( 1) = zzet * zzet      * l1 - zzet * bf / 8d0    !shape( 1) = zzet * zzet       *f1
       shape( 2) = xita * xita      * l1 - xita * bf / 8d0    !shape( 2) = xita * xita       *f1
       shape( 3) = eta  * eta       * l1 - eta  * bf / 8d0    !shape( 3) = eta  * eta        *f1
       shape( 7) = 2d0* xita * zzet * l1                      !shape( 7) = 2d0* xita * zzet  *f1
       shape( 8) = 2d0* eta  * xita * l1                      !shape( 8) = 2d0* eta  * xita  *f1
       shape( 9) = 2d0* zzet * eta  * l1                      !shape( 9) = 2d0* zzet * eta   *f1
       shape( 4) = zzet * zzet      * l2 - zzet * bf / 8d0    !shape( 4) = zzet * zzet       *f2
       shape( 5) = xita * xita      * l2 - xita * bf / 8d0    !shape( 5) = xita * xita       *f2
       shape( 6) = eta  * eta       * l2 - eta  * bf / 8d0    !shape( 6) = eta  * eta        *f2
       shape(13) = 2d0* xita * zzet * l2                      !shape(13) = 2d0* xita * zzet  *f2
       shape(14) = 2d0* eta  * xita * l2                      !shape(14) = 2d0* eta  * xita  *f2
       shape(15) = 2d0* zzet * eta  * l2                      !shape(15) = 2d0* zzet * eta   *f2
       shape(10) =                         zzet * bf / 4d0    !shape(10) = zzet              *fb
       shape(11) =                         xita * bf / 4d0    !shape(11) = xita              *fb
       shape(12) =                         eta  * bf / 4d0    !shape(12) = eta               *fb
                                                              !
       deriv( 1,1) = -2d0*zzet        * l1 + bf/8d0           !deriv( 1,1) = -2d0*zzet       *f1
       deriv( 2,1) =  2d0*xita        * l1 - bf/8d0           !deriv( 2,1) =  2d0*xita       *f1
       deriv( 3,1) =  0d0                                     !deriv( 3,1) =  0d0
       deriv( 7,1) =  2d0*(zzet-xita) * l1                    !deriv( 7,1) =  2d0*(zzet-xita)*f1
       deriv( 8,1) =  2d0*eta         * l1                    !deriv( 8,1) =  2d0*eta        *f1
       deriv( 9,1) = -2d0*eta         * l1                    !deriv( 9,1) = -2d0*eta        *f1
       deriv( 4,1) = -2d0*zzet        * l2 + bf/8d0           !deriv( 4,1) = -2d0*zzet       *f2
       deriv( 5,1) =  2d0*xita        * l2 - bf/8d0           !deriv( 5,1) =  2d0*xita       *f2
       deriv( 6,1) =  0d0                                     !deriv( 6,1) =  0d0
       deriv(13,1) =  2d0*(zzet-xita) * l2                    !deriv(13,1) =  2d0*(zzet-xita)*f2
       deriv(14,1) =  2d0*eta         * l2                    !deriv(14,1) =  2d0*eta        *f2
       deriv(15,1) = -2d0*eta         * l2                    !deriv(15,1) = -2d0*eta        *f2
       deriv(10,1) =                        -bf/4d0           !deriv(10,1) = -                fb
       deriv(11,1) =                         bf/4d0           !deriv(11,1) =                  fb
       deriv(12,1) =                         0d0              !deriv(12,1) =  0d0
                                                              !
       deriv( 1,2) = -2d0*zzet        * l1 + bf/8d0                    !deriv( 1,2) = -2d0*zzet       *f1
       deriv( 2,2) =  0d0                                     !deriv( 2,2) =  0d0
       deriv( 3,2) =  2d0*eta         * l1 - bf/8d0                    !deriv( 3,2) =  2d0*eta        *f1
       deriv( 7,2) = -2d0*xita        * l1                   !deriv( 7,2) = -2d0*xita       *f1
       deriv( 8,2) =  2d0*xita        * l1                   !deriv( 8,2) =  2d0*xita       *f1
       deriv( 9,2) =  2d0*(zzet-eta)  * l1                   !deriv( 9,2) =  2d0*(zzet-eta) *f1
       deriv( 4,2) = -2d0*zzet        * l2 + bf/8d0                    !deriv( 4,2) = -2d0*zzet       *f2
       deriv( 5,2) =  0d0                                     !deriv( 5,2) =  0d0
       deriv( 6,2) =  2d0*eta         * l2  - bf/8d0                    !deriv( 6,2) =  2d0*eta        *f2
       deriv(13,2) = -2d0*xita        * l2                   !deriv(13,2) = -2d0*xita       *f2
       deriv(14,2) =  2d0*xita        * l2                   !deriv(14,2) =  2d0*xita       *f2
       deriv(15,2) =  2d0*(zzet-eta)  * l2                   !deriv(15,2) =  2d0*(zzet-eta) *f2
       deriv(10,2) =                        - bf/4d0         !deriv(10,2) = -                fb
       deriv(11,2) =                          0d0            !deriv(11,2) =  0d0
       deriv(12,2) =                        + bf/4d0         !deriv(12,2) =                  fb
                                                              !
                                                              !f1  =-fq*l1 - fl/2d0
                                                              !f2  = fq*l2 + fl/2d0
                                                              !fb  =-fq*zeta
       deriv( 1,3) = -zzet*zzet/2d0 + zzet*zeta/4d0           !deriv( 1,3) = zzet * zzet      *f1
       deriv( 2,3) = -xita*xita/2d0 + xita*zeta/4d0           !deriv( 2,3) = xita * xita      *f1
       deriv( 3,3) = -eta * eta/2d0 + eta *zeta/4d0           !deriv( 3,3) = eta  * eta       *f1
       deriv( 7,3) = -zzet*xita                               !deriv( 7,3) = 2d0* xita * zzet *f1
       deriv( 8,3) = -xita*eta                                !deriv( 8,3) = 2d0* eta  * xita *f1
       deriv( 9,3) = -eta *zzet                               !deriv( 9,3) = 2d0* zzet * eta  *f1
       deriv( 4,3) =  zzet*zzet/2d0 + zzet*zeta/4d0           !deriv( 4,3) = zzet * zzet      *f2
       deriv( 5,3) =  xita*xita/2d0 + xita*zeta/4d0           !deriv( 5,3) = xita * xita      *f2
       deriv( 6,3) =  eta * eta/2d0 + eta *zeta/4d0           !deriv( 6,3) = eta  * eta       *f2
       deriv(13,3) =  zzet*xita                               !deriv(13,3) = 2d0* xita * zzet *f2
       deriv(14,3) =  xita*eta                                !deriv(14,3) = 2d0* eta  * xita *f2
       deriv(15,3) =  eta *zzet                               !deriv(15,3) = 2d0* zzet * eta  *f2
       deriv(10,3) = -zzet*zeta /2d0                          !deriv(10,3) = zzet             *fb
       deriv(11,3) = -xita*zeta /2d0                          !deriv(11,3) = xita             *fb
       deriv(12,3) = -eta *zeta /2d0                          !deriv(12,3) = eta              *fb

       IF( ord ) THEN
         shape(:)   = shape((/1:3,7:9,4:6,13:15,10:12/))
         deriv(:,1) = deriv((/1:3,7:9,4:6,13:15,10:12/),1)
         deriv(:,2) = deriv((/1:3,7:9,4:6,13:15,10:12/),2)
         deriv(:,3) = deriv((/1:3,7:9,4:6,13:15,10:12/),3)
       END IF

     ELSE  !standard serendipit prism

       bf = 1d0 - zeta**2                      !bubble function
       shape(1) = zzet*(2d0*zzet-1d0) *l1 - zzet*bf/2d0
       shape(2) = xita*(2d0*xita-1d0) *l1 - xita*bf/2d0
       shape(3) = eta *(2d0*eta -1d0) *l1 - eta *bf/2d0
       shape(7) = 4d0*zzet*xita       *l1
       shape(8) = 4d0*xita*eta        *l1
       shape(9) = 4d0*eta *zzet       *l1
       shape(10) = zzet               *bf
       shape(11) = xita               *bf
       shape(12) = eta                *bf
       shape(4) = zzet*(2d0*zzet-1d0) *l2 - zzet*bf/2d0
       shape(5) = xita*(2d0*xita-1d0) *l2 - xita*bf/2d0
       shape(6) = eta *(2d0*eta -1d0) *l2 - eta *bf/2d0
       shape(13) = 4d0*zzet*xita      *l2
       shape(14) = 4d0*xita*eta       *l2
       shape(15) = 4d0*eta *zzet      *l2

       deriv(1,1) = (-4d0*zzet+1d0)   *l1 + bf/2d0
       deriv(2,1) =  (4d0*xita-1d0)   *l1 - bf/2d0
       deriv(3,1) =  0d0
       deriv(7,1) =  4d0*(zzet-xita)  *l1
       deriv(8,1) =  4d0*eta          *l1
       deriv(9,1) = -4d0*eta          *l1
       deriv(10,1) = -                 bf
       deriv(11,1) =                   bf
       deriv(12,1) =  0d0
       deriv(4,1) = (-4d0*zzet+1d0)   *l2 + bf/2d0
       deriv(5,1) =  (4d0*xita-1d0)   *l2 - bf/2d0
       deriv(6,1) =  0d0
       deriv(13,1) =  4d0*(zzet-xita) *l2
       deriv(14,1) =  4d0*eta         *l2
       deriv(15,1) = -4d0*eta         *l2

       deriv(1,2) = (-4d0*zzet+1d0)   *l1 + bf/2d0
       deriv(2,2) =  0d0
       deriv(3,2) =  (4d0*eta -1d0)   *l1 - bf/2d0
       deriv(7,2) = -4d0*xita         *l1
       deriv(8,2) =  4d0*xita         *l1
       deriv(9,2) =  4d0*(-eta+zzet)  *l1
       deriv(10,2) = -                 bf
       deriv(11,2) =  0d0
       deriv(12,2) =                   bf
       deriv(4,2) = (-4d0*zzet+1d0)   *l2 + bf/2d0
       deriv(5,2) =  0d0
       deriv(6,2) =  (4d0*eta -1d0)   *l2 - bf/2d0
       deriv(13,2) = -4d0*xita        *l2
       deriv(14,2) =  4d0*xita        *l2
       deriv(15,2) = 4d0*(-eta+zzet)  *l2

       deriv(1,3) = -zzet*(2d0*zzet-1d0)/2d0 + zzet*zeta
       deriv(2,3) = -xita*(2d0*xita-1d0)/2d0 + xita*zeta
       deriv(3,3) = -eta *(2d0*eta -1d0)/2d0 + eta *zeta
       deriv(7,3) = -2d0*zzet*xita
       deriv(8,3) = -2d0*xita*eta
       deriv(9,3) = -2d0*eta *zzet
       deriv(10,3) = -zzet*2d0*zeta
       deriv(11,3) = -xita*2d0*zeta
       deriv(12,3) = -eta *2d0*zeta
       deriv(4,3) =  zzet*(2d0*zzet-1d0)/2d0 + zzet*zeta
       deriv(5,3) =  xita*(2d0*xita-1d0)/2d0 + xita*zeta
       deriv(6,3) =  eta *(2d0*eta -1d0)/2d0 + eta *zeta
       deriv(13,3) = 2d0*zzet*xita
       deriv(14,3) = 2d0*xita*eta
       deriv(15,3) = 2d0*eta *zzet
     END IF
   CASE (18)
     IF( bez )THEN
       ! bezier quadratic polynomials
       l12= l1*l1
       l22= l2*l2
       bf  = 1d0-l12-l22  ! (1d0 - zeta*zeta)/2d0

       shape( 1) = zzet * zzet          * l12
       shape( 2) = xita * xita          * l12
       shape( 3) = eta  * eta           * l12
       shape( 7) = 2.0d0  * xita * zzet * l12
       shape( 8) = 2.0d0  * eta  * xita * l12
       shape( 9) = 2.0d0  * zzet * eta  * l12
       shape(10) = zzet * zzet          * bf
       shape(11) = xita * xita          * bf
       shape(12) = eta  * eta           * bf
       shape(16) = 2.0d0  * xita * zzet * bf
       shape(17) = 2.0d0  * eta  * xita * bf
       shape(18) = 2.0d0  * zzet * eta  * bf
       shape( 4) = zzet * zzet          * l22
       shape( 5) = xita * xita          * l22
       shape( 6) = eta  * eta           * l22
       shape(13) = 2.0d0  * xita * zzet * l22
       shape(14) = 2.0d0  * eta  * xita * l22
       shape(15) = 2.0d0  * zzet * eta  * l22

       deriv( 1,1) = -2.0d0  * zzet        * l12
       deriv( 2,1) =  2.0d0  * xita        * l12
       deriv( 3,1) =  0.d0
       deriv( 7,1) =  2.0d0  *(zzet - xita)* l12
       deriv( 8,1) =  2.0d0  * eta         * l12
       deriv( 9,1) = -2.0d0  * eta         * l12
       deriv(10,1) = -2.0d0  * zzet        * bf
       deriv(11,1) =  2.0d0  * xita        * bf
       deriv(12,1) =  0.d0
       deriv(16,1) =  2.0d0  *(zzet - xita)* bf
       deriv(17,1) =  2.0d0  * eta         * bf
       deriv(18,1) = -2.0d0  * eta         * bf
       deriv( 4,1) = -2.0d0  * zzet        * l22
       deriv( 5,1) =  2.0d0  * xita        * l22
       deriv( 6,1) =  0.d0
       deriv(13,1) =  2.0d0  *(zzet - xita)* l22
       deriv(14,1) =  2.0d0  * eta         * l22
       deriv(15,1) = -2.0d0  * eta         * l22

       deriv( 1,2) = -2.0d0  * zzet        * l12
       deriv( 2,2) =  0.0d0
       deriv( 3,2) =  2.0d0  * eta         * l12
       deriv( 7,2) = -2.0d0  * xita        * l12
       deriv( 8,2) =  2.0d0  * xita        * l12
       deriv( 9,2) =  2.0d0  *(zzet - eta )* l12
       deriv(10,2) = -2.0d0  * zzet        * bf
       deriv(11,2) =  0.0d0
       deriv(12,2) =  2.0d0  * eta         * bf
       deriv(16,2) = -2.0d0  * xita        * bf
       deriv(17,2) =  2.0d0  * xita        * bf
       deriv(18,2) =  2.0d0  *(zzet - eta )* bf
       deriv( 4,2) = -2.0d0  * zzet        * l22
       deriv( 5,2) =  0.0d0
       deriv( 6,2) =  2.0d0  * eta         * l22
       deriv(13,2) = -2.0d0  * xita        * l22
       deriv(14,2) =  2.0d0  * xita        * l22
       deriv(15,2) =  2.0d0  *(zzet - eta )* l22

       deriv( 1,3) = -zzet * zzet          * l1
       deriv( 2,3) = -xita * xita          * l1
       deriv( 3,3) = -eta  * eta           * l1
       deriv( 7,3) = -2.0d0  * xita * zzet * l1
       deriv( 8,3) = -2.0d0  * eta  * xita * l1
       deriv( 9,3) = -2.0d0  * zzet * eta  * l1
       deriv(10,3) = -zzet * zzet          * zeta
       deriv(11,3) = -xita * xita          * zeta
       deriv(12,3) = -eta  * eta           * zeta
       deriv(16,3) = -2.0d0  * xita * zzet * zeta
       deriv(17,3) = -2.0d0  * eta  * xita * zeta
       deriv(18,3) = -2.0d0  * zzet * eta  * zeta
       deriv( 4,3) =  zzet * zzet          * l2
       deriv( 5,3) =  xita * xita          * l2
       deriv( 6,3) =  eta  * eta           * l2
       deriv(13,3) =  2.0d0  * xita * zzet * l2
       deriv(14,3) =  2.0d0  * eta  * xita * l2
       deriv(15,3) =  2.0d0  * zzet * eta  * l2

     ELSE  !standard quadratic prism

       l12 = -zeta*(1d0-zeta)/2d0  !quadratic functions N1
       bf  = 1d0 - zeta**2         !bubble function (N2)
       l22 =  zeta*(1d0+zeta)/2d0  !quadratic functions N3
       shape( 1) = zzet*(2d0*zzet-1d0) *l12
       shape( 2) = xita*(2d0*xita-1d0) *l12
       shape( 3) = eta *(2d0*eta -1d0) *l12
       shape( 7) = 4d0*zzet*xita       *l12
       shape( 8) = 4d0*xita*eta        *l12
       shape( 9) = 4d0*eta *zzet       *l12
       shape(10) = zzet*(2d0*zzet-1d0) *bf
       shape(11) = xita*(2d0*xita-1d0) *bf
       shape(12) = eta *(2d0*eta -1d0) *bf
       shape(16) = 4d0*zzet*xita       *bf
       shape(17) = 4d0*xita*eta        *bf
       shape(18) = 4d0*eta *zzet       *bf
       shape( 4) = zzet*(2d0*zzet-1d0) *l22
       shape( 5) = xita*(2d0*xita-1d0) *l22
       shape( 6) = eta *(2d0*eta -1d0) *l22
       shape(13) = 4d0*zzet*xita       *l22
       shape(14) = 4d0*xita*eta        *l22
       shape(15) = 4d0*eta *zzet       *l22

       deriv( 1,1) =-(4d0*zzet-1d0)      *l12
       deriv( 2,1) = (4d0*xita-1d0)      *l12
       deriv( 3,1) =  0d0
       deriv( 7,1) =  4d0*(zzet-xita)    *l12
       deriv( 8,1) =  4d0*eta            *l12
       deriv( 9,1) = -4d0*eta            *l12
       deriv(10,1) =-(4d0*zzet-1d0)      *bf
       deriv(11,1) = (4d0*xita-1d0)      *bf
       deriv(12,1) =  0d0
       deriv(16,1) =  4d0*(zzet-xita)    *bf
       deriv(17,1) =  4d0*eta            *bf
       deriv(18,1) = -4d0*eta            *bf
       deriv( 4,1) =-(4d0*zzet-1d0)      *l22
       deriv( 5,1) = (4d0*xita-1d0)      *l22
       deriv( 6,1) =  0d0
       deriv(13,1) =  4d0*(zzet-xita)    *l22
       deriv(14,1) =  4d0*eta            *l22
       deriv(15,1) = -4d0*eta            *l22

       deriv( 1,2) =-(4d0*zzet-1d0)      *l12
       deriv( 2,2) =  0d0
       deriv( 3,2) = (4d0*eta -1d0)      *l12
       deriv( 7,2) = -4d0*xita           *l12
       deriv( 8,2) =  4d0*xita           *l12
       deriv( 9,2) =  4d0*(zzet-eta)     *l12
       deriv(10,2) =-(4d0*zzet-1d0)      *bf
       deriv(11,2) =  0d0
       deriv(12,2) = (4d0*eta -1d0)      *bf
       deriv(16,2) = -4d0*xita           *bf
       deriv(17,2) =  4d0*xita           *bf
       deriv(18,2) =  4d0*(zzet-eta)     *bf
       deriv( 4,2) =-(4d0*zzet-1d0)      *l22
       deriv( 5,2) =  0d0
       deriv( 6,2) = (4d0*eta -1d0)      *l22
       deriv(13,2) = -4d0*xita           *l22
       deriv(14,2) =  4d0*xita           *l22
       deriv(15,2) =  4d0*(zzet-eta)     *l22

       deriv( 1,3) =  zzet*(2d0*zzet-1d0)*(zeta-0.5d0)
       deriv( 2,3) =  xita*(2d0*xita-1d0)*(zeta-0.5d0)
       deriv( 3,3) =  eta *(2d0*eta -1d0)*(zeta-0.5d0)
       deriv( 7,3) =  4d0*zzet*xita      *(zeta-0.5d0)
       deriv( 8,3) =  4d0*xita*eta       *(zeta-0.5d0)
       deriv( 9,3) =  4d0*eta *zzet      *(zeta-0.5d0)
       deriv(10,3) =  zzet*(2d0*zzet-1d0)*(-2d0*zeta)
       deriv(11,3) =  xita*(2d0*xita-1d0)*(-2d0*zeta)
       deriv(12,3) =  eta *(2d0*eta -1d0)*(-2d0*zeta)
       deriv(16,3) =  4d0*zzet*xita      *(-2d0*zeta)
       deriv(17,3) =  4d0*xita*eta       *(-2d0*zeta)
       deriv(18,3) =  4d0*eta *zzet      *(-2d0*zeta)
       deriv( 4,3) =  zzet*(2d0*zzet-1d0)*(zeta+0.5d0)
       deriv( 5,3) =  xita*(2d0*xita-1d0)*(zeta+0.5d0)
       deriv( 6,3) =  eta *(2d0*eta -1d0)*(zeta+0.5d0)
       deriv(13,3) =  4d0*zzet*xita      *(zeta+0.5d0)
       deriv(14,3) =  4d0*xita*eta       *(zeta+0.5d0)
       deriv(15,3) =  4d0*eta *zzet      *(zeta+0.5d0)
     END IF
   END SELECT
   ! verif
   l1 = SUM(shape)
   l2 = SUM(deriv(:,1))
   l12= SUM(deriv(:,2))
   l22= SUM(deriv(:,3))

 RETURN
 END SUBROUTINE shape4
