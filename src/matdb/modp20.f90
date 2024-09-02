 SUBROUTINE modp20(dmatx,ntype,young,poiss,prop,efpst,defps, &
                   bbar,elast,newmt)
 !*********************************************************************
 !
 !**** this SUBROUTINE evaluates the elastic d-matrix (upper part only)
 !
 !
 !**********************************************************************
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: bbar,elast
 LOGICAL, INTENT(IN OUT) :: newmt
 INTEGER (kind=4), INTENT(IN) :: ntype
 REAL (kind=8), INTENT(IN) :: prop(5),young,poiss,defps,efpst
 REAL (kind=8), INTENT(OUT) :: dmatx(4,4)

 ! --------- imodf=0 elastic  matrix
 ! --------- imodf=1 elastic-plastic matrix
 !INTEGER (kind=4), PARAMETER :: imodf = 1  !elastic-plastic matrix

 INTEGER (kind=4), SAVE :: is
 REAL (kind=8), SAVE :: kvol,g1,g2,d(4,4),c0,e0,expn,hard,rfs
 REAL (kind=8) :: g1n,g2n,ddiag,dndia,pmult,yield

 IF( newmt )THEN
   newmt = .FALSE.
   g1 = young/(1-poiss)
   g2 = young/(1d0+poiss)
   IF(bbar) THEN
     kvol = 0d0
   ELSE
     kvol = young/(1d0-2d0*poiss)/3d0
   END IF

   d(1:2,3) = 0d0
   d(3,4) = 0d0
   IF( ntype == 1 ) THEN
     d(1,1) = (g1+g2)/2d0
     d(1,2) = (g1-g2)/2d0
     d(2,2) = (g1+g2)/2d0
     d(3,3) =      g2/2d0
   ELSE
     ddiag = kvol + 2d0 * g2 /3d0
     dndia = kvol - g2/3d0

     d(1,1) = ddiag
     d(1,2) = dndia
     d(1,4) = dndia
     d(2,2) = ddiag
     d(2,4) = dndia
     d(3,3) = g2/2d0
     d(4,4) = ddiag
   END IF
   dmatx = d               ! to avoid multiple calls and assignation in elastic problems
   IF( .NOT.elast )THEN    ! keep constant of isotropic hardening
     is =  INT(prop(5))    ! type of isotropic hardening
     SELECT CASE ( is )
     CASE (1)            !no hardening
       c0    = prop(1)     !Initial yield
     CASE (2)            !linear hardening
       c0    = prop(1)     !Initial yield
       hard  = prop(2)     !hardening modulus
     CASE (3)            !exponentical hardening
       c0   = prop(1)      !c0 constant
       e0   = prop(2)      !offset value
       expn = prop(3)      !exponent for non-linear hardening
     CASE (4)            !linear + saturation law
       c0   = prop(1)      !c0 constant
       hard = prop(2)      !linear hardening
       expn = prop(3)      !exponent for non-linear hardening
       rfs  = prop(4)      !residual flow stress
     END SELECT
   END IF
   RETURN
 END IF

 IF( elast ) RETURN        !do not assign again
 !IF( defps == 0d0)THEN     !if the step is elastic
   dmatx = d               !assign elastic modulus
   RETURN                  ! & return
 !END IF
 !    setup yield FUNCTION radius
 SELECT CASE (is)
 CASE (1)
   yield = c0                  !Initial yield
 CASE (2)
   yield = c0 + hard*efpst     !linear hardening
 CASE (3)
   yield = c0*(e0+efpst)**expn !non-linear (exponential) hardening
 CASE (4)
   yield = c0+hard*efpst+(rfs-c0)*(1d0-1d0/EXP(expn*efpst)) !linear + saturation law hardening
 END SELECT
 ! set pmult (beta factor)
 pmult = g2*defps/yield

 dmatx(1:2,3) = 0d0
 dmatx(3,4) = 0d0

 IF(ntype == 1) THEN            ! D matrix for plane stress case

   !IF(imodf == 1) THEN         ! elastic-plastic tangent matrix
     g1n = g1 / (1d0 + pmult*g1/3d0)
     g2n = g2 / (1d0 + pmult*g2)
   !END IF

   dmatx(1,1) = (g1n+g2n)/2d0
   dmatx(1,2) = (g1n-g2n)/2d0
   dmatx(2,2) = (g1n+g2n)/2d0
   dmatx(3,3) =      g2n/2d0

 ELSE     !ntype = 2 plane strain case or
          !ntype = 3 axisymmetric case

   !IF(imodf == 1) g2n = g2 * (1d0-pmult) ! elastic plastic tangent matrix
   g2n = g2 * (1d0-pmult) ! elastic-plastic tangent matrix

   ddiag = kvol + 2d0 * g2n /3d0
   dndia = kvol - g2n/3d0

   dmatx(1,1) = ddiag
   dmatx(1,2) = dndia
   dmatx(1,4) = dndia
   dmatx(2,2) = ddiag
   dmatx(2,4) = dndia
   dmatx(3,3) = g2n/2d0
   dmatx(4,4) = ddiag

 END IF
 RETURN
 END SUBROUTINE modp20
