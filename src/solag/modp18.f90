 SUBROUTINE modp18(dmatx,young,poiss,prop,efpst,defps, &
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
 REAL (kind=8), INTENT(IN) :: prop(5),young,poiss,defps,efpst
 REAL (kind=8), INTENT(OUT) :: dmatx(6,6)

 ! --------- imodf=0 elastic  matrix
 ! --------- imodf=1 elastic-plastic matrix
 !INTEGER (kind=4), PARAMETER :: imodf = 1  !elastic-plastic matrix

 INTEGER (kind=4), SAVE :: is
 REAL (kind=8), SAVE :: kvol,g,d(6,6),c0,e0,expn,hard,rfs
 REAL (kind=8) :: gn,ddiag,dndia,pmult,yield

 IF( newmt )THEN
   newmt = .FALSE.
   g = young/(1d0+poiss)/2d0
   IF(bbar) THEN
     kvol = 0d0
   ELSE
     kvol = young/(1d0-2d0*poiss)/3d0
   END IF

   d(1:3,4:6) = 0d0
   d(4,5) = 0d0
   d(4,6) = 0d0
   d(5,6) = 0d0

   ddiag = kvol + 4d0 * g /3d0
   dndia = kvol - 2d0 * g /3d0

   d(1,1) = ddiag
   d(1,2) = dndia
   d(1,3) = dndia
   d(2,2) = ddiag
   d(2,3) = dndia
   d(3,3) = ddiag
   d(2,1) = dndia
   d(3,1) = dndia
   d(3,2) = dndia
   d(4,4) = g
   d(5,5) = g
   d(6,6) = g

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
 IF( defps == 0d0)THEN     !if the step is elastic
   dmatx = d               !assign elastic modulus
   RETURN                  ! & return
 END IF
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
 pmult = g*2d0*defps/yield

 !IF(modf) g = g * (1d0-pmult) ! consistent tangent matrix
 gn = g * (1d0-pmult) ! consistent tangent matrix

 dmatx(1:3,4:6) = 0d0
 dmatx(4,5) = 0d0
 dmatx(4,6) = 0d0
 dmatx(5,6) = 0d0

 ddiag = kvol + 4d0 * gn /3d0
 dndia = kvol - 2d0 * gn /3d0

 dmatx(1,1) = ddiag
 dmatx(1,2) = dndia
 dmatx(1,3) = dndia
 dmatx(2,2) = ddiag
 dmatx(2,3) = dndia
 dmatx(3,3) = ddiag
 dmatx(4,4) = gn
 dmatx(5,5) = gn
 dmatx(6,6) = gn

 RETURN
 END SUBROUTINE modp18
