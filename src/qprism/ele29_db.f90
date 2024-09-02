 MODULE ele29_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input  !All?
  IMPLICIT NONE

   ! Total Lagrangial (Log strain) 3-D Bezier prism Solid-Shell element
   INTEGER, PARAMETER :: nvare = 14, &!number of internal variables per integration point
                         nstre =  6, &! number of stress measures
                         nnode = 15, &! number of nodes per element
                         ndofe = 45, &! number of DOFs per element
                         nasmm =  9, &! number of Assumed Strain for membrane
                         nnas1 =  6, &! number of shear assumed strain points (lineal OZTS)
                         nnas2 =  8, &! number of shear assumed strain points (lineal + Quad)
                         nface =  2, &! number of faces
                         ngaup =  3, &! number of in-plane Gauss points
                         ngaud =  3   ! number of Gauss point to define volume

!  In-plane Gauss points (zeta,xita,eta)                   zeta   xita   eta   for membrane interpolation
   REAL (kind=8), PARAMETER :: psg_in(2,ngaup) = RESHAPE ( (/ 0.166666666667D0,   0.166666666667D0, &
                                                              0.666666666666D0,   0.166666666667D0, &
                                                              0.166666666667D0,   0.666666666666D0 /), (/2,3/))
   REAL (kind=8), PARAMETER :: psg_ms(2,ngaup) = RESHAPE ( (/ 0.5d0, 0.0d0, &
                                                              0.5D0, 0.5D0, &
                                                              0.0d0, 0.5D0 /), (/2,3/))


! PARAMETERS FOR MEMBRANE FORMULATION
! For model 2: Sampling points at mid-side of each sub-triangle
!   sampling points
    REAL (kind=8), PARAMETER :: gpmm(2,nasmm) =  RESHAPE ((/  &        !
       0.25D0, 0.00D0, 0.75D0, 0.00D0, 0.25D0, 0.50D0, 0.00D0, 0.25D0, 0.00D0, 0.75D0, &
       0.50D0, 0.25D0, 0.75D0, 0.25D0, 0.25D0, 0.75D0, 0.25D0, 0.25D0 /),(/2,nasmm/))
!  tangent derivatives at sampling points
    REAL (kind=8) :: ntan2(nnode,nasmm,nface)

! interpolates the assumed strain values to Gauss points
   REAL (kind=8), PARAMETER :: pa2_in(3,nasmm,ngaup) =  RESHAPE ((/                     &
+0.8333333333333D0, +0.0000000000000D0, +0.4166666666667D0, -0.1666666666667D0, +0.0000000000000D0, -0.0833333333333D0,  &
+0.3333333333333D0, +0.0000000000000D0, +0.1666666666667D0, +0.0000000000000D0, +0.8333333333333D0, +0.4166666666667D0,  &
+0.0000000000000D0, -0.1666666666667D0, -0.0833333333333D0, +0.0000000000000D0, +0.3333333333333D0, +0.1666666666667D0,  &
+0.0000000000000D0, +0.0000000000000D0, +0.1666666666667D0, +0.0000000000000D0, +0.0000000000000D0, +0.1666666666667D0,  &
+0.0000000000000D0, +0.0000000000000D0, -1.3333333333333D0,                                                              &
-0.1666666666667D0, +0.0000000000000D0, -0.0833333333333D0, +0.8333333333333D0, +0.0000000000000D0, +0.4166666666667D0,  &
+0.3333333333333D0, +0.0000000000000D0, +0.1666666666667D0, +0.0000000000000D0, -0.1666666666667D0, -0.0833333333333D0,  &
+0.0000000000000D0, -0.1666666666667D0, -0.0833333333333D0, +0.0000000000000D0, +1.3333333333333D0, +0.6666666666667D0,  &
+0.0000000000000D0, +0.0000000000000D0, -0.8333333333333D0, +0.0000000000000D0, +0.0000000000000D0, +0.1666666666667D0,  &
+0.0000000000000D0, +0.0000000000000D0, -0.3333333333333D0,                                                              &
-0.1666666666667D0, +0.0000000000000D0, -0.0833333333333D0, -0.1666666666667D0, +0.0000000000000D0, -0.0833333333333D0,  &
+1.3333333333333D0, +0.0000000000000D0, +0.6666666666667D0, +0.0000000000000D0, -0.1666666666667D0, -0.0833333333333D0,  &
+0.0000000000000D0, +0.8333333333333D0, +0.4166666666667D0, +0.0000000000000D0, +0.3333333333333D0, +0.1666666666667D0,  &
+0.0000000000000D0, +0.0000000000000D0, +0.1666666666667D0, +0.0000000000000D0, +0.0000000000000D0, -0.8333333333333D0,  &
+0.0000000000000D0, +0.0000000000000D0, -0.3333333333333D0  /),(/ 3,nasmm,ngaup /))
   REAL (kind=8), PARAMETER :: pa2b_in(3,nasmm,ngaup) =  RESHAPE ((/                    &
 0.8333333333333D0,  0.0000000000000D0,  0.8333333333333D0, -0.1666666666667D0,  0.0000000000000D0, -0.1666666666667D0,  &
 0.3333333333333D0,  0.0000000000000D0,  0.3333333333333D0,  0.0000000000000D0,  0.8333333333333D0,  0.8333333333333D0,  &
 0.0000000000000D0, -0.1666666666667D0, -0.1666666666667D0,  0.0000000000000D0,  0.3333333333333D0,  0.3333333333333D0,  &
 0.0000000000000D0,  0.0000000000000D0,  0.3333333333333D0,  0.0000000000000D0,  0.0000000000000D0,  0.3333333333333D0,  &
 0.0000000000000D0,  0.0000000000000D0, -2.6666666666667D0,                                                              &
-0.1666666666667D0,  0.0000000000000D0, -0.1666666666667D0,  0.8333333333333D0,  0.0000000000000D0,  0.8333333333333D0,  &
 0.3333333333333D0,  0.0000000000000D0,  0.3333333333333D0,  0.0000000000000D0, -0.1666666666667D0, -0.1666666666667D0,  &
 0.0000000000000D0, -0.1666666666667D0, -0.1666666666667D0,  0.0000000000000D0,  1.3333333333333D0,  1.3333333333333D0,  &
 0.0000000000000D0,  0.0000000000000D0, -1.6666666666667D0,  0.0000000000000D0,  0.0000000000000D0,  0.3333333333333D0,  &
 0.0000000000000D0,  0.0000000000000D0, -0.6666666666667D0,                                                              &
-0.1666666666667D0,  0.0000000000000D0, -0.1666666666667D0, -0.1666666666667D0,  0.0000000000000D0, -0.1666666666667D0,  &
 1.3333333333333D0,  0.0000000000000D0,  1.3333333333333D0,  0.0000000000000D0, -0.1666666666667D0, -0.1666666666667D0,  &
 0.0000000000000D0,  0.8333333333333D0,  0.8333333333333D0,  0.0000000000000D0,  0.3333333333333D0,  0.3333333333333D0,  &
 0.0000000000000D0,  0.0000000000000D0,  0.3333333333333D0,  0.0000000000000D0,  0.0000000000000D0, -1.6666666666667D0,  &
 0.0000000000000D0,  0.0000000000000D0, -0.6666666666667D0    /),(/ 3,nasmm,ngaup /))

   REAL (kind=8), PARAMETER :: pa2_ms(3,nasmm,ngaup) =  RESHAPE ((/                     &
    +0.50D0, +0.00D0, +0.25D0, +0.50D0, +0.00D0, +0.25D0, +0.00D0, +0.00D0, +0.00D0, &
    +0.00D0, +0.50D0, +0.25D0, +0.00D0, -0.50D0, -0.25D0, +0.00D0, +1.00D0, +0.50D0, &
    +0.00D0, +0.00D0, -0.50D0, +0.00D0, +0.00D0, +0.50D0, +0.00D0, +0.00D0, -1.00D0, &
    -0.50D0, +0.00D0, -0.25D0, +0.50D0, +0.00D0, +0.25D0, +1.00D0, +0.00D0, +0.50D0, &
    +0.00D0, -0.50D0, -0.25D0, +0.00D0, +0.50D0, +0.25D0, +0.00D0, +1.00D0, +0.50D0, &
    +0.00D0, +0.00D0, -0.50D0, +0.00D0, +0.00D0, -0.50D0, +0.00D0, +0.00D0, +0.00D0, &
    +0.50D0, +0.00D0, +0.25D0, -0.50D0, +0.00D0, -0.25D0, +1.00D0, +0.00D0, +0.50D0, &
    +0.00D0, +0.50D0, +0.25D0, +0.00D0, +0.50D0, +0.25D0, +0.00D0, +0.00D0, +0.00D0, &
    +0.00D0, +0.00D0, +0.50D0, +0.00D0, +0.00D0, -0.50D0, +0.00D0, +0.00D0, -1.00D0  /),(/ 3,nasmm,ngaup /))

   REAL (kind=8), PARAMETER :: pa2b_ms(3,nasmm,ngaup) =  RESHAPE ((/                    &
    +0.50D+0, +0.00D+0, +0.50D+0, +0.50D+0, +0.00D+0, +0.50D+0, +0.00D+0, +0.00D+0, +0.00D+0,  &
    +0.00D+0, +0.50D+0, +0.50D+0, +0.00D+0, -0.50D+0, -0.50D+0, +0.00D+0, +1.00D+0, +1.00D+0,  &
    +0.00D+0, +0.00D+0, -1.00D+0, +0.00D+0, +0.00D+0, +1.00D+0, +0.00D+0, +0.00D+0, -2.00D+0,  &
    -0.50D+0, +0.00D+0, -0.50D+0, +0.50D+0, +0.00D+0, +0.50D+0, +1.00D+0, +0.00D+0, +1.00D+0,  &
    +0.00D+0, -0.50D+0, -0.50D+0, +0.00D+0, +0.50D+0, +0.50D+0, +0.00D+0, +1.00D+0, +1.00D+0,  &
    +0.00D+0, +0.00D+0, -1.00D+0, +0.00D+0, +0.00D+0, -1.00D+0, +0.00D+0, +0.00D+0, +0.00D+0,  &
    +0.50D+0, +0.00D+0, +0.50D+0, -0.50D+0, +0.00D+0, -0.50D+0, +1.00D+0, +0.00D+0, +1.00D+0,  &
    +0.00D+0, +0.50D+0, +0.50D+0, +0.00D+0, +0.50D+0, +0.50D+0, +0.00D+0, +0.00D+0, +0.00D+0,  &
    +0.00D+0, +0.00D+0, +1.00D+0, +0.00D+0, +0.00D+0, -1.00D+0, +0.00D+0, +0.00D+0, -2.00D+0/),(/ 3,nasmm,ngaup /))

! PARAMETERS FOR SHEAR FORMULATION (mixed assumed shear)

   REAL (kind=8), PARAMETER :: am1(nnas1,nnas1) =  RESHAPE ((/  &        ! A-matrix (Lineal)
   1.3660254037844, -1.7320508075689, -1.3660254037844,  0.0000000000000, -0.3660254037844, 0.0000000000000, &
  -0.3660254037844,  1.7320508075689,  0.3660254037844,  0.0000000000000,  1.3660254037844, 0.0000000000000, &
   0.0000000000000,  0.0000000000000,  0.5176380902050,  0.0000000000000,  1.9318516525781, 0.0000000000000, &
   0.0000000000000,  0.0000000000000, -1.9318516525781,  0.0000000000000, -0.5176380902050, 0.0000000000000, &
   0.0000000000000,  0.0000000000000,  1.3660254037844, -0.3660254037844,  0.3660254037844, 1.7320508075689, &
   0.0000000000000,  0.0000000000000, -0.3660254037844,  1.3660254037844, -1.3660254037844,-1.7320508075689 /),(/nnas1,nnas1/))

   REAL (kind=8), PARAMETER :: am2(nnas2,nnas2) =  RESHAPE ((/  &        ! A-matrix (BATHE et al)
 1.3660254037844,-1.7320508075689,-3.7320508075689, 1.7320508075689, 0.0000000000000, 1.3660254037844, 0.0000000000000,-2.3660254037844, &
-0.3660254037844, 1.7320508075689,-0.2679491924311,-1.7320508075689, 0.0000000000000,-0.3660254037844, 0.0000000000000,-0.6339745962156, &
 0.0000000000000, 0.0000000000000, 1.4142135623731,-3.3460652149512, 0.0000000000000,-1.4142135623731, 0.0000000000000, 0.8965754721681, &
 0.0000000000000, 0.0000000000000, 1.4142135623731,-0.8965754721681, 0.0000000000000,-1.4142135623731, 0.0000000000000, 3.3460652149512, &
 0.0000000000000, 0.0000000000000,-0.3660254037844,-0.6339745962156,-0.3660254037844,-0.2679491924311, 1.7320508075689,-1.7320508075689, &
 0.0000000000000, 0.0000000000000, 1.3660254037844,-2.3660254037844, 1.3660254037844,-3.7320508075689,-1.7320508075689, 1.7320508075689, &
 0.0000000000000, 0.0000000000000, 6.0000000000000,-3.0000000000000, 0.0000000000000,-3.0000000000000, 0.0000000000000, 6.0000000000000, &
 0.0000000000000, 0.0000000000000,-3.0000000000000, 6.0000000000000, 0.0000000000000, 6.0000000000000, 0.0000000000000,-3.0000000000000  &
   /),(/nnas2,nnas2/))

!   tangent derivatives at each assumed strain point
   REAL (kind=8) :: ntan(nnode,nnas2,nface)

!   PA Matrix Product at in-plane Gauss Points (incomplete quadratic approach)
   REAL (kind=8), PARAMETER :: pag_in(2,nnas2,ngaup) =  RESHAPE ((/  &
+0.5691772515768,+0.1138354503154, -0.1525105849102,-0.0305021169820,   &
+0.1178511301978,-0.1178511301978, +0.1178511301978,-0.1178511301978,   &
-0.0305021169820,-0.1525105849102, +0.1138354503154,+0.5691772515769,   &
+0.7500000000000,-0.2500000000000, -0.2500000000000,+0.7500000000000,   &
-0.1525105849102,-0.1220084679281, +0.5691772515768,+0.4553418012615,   &
-0.1609876377148,+0.6439505508594, +0.0431365075171,-0.1725460300683,   &
-0.0833333333333,-0.1666666666667, -0.0833333333333,-0.1666666666667,   &
+0.5000000000000,+0.0000000000000, +0.2500000000000,+1.0000000000000,   &
-0.1666666666667,-0.0833333333333, -0.1666666666667,-0.0833333333333,   &
+0.1725460300683,-0.0431365075171, -0.6439505508594,+0.1609876377148,   &
+0.4553418012615,+0.5691772515768, -0.1220084679281,-0.1525105849102,   &
+1.0000000000000,+0.2500000000000, +0.0000000000000,+0.5000000000000    &
/),(/ 2,nnas2,ngaup /))


   REAL (kind=8), PARAMETER :: pag_ms(2,nnas2,ngaup) =  RESHAPE ((/  &
+0.5000000000000,+0.2500000000000, +0.5000000000000,+0.2500000000000,   &
+0.0000000000000,+0.1294095225513, +0.0000000000000,-0.4829629131445,   &
+0.0000000000000,-0.3415063509461, +0.0000000000000,+0.0915063509461,   &
+0.0000000000000,-0.7500000000000, +0.0000000000000,+1.5000000000000,   &
-0.3415063509461,-0.3415063509461, +0.0915063509461,+0.0915063509461,   &
-0.3535533905933,+0.3535533905933, -0.3535533905933,+0.3535533905933,   &
+0.0915063509461,+0.0915063509461, -0.3415063509461,-0.3415063509461,   &
+0.7500000000000,+0.7500000000000, +0.7500000000000,+0.7500000000000,   &
+0.0915063509461,+0.0000000000000, -0.3415063509461,+0.0000000000000,   &
+0.4829629131445,+0.0000000000000, -0.1294095225513,+0.0000000000000,   &
+0.2500000000000,+0.5000000000000, +0.2500000000000,+0.5000000000000,   &
+1.5000000000000,+0.0000000000000, -0.7500000000000,+0.0000000000000    &
/),(/ 2,nnas2,ngaup /))

!   PA Matrix Product at in-plane Gauss Points (linear approach)
   REAL (kind=8), PARAMETER :: pag0_in(2,nnas1,ngaup) =  RESHAPE ((/  &
 0.8496793685589,-0.0610042339641,-0.0163460352256,-0.2440169358563, 0.1666666666667,-0.0610042339641, &
-0.0163460352256, 0.2276709006307, 0.8496793685589, 0.9106836025230, 0.1666666666667, 0.2276709006307, &
 0.0862730150342, 0.3219752754297, 0.0862730150342, 1.2879011017188, 0.3450920601367, 0.3219752754297, &
-0.3219752754297,-0.0862730150342,-0.3219752754297,-0.3450920601367,-1.2879011017188,-0.0862730150342, &
 0.2276709006307,-0.0163460352256, 0.2276709006307, 0.1666666666667, 0.9106836025230, 0.8496793685589, &
-0.0610042339641, 0.8496793685589,-0.0610042339641, 0.1666666666667,-0.2440169358563,-0.0163460352256  &
/),(/ 2,nnas1,ngaup /))

   REAL (kind=8), PARAMETER :: pag0_ms(2,nnas1,ngaup) =  RESHAPE ((/  &
 0.5000000000000,-0.1830127018922, 0.5000000000000, 0.6830127018922, 0.0000000000000, 0.9659258262891, &
 0.0000000000000,-0.2588190451025, 0.0000000000000,-0.1830127018922, 0.0000000000000, 0.6830127018922, &
-0.1830127018922,-0.1830127018922, 0.6830127018922, 0.6830127018922, 0.2588190451025, 0.9659258262891, &
-0.9659258262891,-0.2588190451025, 0.6830127018922, 0.6830127018922,-0.1830127018922,-0.1830127018922, &
 0.6830127018922, 0.0000000000000,-0.1830127018922, 0.0000000000000, 0.2588190451025, 0.0000000000000, &
-0.9659258262891, 0.0000000000000, 0.6830127018922, 0.5000000000000,-0.1830127018922, 0.5000000000000  &
/),(/ 2,nnas1,ngaup /))



! Gauss points position of Assumed Shear Strain Sampling Points
 REAL (kind=8), PARAMETER :: r3_1 = 0.577350269189626D+00, r3 = 1.73205080756888D0, r2= 1.4142135623731D0

 REAL (kind=8), PARAMETER :: gps(2,nnas2,nface) = RESHAPE((/   &
            0.2113248654052D0, 0.0000000000000D0, &      !bottom face
            0.7886751345948D0, 0.0000000000000D0, &
            0.7886751345948D0, 0.2113248654052D0, &
            0.2113248654052D0, 0.7886751345948D0, &
            0.0000000000000D0, 0.7886751345948D0, &
            0.0000000000000D0, 0.2113248654052D0, &
            0.3333333333333D0, 0.3333333333333D0, &
            0.3333333333333D0, 0.3333333333333D0, &
            0.2113248654052D0, 0.0000000000000D0, &      !top face
            0.7886751345948D0, 0.0000000000000D0, &
            0.7886751345948D0, 0.2113248654052D0, &
            0.2113248654052D0, 0.7886751345948D0, &
            0.0000000000000D0, 0.7886751345948D0, &
            0.0000000000000D0, 0.2113248654052D0, &
            0.3333333333333D0, 0.3333333333333D0, &
            0.3333333333333D0, 0.3333333333333D0  /),(/2,nnas2,nface/))
! options for in-plane Gauss Points
 LOGICAL :: midsidegausspoints = .FALSE.  !.FALSE. internal points   .TRUE. mid-side points
 REAL (kind=8) :: psg(2,ngaup),         & ! PoSitions of in-plane Gauss points
                  pa2(3,nasmm,ngaup),   & ! PA matrix for membrane approach
                  pa2b(3,nasmm,ngaup),  & ! (bar) PA matrix for membrane approach
                  pag(2,nnas2,ngaup)      ! PA matrix for shear approach

   ! Derived type for an ELE29 element

   TYPE ele29
     INTEGER (kind=4) :: numel                     ! element label
     INTEGER (kind=4) :: matno                     ! Material number
     INTEGER (kind=4) :: lnods(nnode)              ! Conectivities
     REAL (kind=8) :: angle,                &      ! Euler angle between local t1-t2 and orthotropic
                      dvol(ngaud,ngaup),    &      ! Initial volume
                      cartd(nnode,2,ngaup), &      ! cartesian der (y3) of Shape Functions at in-plane gauss points at faces
                      sem(3,3,2),           &      ! integrated stresses for membrane geometric stiffness
                      ses(nnas2,2),         &      ! integrated stresses for shear geometric stiffness
                      set(ngaup,2),         &      ! and integrated stresses for transvers stress
                      cqi(3,2,ngaup),           &  ! initial metric tensor at each face at each in-plane gauss point
                      ccasi(2,2,ngaup),         &  ! Initial shear stain at each face at each in-plane gauss point
                      c33i(ngaup,2)                ! initial strains for transvers strains
     REAL (kind=8), POINTER :: ipcdm(:,:,:,:),  &  ! InPlaneCartesianDerivativeMembrane (nnode,2,nface,ngaup)
                               nfdas(:,:,:),    &  ! Nodal Function Derivatives at the Assumed Strain points (Shear)(nnode,nassp,nface)
                               jacin(:,:,:,:),  &  ! in-plane inverse jacobian for membrane and shear ANS (2,2,nface,ngaud)
                               stint(:,:,:),    &  ! Actual Kirchhoff stresses (6, ngaus, ngaup)
                               gausv(:,:,:),    &   !Gauss-point internal variables  (7, ngaus, ngaup)
                               gaus0(:,:,:)         !Gauss-point initial strains  (7, ngaus, ngaup)
     TYPE (ele29), POINTER :: next                  !pointer to next element
   END TYPE ele29

   ! Derived type for a set of ELE29 elements
   TYPE ele29_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for History output
                         narch, &  ! number of output unit
                         ngaus, &  ! number of integration points
                         nassp, &  ! number of Assumed Strain for Shear
                         anssh, &  ! shear option 0: displac.  1: linear 2:imcomplete quadratic
                         ansmm     ! membrane option
                                   ! 0:standard
                                   ! 1:ANS at corner nodes (not used)
                                   ! 2:ANS at mid-side points of each subtriangle
                                   ! 3:ANS 6 values at sides and 3 at element center (not implemented)
     LOGICAL :: gauss      ! .FALSE. -> Initial constants not defined or not updated
     LOGICAL :: small      ! .TRUE. -> use Green strain tensor
                           ! .FALSE. -> Use log strains (Default)
     INTEGER :: plstr      ! compute Plastic Strain Flag
                           ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax            ! local x definition option
     REAL (kind=8) :: angdf      ! Default Euler angles between
                                 ! Global X-Y-Z and orthotropic system
     REAL (kind=8) :: psgo(2,2)  ! factors for output at external faces

     INTEGER(kind=4) :: isgo(2,2)! inner gauss points for output

     TYPE (ele29), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele29_set), POINTER :: next      !pointer to next set
   END TYPE ele29_set
   TYPE (ele29_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE ini_ele29 (head, tail)
     !initialize a list of ELE29 sets

     !Dummy arguments
     TYPE (ele29_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele29

   SUBROUTINE add_ele29 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele29_set), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add a segment to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     END IF
   END SUBROUTINE add_ele29

   SUBROUTINE srch_ele29 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele29_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF(TRIM(posic%sname) == TRIM(name)) THEN    !check name
           found = .TRUE.                !set flag to .TRUE.
           EXIT                          !O.K. exit search
         END IF
         IF (ASSOCIATED(posic%next) ) THEN   !there is a next set
           anter => posic                    !point anter to present
           posic => posic%next               !point present to next
         ELSE
           EXIT                              !list exhausted, Exit search
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)   !set not found, null pointers
   END SUBROUTINE srch_ele29

   SUBROUTINE del_ele29 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele29_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele29), POINTER :: ea,ep
     INTEGER (kind=4) :: iel

     IF (.NOT.ASSOCIATED (anter)) THEN  !if anter pointer is null => head
       head => posic%next               !point first to next
     ELSE
       anter%next => posic%next         !skip posic
     END IF

     ! deallocation of the set memory is next done
     NULLIFY( ea )                  !nullify previous element in list
     ep => posic%head               !point present element to first
     DO iel = 1,posic%nelem         !for each element in the set
       CALL del_ele29e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele29

   !----------- Element management routines

   SUBROUTINE ini_ele29e (head, tail)
     !initialize a list of ELE29 elements

     ! dummy arguments
     TYPE (ele29), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele29e

   SUBROUTINE new_ele29e(elm)
   !Create a new element of ELE29 sets

     TYPE(ele29),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods = 0d0      !Initializes connectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvol = 0d0       !     "     elemental volumes
     elm%sem=0d0          !Initializes equivalent stresses for membrane
     elm%ses=0d0          !Initializes equivalent stresses for shear
     elm%set=0d0          !Initializes equivalent stresses for transverse stress
     elm%cartd=0d0        !Initializes cartesyan derivatives (y3)
     elm%c33i = 0d0
     elm%set  = 0d0
     NULLIFY(elm%ipcdm,elm%nfdas,elm%jacin,elm%stint,elm%gausv,elm%gaus0) !initializes pointers to allocatable arrays
     NULLIFY(elm%next)    !Initializes pointer to next element in the list
   RETURN
   END SUBROUTINE new_ele29e

   SUBROUTINE add_ele29e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele29), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN       !list is empty, start it
       head => new                           !first element
       tail => new                           !last element
       NULLIFY (tail%next)                   !last poit to nothing

     ELSE                                    !add a segment to the list
       tail%next => new                      !point to the new element
       NULLIFY (new%next)                    !nothing beyond the last
       tail => new                           !new last element

     END IF
   END SUBROUTINE add_ele29e

   SUBROUTINE srch_ele29e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele29), POINTER :: head, anter, posic

     found = .FALSE.             !initializes flag and pointers
     NULLIFY (posic,anter)

     IF (ASSOCIATED (head)) THEN          !Check if a list is empty
       posic => head                      !begin at top
       DO
         IF(posic%numel == kelem) THEN    !check if label found
           found = .TRUE.                 !Found
           EXIT                           !element found, EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN    !check if there are more els
           anter => posic                     !remember previous element
           posic => posic%next                !new present element
         ELSE
           EXIT                               !no more elements EXIT
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)    !point to nothing
     RETURN
   END SUBROUTINE srch_ele29e

   SUBROUTINE del_ele29e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele29), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele29), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%stint)              !deallocate stress array
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele29e

   SUBROUTINE cut_ele29e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele29), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele29e

!    INCLUDE 'actu29.fi'
   INCLUDE 'acvd29.fi'
   INCLUDE 'bmat29.fi'
   INCLUDE 'bmai29.fi'
   INCLUDE 'bmma29.fi'
   INCLUDE 'btma29.fi'
   INCLUDE 'bsma29.fi'
   INCLUDE 'comm29.fi'
!!   INCLUDE 'dump29.fi'
   INCLUDE 'elmd29.fi'
!!   INCLUDE 'expo29.fi'
   INCLUDE 'gaus29.fi'
!!   INCLUDE 'impo29.fi'
   INCLUDE 'jacob29.fi'
   INCLUDE 'kgmm29.fi'
   INCLUDE 'kgms29.fi'
!!   INCLUDE 'kgmt29.fi'
   INCLUDE 'lcsy29.fi'
   INCLUDE 'load29.fi'
   INCLUDE 'mase29.fi'
   INCLUDE 'masm29.fi'
   INCLUDE 'nodx29.fi'
   INCLUDE 'outp29.fi'
!!   INCLUDE 'rest29.fi'
   INCLUDE 'resv29.fi'
   INCLUDE 'stif29.fi'
   INCLUDE 'stra29.fi'

 END MODULE ele29_db
