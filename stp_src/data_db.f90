MODULE data_db

  ! global data base for post-process interface of SIMPACT/DELTA

  IMPLICIT NONE

  CHARACTER(len=20) :: text           !title

  INTEGER :: ndime, & !Dimension problem
             ndofn, & !Number of Dofs per node
             nrotd, & !Number of Rotation Dofs per node
             nadd,  & !Number of additional DOFs per node
             ndoft, & !Number of Temperature Dofs per node
             neulr, & !Code fot Euler angles
             npoin, & !number of points in the mesh
             narg , & !number of arguments in command line
             numsec,& !number of sections
             last_label, & !last label used for auxiliar points
             ip       !post-process tool

  INTEGER :: ntype, & ! problem type for 2-D
             idyna    ! dynamic or static ploblem


  INTEGER, PARAMETER :: maxs = 100  !maximum number of sets
  INTEGER :: setyp(maxs),  & !set type
             nsets           !number of sets

  LOGICAL :: fin    = .FALSE. , & !.TRUE. end of file detected
             first  = .TRUE.  , & !.TRUE. first time for coordinates
             addof  = .FALSE.     !.TRUE. if additional DOFs exist

  REAL (kind=8) :: ttime       !step time
  CHARACTER (len=10) :: loadstep


  ! flags associated to kinematics

  LOGICAL :: wtdisp  , & ! write nodal total displacements
             wsdisp  , & ! write nodal stage displacements
             widisp  , & ! write nodal initial stage displacements (-)
             wveloc  , & ! write nodal velocities
             wtempe  , & ! write nodal tempertures
             waccel  , & ! write nodal accelerations
             weuler  , & ! write nodal Euler angles
             wangve  , & ! write nodal angular velocites
             wangac  , & ! write nodal angular accelerations
             waddof      ! write nodal additional DOFs


  CHARACTER(len=15) :: tdisp_l(4) ,tdisp_u , & ! label nodal total displacements
                       sdisp_l(4) ,sdisp_u , & ! label nodal stage displacements
                       idisp_l(4) ,idisp_u , & ! label nodal initial stage displacements (-)
                       veloc_l(4) ,veloc_u , & ! label nodal velocities
                       tempe_l(2) ,tempe_u , & ! label nodal temperatures
                       accel_l(4) ,accel_u , & ! label nodal accelerations
                       euler_l(4) ,euler_u , & ! label nodal Euler angles
                       angve_l(4) ,angve_u , & ! label nodal angular velocites
                       angac_l(4) ,angac_u , & ! label nodal angular accelerations
                       addof_l(4) ,addof_u     ! label nodal additional DOFs

  REAL(kind=8)      :: tdisp_f , & ! factor nodal total displacements
                       sdisp_f , & ! factor nodal stage displacements
                       idisp_f , & ! factor nodal initial stage displacements (-)
                       veloc_f , & ! factor nodal velocities
                       tempe_f , & ! factor nodal temperatures
                       accel_f , & ! factor nodal accelerations
                       euler_f , & ! factor nodal Euler angles
                       angve_f , & ! factor nodal angular velocites
                       angac_f , & ! factor nodal angular accelerations
                       addof_f     ! factor nodal additional DOFs


  INTEGER (kind=4), ALLOCATABLE :: label(:)  !node labels
  LOGICAL, ALLOCATABLE :: meshn(:)  !nodes in the mesh

  REAL (kind=8), POINTER :: coord(:,:), & !(ndime,npoin) original coordinates
                            coors(:,:), & !(ndime,npoin) stage coordinates
                            coorf(:,:), & !(ndime,npoin) present coordinates
                            displ(:,:), & !(ndime,npoin) displacements
                            dispi(:,:), & !(ndime,npoin) initial displacements
                            veloc(:,:), & !(ndime,ndime) velocities
                            tempe(:,:), & !(ndoft,npoin) accelerations
                            accel(:,:), & !(ndime,npoin) accelerations
                            euler(:,:), & !(neulr,npoin) rotation matrix
                            anvel(:,:), & !(nrotd,npoin) ang velocities
                            anacc(:,:), & !(nrotd,npoin) ang accelerations
                            psi(:,:)      !(2,npoin) additional DOFs


  REAL (kind=8), POINTER :: coora(:,:), & !(ndime,ngp) coordinates of added nodes
                            dispa(:,:)
  INTEGER(kind=4), POINTER :: labea(:)    !(ngp) label of added nodes
  INTEGER(kind=4) :: ngp                  ! added nodes
! ****  control variables for TecPlot
  CHARACTER (len=512) :: t_vars              !string with all variables
  CHARACTER (len=256) :: sh_var              !string with shared variables (solution)
  CHARACTER (len=256) :: ps_var              !string with passive variables (solution)
  CHARACTER (len=256) :: lc_var              !string with variables location (solution)
  INTEGER :: t_p,  &  !length of variable t_vars
             t_nv, &  !number of variables to print
             t_ngv,&  !number of global variables
             t_nsv, t_npv, & !number of Shared and Passive variables
             t_ev(2,10) !number of elemental variables (1) nodal (2) cell-centered
  INTEGER(kind=4), POINTER :: sh_vara(:), & !array of shared variables
                              ps_vara(:), & !array of passive variables
                              lc_vara(:)    !array of variables location
  INTEGER(kind=4), PARAMETER :: precision=0    ! REAL(type=4)

! **** Names of Result Range Tables
  TYPE rrt
    CHARACTER (len=22) :: name
    INTEGER (kind=4) :: nval
    REAL(kind=8), ALLOCATABLE :: vals(:)
    CHARACTER (len=15), ALLOCATABLE :: lab(:)
  END TYPE rrt

  TYPE(rrt) :: rrt_fr, rrt_sf, rrt_db


! ****  DATA BASE for Each Element Type

!  Flag values : 0 do not print
!              : 1 print at nodal points
!              : 2 print at Gauss points
!              : 3 print at nodal and Gauss points
!

! ******* SPOT data base ******

  ! Derived type for a SPOT element
  TYPE spot
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, &  ! number of element
               matno, &  ! material associated to st
               ngaus, &  ! (1) number of Gauss points per element
               nstre, &  ! (3) number of variables to read in the Gauss point
               set,   &  ! set number (order)
               nnode     ! number of nodes per element
    INTEGER, POINTER :: lnods(:,:)  !(2/6,nelem)
    REAL (kind=8), POINTER :: elvar(:,:,:)  !(nvarg,ngaus=1,nelem)
    !   variables for 6-node fasteners
    INTEGER :: first_l           ! first label for the set
    REAL(kind=8), POINTER :: a_x(:,:,:)   !(3,spot_npoin,3) coordinates (original,stage,final)
    REAL(kind=8), POINTER :: lc(:,:,:)    !(3,2,nelem)         local coordinates (xita,eta,zeta)
    TYPE (spot), POINTER :: next
  END TYPE spot

  !  flags for SPOT type elements
  INTEGER :: spot_force , & !write axial force
             spot_shear , & !write shear force
             spot_eqpst     !write equivalent plastic strain

  !  labels for SPOT type elements (1): Nodal (2): Gauss points
  CHARACTER (len=15) :: spot_force_l(2) , spot_force_u , & !total force label
                        spot_shear_l(2) , spot_shear_u , & !Cauchy stress label
                        spot_eqpst_l(2) , spot_eqpst_u     !equivalent plastic strain label

  REAL (kind=8) ::      spot_force_f , &  !factor to change units
                        spot_shear_f , &
                        spot_eqpst_f

  !  variables and arrays for smoothed spot values
  INTEGER :: spot_sets , & ! number of spot sets
             spot_nelem, & ! number of spot elements
             spot_nvarg, & ! number of variables to print in Gauss points
             spot_nps      ! number of nodes connected to spot element

  INTEGER, POINTER ::  spot_nodes(:,:)  !(npoin) local number for spot nodes

  TYPE (spot), POINTER :: spot_head, spot_tail
  INTEGER :: spot_auxiliar_nodes = 0

! ******* TRUSS data base ******

  ! Derived type for a TRUSS element
  TYPE truss
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, &  ! number of element
               matno, &  ! material associated to st
               ngaus, &  ! (1) number of Gauss points per element
               nstre, &  ! (3) number of variables to read in the Gauss point
               set,   &  ! set number (order)
               nnode     ! number of nodes per element
    INTEGER, POINTER :: lnods(:,:)  !(2,nelem)
    REAL (kind=8), POINTER :: elvar(:,:,:)  !(nvarg,ngaus=1,nelem)
    TYPE (truss), POINTER :: next
  END TYPE truss

  !  flags for TRUSS type elements
  INTEGER :: truss_force , & !write total force
             truss_stres , & !write Cauchy stress
             truss_eqpst     !write equivalent plastic strain

  !  labels for TRUSS type elements (1): Nodal (2): Gauss points
  CHARACTER (len=15) :: truss_force_l(4) , truss_force_u , & !total force label
                        truss_stres_l(4) , truss_stres_u , & !Cauchy stress label
                        truss_eqpst_l(4) , truss_eqpst_u     !equivalent plastic strain label

  REAL (kind=8)  ::     truss_force_f , & !total force factor
                        truss_stres_f , & !Cauchy stress factor
                        truss_eqpst_f     !equivalent plastic strain factor

  !  variables and arrays for smoothed truss values
  INTEGER :: truss_sets , & ! number of truss sets
             truss_nelem, & ! number of truss elements
             truss_nvarg, & ! number of variables to print in Gauss points
             truss_nvarn, & ! number of variables to print at nodal points
             truss_nps      ! number of nodes connected to truss element

  INTEGER, POINTER :: truss_nodes(:,:)  !(npoin) local number for truss nodes

  REAL(kind=8), POINTER :: truss_accpn(:),   & !(nps) smoothing factor
                           truss_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (truss), POINTER :: truss_head, truss_tail


! ******* 2-D SOLID data base ******

  ! Derived type for a 2-D SOLID element
  TYPE sol2d
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of 2-d elements
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre, & ! number of stress variables
               ver      ! Gauss Point position code
    INTEGER, POINTER :: lnods(:,:), & !(nnode,nelem)       Connectivities
                        matno(:)      !(nelem) associated material
    REAL(kind=8), POINTER :: elvar(:,:,:)    !(nvarg,ngaus,nelem) Gauss point variables
    REAL(kind=8), POINTER :: dirt(:,:)       !(ndime,nelem) thickness direction
    TYPE (sol2d), POINTER :: next
  END TYPE sol2d

  !  flags for 2-D SOLID type elements
  INTEGER :: sol2d_stres  , & !write Cauchy stresses
             sol2d_logst  , & !write logarithmic principal strains
             sol2d_shtst  , & !write sheet strains
             sol2d_thrat  , & !write thickness ratio
             sol2d_eqpst  , & !write effective platic strains
             sol2d_vmise  , & !write von Mises stress
             sol2d_fldma      !write FLD MAP
  LOGICAL :: sol2d_wfldFZ , & !write Forming Zone Diagram
             sol2d_wfldSZ     !write Safety Zone Diagram

  REAL(kind=8) :: sol2d_fldva = 0.0   !default value for FLD MAP

  !  Labels for 2-D SOLID type elements
  CHARACTER(len=15) :: sol2d_stres_l(10)  ,  sol2d_stres_u  , & !Cauchy stresses
                       sol2d_logst_l( 8)  ,  sol2d_logst_u  , & !logarithmic principal strains
                       sol2d_shtst_l( 8)  ,  sol2d_shtst_u  , & !sheet principal strains
                       sol2d_thrat_l( 4)  ,  sol2d_thrat_u  , & !thickness ratio
                       sol2d_eqpst_l( 4)  ,  sol2d_eqpst_u  , & !effective platic strains
                       sol2d_vmise_l( 4)  ,  sol2d_vmise_u  , & !von Mises stress
                       sol2d_fldma_l( 4)  ,                   & !FLD MAP
                       sol2d_fldFZ_l( 1)  ,                   & !Forming zone
                       sol2d_fldSZ_l( 1)                        !Safety zone

   REAL (kind=8) ::     sol2d_stres_f  , & !Cauchy stresses
                        sol2d_logst_f  , & !logarithmic principal strains
                        sol2d_shtst_f  , & !sheet principal strains
                        sol2d_thrat_f  , & !thickness ratio
                        sol2d_eqpst_f  , & !effective platic strains
                        sol2d_vmise_f      !von Mises stress

  !  variables and arrays for smoothed sol2d values
  INTEGER :: sol2d_sets , & ! number of sol2d sets
             sol2d_nelem, & ! number of sol2d element
             sol2d_nvarg, & ! number of variables to print in Gauss points
             sol2d_nvarn, & ! number of variables to print at nodal points
             sol2d_nps      ! number of nodes connected to sol2d element

  INTEGER, POINTER :: sol2d_nodes(:,:)  !(npoin) local number for sol2d nodes
  INTEGER, POINTER :: sol2d_cpx(:,:)    !(npoin) local number for sol2d nodes

  REAL(kind=8), POINTER :: sol2d_accpn(:),   & !(nps) smoothing factor
                           sol2d_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (sol2d), POINTER :: sol2d_head, sol2d_tail

! ******* 3-D SOLID data base ******

  ! Derived type for a 3-D SOLID element
  TYPE sol3d
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre, & ! number of stress variables
               etype    ! element type 4-5-16-18-27
    INTEGER, POINTER :: lnods(:,:)  !(nnode,nelem)       Connectivities
    REAL(kind=8), POINTER :: elvar(:,:,:)    !(nvarg,ngaus,nelem) Gauss point variables
    REAL(kind=8), POINTER :: dirt(:,:,:)     !(ndime,ndime,nelem) thickness direction
    INTEGER(kind=4), POINTER :: matno(:)     !(nelem) associated material
    INTEGER :: first_l     ! first label for the set Etype = 27
    TYPE (sol3d), POINTER :: next
  END TYPE sol3d

  !  flags for 3-D SOLID type elements
  INTEGER :: sol3d_stres  , & !write Cauchy stresses
             sol3d_logst  , & !write logarithmic principal strains
             sol3d_thrat  , & !write thickness ratio
             sol3d_eqpst  , & !write effective platic strains
             sol3d_vmise  , & !write von Mises stress
             sol3d_fldma  , & !write FLD MAP
             sol3d_press      !write nodal pressure
 ! flag associated to gauss point positions in Solids
  LOGICAL :: given = .FALSE.
  LOGICAL :: type27gp = .FALSE.    ! in-plane gauss points .F. internal point, .T. mid-side
  LOGICAL :: gppos27  = .FALSE.    ! in-plane gauss points .F. internal point, .T. mid-side
  LOGICAL :: sol3d_prsmo = .FALSE. ! flag associated to smoothing type for 15-node prism

  REAL(kind=8) :: sol3d_fldva = 0.0            !default value for FLD MAP

  !  Labels for 3-D SOLID type elements
  CHARACTER(len=15):: sol3d_stres_l(14) ,  sol3d_stres_u , & ! Cauchy stresses
                      sol3d_logst_l( 8) ,  sol3d_logst_u , & ! logarithmic principal strains
                      sol3d_thrat_l( 4) ,  sol3d_thrat_u , & ! thickness ratio
                      sol3d_eqpst_l( 4) ,  sol3d_eqpst_u , & ! effective platic strains
                      sol3d_vmise_l( 4) ,  sol3d_vmise_u , & ! von Mises stress
                      sol3d_fldma_l( 4) ,                  & ! FLD MAP
                      sol3d_press_l( 2) ,  sol3d_press_u     ! Nodal Pressure

  REAL (kind=8)  ::    sol3d_stres_f, & ! Cauchy stresses
                       sol3d_logst_f, & ! logarithmic principal strains
                       sol3d_thrat_f, & ! thickness ratio
                       sol3d_eqpst_f, & ! effective platic strains
                       sol3d_vmise_f, & ! von Mises stress
                       sol3d_press_f    ! pressure

  !  variables and arrays for smoothed sol3d values
  INTEGER :: sol3d_sets , & ! number of sol3d sets
             sol3d_nelem, & ! number of sol3d element
             sol3d_nvarg, & ! number of variables to print in Gauss points
             sol3d_nvarn, & ! number of variables to print at nodal points
             sol3d_nps,   & ! number of nodes connected to sol3d element
             sol3d_auxiliar_nodes  ! number of auxiliar nodes connected to sol3d element

  INTEGER, POINTER :: sol3d_nodes(:,:)  !(npoin) local number for sol3d nodes

  REAL(kind=8), POINTER :: sol3d_accpn(:),   & !(nps) smoothing factor
                           sol3d_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (sol3d), POINTER :: sol3d_head, sol3d_tail


! ******* 3-D SHELL (Simo theory elements) data base ******

  ! Derived type for a 3-D SHELL element
  TYPE shl3d
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre    ! number of variables per Gauss point
    INTEGER, POINTER :: lnods(:,:), &  !(nnode,nelem)       Connectivities
                        matno(:)       ! nelem  associated material
    REAL(kind=8), POINTER :: elvar(:,:,:)   !(nvarg,ngaus,nelem) Gauss point variables
    TYPE (shl3d), POINTER :: next
  END TYPE shl3d

  !  flags for 3-D SHELL type elements
  INTEGER :: shl3d_force  , & !write integrated stresses
             shl3d_momen  , & !write integrated moments
             shl3d_shear  , & !write integrated shears
             shl3d_logst  , & !write logarithmic principal strains
             shl3d_eqpst  , & !write effective platic strains
             shl3d_vmise  , & !write von Mises stress
             shl3d_thrat      !write thickness ratio

  !  labels for 3-D SHELL type elements
  CHARACTER(len=15) :: shl3d_force_l(8)  , shl3d_force_u  , & ! integrated stresses
                       shl3d_momen_l(8)  , shl3d_momen_u  , & ! integrated moments
                       shl3d_shear_l(6)  , shl3d_shear_u  , & ! integrated shears
                       shl3d_logst_l(8)  , shl3d_logst_u  , & ! logarithmic principal strains
                       shl3d_eqpst_l(4)  , shl3d_eqpst_u  , & ! effective platic strains
                       shl3d_vmise_l(4)  , shl3d_vmise_u  , & ! von Mises stress
                       shl3d_thrat_l(4)  , shl3d_thrat_u      ! thickness ratio

   REAL (kind=8) ::    shl3d_force_f  , & ! integrated stresses
                       shl3d_momen_f  , & ! integrated moments
                       shl3d_shear_f  , & ! integrated shears
                       shl3d_logst_f  , & ! logarithmic principal strains
                       shl3d_eqpst_f  , & ! effective platic strains
                       shl3d_vmise_f  , & ! von Mises stress
                       shl3d_thrat_f      ! thickness ratio

  !  variables and arrays for smoothed shl3d values
  INTEGER :: shl3d_sets , & ! number of shl3d sets
             shl3d_nelem, & ! number of shl3d element
             shl3d_nvarg, & ! number of variables to print in Gauss points
             shl3d_nvarn, & ! number of variables to print at nodal points
             shl3d_nps      ! number of nodes connected to shl3d element

  INTEGER, POINTER :: shl3d_nodes(:,:)  !(npoin) local number for shl3d nodes

  REAL(kind=8), POINTER :: shl3d_accpn(:),   & !(nps) smoothing factor
                           shl3d_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (shl3d), POINTER :: shl3d_head, shl3d_tail

! ******* 3-D BEAM (Simo theory elements) data base ******

  ! Derived type for a 3-D BEAM element
  TYPE beame
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               matno, & ! associated material
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre    ! number of variables per Gauss point
    INTEGER, POINTER :: lnods(:,:)  !(nnode,nelem)       Connectivities
    REAL(kind=8), POINTER :: elvar(:,:,:)    !(nvarg,ngaus,nelem) Gauss point variables
    TYPE (beame), POINTER :: next
  END TYPE beame

  !  flags for 3-D BEAM type elements
  INTEGER :: beame_force  , & !write integrated stresses (normal and shears)
             beame_momen  , & !write integrated moments (flector and torsor)
             beame_eqpst  , & !write effective platic strains
             beame_vmise  , & !write von Mises stress
             beame_arrat      !write area ratio

  !  labels for 3-D BEAM type elements
  CHARACTER(Len=15) :: beame_force_l(8) , beame_force_u , & ! integrated stresses (normal and shears)
                       beame_momen_l(8) , beame_momen_u , & ! integrated moments (flector and torsor)
                       beame_eqpst_l(4) , beame_eqpst_u , & ! effective platic strains
                       beame_vmise_l(4) , beame_vmise_u , & ! von Mises stress
                       beame_arrat_l(4) , beame_arrat_u     ! area ratio

  REAL (kind=8) ::     beame_force_f , & ! integrated stresses (normal and shears)
                       beame_momen_f , & ! integrated moments (flector and torsor)
                       beame_eqpst_f , & ! effective platic strains
                       beame_vmise_f , & ! von Mises stress
                       beame_arrat_f     ! area ratio

  !  variables and arrays for smoothed beame values
  INTEGER :: beame_sets , & ! number of beame sets
             beame_nelem, & ! number of beame element
             beame_nvarg, & ! number of variables to print in Gauss points
             beame_nvarn, & ! number of variables to print at nodal points
             beame_nps      ! number of nodes connected to beame element

  INTEGER, POINTER :: beame_nodes(:,:)  !(npoin) local number for beame nodes

  REAL(kind=8), POINTER :: beame_accpn(:),   & !(nps) smoothing factor
                           beame_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (beame), POINTER :: beame_head, beame_tail


! ******* 2-D SHELL (Simo theory elements) data base ******

  ! Derived type for a 2-D SHELL element
  TYPE shrev
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               matno, & ! associated material
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre    ! number of variables per Gauss point
    INTEGER, POINTER :: lnods(:,:)  !(nnode,nelem)       Connectivities
    REAL(kind=8), POINTER :: elvar(:,:,:)   !(nvarg,ngaus,nelem) Gauss point variables
    TYPE (shrev), POINTER :: next
  END TYPE shrev

  !  flags for 2-D SHELL type elements
  INTEGER :: shrev_force  , & !write integrated stresses
             shrev_shear  , & !write integrated transverse shear stresses
             shrev_momen  , & !write integrated moments
             shrev_eqpst  , & !write effective platic strains
             shrev_vmise  , & !write von Mises stress
             shrev_thrat      !write thickness ratio

  !  Labels for 2-D SHELL type elements
  CHARACTER(len=15) :: shrev_force_l(6) , shrev_force_u , & ! integrated stresses
                       shrev_shear_l(4) , shrev_shear_u , & ! integrated transverse shear stresses
                       shrev_momen_l(6) , shrev_momen_u , & ! integrated moments
                       shrev_eqpst_l(4) , shrev_eqpst_u , & ! effective platic strains
                       shrev_vmise_l(4) , shrev_vmise_u , & ! von Mises stress
                       shrev_thrat_l(4) , shrev_thrat_u     ! thickness ratio

   REAL(kind=8) ::     shrev_force_f , & ! integrated stresses
                       shrev_shear_f , & ! integrated transverse shear stresses
                       shrev_momen_f , & ! integrated moments
                       shrev_eqpst_f , & ! effective platic strains
                       shrev_vmise_f , & ! von Mises stress
                       shrev_thrat_f     ! thickness ratio

  !  variables and arrays for smoothed shrev values
  INTEGER :: shrev_sets , & ! number of shrev sets
             shrev_nelem, & ! number of shrev element
             shrev_nvarg, & ! number of variables to print in Gauss points
             shrev_nvarn, & ! number of variables to print at nodal points
             shrev_nps      ! number of nodes connected to shrev element

  INTEGER, POINTER :: shrev_nodes(:,:)  !(npoin) local number for shrev nodes

  REAL(kind=8), POINTER :: shrev_accpn(:),   & !(nps) smoothing factor
                           shrev_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (shrev), POINTER :: shrev_head, shrev_tail


! ******* BST SHELL data base ******

  ! Derived type for a BST SHELL element
  TYPE bst
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               set  , & ! set number
               ngaus, & ! number of Gauss points per element
               nnode, & ! number of nodes per element
               nstre, & ! number of variables to read per Gauss point
               locax    !local system option
    INTEGER, POINTER :: lnods(:,:), &  !(nnode,nelem)       Connectivities
                        matno(:)       !nelem  associated material
    REAL(kind=8), POINTER :: elvar(:,:,:), & !(nvarg,ngaus,nelem) Gauss point variables
                             angle(:)        !(nelem) Local angle for stress definition

    TYPE (bst), POINTER :: next      !pointer to next set
  END TYPE bst

  !  flags for BST SHELL type elements
  INTEGER :: bst_force  , & !write Cauchy integrated stresses
             bst_momen  , & !write Cauchy integrated stress moments
             bst_shear  , & !write Cauchy integrated transverse stresses
             bst_logst  , & !write logarithmic principal strains
             bst_curva  , & !write principal curvatures
             bst_thrat  , & !write thickness ratio
             bst_eqpst  , & !write effective platic strains
             bst_vmise  , & !write von mises equivalent stress
             bst_fldma      !write FLD MAP
  LOGICAL :: bst_wfldFZ , & !write Forming Zone diagram
             bst_wfldSZ     !write Safety Zone diagram

  REAL(kind=8) :: bst_fldva = 0.0            !default value for FLD MAP

  !  Labels for BST SHELL type elements
  CHARACTER(len=15) :: bst_force_l(8) , bst_force_u , & ! Cauchy integrated stresses
                       bst_momen_l(8) , bst_momen_u , & ! Cauchy integrated stress moments
                       bst_shear_l(6) , bst_shear_u , & ! Cauchy integrated transverse stresses
                       bst_logst_l(8) , bst_logst_u , & ! logarithmic principal strains
                       bst_curva_l(8) , bst_curva_u , & ! principal curvatures
                       bst_thrat_l(4) , bst_thrat_u , & ! thickness ratio
                       bst_eqpst_l(4) , bst_eqpst_u , & ! effective platic strains
                       bst_vmise_l(4) , bst_vmise_u , & ! von mises equivalent stress
                       bst_fldma_l(4) ,               & ! FLD MAP
                       bst_fldFZ_l(1) ,               & ! Forming Zone
                       bst_fldSZ_l(1)                   ! Safety Zone

   REAL(kind=8) :: bst_force_f , & ! Cauchy integrated stresses
                   bst_momen_f , & ! Cauchy integrated stress moments
                   bst_shear_f , & ! Cauchy integrated transverse stresses
                   bst_logst_f , & ! logarithmic principal strains
                   bst_curva_f , & ! principal curvatures
                   bst_thrat_f , & ! thickness ratio
                   bst_eqpst_f , & ! effective platic strains
                   bst_vmise_f     ! von mises equivalent stress

  !  variables and arrays for smoothed bst values
  INTEGER :: bst_sets , & ! number of bst sets
             bst_nelem, & ! number of bst element
             bst_nvarg, & ! number of variables to print in Gauss points
             bst_nvarn, & ! number of variables to print at nodal points
             bst_nps      ! number of nodes connected to bst element

  INTEGER, POINTER ::   bst_nodes(:,:)  !(npoin) local number for bst nodes

  REAL(kind=8), POINTER :: bst_accpn(:),   & !(nps) smoothing factor
                           bst_vargs(:,:)    !(nvarn,nps) smoothed variables

  TYPE (bst), POINTER :: bst_head, bst_tail



! ******* RIGID SOLIDS data base ******

  ! Derived type for a 2/3-D RIGID SOLID element
  TYPE RIGID
    CHARACTER (len=30) :: sname
    INTEGER :: nelem, & ! number of elements
               matno, & ! associated material or type
               set,   & ! set number
               nnode, & ! number of nodes per element
               ntype, & ! 1: nodes, 2:surface 3:solid
               nmast    ! /= 0 :master node

    INTEGER, POINTER :: lnods(:,:)  !(nnode,nelem)       Connectivities
    TYPE (RIGID), POINTER :: next
  END TYPE RIGID

  LOGICAL :: rigid_show     ! write
  INTEGER :: rigid_sets,  & ! number of sets
             rigid_nps      ! number of nodes connected to Rigid elements

  INTEGER, POINTER :: rigid_nodes(:,:)  !(npoin) local number for rigid nodes

  TYPE (RIGID), POINTER :: rigid_head, rigid_tail

! ******* Draw Beads data base ******

  ! Derived type for a 3-D DRAW BEADS affected sheets
  TYPE DRAWB
    CHARACTER (len=30) :: pname  !drawbead name
    CHARACTER (len=30) :: sname  !drawbead sheet
    INTEGER :: nelem, & ! number of elements
               nnode, & ! number of nodes per element
               icode    ! if sheet is associated with an shell element set
    LOGICAL :: eall

    INTEGER, POINTER :: lnods(:,:)  !(nnode,nelem)       Connectivities
    LOGICAL, POINTER :: daz(:)  !(nelem)       values
    INTEGER, POINTER :: ass(:)  !(nelem)       relation
    TYPE (DRAWB), POINTER :: next
  END TYPE DRAWB

  LOGICAL :: drawb_wafz     ! write drawbeads affected zones
  INTEGER :: drawb_sets,  & ! number of sets
             drawb_nps,   & ! number of nodes connected to drawbead sheets
             drawb_shs,   & ! number of different sheets (associated to sets)
             drawb_shn      ! number different sheets (no associated to sets)
  INTEGER, ALLOCATABLE :: drawb_shnel(:) !number of elements in the associated sheets

  TYPE (DRAWB), POINTER :: drawb_head, drawb_tail

! ******* Internal variables of user defined materials

   ! Derived type for post-process variables

  TYPE udmat
     INTEGER (kind=4) :: matno         !associated section (material)
     INTEGER (kind=4) :: nvar          !number of variables
     INTEGER (kind=4) :: nvarv         !number of total scalars
     INTEGER (kind=4), ALLOCATABLE :: type(:)       !type 0:scalar 1:vector 2:tensor
     INTEGER (kind=4), ALLOCATABLE :: dim(:)        !space dimension 1 2 3
     CHARACTER (len=15), ALLOCATABLE :: name(:,:)   !variable name and components
     TYPE(udmat), POINTER :: next      !pointer to next
  END TYPE udmat

  TYPE (udmat), POINTER :: udmat_head, udmat_tail

 CONTAINS

 SUBROUTINE new_post_data(new,n)
 IMPLICIT NONE
 TYPE (udmat), POINTER :: new
 INTEGER (kind=4) :: n
 ALLOCATE (new)
 ALLOCATE (new%type(n), new%dim(n), new%name(7,n))
 new%nvar = n
 new%name = ''
 IF( ASSOCIATED( udmat_head ) )THEN
   udmat_tail%next => new
 ELSE
   udmat_head => new
 END IF
 udmat_tail => new
 RETURN
 END SUBROUTINE new_post_data

END MODULE data_db
