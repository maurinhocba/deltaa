 SUBROUTINE rdconf ()

 ! read configuration file
 !   TASK missing, read configuration data for Additional DOFs

 USE data_db
 USE cont_db
 IMPLICIT NONE

 CHARACTER (len= 1) :: c,flags(2)
 CHARACTER (len=99) :: lin
 CHARACTER (len=5)  :: ctype,var
 CHARACTER (len=15) :: names(2),cname(6,2),units
 CHARACTER (len=256):: exepath
 INTEGER (kind=4) :: pos

 REAL(kind=8) :: rval

 INTEGER :: set_value !(FUNCTION)
 INTEGER :: ios

 INTERFACE
   SUBROUTINE readln (c,lin,word,value,flags,names,cname,units)

   IMPLICIT NONE

   !***  arguments
   CHARACTER(len=1 ) ,INTENT(IN OUT) :: c
   CHARACTER(len=99) ,INTENT(IN OUT) :: lin
   CHARACTER(len=5 ) ,INTENT(OUT) :: word
   CHARACTER(len=1 ) ,INTENT(OUT) :: flags(2)
   CHARACTER(len=15) ,INTENT(OUT) :: names(2)
   CHARACTER(len=15) ,INTENT(OUT) :: cname(6,2),units
   REAL(kind=8) ,INTENT(OUT) :: value

   END SUBROUTINE readln

   SUBROUTINE as_names(nd,nm,flags,names,cname,fd,nn,nc,old_l,rval,units,var_f,var_u)

   IMPLICIT NONE
   INTEGER (kind=4), INTENT(IN) :: nd, & !number of data variables (READ)
                                   nm, & !number of components (READ)
                                   nn, & !number of variables to be included
                                   nc    !number of components to be includes
   CHARACTER (len=1), INTENT(IN) :: flags(nd), & !flags read
                                    fd(nn)       !flags to be considered
   CHARACTER (len=15), INTENT(IN) :: names(nd), & !variable names (READ)
                                     cname(nm,nd) !variable components (READ)
   CHARACTER (len=15), INTENT(IN OUT) :: old_l(nn*(nc+1)) !variables and components names
   REAL(kind=8), OPTIONAL :: rval, var_f        !factor to change units
   CHARACTER (len=15), OPTIONAL :: units,var_u  !units
   END SUBROUTINE as_names

 END INTERFACE

 ! *************** set DEFAULT values ****************

  ! post-process tool

  ip = 4    !GiD binary

  !    flags associated to kinematics

  wtdisp = .TRUE.  ! write total displacements
  wsdisp = .FALSE. ! write stage displacements
  widisp = .FALSE. ! write initial stage displacements
  wveloc = .FALSE. ! write nodal velocities
  wtempe = .TRUE.  ! write nodal temperatures
  waccel = .FALSE. ! write nodal accelerations
  weuler = .FALSE. ! write nodal Euler angles
  wangve = .FALSE. ! write nodal angular velocites
  wangac = .FALSE. ! write nodal angular accelerations
  waddof = .TRUE.  ! write nodal additional DOFs

  tdisp_l = (/ 'Total_Disp',    'X_TDisp','Y_TDisp','Z_Tdisp' /) ; tdisp_u = ''  ; tdisp_f = 1D0
  sdisp_l = (/ 'Stage_Disp',    'X_SDisp','Y_SDisp','Z_SDisp' /) ; sdisp_u = ''  ; sdisp_f = 1D0
  idisp_l = (/ 'Initial   ',    'X_IDisp','Y_IDisp','Z_IDisp' /) ; idisp_u = ''  ; idisp_f = 1D0
  veloc_l = (/ 'Velocities',    'X_Vel  ','Y_Vel  ','Z_Vel'   /) ; veloc_u = ''  ; veloc_f = 1D0
  tempe_l = (/ 'NodalTemper',   'Temper'/)                       ; tempe_u = ''  ; tempe_f = 1D0
  accel_l = (/ 'Accelerations', 'X_Accl ','Y_Accl ','Z_Accl'  /) ; accel_u = ''  ; accel_f = 1D0
  euler_l = (/ 'Local_Axes',    'Alpha  ','Beta   ','Gamma'   /) ; euler_u = ''  ; euler_f = 1D0
  angve_l = (/ 'Angular_Vel',   'X_Omega','Y_Omega','Z_Omega' /) ; angve_u = ''  ; angve_f = 1D0
  angac_l = (/ 'Angular_Accel', 'X_AA   ','Y_AA   ','Z_AA'    /) ; angac_u = ''  ; angac_f = 1D0
  addof_l = (/ 'Add_Displac.',  'Psi_x  ','Psi_y  ','    '    /) ; addof_u = ''  ; addof_f = 1D0

  !    flags for SPOT type elements
  spot_force =  2   !write axial force
  spot_shear =  1   !write shear force
  spot_eqpst =  0   !write equivalent plastic strain

  spot_force_l = (/'Spot_Force', 'Force'    /) ; spot_force_u = ''  ; spot_force_f = 1d0
  spot_shear_l = (/'Spot_Shear', 'Shear'    /) ; spot_shear_u = ''  ; spot_shear_f = 1d0
  spot_eqpst_l = (/'Spot_EqPlSt','Eq_Pl_St' /) ; spot_eqpst_u = ''  ; spot_eqpst_f = 1d0

  !    flags for TRUSS type elements
  truss_force =  1   !write total force
  truss_stres =  0   !write Cauchy stress
  truss_eqpst =  1   !write equivalent plastic strain

  truss_force_l = (/' N_Tr_Force' , 'Force'          , 'Truss_Force' , 'Force'    /)
  truss_stres_l = (/' N_Tr_Stress', 'Stress'         , 'Truss_Stress', 'Stress'   /)
  truss_eqpst_l = (/' N_Tr_EqPlSt', 'Eq. Plastic St.', 'Truss_EqPlSt', 'Eq_Pl_St' /)

  truss_force_u = ''  ;  truss_force_f = 1d0
  truss_stres_u = ''  ;  truss_stres_f = 1d0
  truss_eqpst_u = ''  ;  truss_eqpst_f = 1d0

   !  flags for 2-D SOLID type elements
  sol2d_stres =  0   !write Cauchy stresses
  sol2d_logst =  0   !write logarithmic principal strains
  sol2d_shtst =  0   !write shett principal strains
  sol2d_thrat =  0   !write thickness ratio
  sol2d_eqpst =  1   !write effective platic strains
  sol2d_vmise =  1   !write von Mises stress
  sol2d_fldma =  0   !write FLD MAP

  sol2d_fldva = 0.0  !default value for FLD MAP

  sol2d_stres_l = (/ 'N_Stress'  ,'S_xx','S_yy','S_xy','S_zz', &
                &    'G_Stress'  ,'S_xx','S_yy','S_xy','S_zz' /)
  sol2d_logst_l = (/ 'N_Log_St'  ,'E1','E2','E3(circ)',        &
                &    'G_Log_St'  ,'E1','E2','E3(circ)'        /)
  sol2d_shtst_l = (/ 'N_Sht_St'  ,'E1','E2','E3(th)',        &
                &    'G_Sht_St'  ,'E1','E2','E3(th)'        /)
  sol2d_thrat_l = (/ 'N_Th_Rat'  ,'Thickness Ratio','G_Th_Rat'  ,'Thickness Ratio' /)
  sol2d_eqpst_l = (/ 'N_Eq_Pl_St','Eq. Plastic St.','G_Eq_Pl_St','Eq. Plastic St.' /)
  sol2d_vmise_l = (/ 'N_Eq_VM_St','V Mises Stress' ,'G_Eq_VM_St','V Mises Stress'  /)
  sol2d_fldma_l = (/ 'N_FLD_Map' ,'FLD_Map'        ,'G_FLD_Map' ,'FLD_Map'         /)
  sol2d_fldFZ_l = (/ 'FORMING ZONE'/)
  sol2d_fldSZ_l = (/ 'SAFETY ZONE'/)

  sol2d_stres_u = ''  ; sol2d_stres_f =   1d0
  sol2d_logst_u = ''  ; sol2d_logst_f =   1d0
  sol2d_shtst_u = ''  ; sol2d_shtst_f =   1d0
  sol2d_thrat_u = ''  ; sol2d_thrat_f =   1d0
  sol2d_eqpst_u = ''  ; sol2d_eqpst_f =   1d0
  sol2d_vmise_u = ''  ; sol2d_vmise_f =   1d0

  !  flags for 3-D SOLID type elements
  sol3d_stres =  0   !write Cauchy stresses
  sol3d_logst =  0   !write logarithmic principal strains
  sol3d_thrat =  0   !write thickness ratio
  sol3d_eqpst =  1   !write effective platic strains
  sol3d_vmise =  1   !write von Mises stress
  sol3d_fldma =  0   !write FLD MAP
  sol3d_prsmo = .FALSE. ! flag associated to smoothing type for 15-node prism
  sol3d_press =  0   !write nodal press

  sol3d_fldva = 0.0            !default value for FLD MAP

  sol3d_stres_l = (/ 'N_Stress'  ,'S_xx','S_yy','S_zz','S_xy','S_yz','S_xz', &
                &    'G_Stress'  ,'S_xx','S_yy','S_zz','S_xy','S_yz','S_xz' /)
  sol3d_logst_l = (/ 'N_Log_St'  ,'E1','E2','E3(circ)',        &
                &    'G_Log_St'  ,'E1','E2','E3(circ)'        /)
  sol3d_thrat_l = (/ 'N_Th_Rat'  ,'Thickness Ratio','G_Th_Rat'  ,'Thickness Ratio' /)
  sol3d_eqpst_l = (/ 'N_Eq_Pl_St','Eq. Plastic St.','G_Eq_Pl_St','Eq. Plastic St.' /)
  sol3d_vmise_l = (/ 'N_Eq_VM_St','V Mises Stress' ,'G_Eq_VM_St','V Mises Stress'  /)
  sol3d_fldma_l = (/ 'N_FLD_Map' ,'FLD_Map'        ,'G_FLD_Map' ,'FLD_Map'         /)
  sol3d_press_l = (/ 'N_Pressure','Pressure' /)

  sol3d_stres_u = '' ;  sol3d_stres_f = 1D0
  sol3d_logst_u = '' ;  sol3d_logst_f = 1D0
  sol3d_thrat_u = '' ;  sol3d_thrat_f = 1D0
  sol3d_eqpst_u = '' ;  sol3d_eqpst_f = 1D0
  sol3d_vmise_u = '' ;  sol3d_vmise_f = 1D0
  sol3d_press_u = '' ;  sol3d_press_f = 1D0

  !  flags for 3-D SHELL type elements
  shl3d_force =  1   !write integrated stresses
  shl3d_momen =  1   !write integrated moments
  shl3d_shear =  0   !write integrated transverse shear stresses
  shl3d_logst =  0   !write log strains
  shl3d_eqpst =  0   !write effective platic strains
  shl3d_vmise =  0   !write von Mises stress
  shl3d_thrat =  0   !write thickness ratio

  shl3d_force_l = (/ 'N_Sh_Forces' ,'N_xx','N_yy','N_xy', &
                &    'G_Sh_Forces' ,'N_xx','N_yy','N_xy' /)
  shl3d_momen_l = (/ 'N_Sh_Moments','M_xx','M_yy','M_xy', &
                &    'G_Sh_Moments','M_xx','M_yy','M_xy' /)
  shl3d_shear_l = (/ 'N_Sh_Shears' ,'NQ_x','NQ_y',          &
                &    'G_Sh_Shears' ,'GQ_x','GQ_y' /)
  shl3d_logst_l = (/ 'N_Sh_Log_St' ,'E1','E2','E3(thick)',&
                &    'G_Sh_Log_St' ,'E1','E2','E3(thick)'/)
  shl3d_eqpst_l = (/ 'N_Sh_EPS','Eq. Plastic St.','G_Sh_EPS','Eq. Plastic St.' /)
  shl3d_vmise_l = (/ 'N_Sh_VMS','V Mises Stress' ,'G_Sh_VMS','V Mises Stress'  /)
  shl3d_thrat_l = (/ 'N_Sh_ThR','Thickness Ratio','G_Sh_ThR','Thickness Ratio' /)

  shl3d_momen_u = ''   ;  shl3d_momen_f = 1D0
  shl3d_shear_u = ''   ;  shl3d_shear_f = 1D0
  shl3d_logst_u = ''   ;  shl3d_logst_f = 1D0
  shl3d_eqpst_u = ''   ;  shl3d_eqpst_f = 1D0
  shl3d_vmise_u = ''   ;  shl3d_vmise_f = 1D0
  shl3d_thrat_u = ''   ;  shl3d_thrat_f = 1D0

  !  flags for 3-D BEAM (Simo theory) type elements
  beame_force =  1   !write integrated stresses
  beame_momen =  1   !write integrated moments
  beame_eqpst =  0   !write effective platic strains
  beame_vmise =  0   !write von Mises stress
  beame_arrat =  0   !write thickness ratio

  beame_force_l = (/ 'N_Beam_Forces' ,'Axial' ,'Shear_Y'  ,'Shear_Z'  , &
                &    'G_Beam_Forces' ,'Axial' ,'Shear_Y'  ,'Shear_Z'    /)
  beame_momen_l = (/ 'N_Beam_Moments','Torsor','Bending_Y','Bending_Z', &
                &    'G_Beam_Moments','Torsor','Bending_Y','Bending_Z'  /)
  beame_eqpst_l = (/ 'N_Beam_EPS','Eq. Plastic St.','G_Beam_EPS','Eq. Plastic St.' /)
  beame_vmise_l = (/ 'N_Beam_VMS','V Mises Stress' ,'G_Beam_VMS','V Mises Stress'  /)
  beame_arrat_l = (/ 'N_Beam_ArR','Area Ratio'     ,'G_Beam_ArR','Area Ratio'      /)

  beame_force_u = '' ;  beame_force_f = 1D0
  beame_momen_u = '' ;  beame_momen_f = 1D0
  beame_eqpst_u = '' ;  beame_eqpst_f = 1D0
  beame_vmise_u = '' ;  beame_vmise_f = 1D0
  beame_arrat_u = '' ;  beame_arrat_f = 1D0


  !  flags for 2-D SHELL (Simo theory) type elements
  shrev_force =  1   !write integrated stresses
  shrev_shear =  0   !write integrated shear stresses
  shrev_momen =  1   !write integrated moments
  shrev_eqpst =  0   !write effective platic strains
  shrev_vmise =  0   !write von Mises stress
  shrev_thrat =  0   !write thickness ratio

  shrev_force_l = (/ 'N_Sh_Force','N11','N22',  &
                &    'G_Sh_Force','N11','N22'   /)
  shrev_shear_l = (/ 'N_Sh_Shear','N_Sh_Shear','G_Sh_Shear','G_Sh_Shear'  /)
  shrev_momen_l = (/ 'N_Sh_Moment','M11','M22', &
                &    'G_Sh_Moment','M11','M22'  /)
  shrev_eqpst_l = (/ 'N_Sh_EPS','Eq. Plastic St.', 'G_Sh_EPS','Eq. Plastic St.' /)
  shrev_vmise_l = (/ 'N_Sh_VMS','V Mises Stress' , 'G_Sh_VMS','V Mises Stress'  /)
  shrev_thrat_l = (/ 'N_Sh_ThR','Thickness Ratio', 'G_Sh_ThR','Thickness Ratio' /)

  shrev_force_l = ''  ;  shrev_force_f =  1d0
  shrev_shear_l = ''  ;  shrev_shear_f =  1d0
  shrev_momen_l = ''  ;  shrev_momen_f =  1d0
  shrev_eqpst_l = ''  ;  shrev_eqpst_f =  1d0
  shrev_vmise_l = ''  ;  shrev_vmise_f =  1d0
  shrev_thrat_l = ''  ;  shrev_thrat_f =  1d0

   !  flags for BST SHELL type elements
  bst_force =  0      !write Cauchy integrated stresses
  bst_momen =  0      !write Cauchy integrated stress moments
  bst_shear =  0      !write integrated transverse shear stresses
  bst_logst =  1      !write logarithmic principal strains
  bst_curva =  0      !write principal curvatures
  bst_thrat =  1      !write thickness ratio
  bst_vmise =  1      !write von mises equivalent stress
  bst_eqpst =  1      !write effective platic strains
  bst_fldma =  0      !write FLD MAP
  bst_wfldFZ = .FALSE. !write FLD diagrams
  bst_wfldSZ = .FALSE. !write FLD diagrams

  bst_fldva = 0.0            !default value for FLD MAP

  bst_force_l = (/ 'N_BST_Forces', 'N_xx','N_yy','N_xy',  &
              &    'G_BST_Forces', 'N_xx','N_yy','N_xy'   /)
  bst_momen_l = (/ 'N_BST_Moments','M_xx','M_yy','M_xy',  &
              &    'G_BST_Moments','M_xx','M_yy','M_xy'   /)
  bst_shear_l = (/ 'N_BST_Shears' ,'Q_x','Q_y',           &
              &    'G_BST_Shears' ,'Q_x','Q_y' /)
  bst_logst_l = (/ 'N_BST_Log_St', 'E1','E2','E3(thick)', &
              &    'G_BST_Log_St', 'E1','E2','E3(thick)'  /)
  bst_curva_l = (/ 'N_BST_Curvat', 'C11','C22','C12',     &
              &    'G_BST_Curvat', 'C11','C22','C12'      /)
  bst_thrat_l = (/ 'N_BST_Th_Rat','Thickness Ratio', 'G_BST_Th_Rat','Thickness Ratio' /)
  bst_eqpst_l = (/ 'N_BST_EPS'   ,'Eq. Plastic St.', 'G_BST_EPS'   ,'Eq. Plastic St.' /)
  bst_vmise_l = (/ 'N_BST_VMS'   ,'V Mises Stress' , 'G_BST_VMS'   ,'V Mises Stress'  /)
  bst_fldma_l = (/ 'N_BST_FLD'   ,'FLD Map'        , 'G_BST_FLD'   ,'FLD Map'         /)
  bst_fldFZ_l = (/ 'FORMING ZONE'/)
  bst_fldSZ_l = (/ 'SAFETY ZONE'/)

  bst_force_u = ''   ;   bst_force_f = 1d0
  bst_momen_u = ''   ;   bst_momen_f = 1d0
  bst_shear_u = ''   ;   bst_shear_f = 1d0
  bst_logst_u = ''   ;   bst_logst_f = 1d0
  bst_curva_u = ''   ;   bst_curva_f = 1d0
  bst_thrat_u = ''   ;   bst_thrat_f = 1d0
  bst_eqpst_u = ''   ;   bst_eqpst_f = 1d0
  bst_vmise_u = ''   ;   bst_vmise_f = 1d0

  ! Flags for RIGID SOLIDS

  rigid_show = .TRUE.

  ! Flags for Drawbead affected zones

  drawb_wafz = .TRUE.


 ! ******* CONTACT data base ******

  wpress  = .FALSE.  !write binder pressure
  wwrink  = .FALSE.  !write gap information
  wwear   = .FALSE.  !write wearing information

  wear_l  = (/ 'Wearing', 'Blank_Marking','Tool Wearing' /)   ; wear_u  = ''  ; wear_f  = 1d0
  wrink_l = (/ 'Gaps','Actual','Minimum','Difference' /)      ; wrink_u = ''  ; wrink_f = 1d0
  press_l = (/ 'C_Press', 'Binder_Press','Contact_Press' /)   ; press_u = ''  ; press_f = 1d0

 ! generate Result Range Tables

 ! forming
 rrt_fr%name = 'FORMING'
 rrt_fr%nval = 7
 ALLOCATE(rrt_fr%vals(rrt_fr%nval+1),rrt_fr%lab(rrt_fr%nval))
 rrt_fr%vals = (/ 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5D0 /)
 rrt_fr%lab  = (/'TIGHT','PLANE STRAIN','SEMITIGHT','LOOSE','LOW STRAIN' , &
                 'WRINKLING' ,'STRONG WRKL.' /)

 ! safety
 rrt_sf%name = 'SAFETY'
 rrt_sf%nval = 7
 ALLOCATE(rrt_sf%vals(rrt_sf%nval+1),rrt_sf%lab(rrt_sf%nval))
 rrt_sf%vals = (/ 0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5D0 /)
 rrt_sf%lab  = (/'STRONG WRKL.','WRINKLING','LOW STRAIN','SAFE','MARGINAL', &
                 'THINNING', 'FAIL' /)

 ! Drawbead Affected Zone
 rrt_db%name = 'Drawbead Affected Zone'
 rrt_db%nval = 2
 ALLOCATE(rrt_db%vals(rrt_db%nval+1),rrt_db%lab(rrt_db%nval))
 rrt_db%vals = (/ -0.5d0, 0.5d0, 1.5d0 /)
 rrt_db%lab  = (/'Free Zone','Affected Zone'/)


 ! Check if configuration file exists

 ! Search in current path (path of the example)
 OPEN ( 7,FILE='stp.cfg',STATUS='OLD',IOSTAT=ios, FORM='FORMATTED',READONLY)
 ! Search in the path of the executable
 IF( ios /= 0 ) THEN
   CALL get_command_argument(0, exepath)
   pos= INDEX(exepath,'/',.TRUE.)
   IF (pos==0) pos= INDEX(exepath,'\',.TRUE.)
   exepath=exepath(1:pos)
   OPEN ( 7,FILE=TRIM(exepath)//'stp.cfg',STATUS='OLD',IOSTAT=ios,  &
                      FORM='FORMATTED',READONLY)
 END IF
 ! Search in developers' path
 IF( ios /= 0 )OPEN ( 7,FILE='c:\uti\stp.cfg',STATUS='OLD',IOSTAT=ios,  &
                      FORM='FORMATTED',READONLY)
 IF( ios /= 0 )THEN
   WRITE(6,"('    No configuration file found in present directory',/, &
 &           '    DEFAULT values will be used')")
   RETURN
 END IF

 DO

   READ(7,'(a1,a99)',IOSTAT=ios)c,lin   !read a line
   IF( ios /= 0 )THEN
     WRITE(6,"('    Error reading configuration file last line read',/, &
              & a1,a99)")c,lin
     EXIT
   END IF

   SELECT CASE (c)                   !according to first character
   CASE ('$')                        !comment
     CYCLE                           !do nothing

   CASE ('#')  !New element type, set defaults all to false
     c = ' '
     CALL readln( c,lin,var,rval,flags,names,cname,units )  !nullify first character to use readln
     ctype = var                  !get data set type

     SELECT CASE (ctype)             !according to key-word

     CASE ('TOOL ')                  !modify post-process tool
       SELECT CASE (flags(1))
       CASE ('F')              !TecPlot Binary (360)
         ip = 1
       CASE ('G')              !GiD
         ip = 2
       CASE ('T')              !TecPlot ASCII (10)
         ip = 3
       CASE ('B')              !Binary GiD
         ip = 4
       CASE ('A')              !ASCII GiD using GiD routines
         ip = 5
       END SELECT

     CASE ('NODAL')                  !Initializes nodal kinematics
       wtdisp =  .TRUE.   ! write total displacements
       wsdisp =  .FALSE.  ! write stage displacements
       widisp =  .FALSE.  ! write initial stage displacements
       wveloc =  .FALSE.  ! write nodal velocities
       wtempe =  .TRUE.   ! write nodal temperatures
       waccel =  .FALSE.  ! write nodal accelerations
       weuler =  .FALSE.  ! write nodal Euler angles
       wangve =  .FALSE.  ! write nodal angular velocites
       wangac =  .FALSE.  ! write nodal angular accelerations

     CASE ('SPOT')                  !Initializest SPOT element
       spot_force =  0  !write axial force
       spot_shear =  0  !write shear force
       spot_eqpst =  0  !write equivalent plastic strain

     CASE ('TRUSS')                  !Initializest TRUSS element
       truss_force =  0  !write total force
       truss_stres =  0  !write Cachy stress
       truss_eqpst =  0  !write equivalent plastic strain

     CASE ('SOL2D')      !Initializes 2-D SOLID element
       sol2d_stres =  0   !write Cauchy stresses
       sol2d_logst =  0   !write logarithmic principal strains
       sol2d_shtst =  0   !write sheet strains
       sol2d_thrat =  0   !write thickness ratio
       sol2d_eqpst =  0   !write effective platic strains
       sol2d_vmise =  0   !write von Mises stress
       sol2d_fldma =  0     !write FLD MAP
       sol2d_wfldSZ = .FALSE.  !write Safety Zone
       sol2d_wfldFZ = .FALSE.  !write Forming Zone
       sol2d_fldva = 0d0    !default value for FLD MAP

     CASE ('SOL3D')                  !Initializes 3-D SOLID element
       sol3d_stres =  0   !write Cauchy stresses
       sol3d_logst =  0   !write logarithmic principal strains
       sol3d_thrat =  0   !write thickness ratio
       sol3d_eqpst =  1   !write effective platic strains
       sol3d_vmise =  0   !write von Mises stress
       sol3d_fldma =  0   !write FLD MAP
       sol3d_fldva = 0.0  !default value for FLD MAP
       sol3d_prsmo = .FALSE. ! standard smoothing type for 15-node prism
       IF( flags(2) == 'S' ) sol3d_prsmo = .TRUE.
       sol3d_press =  0   !write nodal pressure

     CASE ('BEAME')                  !Initializes 3-D BEAM element
       beame_force =  0   !write integrated stresses
       beame_momen =  0   !write integrated moments
       beame_eqpst =  0   !write effective platic strains
       beame_vmise =  0   !write von Mises stress
       beame_arrat =  0   !write thickness ratio

     CASE ('SHREV')                  !Initializes 2-D SHELL/BEAM element
       shrev_force =  0   !write integrated stresses
       shrev_shear =  0   !write integrated shear
       shrev_momen =  0   !write integrated moments
       shrev_eqpst =  0   !write effective platic strains
       shrev_vmise =  0   !write von Mises stress
       shrev_thrat =  0   !write thickness ratio

     CASE ('SHL3D')                  !Initializes 3-D SHELL(Simo) element
       shl3d_force =  0   !write integrated stresses
       shl3d_momen =  0   !write integrated moments
       shl3d_shear =  0   !write integrated transvers shear stresses
       shl3d_logst =  0   !write log strains
       shl3d_eqpst =  0   !write effective platic strains
       shl3d_vmise =  0   !write von Mises stress
       shl3d_thrat =  0   !write thickness ratio

     CASE ('BST  ')                  !Initializes 3-D BST Shell element
       bst_force =  0     !write Cauchy integrated stresses
       bst_momen =  0     !write Cauchy integrated stress moments
       bst_shear =  0     !write integrated transvers shear stresses
       bst_logst =  0     !write logarithmic principal strains
       bst_curva =  0     !write principal curvatures
       bst_thrat =  0     !write thickness ratio
       bst_vmise =  0     !write von mises equivalent stress
       bst_eqpst =  0     !write effective platic strains
       bst_fldma =  0     !write FLD MAP
       bst_wfldSZ = .FALSE.  !write Safety Zone
       bst_wfldFZ = .FALSE.  !write Forming Zone
       bst_fldva = 0d0    !default value for FLD MAP

     CASE ('RIGID')                  !Initializes RIGID element
       rigid_show = .FALSE.

     CASE ('CONTA')                  !Initializes CONTACT data
       wpress  = .FALSE.
       wwrink  = .FALSE.
       wwear   = .FALSE.
       drawb_wafz = .TRUE.
     CASE ('END  ')                  !end of configuration file
       EXIT               ! exit of configuration file

     END SELECT

   CASE DEFAULT           ! read data for a particular value to process

     CALL readln( c ,lin,var,rval,flags,names,cname,units )  !process words in line read

     SELECT CASE (ctype)            !according to set type

     CASE ('NODAL')                 !NODAL kinematics
       SELECT CASE (var)
       CASE ('TDISP')                    !total displacement
         wtdisp = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,tdisp_l,rval,units,tdisp_f,tdisp_u)
       CASE ('SDISP')                    !stage displacement
         wsdisp = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,sdisp_l,rval,units,sdisp_f,sdisp_u)
       CASE ('IDISP')                    !initial displacement
         widisp = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,idisp_l,rval,units,idisp_f,idisp_u)
       CASE ('VELOC')                    !translational velocities
         wveloc = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,veloc_l,rval,units,veloc_f,veloc_u)
       CASE ('TEMPE')                    !nodal temperatures
         wtempe = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,tempe_l,rval,units,tempe_f,tempe_u)
       CASE ('ACCEL')                    !translational accelerations
         waccel = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,accel_l,rval,units,accel_f,accel_u)
       CASE ('EULER')                    !Euler angles
         weuler = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,euler_l,rval,units,euler_f,euler_u)
       CASE ('ANGVE')                    !angular velocities
         wangve = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,angve_l,rval,units,angve_f,angve_u)
       CASE ('ANGAC')                    !angular accelerations
         wangac = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,angac_l,rval,units,angac_f,angac_u)
       END SELECT

     CASE ('SPOT')                 ! SPOT element
       SELECT CASE (var)
       CASE ('FORCE')                    !section force
         spot_force = set_value(flags)
         CALL as_names( 1, 1,flags,names,cname,'G',1,1,spot_force_l,rval,units,spot_force_f,spot_force_u)
       CASE ('SHEAR')                    !normal stress
         spot_shear = set_value(flags)
         CALL as_names( 1, 1,flags,names,cname,'G',1,1,spot_shear_l,rval,units,spot_shear_f,spot_shear_u)
       CASE ('EQPST')                    !equivalent plastic strain
         spot_eqpst = set_value(flags)
         CALL as_names( 1, 1,flags,names,cname,'G',1,1,spot_eqpst_l,rval,units,spot_eqpst_f,spot_eqpst_u)
       END SELECT

     CASE ('TRUSS')                 ! TRUSS element
       SELECT CASE (var)
       CASE ('FORCE')                    !section force
         truss_force = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,truss_force_l,rval,units,truss_force_f,truss_force_u)
       CASE ('STRES')                    !normal stress
         truss_stres = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,truss_stres_l,rval,units,truss_stres_f,truss_stres_u)
       CASE ('EQPST')                    !equivalent plastic strain
         truss_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,truss_eqpst_l,rval,units,truss_eqpst_f,truss_eqpst_u)
       END SELECT

     CASE ('SOL2D')     ! 2D-SOLID elements
       SELECT CASE (var)
       CASE ('STRES')                    !Cauchy stresses
         sol2d_stres = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,4,sol2d_stres_l,rval,units,sol2d_stres_f,sol2d_stres_u)
       CASE ('LOGST')                    !logarithmic strains
         sol2d_logst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,sol2d_logst_l,rval,units,sol2d_logst_f,sol2d_logst_u)
       CASE ('SHTST')                    !logarithmic strains
         sol2d_shtst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,sol2d_shtst_l,rval,units,sol2d_shtst_f,sol2d_shtst_u)
       CASE ('THRAT')                    !thickness ratio
         sol2d_thrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol2d_thrat_l,rval,units,sol2d_thrat_f,sol2d_thrat_u)
       CASE ('EQPST')                    !equivalent plastic strain
         sol2d_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol2d_eqpst_l,rval,units,sol2d_eqpst_f,sol2d_eqpst_u)
       CASE ('VMISE')                    !von Mises equivalent stress
         sol2d_vmise = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol2d_vmise_l,rval,units,sol2d_vmise_f,sol2d_vmise_u)
       CASE ('FLDMA')                    !forming limit diagram
         sol2d_fldma = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol2d_fldma_l)
       CASE ('FORMI')                    !Safety Zone diagram
         sol2d_wfldFZ = flags(1) == 'Y'
         CALL as_names( 1, 1,flags,names,cname,'Y',1,0,sol2d_fldFZ_l)
       CASE ('SAFET')                    !Safety Zone diagram
         sol2d_wfldSZ = flags(1) == 'Y'
         CALL as_names( 1, 1,flags,names,cname,'Y',1,0,sol2d_fldSZ_l)
       CASE ('FLDVA')                    !forming limit diagram default value
         sol2d_fldva = rval
       END SELECT

     CASE ('SOL3D')                    ! 3D-SOLID elements
       SELECT CASE (var)
       CASE ('STRES')                    !Cauchy stresses
         sol3d_stres = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,6,sol3d_stres_l,rval,units,sol3d_stres_f,sol3d_stres_u)
       CASE ('LOGST')                    !logarithmic strains
         sol3d_logst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,sol3d_logst_l,rval,units,sol3d_logst_f,sol3d_logst_u)
       CASE ('THRAT')                    !thickness ratio
         sol3d_thrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol3d_thrat_l,rval,units,sol3d_thrat_f,sol3d_thrat_u)
       CASE ('EQPST')                    !equivalent plastic strain
         sol3d_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol3d_eqpst_l,rval,units,sol3d_eqpst_f,sol3d_eqpst_u)
       CASE ('VMISE')                    !von Mises equivalent stress
         sol3d_vmise = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol3d_vmise_l,rval,units,sol3d_vmise_f,sol3d_vmise_u)
       CASE ('FLDMA')                    !forming limit diagram
         sol3d_fldma = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,sol3d_fldma_l)
       CASE ('FLDVA')                    !forming limit diagram default value
         sol3d_fldva = rval
       CASE ('PRESS')                    !Cauchy stresses
         IF( MOD(sol3d_stres,2) == 1 )THEN
           sol3d_press = set_value(flags)
           CALL as_names( 2, 6,flags,names,cname,'N',1,1,sol3d_press_l,rval,units,sol3d_press_f,sol3d_press_u)
         END IF
       END SELECT

     CASE ('BEAME')                    !3-D beam element (Simo theory)
       SELECT CASE (var)
       CASE ('FORCE')                    !integrated stresses
         beame_force = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,beame_force_l,rval,units,beame_force_f,beame_force_u)
       CASE ('MOMEN')                    !integrated moments
         beame_momen = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,beame_momen_l,rval,units,beame_momen_f,beame_momen_u)
       CASE ('EQPST')                    !equivalent plastic strain
         beame_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,beame_eqpst_l,rval,units,beame_eqpst_f,beame_eqpst_u)
       !CASE ('VMISE')                    !von Mises equivalent stress
       !  beame_vmise = set_value(flags)
       CASE ('ARRAT')                    !area ratio
         beame_arrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,beame_arrat_l,rval,units,beame_arrat_f,beame_arrat_u)
       END SELECT

     CASE ('SHREV')                    !2-D shell/beam element (Simo theory)
       SELECT CASE (var)
       CASE ('FORCE')                    !integrated stresses
         shrev_force = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,2,shrev_force_l,rval,units,shrev_force_f,shrev_force_u)
       CASE ('SHEAR')                    !integrated stresses
         shrev_shear = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shrev_shear_l,rval,units,shrev_shear_f,shrev_shear_u)
       CASE ('MOMEN')                    !integrated moments
         shrev_momen = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,2,shrev_momen_l,rval,units,shrev_momen_f,shrev_momen_u)
       CASE ('EQPST')                    !equivalent plastic strain
         shrev_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shrev_eqpst_l,rval,units,shrev_eqpst_f,shrev_eqpst_u)
       CASE ('VMISE')                    !von Mises equivalent stress
         shrev_vmise = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shrev_vmise_l,rval,units,shrev_vmise_f,shrev_vmise_u)
       CASE ('THRAT')                    !thickness ratio
         shrev_thrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shrev_thrat_l,rval,units,shrev_thrat_f,shrev_thrat_u)
       END SELECT

     CASE ('SHL3D')                    !3-D shell element (Simo theory)
       SELECT CASE (var)
       CASE ('FORCE')                    !integrated stresses
         shl3d_force = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,shl3d_force_l,rval,units,shl3d_force_f,shl3d_force_u)
       CASE ('MOMEN')                    !integrated moments
         shl3d_momen = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,shl3d_momen_l,rval,units,shl3d_momen_f,shl3d_momen_u)
       CASE ('SHEAR')                    !integrated transverse shear stresses
         shl3d_shear = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,2,shl3d_shear_l,rval,units,shl3d_shear_f,shl3d_shear_u)
       CASE ('LOGST')                    !log strains
         shl3d_logst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,shl3d_logst_l,rval,units,shl3d_logst_f,shl3d_logst_u)
       CASE ('EQPST')                    !equivalent plastic strain
         shl3d_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shl3d_eqpst_l,rval,units,shl3d_eqpst_f,shl3d_eqpst_u)
       CASE ('VMISE')                    !von Mises equivalent stress
         shl3d_vmise = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shl3d_vmise_l,rval,units,shl3d_vmise_f,shl3d_vmise_u)
       CASE ('THRAT')                    !thickness ratio
         shl3d_thrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,shl3d_thrat_l,rval,units,shl3d_thrat_f,shl3d_thrat_u)
       END SELECT

     CASE ('BST  ')                      !3-D shell (BST type element)
       SELECT CASE (var)
       CASE ('FORCE')                    !integrated stresses
         bst_force = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,bst_force_l,rval,units,bst_force_f,bst_force_u)
       CASE ('MOMEN')                    !integrated moments
         bst_momen = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,bst_momen_l,rval,units,bst_momen_f,bst_momen_u)
       CASE ('SHEAR')                    !integrated transverse shear stresses
         bst_shear = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,2,bst_shear_l,rval,units,bst_shear_f,bst_shear_u)
       CASE ('LOGST')                    !logarithmic strains
         bst_logst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,bst_logst_l,rval,units,bst_logst_f,bst_logst_u)
       CASE ('CURVA')                    !curvatures
         bst_curva = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,3,bst_curva_l,rval,units,bst_curva_f,bst_curva_u)
       CASE ('THRAT')                    !thickness ratio
         bst_thrat = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,bst_thrat_l,rval,units,bst_thrat_f,bst_thrat_u)
       CASE ('EQPST')                    !equivalent plastic strain
         bst_eqpst = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,bst_eqpst_l,rval,units,bst_eqpst_f,bst_eqpst_u)
       CASE ('VMISE')                    !von Mises equivalent stress
         bst_vmise = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,bst_vmise_l,rval,units,bst_vmise_f,bst_vmise_u)
       CASE ('FLDMA')                    !forming limit diagram
         bst_fldma = set_value(flags)
         CALL as_names( 2, 6,flags,names,cname,'NG',2,1,bst_fldma_l)
       CASE ('FORMI')                    !Safety Zone diagram
         bst_wfldFZ = flags(1) == 'Y'
         CALL as_names( 1, 1,flags,names,cname,'Y',1,0,bst_fldFZ_l)
       CASE ('SAFET')                    !Safety Zone diagram
         bst_wfldSZ = flags(1) == 'Y'
         CALL as_names( 1, 1,flags,names,cname,'Y',1,0,bst_fldSZ_l)
       CASE ('FLDVA')                    !forming limit diagram default value
         bst_fldva = rval
       END SELECT

     CASE ('RIGID')
       SELECT CASE (var)
       CASE ('SHOW ')                    !show rigid bodies and contact surfaces
         rigid_show = flags(1) == 'Y'
       END SELECT

     CASE ('CONTA')
       SELECT CASE (var)
       CASE ('PRESS')                    !print binder pressure
         wpress = flags(1) == 'Y'
         CALL as_names( 1, 2,flags,names,cname,'Y',1,2,press_l,rval,units,press_f,press_u)
       CASE ('GAPS ')                    !print gaps information
         wwrink = flags(1) == 'Y'
         CALL as_names( 1, 3,flags,names,cname,'Y',1,3,wrink_l,rval,units,wrink_f,wrink_u)
       CASE ('WEAR ')                    !print wearing information
         wwear  = flags(1) == 'Y'
         CALL as_names( 1, 2,flags,names,cname,'Y',1,2,wear_l,rval,units,wear_f,wear_u)
       CASE ('DBAFZ ')                   !Write DrawBeads AFectted Zones
         drawb_wafz = flags(1) == 'Y'
       END SELECT
     END SELECT

   END SELECT
 END DO

 CLOSE(7)
 RETURN

 END SUBROUTINE rdconf

 FUNCTION set_value (flags)
 !
 !  Compute printing flag
 !
 !  Flag values : 0 do not print
 !              : 1 print at nodal points
 !              : 2 print at Gauss points
 !              : 3 print at nodal and Gauss points
 !
 IMPLICIT NONE
 CHARACTER(len=1), INTENT(IN) :: flags(2)
 INTEGER :: set_value

 set_value = 0                                                       !initializes
 IF(flags(1) == 'G' .OR. flags(2) == 'G' ) set_value = 2             !at integration points
 IF(flags(1) == 'N' .OR. flags(2) == 'N' ) set_value = set_value + 1 !at nodal points

 RETURN
 END FUNCTION set_value

 SUBROUTINE as_names(nd,nm,flags,names,cname,fd,nn,nc,old_l,rval,units,var_f,var_u)
 ! puts labels read into data base (ASsign NAMES)
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nd, & !number of data variables (READ)
                                 nm, & !number of components (READ)
                                 nn, & !number of variables to be included
                                 nc    !number of components to be included
 CHARACTER (len=1), INTENT(IN) :: flags(nd), & !flags read
                                  fd(nn)       !flags to be considered
 CHARACTER (len=15), INTENT(IN) :: names(nd), & !variable names (READ)
                                   cname(nm,nd) !variable components (READ)
 CHARACTER (len=15), INTENT(IN OUT) :: old_l(nn*(nc+1)) !variables and components names
 !
 REAL(kind=8), OPTIONAL :: rval, var_f        !factor to change units
 CHARACTER (len=15), OPTIONAL :: units,var_u  !units
 ! local variables
 INTEGER (kind=4) :: i,j,k,l
 CHARACTER (len=1) :: letra

 DO j=1,nn   !for each searched variable
   letra = fd(j)      !possible flag
   l = j*nc + j - nc  !first position in OLD_L
   DO i=1,nd !for each data variable (NAMES and CNAME)
     IF(letra == flags(i))THEN !if searched flag found
       IF( LEN_TRIM(names(i)) > 0) THEN !see if a name was input
         old_l(l) = names(i)            !new name
         DO k=1,nc                      !check and modify each component
           IF( LEN_TRIM(cname(k,i)) > 0 ) old_l(l+k) = cname(k,i)
         END DO
       END IF
       EXIT
     END IF
   END DO
 END DO
 ! if present consider factor and units
 IF( PRESENT(rval) )THEN
   IF( rval /= 0D0 )THEN
     var_f = rval
   ELSE
     var_f = 1d0
   END IF
   var_u = units
 END IF
 RETURN
 END SUBROUTINE as_names
