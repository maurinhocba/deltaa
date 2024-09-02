 SUBROUTINE tecvar ( )

 USE data_db
 USE cont_db
 IMPLICIT NONE
 ! local variables
 CHARACTER(len=3 ) :: ch,comp  !function to convert a 3 digit integer into a string
 INTEGER :: nnode,nelem,i,iv,nd,node,pv,j,k


 ! INITIALIZES LIST OF VARIABLES

 t_p  = 0  !initializes pointer
 t_nv = 0  !initializes number of variables
 !          initializes strings
 t_vars = ''
 sh_var = ''

 !*************************************************************
 ! Generates list of variables (NOT using STP.CFG data yet)
 ! TO BE modified

 ! FIRST compute number of global nodal variables

 ! deformed configuration
 t_nv = t_nv + ndime
 pv  = t_p + 1
 t_p = t_p + 10
 t_vars(pv:t_p) = '"XX" "YY" '
 IF( ndime == 3 )THEN
   pv = t_p + 1
   t_p  = t_p + 5
   t_vars(pv:t_p) = '"ZZ" '
 END IF

 !Total displacements
 IF( wtdisp )THEN
   DO i=1,ndime
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(tdisp_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(tdisp_l(i+1))//'" '
   END DO
   t_nv = t_nv + ndime
 END IF

 !Stage displacements
 IF( wsdisp )THEN
   DO i=1,ndime
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(sdisp_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(sdisp_l(i+1))//'" '
   END DO
   t_nv = t_nv + ndime
 END IF

 !Initial displacements
 IF( widisp )THEN
   DO i=1,ndime
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(idisp_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(idisp_l(i+1))//'" '
   END DO
   t_nv = t_nv + ndime
 END IF

 !additional displacements
 IF( addof .AND. waddof )THEN
   DO i=1,2
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(addof_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(addof_l(i+1))//'" '
   END DO
   t_nv = t_nv + 2
 END IF

 !Velocities
 IF( wveloc )THEN
   DO i=1,ndime
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(veloc_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(veloc_l(i+1))//'" '
   END DO
   t_nv = t_nv + ndime
 END IF

 !Accelerations
 IF( waccel )THEN
   DO i=1,ndime
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(accel_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(accel_l(i+1))//'" '
   END DO
   t_nv = t_nv + ndime
 END IF

 !Euler Angles
 IF( weuler )THEN
   pv = t_p + 1
   t_p  = t_p + 3+ LEN_TRIM(euler_l(2))
   t_vars(pv:t_p) = '"'//TRIM(euler_l(2))//'" '
   t_nv = t_nv + 1
   IF( ndime == 3 )THEN
     DO i=2,ndime
       pv = t_p + 1
       t_p  = t_p + 3+ LEN_TRIM(euler_l(i+1))
       t_vars(pv:t_p) = '"'//TRIM(euler_l(i+1))//'" '
     END DO
     t_nv = t_nv + 2
   END IF
 END IF

 !Angular velocities
 IF( wangve )THEN
   nd = ndofn - ndime
   DO i=1,nd
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(angve_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(angve_l(i+1))//'" '
   END DO
   t_nv = t_nv + nd
 END IF

 !Angular accelerations
 IF( wangac )THEN
   nd = ndofn - ndime
   DO i=1,nd
     pv = t_p + 1
     t_p  = t_p + 3+ LEN_TRIM(angac_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(angac_l(i+1))//'" '
   END DO
   t_nv = t_nv + nd
 END IF

 !Temperature
 IF( wtempe )THEN
   DO i=1,ndoft
     pv = t_p + 1
     comp=ch(i)
     t_p  = t_p + 3 + LEN_TRIM(tempe_l(2)) + 1
     t_vars(pv:t_p) = '"'//TRIM(tempe_l(2))//comp(3:3)//'" '
   END DO
   t_nv = t_nv + ndoft
 END IF

 ! Wearing Work
 IF( wwear  )THEN
   pv = t_p + 1
   t_p  = t_p + 3 + LEN_TRIM(wear_l(2))
   t_vars(pv:t_p) = '"'//TRIM(wear_l(2))//'" '
   t_nv = t_nv + 1
 END IF

 ! Gap-information
 IF( wwrink )THEN
   DO i=1,3
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(wrink_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(wrink_l(i+1))//'" '
   END DO
   t_nv = t_nv + 3
 END IF

 ! binder-Pressure
 IF( wpress )THEN
   DO i=1,2
     pv = t_p + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(press_l(i+1))
     t_vars(pv:t_p) = '"'//TRIM(press_l(i+1))//'" '
   END DO
   t_nv = t_nv + 2
 END IF

 t_ngv = t_nv  !global nodal variables

 ! SECOND compute elemental nodal variables
 t_ev    = 0 !initializes no elemental variables

 !  spot not included yet

 IF( truss_sets > 0 )THEN  !          T R U S S
   !       smoothed variables
   j = t_nv
   IF( MOD(truss_force,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_force_l(2))
     t_vars(pv:t_p) = '"'//TRIM(truss_force_l(2))//'" '
   END IF
   IF( MOD(truss_stres,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_stres_l(2))
     t_vars(pv:t_p) = '"'//TRIM(truss_stres_l(2))//'" '
   END IF
   IF( MOD(truss_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(truss_eqpst_l(2))//'" '
   END IF
   t_ev(1,2) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF(  truss_force > 1 ) THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_force_l(4))
     t_vars(pv:t_p) = '"'//TRIM(truss_force_l(4))//'" '
   END IF
   IF( truss_stres > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_stres_l(4))
     t_vars(pv:t_p) = '"'//TRIM(truss_stres_l(4))//'" '
   END IF
   IF(truss_eqpst > 1 )THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(truss_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(truss_eqpst_l(4))//'" '
   END IF
   t_ev(2,2) = t_nv - j
 END IF

 IF( sol2d_sets > 0 )THEN        ! S O L I D   2 - D
   j = t_nv
   IF( MOD(sol2d_stres,2) == 1)THEN
     DO i=1,4
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_stres_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_stres_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(sol2d_logst,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_logst_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_logst_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(sol2d_shtst,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_shtst_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_shtst_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(sol2d_thrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_thrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_thrat_l(2))//'" '
   END IF
   IF( MOD(sol2d_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_eqpst_l(2))//'" '
   END IF
   IF( MOD(sol2d_vmise,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_vmise_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_vmise_l(2))//'" '
   END IF
   IF( MOD(sol2d_fldma,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_fldma_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_fldma_l(2))//'" '
   END IF
   !  Forming Zones Diagram  &  Safety Zone Diagram only as CELL-CENTERED

   t_ev(1,3) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( sol2d_stres > 1 ) THEN
     DO i=1,4
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_stres_l(6+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_stres_l(6+i))//'" '
     END DO
   END IF
   IF( sol2d_logst > 1 ) THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_logst_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_logst_l(5+i))//'" '
     END DO
   END IF
   IF( sol2d_shtst > 1 ) THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_shtst_l(4+i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_shtst_l(4+i))//'" '
     END DO
   END IF
   IF( sol2d_thrat > 1 ) THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_thrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_thrat_l(4))//'" '
   END IF
   IF( sol2d_eqpst > 1 ) THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_eqpst_l(4))//'" '
   END IF
   IF( sol2d_vmise > 1 ) THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_vmise_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_vmise_l(4))//'" '
   END IF
   IF( sol2d_fldma > 1 ) THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol2d_fldma_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol2d_fldma_l(4))//'" '
   END IF
   IF( sol2d_wfldFZ .OR. sol2d_wfldSZ ) THEN
     DO i=1,2
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol2d_fldfz_l(i))
       t_vars(pv:t_p) = '"'//TRIM(sol2d_fldfz_l(i))//'" '
     END DO
   END IF
   t_ev(2,3) = t_nv - j
 END IF

 IF( sol3d_sets > 0 )THEN  ! S O L I D   3 - D
   j = t_nv
   IF( MOD(sol3d_stres,2) == 1)THEN
     DO i=1,6
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol3d_stres_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(sol3d_stres_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(sol3d_logst,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol3d_logst_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(sol3d_logst_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(sol3d_thrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_thrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_thrat_l(2))//'" '
   END IF
   IF( MOD(sol3d_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_eqpst_l(2))//'" '
   END IF
   IF( MOD(sol3d_vmise,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_vmise_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_vmise_l(2))//'" '
   END IF
   IF( MOD(sol3d_fldma,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_fldma_l(2))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_fldma_l(2))//'" '
   END IF
   t_ev(1,4) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( sol3d_stres > 1)THEN
     DO i=1,6
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol3d_stres_l(8+i))
       t_vars(pv:t_p) = '"'//TRIM(sol3d_stres_l(8+i))//'" '
     END DO
   END IF
   IF( sol3d_logst > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(sol3d_logst_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(sol3d_logst_l(5+i))//'" '
     END DO
   END IF
   IF( sol3d_thrat > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_thrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_thrat_l(4))//'" '
   END IF
   IF( sol3d_eqpst > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_eqpst_l(4))//'" '
   END IF
   IF( sol3d_vmise > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_vmise_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_vmise_l(4))//'" '
   END IF
   IF( sol3d_fldma > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(sol3d_fldma_l(4))
     t_vars(pv:t_p) = '"'//TRIM(sol3d_fldma_l(4))//'" '
   END IF
   t_ev(2,4) = t_nv - j
 END IF

 IF( shl3d_sets > 0 )THEN    !Shear deformable shells 3D
   j = t_nv
   IF( MOD(shl3d_force,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_force_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_force_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shl3d_momen,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_momen_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_momen_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shl3d_shear,2) == 1)THEN
     DO i=1,2
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_shear_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_shear_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shl3d_logst,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_logst_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_logst_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shl3d_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_eqpst_l(2))//'" '
   END IF
   IF( MOD(shl3d_vmise,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_vmise_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_vmise_l(2))//'" '
   END IF
   IF( MOD(shl3d_thrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_thrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_thrat_l(2))//'" '
   END IF
   t_ev(1,5) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( shl3d_force > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_force_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_force_l(5+i))//'" '
     END DO
   END IF
   IF( shl3d_momen > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_momen_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_momen_l(5+i))//'" '
     END DO
   END IF
   IF( shl3d_shear > 1)THEN
     DO i=1,2
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_shear_l(4+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_shear_l(4+i))//'" '
     END DO
   END IF
   IF( shl3d_logst > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(shl3d_logst_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(shl3d_logst_l(5+i))//'" '
     END DO
   END IF
   IF( shl3d_eqpst > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_eqpst_l(4))//'" '
   END IF
   IF( shl3d_vmise > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_vmise_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_vmise_l(4))//'" '
   END IF
   IF( shl3d_thrat > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shl3d_thrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shl3d_thrat_l(4))//'" '
   END IF
   t_ev(2,5) = t_nv - j
 END IF

 IF( beame_sets > 0 )THEN  !3D beams
   j = t_nv
   IF( MOD(beame_force,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(beame_force_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(beame_force_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(beame_momen,2) == 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(beame_momen_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(beame_momen_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(beame_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(beame_eqpst_l(2))//'" '
   END IF
   IF( MOD(beame_vmise,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_vmise_l(2))
     t_vars(pv:t_p) = '"'//TRIM(beame_vmise_l(2))//'" '
   END IF
   IF( MOD(beame_arrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_arrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(beame_arrat_l(2))//'" '
   END IF
   t_ev(1,6) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( beame_force > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(beame_force_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(beame_force_l(5+i))//'" '
     END DO
   END IF
   IF( beame_momen > 1)THEN
     DO i=1,3
       t_nv = t_nv + 1
       pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(beame_momen_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(beame_momen_l(5+i))//'" '
     END DO
   END IF
   IF( beame_eqpst > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(beame_eqpst_l(4))//'" '
   END IF
   IF( beame_vmise > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_vmise_l(4))
     t_vars(pv:t_p) = '"'//TRIM(beame_vmise_l(4))//'" '
   END IF
   IF( beame_arrat > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(beame_arrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(beame_arrat_l(4))//'" '
   END IF
   t_ev(2,6) = t_nv - j
 END IF

 IF( shrev_sets > 0 )THEN   !2D beams and Shells of revolution
   j = t_nv
   k = MIN(ntype,2)
   IF( MOD(shrev_force,2) == 1)THEN
     DO i=1,k
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(shrev_force_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shrev_force_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shrev_shear,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_shear_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shrev_shear_l(2))//'" '
   END IF
   IF( MOD(shrev_momen,2) == 1)THEN
     DO i=1,k
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(shrev_momen_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(shrev_momen_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(shrev_eqpst,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_eqpst_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shrev_eqpst_l(2))//'" '
   END IF
   IF( MOD(shrev_vmise,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_vmise_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shrev_vmise_l(2))//'" '
   END IF
   IF( MOD(shrev_thrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_thrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(shrev_thrat_l(2))//'" '
   END IF
   t_ev(1,7) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( shrev_force > 1)THEN
     DO i=1,k
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(shrev_force_l(4+i))
       t_vars(pv:t_p) = '"'//TRIM(shrev_force_l(4+i))//'" '
     END DO
   END IF
   IF( shrev_shear > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_shear_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shrev_shear_l(4))//'" '
   END IF
   IF( shrev_momen > 1)THEN
     DO i=1,k
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(shrev_momen_l(4+i))
       t_vars(pv:t_p) = '"'//TRIM(shrev_momen_l(4+i))//'" '
     END DO
   END IF
   IF( shrev_eqpst > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_eqpst_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shrev_eqpst_l(4))//'" '
   END IF
   IF( shrev_vmise > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_vmise_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shrev_vmise_l(4))//'" '
   END IF
   IF( shrev_thrat > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(shrev_thrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(shrev_thrat_l(4))//'" '
   END IF
   t_ev(2,7) = t_nv - j
 END IF

 IF( bst_sets > 0 )THEN  !Rotation Free Shells
   j = t_nv
   IF( MOD(bst_force,2) == 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_force_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_force_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(bst_momen,2) == 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_momen_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_momen_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(bst_shear,2) == 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_shear_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_shear_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(bst_logst,2) == 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_logst_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_logst_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(bst_curva,2) == 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_curva_l(1+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_curva_l(1+i))//'" '
     END DO
   END IF
   IF( MOD(bst_thrat,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(bst_thrat_l(2))
     t_vars(pv:t_p) = '"'//TRIM(bst_thrat_l(2))//'" '
   END IF
   IF( MOD(bst_eqpst,2) == 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_eqpst_l(2)) + 1
       t_vars(pv:t_p) = '"'//TRIM(bst_eqpst_l(2))//'1" '
     END DO
     t_vars(t_p-2:t_p-2) = 'N' !correct
   END IF
   IF( MOD(bst_vmise,2) == 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_vmise_l(2)) + 1
       t_vars(pv:t_p) = '"'//TRIM(bst_vmise_l(2))//'1" '
     END DO
     t_vars(t_p-2:t_p-2) = 'N' !correct
   END IF
   IF( MOD(bst_fldma,2) == 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(bst_fldma_l(2))
     t_vars(pv:t_p) = '"'//TRIM(bst_fldma_l(2))//'" '
   END IF
   t_ev(1,9) = t_nv - j
   ! Gaussian variables
   j = t_nv
   IF( bst_force > 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_force_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_force_l(5+i))//'" '
     END DO
   END IF
   IF( bst_shear > 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_shear_l(4+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_shear_l(4+i))//'" '
     END DO
   END IF
   IF( bst_momen > 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_momen_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_momen_l(5+i))//'" '
     END DO
   END IF
   IF( bst_logst > 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_logst_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_logst_l(5+i))//'" '
     END DO
   END IF
   IF( bst_curva > 1)THEN
     DO i=1,3
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_curva_l(5+i))
       t_vars(pv:t_p) = '"'//TRIM(bst_curva_l(5+i))//'" '
     END DO
   END IF
   IF( bst_thrat > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(bst_thrat_l(4))
     t_vars(pv:t_p) = '"'//TRIM(bst_thrat_l(4))//'" '
   END IF
   IF( bst_eqpst > 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_eqpst_l(4)) + 1
       t_vars(pv:t_p) = '"'//TRIM(bst_eqpst_l(4))//'1" '
     END DO
     t_vars(t_p-2:t_p-2) = 'N' !correct
   END IF
   IF( bst_vmise > 1)THEN
     DO i=1,2
       pv = t_p + 1
       t_nv = t_nv + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_vmise_l(4)) + 1
       t_vars(pv:t_p) = '"'//TRIM(bst_vmise_l(4))//'1" '
     END DO
     t_vars(t_p-2:t_p-2) = 'N' !correct
   END IF
   IF( bst_fldma > 1)THEN
     t_nv = t_nv + 1
     pv = t_p + 1
       t_p  = t_p + 3 + LEN_TRIM(bst_fldma_l(4))
       t_vars(pv:t_p) = '"'//TRIM(bst_fldma_l(4))//'" '
   END IF
   IF( bst_wfldFZ )THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(bst_fldFZ_l(1))
     t_vars(pv:t_p) = '"'//TRIM(bst_fldFZ_l(1))//'" '
   END IF
   IF( bst_wfldSZ )THEN
     t_nv = t_nv + 1
     pv = t_p + 1
     t_p  = t_p + 3 + LEN_TRIM(bst_fldSZ_l(1))
     t_vars(pv:t_p) = '"'//TRIM(bst_fldSZ_l(1))//'" '
   END IF
   t_ev(2,9) = t_nv - j
 END IF

 IF( ASSOCIATED(sh_vara) )DEALLOCATE(sh_vara,ps_vara,lc_vara)
 ALLOCATE(sh_vara(t_nv),ps_vara(t_nv),lc_vara(t_nv))

 RETURN
 END SUBROUTINE tecvar

