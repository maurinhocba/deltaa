 SUBROUTINE restar()

 ! reads data for a restart

 USE c_input, ONLY : openfi
       ! Global data bases
 USE cms_db
 USE curv_db    !Variables for curve information
 USE dyna_db    !equation information
 USE esets_db
 USE ctrl_db, ONLY : npoin,ndime,neq,maxa,nrotd,neulr,numct,text,nload,ndyna
 USE kinc_db, ONLY :nvelr,rest_kinc
 USE loa_db    !load sets
 USE mat_dba     !material data
 USE npo_db    !point information
 USE nsets_db, ONLY : rest_nsdb
 USE outp_db    !output asked, frequencies
 USE solv_db    !stiffness matrices, iterative vectors & control parameters
 USE name_db,ONLY : namr  !output
 USE nsld_db,ONLY: nnsld

 IMPLICIT NONE

 !     Other auxiliar variables

 !     local variables
 REAL (kind=8) :: cpuol, cpuac, cpuio

 INTERFACE
   INCLUDE 'timuse.h'
   INCLUDE 'elemnt.h'
   INCLUDE 'contac.h'
   !INCLUDE 'restda.h'
 END INTERFACE

 CALL timuse(cpuol)                      !present CPU time

 CALL openfi(nfile=51,rest=TRIM(namr),flag=0)         !open RESTART file

 READ(51) text,time,cpuio   !first set of variables
 ! READ curve information
 CALL rest_curv( )                       ! reads curves information
 CALL rest_nsdb( )                       ! reads nodal sets database
 CALL rest_npo( )                       !point information
 CALL rest_mate( )                       !material data
 CALL rest_esets( )
 CALL elemnt('RESTAR')
 CALL rest_kinc(nnsld)    !restores kinematic constrains parameters
 CALL rest_dyna(neq,maxa)              !dynamic parameters
 CALL rest_outp( )                     !output asked, frequencies
 CALL rest_lo(ndime,nrotd,nload) !load sets
 CALL rest_solv(maxa,neq)              !stiffness matrices, iterative vectors &
                                       !control parameters
 IF(numct > 0) CALL contac('RESTAR')   !read contact data base
 CALL rest_cm( nrotd )                 !concentrated mass information

 READ (51) cpuac                         !read CPU time at dumping

 CLOSE(51,STATUS='keep')                 !close RESTART file

 IF(ndyna == 0)THEN                  !for static problems allocate auxiliar arrays
   ALLOCATE( coor1(ndime,npoin) )        !coordinates
   IF(neulr > 0)THEN                     !local coordinate systems
     ALLOCATE( locs1(neulr,npoin), eule1(neulr,npoin) )
   ELSE
     ALLOCATE( locs1(1,1), eule1(1,1) )
   END IF
 END IF
 ! this lines must be rewritten
 cpui = cpui + cpuio - cpuac             !corrected initial CPU time
 CALL timuse(cpuac)                      !present CPU time
 time(1) = time(1) - cpuol + cpuac !add time devoted to restart

 WRITE(*,"('+ PROGRAM SUCCESFULLY RESTARTED                  ')")

 RETURN
 END SUBROUTINE restar
