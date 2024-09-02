 SUBROUTINE dumpin(ircon,flag)
 !
 !  dumping of problem data for restart
 !
 USE param_db,ONLY: mlin
 USE c_input,ONLY: openfi
       ! Global data bases
 USE cms_db
 USE curv_db    !Variables for curve information
 USE dyna_db    !equation information
 USE esets_db
 USE ctrl_db, ONLY : ndime,nrotd,npoin,numct,neq,maxa,text,nload
 USE kinc_db, ONLY : nvelr,dump_kinc
 USE loa_db     !load sets
 USE mat_dba    !material data
 USE npo_db     !point information
 USE nsets_db, ONLY : dump_nsdb
 USE outp_db    !output asked, frequencies
 USE solv_db    !stiffness matrices, iterative vectors & control parameters
 USE nsld_db,ONLY: nnsld

 IMPLICIT NONE

 !     Variables for DUMPING control

 INTEGER (kind=4), INTENT(IN OUT) :: ircon,flag

 !     local variables
 CHARACTER(len=mlin) :: namei
 CHARACTER(len=1),PARAMETER :: ni(0:9) = &
                               (/'0','1','2','3','4','5','6','7','8','9'/)
 INTEGER (kind=4) :: ircod,j2,j3,j4
 REAL (kind=8) :: cpuac

 INTERFACE
   INCLUDE 'elemnt.h'
   INCLUDE 'timuse.h'
   INCLUDE 'contac.h'
 END INTERFACE

 !       ircon defines the file where to WRITE
 ircon = ircon+1
 !       ircod defines the old file to delete (if exist)
 ircod = ircon-3  !keeps latest three files for restart
 IF(ircod > 0) THEN
   j2 = MOD(ircod,1000)/100
   j3 = MOD(ircod,100)/10
   j4 = MOD(ircod,10)
   namei = '.'//ni(j2)//ni(j3)//ni(j4)
   CALL openfi(nfile=50,rest=namei)
   CLOSE(50,STATUS='delete')
 END IF
 j2 = MOD(ircon,1000)/100
 j3 = MOD(ircon,100)/10
 j4 = MOD(ircon,10)
 namei = '.'//ni(j2)//ni(j3)//ni(j4)
 CALL openfi(nfile=50,rest=namei,flag=flag+1)

 WRITE(50,ERR=9999) text,time,cpui

 ! writes curve information
 CALL dump_curv( )    ! dumps curve database
 CALL dump_nsdb( )    ! dumps nodal sets database
 CALL dump_npo( )                     !point information
 CALL dump_mate( )                     !material data
 CALL dump_esets( )                    !element sets
 CALL elemnt('DUMPIN')                 !element data base
 WRITE(50,ERR=9999)'ENDSET'            !flag to indicate end of element sets
 CALL dump_kinc(nnsld)  ! dumps kinematic constrains parameters
 CALL dump_dyna(neq,maxa)              !dynamic parameters
 CALL dump_outp( )                     !output asked, frequencies
 CALL dump_lo(ndime,nrotd,nload)       !load sets
 CALL dump_solv(neq)                   !stiffness matrices, iterative vectors &
                                       !control parameters
 IF( numct > 0 )CALL contac('DUMPIN')
 CALL dump_cm( nrotd )                 !concentrated mass information

 CALL timuse(cpuac)
 WRITE(50,ERR=9999) cpuac

 CLOSE(50,STATUS='keep')

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dumpin
