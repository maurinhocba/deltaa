 PROGRAM alpha
 !************************************************************************
 !                                                                       *
 !                                                                       *
 !                                                                       *
 !          GRUPO DE METODOS NUMERICOS EN MECANICA ESTRUCTURAL           *
 !                     DEPARTAMENTO DE ESTRUCTURAS                       *
 !           FACULTAD DE CIENCIAS EXACTAS FISICAS Y NATURALES            *
 !                   UNIVERSIDAD NACIONAL DE CORDOBA                     *
 !                                  &                                    *
 !                            C O N I C E T                              *
 !                                                                       *
 !.......................................................................*
 !                                                                       *
 !                             A  L  P  H  A                             *
 !                                                                       *
 !                                                                       *
 !     P R O G R A M :    for elastic (static-dynamic)                   *
 !                        analysis of structures                         *
 !                        by the finite element method                   *
 !                                                                       *
 !     author:   Fernando  G.   F L O R E S                              *
 !                                                                       *
 !     version  '  6.1 '        date: May    2016                        *
 !                                                                       *
 !************************************************************************
       ! Global data bases
 USE DFLIB, ONLY : BEEPQQ  !call atention in interactive mode
 USE param_db,ONLY : milb, mnam     !Global parameters
 USE ctrl_db
 USE npo_db, ONLY : emass,ifpre,mass  !point information
 USE outp_db, ONLY : nreqc,nprqc,time,cpui             !required output frequencies
 USE export_db, ONLY : exp_data, export, neset              !export data
 USE c_input                        !input managment routines
 USE solv_db, ONLY : ecdis,ncdis,ddisp,nbuck,lambd,dlamb,buckl   !stiffness matrices, iterative vectors & control parameters

 USE name_db, ONLY : gen_data_file
 IMPLICIT NONE

 !     Variables for DUMPING control
 INTEGER  (kind=4) :: ircon,   &  ! counter of restart files
                      nrest,   &  ! restart dumping period
                      kstep


 !     Other auxiliar variables

 CHARACTER (len=milb) :: actio      !new problem or restart
 INTEGER (kind=4) :: istop, &  !flag to stop execution
                     chnode

 INTERFACE
    INCLUDE 'byebye.h'
    INCLUDE 'conini.h'
    INCLUDE 'contac.h'
    INCLUDE 'coract.h'
    INCLUDE 'dynast.h'
    INCLUDE 'elemnt.h'
    !INCLUDE 'exten1.h'
    INCLUDE 'extend.h'
    INCLUDE 'iterat.h'
    !INCLUDE 'outdyn.h'
    !INCLUDE 'restar.h'
    INCLUDE 'switch.h'
    INCLUDE 'timing.h'
 END INTERFACE

 !     INITIALIZATION PHASE

 !  actio
 !     assigned in beginp:
 !     NEW     - new problem
 !     RESTAR  - restart an old problem (no changes or minor)
 !     NSTRA0  - restart an old problem with a new strategy
 !     can be changed later to:
 !     NSTRA1  - a new strategy will involve change of boundary conditions without changing nodes
 !     NSTRA2  - a new strategy will involve change of both boundary conditions and nodes

 CALL beginp(actio,nrest,memo,cpui,kstep)  !read

 IF (TRIM(actio) == 'RESTAR' .OR. TRIM(actio) == 'NSTRA0') THEN !OLD problem
   CALL timing(3,.TRUE.)        !restart input/output
   CALL restar( )
   CALL timing(3,.FALSE.)
 ELSE                                                           !NEW problem
   time = 0d0   !initializes process time
   CALL timing(1,.TRUE.)  !initializes CPU time for the current model
 END IF

 strategies : DO ! loop over strategies  =====================================

   IF(TRIM(actio) == 'RESTAR' ) THEN
     ! re-open post-process file for restart
     CALL wrtpos( 1 )                   ! write data for post-processing
   ELSE
     CALL timing(2,.TRUE.) !input data
     istop = 0; ircon = 0 !

     CALL contol(actio)
     !WRITE(*,"(/,5X,'***  BEGIN STAGE ',A,': ',A,/)") TRIM(inttoch(nstra,0)),text

     CALL inpdat( actio)

     CALL histor( )                       !read data for history curves

     CALL loadpl( )

     CALL conini (actio,memo)

     CALL timing(8,.TRUE.)
     CALL intime( )  !
     CALL timing(8,.FALSE.)

     CALL timing(2,.FALSE.)

     CALL timing(4,.TRUE.)

     !          Export data input
     CALL timing(3,.FALSE.)
     CALL exp_data ( )                    !read EXPort DATA for future use
     CALL timing(3,.TRUE.)
     CALL endstr( )      ! end of input data for one strategy

     CALL eqnums (actio)   !compute active DOF numbers

     IF( ncdis > 10 )THEN    !update control equation for output
       ecdis = MOD(ncdis,10)
       ncdis = chnode(ncdis/10)*10 + ecdis
       ecdis = MAX(ifpre(ecdis,ncdis/10),1) !compute control DOF
     ELSE
       ecdis = 1
     END IF

     CALL elemnt('GAUSSC', istop=istop, flag1=.TRUE. )
     IF(istop /= 0) EXIT strategies

     !           Initial conditions and time history data input

     CALL loaini ( )

     CALL timing(4,.FALSE.)
     istep = 0  !initializes strategy step

     CALL timing(11,.TRUE.)             !post-process output
     CALL newprb( )                     !set output file name root
     CALL wrtpos( 0 )                   ! write data for post-processing
     CALL output(.TRUE.,.FALSE.,ecdis)  ! initial state output for post-process
     CALL timing(11,.FALSE.)
   END IF
   ! Mass matrix is considered constant even if slave DOFs exist !!!
   CALL timing(5,.TRUE.)       !mass matrix computation
   CALL masmtx(istop)
   CALL timing(5,.FALSE.)      !mass matrix computation
   IF(istop /= 0) EXIT strategies

   IF(ndyna > 0) THEN
     CALL timing(5,.TRUE.)       !mass matrix computation
     IF(MOD(ndyna,2) == 0)CALL ensmal(ndofn*npoin,ifpre(1,1),emass(1,1),mass(1))
     CALL timing(5,.FALSE.)      !mass matrix computation
   END IF

   DO
     istep = istep+1
     CALL timing(10,.TRUE.)  !iterative process
     IF(ndyna == 0) THEN
       IF(nbuck > 0 .AND. istep > MOD(nbuck,1000)                    &
                    .AND. ABS((buckl-lambd)/dlamb) < 1.0d0 )THEN
         !IF(nbuck < 1000)THEN
           CALL extend(istop,actio)
         !ELSE
         !  CALL exten1(istop)
         !END IF
     !!  ELSE IF (nbuck < 0) THEN
     !!    CALL switch(istop,actio)
       ELSE
         CALL iterat(istop,actio)
       END IF
     ELSE
       CALL dynast(istop,actio)
     END IF
     IF( istep == nstep )istop = -1

     CALL timing(10,.FALSE.) !iterative process
     CALL timing(11,.TRUE.)

     CALL output (.FALSE.,istop,ecdis)

     CALL timing(11,.FALSE.)

     IF(nrest >0 .AND.( MOD(istep,nrest) == 0 .OR. istop /= 0)) THEN
       CALL timing(3,.TRUE.)
       CALL dumpin (ircon,istop)
       CALL timing(3,.FALSE.)
     END IF

     IF(istop /= 0 )EXIT
     actio = 'STEP'
   END DO
   IF( neset > 0 ) CALL export( )
   CALL timing(1,.FALSE.) !closes CPU time for the current model
   CALL byebye (time)
   CALL listen('ALPHAP')
   IF( exists('ENDDAT') )EXIT
   !IF(.NOT.exists('NEWSTR') ) WRITE(*,"(' key-word NEW_STRATEGY missing')")
   actio = 'NSTRA0'
   time = 0d0   !initializes process time
   CALL timing( 1,.TRUE.)  !initializes CPU time for the current model
 END DO strategies

 !CALL outdyn( ifpre )

 IF( istop == -1 )THEN
   WRITE(6,"(T53,' NORMAL END OF EXECUTION   ',/)")
 ELSE
   WRITE(6,"(T20,' ABNORMAL END OF EXECUTION   ',/,T20,' CHECK RESULTS')")
 END IF

 IF( gen_data_file )THEN
   CLOSE(ludat,STATUS='delete')  !PROGRAM.dat file
 ELSE
   CLOSE(ludat,STATUS='keep')    !Users data file
 END IF
 CLOSE(lures,STATUS='KEEP') !
 CLOSE (55)                 !extended report file
 CLOSE (58)                 !debug file for developers

 CALL BEEPQQ(5000,50)  !call atention

 STOP
 END PROGRAM alpha
