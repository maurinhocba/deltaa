 PROGRAM st_2_pos
 !******************************************************************
 !
 !   GID                   post-process interface for
 !    &                    Simpact
 !  TECPLOT                Delta
 !
 !   Developed by
 !   Fernando Flores       GMNME-DEC U. N. Córdoba
 !
 !
 !******************************************************************
  USE param_db,ONLY: mlen
  USE data_db
  USE flc_db, ONLY : nflc,lblflc
  USE cont_db
  IMPLICIT NONE
  INTEGER (kind=4) :: freq,   &   !frequency to read from files
                      ifreq,  &   !step read to compare with FREQ
                      iselec, &   !step print option
                      istep,  &   !Stampack/Simpact step
                      leng,bleng,sleng,  &   !length of problem name
                      istat,  &   !flag for consoler input error
                      number, &   !function to convert string into integer
                      i,j,k,l,&   !auxiliar counters
                      nstep,  &   !post-process number at each stage
                      ns,selst,selec(86)   !no more than 86 steps are possible

 REAL (kind=8) :: oldtm        !previous printed step time

 INTEGER (kind=4) :: nargs     !number of arguments in command line
 CHARACTER(len=mlen) :: input  !string for console input
 CHARACTER(len=mlen) :: inpu1  !auxiliar array for console input
 CHARACTER(len=1 ) :: letra    !character for console input
 CHARACTER(len=1 ) :: postype  !postprocess control type
 LOGICAL ::  allst, &  !flag to print all steps
             same,  &  !flag to check if steps are different
             stage, &  !flag to check if only one stage is processed
             rem,   &  !flag to indicate that remeshing/refinement are presentsta
             jone      !flag to indicate that just one post is desired

 CHARACTER(len=3):: inttoch
 INTERFACE
   INCLUDE 'readcw.h'
 END INTERFACE

 CHARACTER (len=8) :: prognm
 INCLUDE 'prog_data.fpp'

  !------------header-----
  j   = LEN_TRIM(prognm)
  WRITE(6,"(///////33X,<j>(A,1X),/)") (prognm(i:i),i=1,j)
  WRITE( 6,"(33X,'version: ',A/33X,'date   : ',A,/)")TRIM(version),TRIM(verdate)
  DO i=1,7
    WRITE( 6,"(t15,a50)")titles(i)
  END DO
  WRITE (6,"(/)")
  WRITE (6,"(16X,'*** Post Processing ....')")
  !------------------------

 CALL rdconf( )         ! read configuration file

 ! read arguments from comand line (NAME and FREQ)

 ttime = 0d0
 narg = nargs()-1       ! number of arguments in command line
 IF(narg == 0)THEN      ! no arguments, read problem name
   ! read data file if exists
   OPEN(1,FILE='stp.dat',FORM='formatted',STATUS='old',IOSTAT=istat)
   IF (istat == 0 ) READ(1,'(A)',IOSTAT=istat) input   !read first line
   IF( istat == 0 ) THEN
     narg = 1          !update number of arguments read successfully
     READ(1,'(a)',IOSTAT=istat) inpu1   !read a line
     IF (istat == 0) narg = -1
     CLOSE(1)            !close Data file
   ELSE
     PRINT *,' Enter input FILE name (WithOUT extension) '
     READ(5,"(A)") input
   END IF
 ELSE                   ! read problem name from command line
   CALL getarg(1,input,istat)
 END IF
 ! INPUT is checked to see if it is a proyect or a stage
 k = INDEX (input, '_rm', BACK = .TRUE.)            !remeshing/refinement
 IF( k == 0)THEN
   i = INDEX (input, '.@', BACK = .TRUE.)
   IF( i == 0)THEN
     stage = .FALSE.
     bleng = LEN_TRIM(input)
     j = 1    !begin with first state
   ELSE
     stage = .TRUE.
     sleng = LEN_TRIM(input) !length of problem name
   END IF
 ELSE
   jone = .TRUE.
 END IF
 ! read frequency or steps to load
 SELECT CASE (narg)
 CASE(-1)
   !nothing, already read
 CASE(0:1)
   PRINT *,' Enter Frequency to read from files or steps to load'
   READ(5,"(a)") inpu1
 CASE(2:)
   CALL getarg(2,inpu1,istat)
 END SELECT
 IF( stage )THEN   !only if stage is defined
   l = INDEX (inpu1, '{')            !see if specific steps are asked
 ELSE
   l = 0
 END IF
 IF( l == 0 )THEN
   i = LEN_TRIM(inpu1)
   freq = number(inpu1(1:i+1))
   IF( freq <= 0) freq = 1 !check
 ELSE
   freq = 1
   selec = 0
   ns = 0
   ns = 0
   DO
     inpu1(1:l) =''                      !erase characters already read
     i = l+1
     l = INDEX (inpu1, ',')              !see if , delimiter exist
     IF( l == 0) THEN                    !if , is not present
       l = INDEX (inpu1, '}')            !find closing }
       freq = -1                         !set exit flag
     END IF
     ns = ns+1                           !update number of steps asked
     selec(ns) = number(inpu1(i:l-1))      !read step
     IF( freq == -1)EXIT
   END DO
 END IF


 DO !until all stages processed
   IF( .NOT.jone)THEN    !not an specific ref/rem file
     IF(.NOT.stage)THEN  !not an specific stage
       input = input(1:bleng)//'.@'//TRIM(inttoch(j,0))
       sleng = LEN_TRIM(input) !length of problem name
     END IF
     IF( k > 0 ) input = input(1:sleng)//'_rm'//TRIM(inttoch(k,0))
   END IF
   !                         OPEN binary files to be used
   leng = LEN_TRIM(input) !length of problem name
   i = leng+5             !position next to problem name + .??? (extension)
   OPEN (16,FILE=input(1:leng)//'.f16',STATUS='OLD',FORM='UNFORMATTED', &
         SHARE='DENYNONE',SHARED,READONLY,iostat=istat)
   IF( istat /= 0 )THEN  !file not found
     IF(jone)THEN  !if just one file
       WRITE(*,*) 'file '//TRIM(input)//'.f16  not found'
       STOP       !stop with an error message
     ELSE
       IF( k == 0) THEN        !for first remesh/refine
         IF( j == 1)THEN          !for first stage
           WRITE(*,*) 'file '//TRIM(input)//'.@1.f16  not found'
           STOP   !stop with an error message
         ELSE
           EXIT !normal exit
         END IF
       ELSE
         j = j + 1  !increase stage counter
         k = 0      !initializes remesh/refine counter
         CYCLE      !go to next stage
       END IF
     END IF
   ELSE
     WRITE(*,*)TRIM(input),' opened to be post-processed'
     k = k + 1      !update rem/ref counter
   END IF

   OPEN (17,FILE=input(1:leng)//'.f17',STATUS='OLD',FORM='UNFORMATTED', &
         SHARE='DENYNONE',SHARED,READONLY)

   READ(17,ERR=9999) ndime,ndofn,ndoft  ! read first card from post-process files

   IF (ndofn > ndime) THEN           !decide about reading rotation matrices
     IF (ndime == 2)THEN
       neulr = 1       !one value per node
       nrotd = 1       !one value per node
       nadd  = 1
       IF (ndofn > 3) addof = .TRUE.   !additional DOFs
     ELSE IF (ndime == 3) THEN
       neulr = 3                  !3x3 matrix at each node
       nrotd = MIN(ndofn-ndime,3) !
       nadd = 2
       IF (ndofn > 6) addof = .TRUE.   !additional DOFs
     END IF
   ELSE
     neulr = 0                       !no rotation matrices
     nrotd = 0
     addof = .FALSE.                 !no additional DOFs
   END IF

   READ(17)npoin,idyna,text          !number of points, static/dynamic, problem title
   idyna = MIN(idyna,1)              !static or dynamic problem, Time analysis Gid

   IF(ALLOCATED(label))DEALLOCATE ( label, coord, meshn, coors, displ, coorf )

   ALLOCATE ( label(npoin), coord(ndime,npoin), meshn(npoin), & !reserve memory
              coors(ndime,npoin), displ(ndime,npoin), coorf(ndime,npoin)  )

   wtempe = wtempe .AND. (ndoft > 0)
   IF( ASSOCIATED(tempe)) DEALLOCATE(tempe)
   IF( wtempe)ALLOCATE(tempe(ndoft,npoin))

   DO i=1,npoin  !read nodal labels and coordinates (original and stage)
     READ(17) label(i),coord(1:ndime,i),coors(1:ndime,i)
     IF(ip == 1 .OR. ip == 3) label(i) = i     !consecutive order (TECPLOT only)
   END DO
   last_label = label(npoin)
   ! read flc data
   READ(17) numsec
   IF( ALLOCATED(lblflc))DEALLOCATE(lblflc)
   ALLOCATE(lblflc(numsec))
   READ(17)lblflc,nflc
   CALL rdflc ( )

   CALL elmdat(input,leng)           !read element data

   IF( widisp ) THEN
     IF( ASSOCIATED(dispi))DEALLOCATE(dispi)
     ALLOCATE(dispi(ndime,npoin))
     DO i=1,npoin
       dispi(:,i) = coord(:,i) - coors(:,i)  !negative initial displacements
     END DO
   END IF

   IF(idyna == 1)THEN                !reserve space for dynamic problems
     IF( ASSOCIATED(veloc) ) DEALLOCATE( veloc )
     IF( wveloc ) ALLOCATE( veloc(ndime,npoin) )  !Write velocities ?
     IF( ASSOCIATED(accel) ) DEALLOCATE( accel )
     IF( waccel ) ALLOCATE( accel(ndime,npoin) )  !Write accelerations ?
   ELSE
     wveloc = .FALSE.                !Set flags to No Write veloc. nor accelr.
     waccel = .FALSE.
     wangve = .FALSE.
     wangac = .FALSE.
   END IF

   IF (neulr >= 1) THEN          !reserve space if nodal local systems present
     IF ( ASSOCIATED(euler) ) DEALLOCATE ( euler )
     IF ( weuler ) ALLOCATE ( euler(neulr,npoin) )
     IF ( ASSOCIATED(anvel) ) DEALLOCATE ( anvel )
     IF ( wangve .AND. idyna == 1 ) ALLOCATE ( anvel(nrotd,npoin) )
     IF ( ASSOCIATED(anacc) ) DEALLOCATE ( anacc )
     IF ( wangac .AND. idyna == 1 ) ALLOCATE ( anacc(nrotd,npoin) )
     IF ( ASSOCIATED(psi) ) DEALLOCATE ( psi )
     IF ( addof .AND. waddof ) ALLOCATE ( psi(nadd,npoin) )
   ELSE
     weuler = .FALSE.
     wangve = .FALSE.
     wangac = .FALSE.
   END IF
   !                  contact data
   press = .FALSE.        !initializes flags
   wrink = .FALSE.

   IF( wpress .OR. wwrink )THEN  !Write_press or Write_wrink ?
     ! check if information exists
     OPEN (44,FILE=input(1:leng)//'.fc4',STATUS='OLD', ERR = 4,           &
           FORM='UNFORMATTED',SHARE='DENYNONE',SHARED,READONLY)
     CALL readcs ( ) ! read surface data to compute surface contact press
     4 CONTINUE
     wpress = press .AND. wpress
     wwrink = wrink .AND. wwrink
     IF ( ALLOCATED(presn) ) DEALLOCATE ( presn )
     IF(wpress) ALLOCATE ( presn(2,npoin) )    !reserve memory
     IF ( ALLOCATED(wrinkn) ) DEALLOCATE ( wrinkn )
     IF(wwrink) ALLOCATE ( wrinkn(2,npoin) )
   END IF

   IF( wwear )THEN      !write tools wear information ?
     ! check if information exists
     wear  = .FALSE.
     OPEN (45,FILE=input(1:leng)//'.fw4',STATUS='OLD', ERR = 5,           &
           FORM='UNFORMATTED',SHARE='DENYNONE',SHARED,READONLY)
     wear  = .TRUE.
     IF ( ALLOCATED(wwork) ) DEALLOCATE ( wwork, areas )
     ALLOCATE( wwork(npoin), areas(npoin) )
     IF( wwear )CALL readcw (areas,coord,ndime,ntype)
     5 wwear = wear
   END IF


   ! writes geometry definition
   meshn = .FALSE.  !initializes
   SELECT CASE(ip)
   CASE (1)         !TecPlot in Block Format (Binary)
     !CALL tecbinmesh(input,leng)
   CASE (2)         !GiD ASCII
     CALL gidmesh(input,leng)
   CASE (3)         !TecPlot in Point Format (ASCII)
     CALL tecmesh(input,leng)
   CASE (4,5)       !GiD using Library  4:Binary 5:ASCII
     CALL gidbinmesh(input,leng)
   END SELECT

   nstep  = 0         !initializes step for the stage
   iselec = 0         !initialize step selection (write it or not)
   ifreq  = 0         !counts number of read steps
   allst  = narg >= 2 !print all steps if nargs >= 2

   !     loop over each time step
   oldtm = 0d0        !initializes previous printed step time
   same  = .FALSE.    !initializes
   selst = 1

   DO
     READ(17,END=10) istep,ttime,postype  ! step and model time (stampack/simpact)
     SELECT CASE (postype)
     CASE ('T')
       loadstep = 'Time Step '
     CASE ('D')
       loadstep = 'Disp. Step'
     CASE ('C')
       loadstep = 'Curve Val.'
     CASE ('L')
       loadstep = 'Load Anal.'
       ttime = REAL(istep,8)
     CASE DEFAULT
       WRITE(*,"('ERROR: Unknown load step condition.')")
       WRITE(55,"('ERROR: Unknown load step condition.')")
       STOP
     END SELECT

     IF( ttime > 0d0 ) same =  ABS(oldtm - ttime)/ttime < 1e-4  !same step as before ?

     oldtm = ttime              !update last printed step
     istep = MOD(istep,1000000) !limit step number to less than 10^6
     IF( .NOT. same )ifreq = ifreq+1  !increase read steps counter
     nstep = nstep + 1          !number of step read
     IF( freq == -1) THEN
       iselec = 2  !default
       IF( selec(selst) == nstep )THEN
         iselec = 1           !load step
         selst = selst + 1    !update counter
       END IF
     ELSE IF(MOD(ifreq,freq) /= 0 .OR. same )THEN   !don't print the step
       iselec = 2
     ELSE IF( allst )THEN                      !print the step
       iselec = 1
     ELSE                                      !ask the user
       IF(idyna == 1) THEN      !for dynamic problems, use TIME
         WRITE(*,"(///,5x,'Results for Step Number:',i8,          &
               &  '  at ',a10,e14.5,//)") istep,loadstep,ttime
       ELSE                     !for static problems use LOAD FACTOR
         WRITE(*,"(///,5x,'Results for Step Number:',i8,          &
               &  '  for Load Factor:',e12.5,//)")istep,ttime
       END IF

       iselec = 0               !initializes

       DO
         !        print options
         WRITE(*,"(15x,'  load this step     ==> 1 (Yes)',/,  &
               &   15x,'do NOT load the step ==> 2 (No) ',/,  &
               &   15x,'load all the steps   ==> 3 (All)',/,  &
               &   15x,'end                  ==> 4 (End)')")
         READ(*,'(a1)')letra           !read an answer
         SELECT CASE (letra)           !process answer
         CASE ( 'y', 'Y', '1' )        !step must be printed
           iselec= 1
         CASE ( 'n', 'N', '2' )        !step must be skipped
           iselec= 2
         CASE ( 'a', 'A', '3' )        !print all succesive steps
           allst = .TRUE.              !modify ALL STEP flag
           iselec= 1
         CASE ( 'e', 'E', '4' )        !end post-processing interface
           iselec= 4
         CASE DEFAULT
           CYCLE                       !a correct answer is necessary, ask again
         END SELECT
         EXIT
       END DO
     END IF

     IF(iselec == 4) EXIT              ! END of job, exit reading loop

     !       READ displacements, velocities, accelerations and pressures
     CALL readva ( )

     IF (fin) EXIT                     ! end of files found, exit loop
     IF(iselec == 2) THEN              ! does not load the step
       ! reading without writing on the neutral file

       PRINT *,' Step over results for ' ,loadstep,ttime,nstep
       IF(wpress .OR. wwrink)CALL readn (nsurf)
       IF(wwear)READ(45)

       CALL gaussn( )

     ELSE       ! IF(iselec == 1 .OR. iselec == 3)    ! loads the step
       PRINT *,' Processing results for ' ,loadstep,ttime,nstep


       ! prints displacements, local systems, velocities, accelerations and pressure
       CALL wridva( )  !Prints nodal values

       ! prints binder pressure
       IF(wpress .OR. wwrink )CALL presno ( )

       ! prints wearing
       IF(wwear)CALL wearin (npoin,ttime,label(1),ip)

       ! from now on reads and prints gaussian values

       CALL prgaus (input,leng)

     END IF
     IF (fin) EXIT
   END DO

   10 WRITE(*,*)' End of Results-File Found for this stage found'
   SELECT CASE (ip)
   CASE (1)

   CASE (2)
     CLOSE(11)
     CLOSE(13)
   CASE (4,5)
     CALL GID_CLOSEPOSTRESULTFILE()
   END SELECT
   IF( jone ) EXIT
 END DO
 WRITE(*,"(//,5x,'           *** End of Post Processing ! ***',/)")
 CLOSE(16)
 CLOSE(17)
 STOP
 9999 WRITE(*,"(//,7x,'         *** Post Processing interrupted!!!',/)")
 STOP

 END PROGRAM st_2_pos
