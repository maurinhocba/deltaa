 SUBROUTINE contol(actio)
 !******************************************************************
 !
 !*** READ control DATA
 !
 !******************************************************************
 USE param_db,ONLY:  mttl
 USE c_input   !INTENT(IN)  : input managment routines and data
 USE dyna_db, ONLY: alpha, beta, gamma, ccm, cck, eigen
 USE outp_db, ONLY: iener, iwrit, ireac, noutd, noutp, energ, postype
 USE solv_db, ONLY: nbuck, ncdis, scdis, nralg, nsymm, arcln, &
                    dlamb, toler, alph1, gamm1, tol1, linear, arcl1, lsearch, neigen, mbite
 USE ctrl_db, ONLY: ndyna, ndime, ndofn, nrotd, ntype, neulr, nstep, nstra, text, dtime, ndoft, &
                    inverse
 USE esets_db, ONLY : nrenu
 USE gvar_db, ONLY : actchk

 IMPLICIT NONE

 !--- Dummy variables
 CHARACTER(len=*),INTENT(IN):: actio
  !--- Local variables
 INTEGER (kind=4) :: diter,every,first,iaux,inode,ittol,karcl, &
                     lauto,nposn,niter
 LOGICAL done(3), lumpd, dampg

 !*** READ the first DATA card, and echo it immediately.
 CALL listen('CONTOL')                           !reads TITLE
 WRITE(lures,"(//,10x,a60,//)") card(1:60)
 text = card(1:MIN0(mttl,LEN_TRIM(card)))         !store title in TEXT

 ! initializes variables

 IF ( TRIM(actio) == 'NEW' )THEN  !for an new problem set default values
   ! main control variables
   ndime = 3 ; neulr = 0 ; ndofn = 3 ; nrotd= ndofn ; nsymm = 0 ; nrenu = 0
   nstra = 1 ; nstep = 0 ; done = .FALSE.
   ! Newton Rapshon type
   first  = 1 ; every  = 1 ; ittol  = 1 ; niter  = 20 ; toler  = 1d-3 ; lsearch = .FALSE.
   ! static parameters
   lauto  = 0 ; karcl  = 0 ; diter  = 4
   ! buckling parameters
   nbuck  = 0  ; alph1 = 1d-3 ; gamm1 = 1d+0 ; tol1 = toler**(1.5d0) ; linear = .FALSE.
   neigen = 1; mbite = 500
   ! dynamic parameters
   ndyna  = 3 ; alpha  = 0.50d0 ; beta  = 0.25d0 ; gamma  = 0.50d0; eigen = .FALSE. ; neigen = 1
   ! output parameters
   iwrit  = 1 ; iener  = 0 ; ireac  = 0  ; noutd = 1 ; noutp = 1 ; inode  = 0; nposn = 0
   ! temperature data
   ndoft  = 0          !no temperature DOFs
   !
   inverse = .FALSE.
 ELSE
   done = .FALSE. ; done(1) = .TRUE.
   first = MOD(nralg,10)
   every = MOD(nralg,100)/10
   ittol = MOD(nralg,1000)/100
   lauto = MOD(nralg,10000)/1000
   karcl = MOD(nralg,100000)/10000
   diter = REAL(INT(MOD(nralg,10000000)/100000),kind=8)    !two digits
   niter = nralg/10000000                                  !two digits
   nposn = MOD(ncdis,10)
   inode = (ncdis-nposn)/10
   nstra = nstra + 1
 END IF

 DO
   CALL listen('CONTOL')
   IF (exists('CONTRO')) THEN
     IF ( done(1) ) CALL runend('Control variables can not be changed')
     WRITE (lures,"(/,'  C O N T R O L   P A R A M E T E R S ',/)")
     ndime=getint('NDIME ',3,' Number of coordinate Dimension ...')
     IF( ndime == 2 )ntype=getint('NTYPE ',3,' Problem type for 2 Dimensional ...')
     neulr=getint('NEULR ',0,' Nodal system code 0:NO 1:YES .....')
     ndofn=getint('NDOFN ',ndime,' Number of Degrees Of Freedom .....')
     nrenu=getint('NRENU ',0,' Renumbering flag .................')
     IF(ndime == 3) THEN
       neulr=9*neulr
       nrotd = MIN(ndofn,6)
     ELSE
       nrotd = MIN(ndofn,3)
     END IF
     ndoft=getint('NDOFT ',1,' Number of Temperature DOFs..........')
     inverse=exists('INVERS')
     IF( inverse )THEN
       WRITE(lures,"('   Inverse problem is active ')")
       nsymm=getint('NSYMM ',1,' Symmetric or Non-Symmetric matrix.')
     END IF
     done(1)=.TRUE.

   ELSE IF (exists('ADVANC')) THEN
     IF( TRIM(actio) /= 'NEW' ) nralg = nralg-first-10*every-100*ittol-10000000*niter
     IF(done(3))CALL runend('CONTOL: both TIME & STATIC not all.')
     WRITE (lures,"(/,'  A D V A N C E   P A R A M E T E R S ',/)")
     nstep=getint('NSTEP ',0,'!Number of steps in the analysis...')
     first=getint('FIRST ',first,' First iteration for K evaluation..')
     every=getint('EVERY ',every,' Frequency of K evaluation ........')
     niter=getint('NITER ',niter,' Maximum No of iteration in a step.')
     ittol=getint('ITTOL ',ittol,' Convergence Norm to be used ......')
     toler=getrea('TOLER ',toler,' Tolerance in convergence norm ....')
     lsearch=exists('LINESE')
     IF( lsearch ) WRITE(lures,"('   Line search is active ')")
     IF(niter <= 0) niter = 3
     IF(niter > 98)THEN
       niter = 98
       WRITE(lures,"(' NITER modified to Maximum Niter in a step. (98)')")
     END IF
     IF(ittol <= 0 .OR. ittol > 3) ittol = 1
     nralg = nralg+first+10*every+100*ittol+10000000*niter
     inode=getint('NODE  ',inode,' Node for displacement control.....')
     nposn=getint('DOF   ',nposn,' DOF for displacement control......')
     IF( nposn > ndime )THEN
       WRITE(lures,"(' control DOF must be translational ')")
       CALL runend('CONTOL: Set Control DOF <= NDIME   ')
     END IF
     postype = 'L'
     ncdis = nposn + 10*inode
     !IF( ncdis > 0 ) postype = 'D'
     done(2)=.TRUE.

   ELSE IF (exists('TIME  ')) THEN
     lumpd = .FALSE.
     IF(done(3))CALL runend('CONTOL: both TIME & STATIC not all.')
     WRITE (lures,"(/,'  T I M E   P A R A M E T E R S ',/)")
     dtime =getrea('DTIME ',0.d0,  '!Time Increment ...................')
     alpha =getrea('ALPHA ',alpha, ' Alpha coeff. in Newmark Scheme ...')
     beta  =getrea('BETA  ',beta  ,' Beta coeff in Newmark Scheme .....')
     gamma =getrea('GAMMA ',gamma, ' Gamma coeff in Newmark Scheme ....')

     ! check and modify if necessary
     IF(alpha < 0.25d0 .OR. alpha > 1d0  ) alpha = 0.5d0
     IF(beta  < 0.00d0 .OR. beta  > 0.5d0) beta  = 0.25d0
     IF(gamma < 0.00d0 .OR. gamma > 1.0d0) gamma = 0.5d0

     IF( exists('LUMPED') .OR. exists('DAMPIN'))ndyna = 3
     lumpd = exists ('LUMPED')
     IF(lumpd) THEN
       ndyna = ndyna+1
       WRITE(lures,"('  Lumped Mass matrix will be used ')")
     ELSE
       WRITE(lures,"('  Consistent Mass matrix will be used ')")
     END IF
     dampg = exists ('DAMPIN')
     IF(dampg)THEN
       WRITE(lures,"('  Structural Damping will be used ')")
       ndyna = ndyna-2
       ccm =getrea('CMASS ',0d0,'!Mass coeff for Rayleigh Damping...')
       cck =getrea('CSTIF ',0d0,'!Stif coeff for Rayleigh Damping...')
     ELSE
       WRITE(lures,"('  No Damping will be used ')")
     END IF
     IF (exists('EIGEN'))THEN
       eigen = .TRUE.
       neigen =getint('NEIGEN',neigen,' Number of eigenvectors ...........')
       CALL openfi(10,rest='_eigen_mode.res')
       WRITE(10,"('GiD Post Results File 1.0',//)")
     END IF
     IF( inverse ) WRITE(lures,"('   Inverse problem De-activated for dynamic problems ')")
     inverse = .FALSE.
     done(3) = .TRUE.

   ELSE IF (exists('STATIC')) THEN
     IF( TRIM(actio) /= 'NEW' ) nralg = nralg-1000*lauto-10000*karcl-100000*diter
     IF(done(3))CALL runend('CONTOL: both TIME & STATIC not all.')
     WRITE (lures,"(/,'  S T A T I C   P A R A M E T E R S ',/)")
     lauto=getint('LAUTO ',lauto,' Automatic modification of params..')
     diter=getint('DITER ',diter,' Desired number of iterations......')
     nbuck=getint('NBUCK ',nbuck,' Frequency of linear buckling anal.')
     IF( nbuck /= 0 )THEN
       alph1 =getrea('ALPHA ',alph1,' Precision Parameter ..............')
       gamm1 =getrea('GAMMA ',gamm1,' Penalty Parameter ................')
       tol1  =getrea('TOLER ',tol1,' Convergence Tolerance ............')
       linear = exists('LINEAR')
       arcl1 = getrea('SWITCH',0d0,'!Increment in most relevant. DOF ..')
       neigen =getint('NEIGEN',neigen,' Number of bukling loads ..........')
       mbite  =getint('MBITE ',Mbite ,' Max No of interation for mode ....')
       IF( exists('NODE  ') .AND. exists('DOF   ') )THEN
         inode =getint('NODE  ',1,'!Node for displacement control.....')
         nposn =getint('DOF   ',1,'!DOF for displacement control......')
         scdis = nposn + 10*inode
       END IF
       CALL openfi(10,rest='_buck_mode.res')
     END IF
     IF(exists('KARCL ',iaux))THEN
       IF( INT(param(iaux)) == 0)THEN
         SELECT CASE (words(iaux+1))
         CASE ('LOAD  ')
           karcl = 0
           WRITE(lures,"(10x,'Load prescribed')")
         CASE ('TANGEN')
           karcl = 1
           WRITE(lures,"(10x,'Tangent plane method')")
         CASE ('UPDATE')
           karcl = 2
           WRITE(lures,"(10x,'Updated Tangent plane method')")
         CASE ('ARCLEN')
           karcl = 3
           WRITE(lures,"(10x,'Arc-length prescribed')")
         CASE ('DISPLC')
           karcl = 4
           WRITE(lures,"(10x,'Displacement prescribed')")
         CASE ('TOOLS ')
           karcl = 5
           WRITE(lures,"(10x,'Tools movement prescribed')")
         CASE DEFAULT
           karcl = 0
           WRITE(lures,"(10x,'Load prescribed')")
         END SELECT
       ELSE
         karcl=getint('KARCL ',0,'!Type of continuation method ......')
       END IF
     ELSE
       karcl=getint('KARCL ',0,'!Type of continuation method ......')
     END IF
     IF(karcl == 0)THEN
       dlamb =getrea('DLAMB ',0d0,'!Initial Increment in load paramet.')
     ELSE
       arcln =getrea('ARCLN ',0d0,'!Initial increment in arc length...')
       IF(arcln == 0d0 .OR. NBUCK /= 0)THEN
         dlamb =getrea('DLAMB ',0d0,'!Initial Increment in load paramet.')
       ELSE
         dlamb =getrea('DLAMB ',0d0,' Initial Increment in load paramet.')
       END IF
     END IF
     IF(lauto < 0 .OR. lauto > 3) lauto = 0
     IF(karcl < 0 .OR. karcl > 5) karcl = 0
     IF(diter <= 0) diter = 3
     IF(diter > 99)THEN
       niter = 99
       WRITE(lures,"(' DITER modified to Maximum Diter in a step. (99)')")
     END IF
     nralg = nralg+1000*lauto+10000*karcl+100000*diter
     ndyna = 0
     done(3)=.TRUE.

   ELSE IF (exists('OUTPUT')) THEN
     WRITE (lures,"(/,'  O U T P U T   P A R A M E T E R S ',/)")
     noutd=getint('NOUTD ',noutd,' History Output Interval  .........')
     noutp=getint('NOUTP ',noutp,' Complete Output IntervaL .........')
     iwrit=getint('IWRIT ',iwrit,' Print-out code -      0:NO 1:Yes .')
     iener=getint('IENER ',iener,' Calculation of Energy 0:NO 1:Yes .')
     ireac=getint('REACT ',ireac,' Print nodal Reactions 0:NO 1:Yes .')
     ALLOCATE( energ(6*iener) )

   ELSE IF (exists('ENDCON') ) THEN
     IF(.NOT.done(1))CALL runend('CONTOL:  Specify CONTROL Parameters')
     IF(.NOT.done(2))CALL runend('CONTOL:  Specify ADVANCE Parameters')
     IF(.NOT.done(3))CALL runend('CONTOL:  Specify TIME or STATIC Par')
     EXIT
   ELSE IF(ALL(done(1:3)))THEN
     IF(.NOT.exists('ENDCON')) backs = .TRUE.
     EXIT
   ELSE
     CALL runend('CONTOL: Error in CONTROL Parameters')
   END IF
 END DO
 IF(actchk) iwrit = 0

 RETURN
 END SUBROUTINE contol
