 SUBROUTINE inpd27 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 27: Bezier PRISM solid-shell element

 USE ele27_db
 USE gvar_db,ONLY : fimpo,overw
 USE npo_db, ONLY : cpx,coord

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset,gpint,mmsur,shsur
 INTEGER (kind=4) :: nreqs, narch, ngaus, nelem, ansmm, anssh, easts, nassp, i

 TYPE (ele27_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN

   ! check if list of sets and set exists and initializes
   CALL srch_ele27 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm27 (1, nelem,  nreqs, narch, elsnam, elset, ngaus, ansmm, anssh, easts, nassp)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD27')  !read a line
     ngaus = getint('NGAUS ',2,' NUMBER OF GAUSS POINTS TTT (>= 2) .')
     IF( ngaus < 2 )THEN
       WRITE(lures,"(' Minimumn number of gauss points is 2 ........')",ERR=9999)
       elset%ngaus = 2
     END IF

     elset%bezier =  .TRUE.
     IF(.NOT.exists('BEZIER').AND.exists('STANDA')) elset%bezier = .FALSE.    !use BEZIER prism?
     WRITE(lures,"(/,5X,'Use of Bezier Prism:',l6)") elset%bezier

     elset%locax = 3       !default
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( INT(param(i)) > 0 .AND. INT(param(i)) < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     elset%angdf = getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')

     gpint = exists('GPINT') .OR. .NOT.exists('MSGP') !position of integration points
     IF( gpint ) THEN
       WRITE(lures,"(/,5X,'Internal Gauss points will be used')")
       gpcoo = RESHAPE ( (/ 0.166666666666667, 0.166666666666667, &
                            0.666666666666667, 0.166666666666667, &
                            0.166666666666667, 0.666666666666667  /), (/2,ngaup/))
     ELSE
       WRITE(lures,"(/,5X,'Mid-side Gauss points will be used')")
       gpcoo = RESHAPE ( (/ 0.5D0, 0.0D0, &
                            0.5D0, 0.5D0, &
                            0.0D0, 0.5D0  /), (/2,ngaup/))
     END IF
     wp = 1D0/6D0

     ansmm =  2 !Use mid-side points of each subtriangle
     IF( exists('ANSMM ',i)) THEN   !Assumed Natural Strain model for Mem
       IF( INT(param(i)) == 0 .OR. INT(param(i)) == 2 )THEN  !if a valid code
         ansmm =  INT(param(i))
         WRITE(lures,"(/,5X,'ANS Membrane Model',i3)") ansmm
       ELSE  !keep default value
         WRITE(lures,"(/,5X,'Error in the option for ANS membrane model', / &
                         5X,'Invalid model: ',i3,' Default value used (2)')") INT(param(i))
       END IF
     END IF
     mmsur = EXISTS('MMSURF')   !position on sampling points
     IF( mmsur ) THEN
       WRITE(lures,"(/,5X,'Internal surfaces points will be used for membrane strains')")
       gpzm = (/ -5.77350269189626D-01, 5.77350269189626D-01 /)
     ELSE
       WRITE(lures,"(/,5X,'External surfaces points will be used for membrane strains')")
       gpzm = (/ -1.D0, 1.D0 /)
     END IF

     SELECT CASE (ANSMM)  !membrane formulation
     !CASE (0) ! standard displacement formulation
     !  ansmp(:,1:ngaup) = gpcoo
     CASE (1) ! Sampling points as suggested by Bathe
       ansmp =  RESHAPE ((/  0.25d0, 0.00D0, 0.75d0, 0.00d0, 0.25d0, 0.50d0,                &
                             0.00d0, 0.25d0, 0.00d0, 0.75d0, 0.50d0, 0.25d0,                &
                             0.75d0, 0.25d0, 0.25d0, 0.75d0, 0.25d0, 0.25d0 /),(/ 2,nasmm /) )
     CASE (2) ! Sampling points at mid-side of each sub-triangle
       ansmp =  RESHAPE ((/  2.1132486540519D-01, 0.0000000000000D+00, 7.8867513459481D-01, 0.0000000000000D+00, &
                             2.1132486540519D-01, 5.7735026918963D-01, 0.0000000000000D+00, 2.1132486540519D-01, &
                             0.0000000000000D+00, 7.8867513459481D-01, 5.7735026918963D-01, 2.1132486540519D-01, &
                             7.8867513459481D-01, 2.1132486540519D-01, 2.1132486540519D-01, 7.8867513459481D-01, &
                             2.1132486540519D-01, 2.1132486540519D-01  /),(/ 2,nasmm /) )
     END SELECT

     anssh =  2 !Use incomplete quadratic
     nassp = nnas2   !number of sampling points
     IF( exists('ANSSH ',i)) THEN   !Assumed Natural Strain model for shear
       IF( INT(param(i)) >= 0 .AND. INT(param(i)) < 3 )THEN  !if a valid code
         anssh =  INT(param(i))
         WRITE(lures,"(/,5X,'ANS Shear Model',i3)") ansmm
         IF( anssh == 0 )nassp = ngaup*2 !????
         IF( anssh == 1 )nassp = nnas1
       ELSE  !keep default value
         WRITE(lures,"(/,5X,'Error in the option for ANS Shear model', / &
                         5X,'Invalid model: ',i3,' Default value used (2)')") INT(param(i))
       END IF
     END IF
     IF( anssh > 0 )ALLOCATE( ams(nassp,nassp) )

     shsur = EXISTS('SHSURF')   !position on sampling points
     IF( shsur ) THEN
       WRITE(lures,"(/,5X,'Internal surfaces points will be used for shear strains')")
       gpzs = (/ -5.77350269189626D-01, 5.77350269189626D-01 /)
     ELSE
       WRITE(lures,"(/,5X,'External surfaces points will be used for shea strains')")
       gpzs = (/ -1.D0, 1.D0 /)
     END IF

     ! Gauss points position of Assumed Shear Strain Sampling Points
     anssp= RESHAPE((/  0.2113248654052D0, 0.0000000000000D0, &
                        0.7886751345948D0, 0.0000000000000D0, &
                        0.7886751345948D0, 0.2113248654052D0, &
                        0.2113248654052D0, 0.7886751345948D0, &
                        0.0000000000000D0, 0.7886751345948D0, &
                        0.0000000000000D0, 0.2113248654052D0, &
                        0.3333333333333D0, 0.3333333333333D0, &
                        0.3333333333333D0, 0.3333333333333D0  /),(/2,nnas2/))


     easts =  2 !Use EAS model
     IF( exists('EASTS ',i)) THEN   !Enhanced Assumed Strain model for Transverse Strain
       IF( INT(param(i)) >= 0 .AND. INT(param(i)) <= 3 )THEN  !if a valid code
         easts =  INT(param(i))
         WRITE(lures,"(/,5X,'EAS Transverse strain Model',i3)") easts
       ELSE  !keep default value
         WRITE(lures,"(/,5X,'Error in the option for EAS transverse strain model', / &
                         5X,'Invalid model: ',i3,' Default value used (1)')") INT(param(i))
       END IF
     END IF
     IF( easts == 0 )THEN
       gpzt = (/ -1.D0, 1.D0 /)
     ELSE IF (easts == 3 )THEN
       gpzt = (/ -5.77350269189626D-01, 5.77350269189626D-01 /)
       easts = 0
     ELSE
       gpzt = 0d0
     END IF


     IF( exists('QUAD   ') .OR. .NOT.exists('6NODE')) THEN   !Post-process option
       elset%quad    =  .TRUE.
       WRITE(lures,"(/,5X,'Set will be post-precessed as 15-node Prism')")
     ELSE  !keep default
       elset%quad    =  .FALSE.
       WRITE(lures,"(/,5X,'Set will be post-precessed as a 6-node prism')")
     END IF

     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%small = exists('SMALL') .AND. .NOT.exists('LARGE')
     IF(elset%small) WRITE(lures,"(' Green strains will be used if possible')")
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele27e (elset%head, elset%tail)
   END IF
   !  read new data or add data to previous data
   CALL elmd27(nelem, elset%head, elset%tail, iwrit, ngaus, ansmm, anssh, easts, nassp )
   ! compute coordinates of control points
   CALL nodx27(nelem,elset%head,coord,cpx,elset%bezier)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele27 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

! ELSE IF (TRIM(task) == 'RESTAR') THEN
!
!   ALLOCATE (elset)          !initializes a list
!   NULLIFY(elset%head)       !nullify head pointer
!   ! read control parameters
!   READ (51) nelem, nreqs, narch, ngaus, &
!             elset%gauss, elset%small, elset%plstr, elset%angdf, &
!             elset%locax
!   ! restore list of elements
!   CALL rest27 (nelem, nreqs, elset%head, elset%tail, &
!                elset%ngrqs, ngaus )
!   ! add to list of elements
!   CALL add_ele27 (elset, head, tail)
!
! ELSE IF (TRIM(task) == 'IMPORT') THEN
!   ! check if list of sets and set exists and initializes
!   IF( overw )THEN
!     CALL srch_ele27 (head, anter, elset, elsnam, oldset)
!     IF (oldset) THEN    !if set exists
!       CALL comm27(1, nelem,  nreqs, narch, elsnam, elset, ngaus)
!     ELSE
!       CALL runen2(' Old set to overwrite does not exist')
!     END IF
!     READ (fimpo)
!   ELSE
!     ALLOCATE (elset)       !reserve memory for set
!     CALL ini_ele27e (elset%head, elset%tail)
!     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
!     READ (fimpo) nelem,ngaus,elset%angdf,elset%small, &
!                  elset%locax
!     nreqs = 0
!     narch = 0
!   END IF
!   ! restore list of elements
!   CALL impo27 ( nelem, ngaus, elset%head, elset%tail)
!   ! add to list of elements
!   IF( .NOT.overw )THEN
!     CALL add_ele27 (elset, head, tail)
!     nelms = nelms + 1 ! increased set counter for this element type
!   END IF
!
 ELSE
   CALL runend('INPD27: NON-EXISTENT TASK .        ')
 END IF
 CALL comm27(0, nelem, nreqs, narch, elsnam, elset, ngaus, ansmm, anssh, easts, nassp)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd27
