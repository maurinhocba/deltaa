     ! Get material and section parameters
     !
     ! Uses the following variables

     sec => psecs(isec)%p             !point to section
     nlayr = sec%iprop(1)             !number of layers
     nvar  = sec%iprop(2)             !number of internal variables per layer
     IF( nvar > 0 )THEN
       ALLOCATE(varin(nvar))            !auxiliar internal variables
       varin = 0d0                      !initializes
     END IF
     secty = sec%secty                !section type
     thick = sec%rprop(1)             !original thickness
     minstr= sec%rprop(4)             !strain threshold to use TTTI
     osec = isec                      !keep present section
     IF( secty == 12 )THEN            !standard solid section
       mat => sec%mtbas                 !point to associated material
       shell = nlayr > 1                !membrane or shell
       mtype = mat%mtype                !type of base material
       natst = logst .OR. mtype == 6    !use log strains
       elast = mat%matdef(3) == 1       !elastic
       CALL gaussq(nlayr,thf(1),wei(1)) !integration points through the thickness
       thf(1:nlayr) = thf(1:nlayr)/2d0  !positions
       wei(1:nlayr) = wei(1:nlayr)/2d0  !weights
       alpha = mat%prope(6)             !Thermical Dilatation Coeff
       !IF( .NOT.itemp ) alpha = 0d0
       dm      = 0d0                    !initializes integrated elasticity matrix
       dm((/1,2,7,12,16,17,19,21/)) = sec%rprop(5:12)        !integrated elasticity matrix

       IF( mtype == 1)THEN              !standard isotropic material
         IF( elast )THEN
           plast = .FALSE.
         ELSE   !for elasto-plastic mats
           ! e1, nu1, uniaxial, efren, consn, r, exponent m, hill 79
           propi(1:4) = mat%propp(1:4)       ! isotropic hardening parameters
           propi(5) = REAL( mat%matdef(4),8) ! isotropic hardening model
           chi    = mat%propp(16:27)         ! hill coefficients
           !deatht = mat%propp(5)             !end of plasticity
           !IF (mat%matdef(4) == 5 )THEN
           !  val => mat%chead%val
           !  numpt = mat%chead%np
           !ELSE
           !  NULLIFY (val)
           !  numpt = 0
           !END IF
           plast = propi(1) > 0 ! .AND. death > ttime !consider plasticity ?
         END IF
         c(1:4) = mat%prope(7:10)          ! plane stress elasticity matrix
         CALL dmat14(c(1),propi(1),chi(1),dummy,dummy,dummy,dummy,newmt)
       ELSE IF( mtype == 5)THEN            ! orhthotropic material
         c(1:4) = mat%prope(16:19)         ! plane stress elasticity matrix
         elast = .TRUE.
         plast = .FALSE.                        !consider plasticity ?
         !IF( .NOT.elast )THEN              ! if elasto-plastic
         !  propi(1:13) = mat%propp(17:29)    !orthotropic hardening parameters
         !  deatht = mat%propp(5)             !end of plasticity
         !  plast = deatht > ttime .AND. propi(1) > 0  !consider plasticity ?
         !END IF
         d = 0d0; d(1,1) = c(1); d(1,2) = c(2); d(2,2) = c(3); d(3,3) = c(4)
         newmt = .FALSE.                        !same material than previous ?
       ELSE IF( mtype == 6)THEN            !hyperelastic isotropic
         chi(1:12) = mat%prope(7:18)       !material properties
         elast = .TRUE.                    !elastic only
         plast = .FALSE.                   !consider plasticity ?
       !ELSE IF( mtype == 30)THEN             !user defined material
       !  plast = .TRUE.   !default
       END IF

     ELSE  !secty == 13  composite laminae
       elast = .TRUE. !elast = sec%iprop(3) == 0               !elastic problem
       plast = .FALSE.!plast = .NOT.elast                      !consider plasticity ?
       shell = .TRUE.
       dm(1:15) = sec%rprop((/ 6:8, 18:20, 9:10, 21:23, 11, 24:26 /))    !linear membrane elastic integrated matrix
       dm(16:21) = sec%rprop(12:17)            !linear bending elastic integrated matrix
       !dm = sec%rprop(6:26)                 !linear elastic integrated matrix
       !IF( plast )THEN  !plastic lamina
       !  IF( ALLOCATED ( cm ) )DEALLOCATE( volfr,cm,prop,rr )
       !  ALLOCATE( volfr(nlayr),cm(4,nlayr),prop(17,nlayr),rr(5,nlayr))
       !  i = 27                                !pointer to volume fractions
       !  CALL vecasi(nlayr,sec%rprop(i),volfr)        !volume fractions
       !  i = i+2*nlayr                         !pointer to C-matrices
       !  CALL vecasi(nlayr*4,sec%rprop(i),cm(1,1))    !elasticity matrix
       !  i = i+nlayr*4                         !pointer to R-matrices
       !  CALL vecasi(nlayr*5,sec%rprop(i),rr(1,1))    !rotation matrix
       !  i = i+nlayr*5                         !pointer to plastic properties
       !  CALL vecasi(nlayr*17,sec%rprop(i),prop(1,1)) !plastic properties
       !  oldm = -1                             !recompute constants
       !END IF
       newmt = .FALSE.                        !same material than previous ?
       write(55,"(6f15.8,/,15x,5f15.8,/,30x,4f15.8,/,45x,3f15.8,/,  &
                &          60x,2f15.8,/,75x,f15.8,/)")dm
     END IF
