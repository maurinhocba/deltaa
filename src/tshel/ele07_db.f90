  MODULE ele07_db
     USE param_db,ONLY: mnam,midn,mlin
     USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv
     USE c_input
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      nnode = 6, &  ! number of nodes per element
      ngaus = 3, &  ! number of GP
      nnas1 = 6, &  ! number of shear assumed strain points (lineal OZTS)
      nnas2 = 8, &  ! number of shear assumed strain points (lineal + Quad)
      nasmm = 9     ! number of membrane assumed strain values and points

!  parameters for Gauss point weights
   REAL (kind=8), PARAMETER :: osix =  0.1666666666667D0
   REAL (kind=8), PARAMETER :: weigp(ngaus) = osix   !Gauss point weights

!  PARAMETERS FOR TRANSVERSE SHEAR FORMULATION
   REAL (kind=8), PARAMETER :: am1(nnas1,nnas1) =  RESHAPE ((/  &        ! A-matrix (Lineal)
   1.3660254037844, -1.7320508075689, -1.3660254037844,  0.0000000000000, -0.3660254037844, 0.0000000000000, &
  -0.3660254037844,  1.7320508075689,  0.3660254037844,  0.0000000000000,  1.3660254037844, 0.0000000000000, &
   0.0000000000000,  0.0000000000000,  0.5176380902050,  0.0000000000000,  1.9318516525781, 0.0000000000000, &
   0.0000000000000,  0.0000000000000, -1.9318516525781,  0.0000000000000, -0.5176380902050, 0.0000000000000, &
   0.0000000000000,  0.0000000000000,  1.3660254037844, -0.3660254037844,  0.3660254037844, 1.7320508075689, &
   0.0000000000000,  0.0000000000000, -0.3660254037844,  1.3660254037844, -1.3660254037844,-1.7320508075689 /),(/nnas1,nnas1/))

   REAL (kind=8), PARAMETER :: am2(nnas2,nnas2) =  RESHAPE ((/  &        ! A-matrix (BATHE et al)
 1.3660254037844,-1.7320508075689,-3.7320508075689, 1.7320508075689, 0.0000000000000, 1.3660254037844, 0.0000000000000,-2.3660254037844, &
-0.3660254037844, 1.7320508075689,-0.2679491924311,-1.7320508075689, 0.0000000000000,-0.3660254037844, 0.0000000000000,-0.6339745962156, &
 0.0000000000000, 0.0000000000000, 1.4142135623731,-3.3460652149512, 0.0000000000000,-1.4142135623731, 0.0000000000000, 0.8965754721681, &
 0.0000000000000, 0.0000000000000, 1.4142135623731,-0.8965754721681, 0.0000000000000,-1.4142135623731, 0.0000000000000, 3.3460652149512, &
 0.0000000000000, 0.0000000000000,-0.3660254037844,-0.6339745962156,-0.3660254037844,-0.2679491924311, 1.7320508075689,-1.7320508075689, &
 0.0000000000000, 0.0000000000000, 1.3660254037844,-2.3660254037844, 1.3660254037844,-3.7320508075689,-1.7320508075689, 1.7320508075689, &
 0.0000000000000, 0.0000000000000, 6.0000000000000,-3.0000000000000, 0.0000000000000,-3.0000000000000, 0.0000000000000, 6.0000000000000, &
 0.0000000000000, 0.0000000000000,-3.0000000000000, 6.0000000000000, 0.0000000000000, 6.0000000000000, 0.0000000000000,-3.0000000000000  &
   /),(/nnas2,nnas2/))

   REAL (kind=8), PARAMETER :: ntan(nnode,nnas2) =  RESHAPE ((/  &        ! Nodal function derivatives at the assumed strain points
   -1.5773502691896, 0.4226497308104, 0.0000000000000, 1.1547005383793, 0.0000000000000, 0.0000000000000, &
   -0.4226497308104, 1.5773502691896, 0.0000000000000,-1.1547005383793, 0.0000000000000, 0.0000000000000, &
    0.0000000000000,-1.1153550716504, 0.2988584907227, 0.0000000000000, 0.8164965809277, 0.0000000000000, &
    0.0000000000000,-0.2988584907227, 1.1153550716504, 0.0000000000000,-0.8164965809277, 0.0000000000000, &
   -0.4226497308104, 0.0000000000000, 1.5773502691896, 0.0000000000000, 0.0000000000000,-1.1547005383793, &
   -1.5773502691896, 0.0000000000000, 0.4226497308104, 0.0000000000000, 0.0000000000000, 1.1547005383793, &
   -0.6666666666667, 0.6666666666667, 0.0000000000000, 0.0000000000000, 0.6666666666667,-0.6666666666667, &
   -0.6666666666667, 0.0000000000000, 0.6666666666667,-0.6666666666667, 0.6666666666667, 0.0000000000000  /),(/nnode,nnas2/))

   REAL (kind=8), PARAMETER :: nsha(nnode,nnas2) =  RESHAPE ((/  &        ! Nodal functions at the assumed strain points
    0.6220084679281, 0.0446581987385, 0.0000000000000, 0.3333333333333, 0.0000000000000, 0.0000000000000, &
    0.0446581987385, 0.6220084679281, 0.0000000000000, 0.3333333333333, 0.0000000000000, 0.0000000000000, &
    0.0000000000000, 0.6220084679281, 0.0446581987385, 0.0000000000000, 0.3333333333333, 0.0000000000000, &
    0.0000000000000, 0.0446581987385, 0.6220084679281, 0.0000000000000, 0.3333333333333, 0.0000000000000, &
    0.0446581987385, 0.0000000000000, 0.6220084679281, 0.0000000000000, 0.0000000000000, 0.3333333333333, &
    0.6220084679281, 0.0000000000000, 0.0446581987385, 0.0000000000000, 0.0000000000000, 0.3333333333333, &
    0.1111111111111, 0.1111111111111, 0.1111111111111, 0.2222222222222, 0.2222222222222, 0.2222222222222, &
    0.1111111111111, 0.1111111111111, 0.1111111111111, 0.2222222222222, 0.2222222222222, 0.2222222222222  /),(/nnode,nnas2/))


! PARAMETERS FOR MEMBRANE FORMULATION
! For model 1: Sampling points at vertex nodes
      REAL (kind=8), PARAMETER :: dn(nnode,2,3) = RESHAPE( (/    &         !nodal shape functions derivatives (ANS for membrane)
          -2.0D0,  0.0D0,  0.0D0,   2.0D0,  0.0D0,   0.0D0,      &
          -2.0D0,  0.0D0,  0.0D0,   0.0D0,  0.0D0,   2.0D0,      &
           0.0D0,  2.0D0,  0.0D0,  -2.0D0,  0.0D0,   0.0D0,      &
           0.0D0,  0.0D0,  0.0D0,  -2.0D0,  2.0D0,   0.0D0,      &
           0.0D0,  0.0D0,  0.0D0,   0.0D0,  2.0D0,  -2.0D0,      &
           0.0D0,  0.0D0,  2.0D0,   0.0D0,  0.0D0,  -2.0D0 /),(/nnode,2,3/) )

! For model 2: Sampling points at mid-side of each sub-triangle
!  Gradient computation at sampling points
   REAL (kind=8), PARAMETER :: ntan2(nnode,nasmm) =  RESHAPE ((/  &        !
   -1.5000000000000,  0.5000000000000, 0.0000000000000,  1.0000000000000,  0.0000000000000,  0.0000000000000,   &
   -0.5000000000000,  1.5000000000000, 0.0000000000000, -1.0000000000000,  0.0000000000000,  0.0000000000000,   &
   -0.5000000000000,  0.5000000000000, 0.0000000000000,  0.0000000000000,  1.0000000000000, -1.0000000000000,   &
   -1.5000000000000,  0.0000000000000, 0.5000000000000,  0.0000000000000,  0.0000000000000,  1.0000000000000,   &
   -0.5000000000000,  0.0000000000000, 1.5000000000000,  0.0000000000000,  0.0000000000000, -1.0000000000000,   &
   -0.5000000000000,  0.0000000000000, 0.5000000000000, -1.0000000000000,  1.0000000000000,  0.0000000000000,   &
    0.0000000000000, -1.0606601717798, 0.3535533905933,  0.0000000000000,  0.7071067811865,  0.0000000000000,   &
    0.0000000000000, -0.3535533905933, 1.0606601717798,  0.0000000000000, -0.7071067811865,  0.0000000000000,   &
    0.0000000000000, -0.3535533905933, 0.3535533905933, -0.7071067811865,  0.0000000000000,  0.7071067811865 /),(/nnode,nasmm/))
! this matrix express Beta parameters in terms of the natural (tangent) strains A-matrix
   REAL(kind=8), PARAMETER :: amat2(nasmm,nasmm) = RESHAPE ((/                                  &
     1.50D0,  -2.00D0,   -2.00D0,    0.00D0,    0.00D0,    0.00D0,    0.75D0,   -1.00D0,   -1.00D0, &
    -0.50D0,   2.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,   -0.25D0,    1.00D0,    0.00D0, &
     0.00D0,   0.00D0,    2.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,    1.00D0, &
     0.00D0,   0.00D0,    0.00D0,    1.50D0,   -2.00D0,   -2.00D0,    0.75D0,   -1.00D0,   -1.00D0, &
     0.00D0,   0.00D0,    0.00D0,   -0.50D0,    0.00D0,    2.00D0,   -0.25D0,    0.00D0,    1.00D0, &
     0.00D0,   0.00D0,    0.00D0,    0.00D0,    2.00D0,    0.00D0,    0.00D0,    1.00D0,    0.00D0, &
     0.00D0,   0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.50D0,   -2.00D0,    0.00D0, &
     0.00D0,   0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.50D0,    0.00D0,   -2.00D0, &
     0.00D0,   0.00D0,    0.00D0,    0.00D0,    0.00D0,    0.00D0,   -2.00D0,    2.00D0,    2.00D0 /),(/nasmm,nasmm/))


   TYPE ele07
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8) :: angle                  ! angle of local system
     REAL (kind=8) :: dvolu(ngaus),        & ! Gauss associated area
                      cartd(nnode,2,ngaus),& ! Cartesian derivatives
                      jacin(2,2,ngaus),    & ! Jacobian inverse
                      stra0(8,ngaus),      & ! initial strains
                      ambda(2,ngaus)         ! thickness ratio
     REAL (kind=8), POINTER :: beta(:)
     REAL (kind=8), POINTER :: & !output & plasticity variables
                      qbar(:),         & ! equivalent shears
                      nangl(:,:),      & !(2,nnode)   Cos - Sin of nodal angles
                      strsg(:,:),      & !(nstre,ngaus)   present forces and moments
                      ehist(:,:),      & !(5,ngaus)       first set
                      strap(:,:,:),    & !(nstre,ngaus,2) second set
                      stres(:,:,:)       !(5,nlayr,ngaus) third set
     !ehist = (/ epstr-actual,eqstr,epstr-lastconverged,consist-param(2) /)
     TYPE (ele07), POINTER :: next              !pointer to next element
   END TYPE ele07

    ! Derived type for a set of SHELT elements
    TYPE ele07_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) :: nstre = 8, &  ! number of stresses/GP
                          nelem, & ! number of elements
                          nnass, & ! formulation type for shear 1:2
                          ansmm, & ! formulation type 0:2
                          nreqs, & ! number of GP for hist. output
                          narch    ! number of output unit
      LOGICAL :: gpint           ! .FALSE. : mid-side points  .TRUE. interior points
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      INTEGER :: locax           ! local x definition option
      REAL (kind=8) :: angdf     ! angle between X_1 and orthotropic dir 1
      REAL (kind=8) :: posgp(2,ngaus)        !Gauss points positions
      REAL (kind=8), POINTER :: ap1(:,:,:)   !(2,nnass,ngaus)shear interpolation matrix
      REAL (kind=8) :: shape(nnode,ngaus)    !nodal shape functions at Gauss points
      REAL (kind=8) :: omat(3,nasmm,2,ngaus) !interpolation matrix for assumed membrane model
      TYPE (ele07), POINTER    :: head, tail !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      TYPE (ele07_set), POINTER :: next
    END TYPE ele07_set
    TYPE (ele07_set), POINTER :: head
    TYPE (ele07_set), POINTER :: tail

  CONTAINS
    SUBROUTINE ini_ele07 (head, tail)
      !initialize a list of SHELT sets

      TYPE (ele07_set), POINTER :: head, tail

      NULLIFY (head, tail)

    END SUBROUTINE ini_ele07

    SUBROUTINE add_ele07 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele07_set), POINTER :: new, head, tail

      !Check if a list is empty
      IF (.NOT. ASSOCIATED (head)) THEN
        !list is empty, start it
        head => new
        tail => new
        NULLIFY (tail%next)

      ELSE
        !add a segment to the list
        tail%next => new
        NULLIFY (new%next)
        tail => new

      ENDIF
    END SUBROUTINE add_ele07

    SUBROUTINE srch_ele07 (head, anter, posic, name, found)
      !This subroutine searches for a set named "name"
      !Dummy arguments
      LOGICAL :: found
      CHARACTER (len=*) :: name ! set name
      TYPE (ele07_set), POINTER :: head, anter, posic

      found = .FALSE.
      NULLIFY (posic,anter)
      !Check if a list is empty
      IF (ASSOCIATED (head)) THEN
        posic => head
        DO
          IF (TRIM(posic%sname) == TRIM(name)) THEN
            found = .TRUE.
            EXIT
          END IF
          IF (ASSOCIATED(posic%next) ) THEN
            anter => posic
            posic => posic%next
          ELSE
            EXIT
          END IF
        END DO
      ENDIF
      IF (.NOT.found) NULLIFY (posic,anter)
    END SUBROUTINE srch_ele07

    SUBROUTINE del_ele07 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele07_set), POINTER :: head, anter, posic

      TYPE (ele07), POINTER :: ea,ep
      INTEGER (kind=4) :: iel

      IF (.NOT.ASSOCIATED (anter)) THEN  !if anter pointer is null => head
        head => posic%next               !point first to next
      ELSE
        anter%next => posic%next         !skip posic
      END IF

      ! deallocation of the set memory is next done
      NULLIFY( ea )                  !nullify previous element in list
      ep => posic%head               !point present element to first
      DO iel = 1,posic%nelem         !for each element in the set
        CALL del_ele07e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO

     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele07


   SUBROUTINE ini_ele07e (head, tail)
     !initialize a list of ELE07 elements

     TYPE (ele07), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele07e

   SUBROUTINE add_ele07e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele07), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN       !list is empty, start it
       head => new                           !first element
       tail => new                           !last element
       NULLIFY (tail%next)                   !last poit to nothing

     ELSE                                    !add a segment to the list
       tail%next => new                      !point to the new element
       NULLIFY (new%next)                    !nothing beyond the last
       tail => new                           !new last element

     END IF
   END SUBROUTINE add_ele07e

   SUBROUTINE srch_ele07e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele07), POINTER :: head, anter, posic

     found = .FALSE.             !initializes flag and pointers
     NULLIFY (posic,anter)

     IF (ASSOCIATED (head)) THEN          !Check if a list is empty
       posic => head                      !begin at top
       DO
         IF(posic%numel == kelem) THEN    !check if label found
           found = .TRUE.                 !Found
           EXIT                           !element found, EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN    !check if there are more els
           anter => posic                     !remember previous element
           posic => posic%next                !new present element
         ELSE
           EXIT                               !no more elements EXIT
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)    !point to nothing
     RETURN
   END SUBROUTINE srch_ele07e

   SUBROUTINE del_ele07e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele07), POINTER :: head, tail, anter, posic
     TYPE (ele07), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     IF(ASSOCIATED(posic%nangl))DEALLOCATE (posic%nangl)
     IF(ASSOCIATED(posic%strsg))DEALLOCATE (posic%strsg)
     IF(ASSOCIATED(posic%ehist))DEALLOCATE (posic%ehist,posic%strap)
     IF(ASSOCIATED(posic%stres))DEALLOCATE (posic%stres)
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele07e

   SUBROUTINE cut_ele07e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele07), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele07e

   SUBROUTINE new_ele07e(elm)
   !Create a new element of ELE07 sets

     TYPE(ele07),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods = 0        !     "     conectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvolu = 0d0      !     "     initial area
     elm%cartd = 0d0      !     "     cartesyan derivatives
     elm%jacin = 0d0      !     "     inverse jacobians
     elm%stra0 = 0d0      !     "     initial strains
     elm%ambda = 0d0      !     "     thickness ratio
     NULLIFY(elm%ehist,elm%strap,elm%stres,elm%strsg,elm%nangl,elm%qbar)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele07e

   INCLUDE 'actua7.fi'
   INCLUDE 'acvdf7.fi'
   INCLUDE 'bmatx7.fi'
   INCLUDE 'bmmt07.fi'
   INCLUDE 'bshem7.fi'
   INCLUDE 'commv7.fi'
   INCLUDE 'comp_ang07.fi'
   !INCLUDE 'dump07.fi'
   INCLUDE 'elmda7.fi'
   INCLUDE 'gauss7.fi'
   INCLUDE 'istgp7.fi'
   INCLUDE 'kgeom7.fi'
   INCLUDE 'kgmt07.fi'
   INCLUDE 'kgshm7.fi'
   INCLUDE 'loadp7.fi'
   INCLUDE 'masel7.fi'
   INCLUDE 'masmt7.fi'
   INCLUDE 'nodxy7.fi'
   INCLUDE 'omat07.fi'
   INCLUDE 'outpu7.fi'
   !INCLUDE 'rest07.fi'
   INCLUDE 'resv07.fi'
   INCLUDE 'shap07.fi'
   INCLUDE 'stiff7.fi'
   INCLUDE 'stran7.fi'
  END MODULE ele07_db
