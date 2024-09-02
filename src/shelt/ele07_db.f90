  MODULE ele07_db
     USE param_db,ONLY: mnam,midn,mlin
     USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv
     USE c_input
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      nnode = 6, &  ! number of nodes per element
      ngaus = 3, &  ! number of GP
      ngamm = 6, &  ! number of shear strain values
      nnass = 6, &  ! number of assumed strain points
      nbeta = 9     !

              !                             r3   = SQRT(3)
              !    r3p1 =(r3+1d0)/2d0       r3m1 = (r3-1d0)/2d0
              !    r3p3 =(r3+3d0)/2d0       r3m3 = (r3-3d0)/2d0
              !    fxa  = 1/6+1/SQRT(12)    fxb  = 1/6-1/SQRT(12)    fxc = 2/3
              !    fxd  = 1+2/SQRT(3)       fxe  =-1+2/SQRT(3)       fxf = 4/SQRT(3)
     REAL (kind=8), PARAMETER ::            r3   = 1.73205080756890d0, &
                 r3p1 = 1.36602540378445d0, r3m1 = 0.36602540378445d0, &
                 r3p3 = 2.36602540378445d0, r3m3 =-0.63397459621555d0, &
                 fxa=  0.455341801261480d0, fxb= -0.122008467928146d0, &
                 fxc=  0.666666666666667d0, fxd=  2.154700538379250d0, &
                 fxe=  0.154700538379252d0, fxf=  2.309401076758500d0
     INTEGER( kind=4), PARAMETER :: kk(3,6) =(/ 1,4,2, 2,4,1, 2,5,3, 3,5,2, 3,6,1, 1,6,3 /)
     REAL( kind=8), PARAMETER :: nf(3) = (/ fxa, fxc, fxb /), nd(3) = (/ -fxd, fxf,-fxe /)

    SAVE

   REAL (kind=8) :: posgp(2,ngaus),     & !Gauss point position
                    weigp(ngaus),       & !Gauss point weigths
                    shape(nnode,ngaus), & !nodal shape functions at Gauss points
                    ap1(2,ngamm,ngaus), & !shear interpolation matrix
                    nfdas(nnode,2,nnass), & !nodal shape functions derivatives (ANS for Shear)
                    dn(nnode,2,3,2)         !nodal shape functions derivatives   (ANS for membrane)

   TYPE ele07
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8) :: angle                  ! angle of local system
     REAL (kind=8) :: dvolu(ngaus),        & ! Gauss associated area
                      cartd(nnode,2,ngaus),& ! Cartesian derivatives
                      jacin(2,2,ngaus),    & ! Jacobian inverse
                      stra0(6,ngaus),      & ! initial strains
                      gamm0(ngamm),        & ! initial shear strains
                      qbar(ngamm),         & ! present equivalent shears
                      ambda(2,ngaus)         ! thickness ratio
     REAL (kind=8), POINTER :: beta(:)
     REAL (kind=8), POINTER :: & !output & plasticity variables
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
                          stype, & ! formulation type 0:3
                          nreqs, & ! number of GP for hist. output
                          narch    ! number of output unit
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      INTEGER :: locax           ! local x definition option
      REAL (kind=8) ::  angdf    ! angle between X_1 and orthotropic dir 1
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
     elm%gamm0 = 0d0      !     "     initial distorsions
     elm%qbar  = 0d0      !     "     equivalent shears
     elm%ambda = 0d0      !     "     thickness ratio
     NULLIFY(elm%ehist,elm%strap,elm%stres,elm%strsg,elm%nangl)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele07e

   INCLUDE 'actua7.fi'
   INCLUDE 'acvdf7.fi'
   INCLUDE 'ap1tm7.fi'
   INCLUDE 'asstr7.fi'
   INCLUDE 'bmatx7.fi'
   INCLUDE 'bmmt27.fi'
   INCLUDE 'bshem7.fi'
   INCLUDE 'commv7.fi'
   INCLUDE 'comp_ang07.fi'
   INCLUDE 'concar.fi'
   INCLUDE 'dump07.fi'
   INCLUDE 'elmda7.fi'
   INCLUDE 'gauss7.fi'
   INCLUDE 'intem7.fi'
   INCLUDE 'intrf7.fi'
   INCLUDE 'kgeom7.fi'
   INCLUDE 'kgmt27.fi'
   INCLUDE 'kgshm7.fi'
   INCLUDE 'loadp7.fi'
   INCLUDE 'masel7.fi'
   INCLUDE 'masmt7.fi'
   INCLUDE 'mbmgp7.fi'
   INCLUDE 'nodxy7.fi'
   INCLUDE 'outpu7.fi'
   INCLUDE 'rest07.fi'
   INCLUDE 'resvp7.fi'
   INCLUDE 'setga7.fi'
   !INCLUDE 'shape7.fi'
   INCLUDE 'stiff7.fi'
   INCLUDE 'stran7.fi'
   INCLUDE 'tanma7.fi'
  END MODULE ele07_db
