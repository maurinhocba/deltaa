  MODULE ele03_db
     ! TLLL shell triangle based on Simo's theory
     USE param_db,ONLY: mnam,midn,mlin
     USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,snn
     USE c_input
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      ngaus = 1, &  ! number of GP
      ngamm = 3     ! number of shear strain values

   REAL (kind=8), PARAMETER :: ap1t(2,ngamm) =     &   !shear interpolation matrix
                  RESHAPE((/ -1d0/3d0,  1d0/3d0,   &
                             -1d0/3d0, -2d0/3d0,   &
                              2d0/3d0,  1d0/3d0 /), (/2,ngamm/))


   !   hh= side element connectivities
   INTEGER(kind=4), PARAMETER :: hh(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )
   REAL(kind=8), PARAMETER :: alp1 = 0.70d0, & !Max angle to apply membrane cuadratic approach
                              alp2 = 1.00d0, & !angle to use standard CST for membrane
                              alp3 = 0.3d0     !alp2-alp1 auxiliar value

   TYPE ele03
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4), POINTER :: lnods(:)   ! Conectivities
     REAL (kind=8) :: angle                  ! angle of local system
     REAL (kind=8) :: dvolu,               & ! Gauss associated area
                      cartd(3,2),          & ! Cartesian derivatives
                      jacin(2,2),          & ! Jacobian inverse
                      stra0(6),            & ! initial strains
                      gamm0(ngamm),        & ! initial shear strains
                      qbar(ngamm),         & ! present equivalent shears
                      ambda(2)               ! thickness ratio
     REAL (kind=8), POINTER :: & !plastictiy variables et al
                      nangl(:,:),    & !(2,nside)       Cos - Sin of nodal angles
                      jas(:,:),      & !(2,ngamm)   local direction at assumed shear strain points
                      strsg(:),      & ! present forces and moments
                      cd(:,:,:),     & !(4,2,3) cartesian derivatives at mid-side points
                      ehist(:),      & !(5)       first set
                      strap(:,:),    & !(nstre,2) second set
                      stres(:,:)       !(5,nlayr) third set
     LOGICAL :: mems(3)                         !use quadratic approach for that side
     !ehist = (/ epstr-actual,eqstr,epstr-lastconverged,consist-param(2) /)
     TYPE (ele03), POINTER :: next              !pointer to next element
   END TYPE ele03

    ! Derived type for a set of SHELT elements
    TYPE ele03_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) :: nelem, & ! number of elements
                          nstre = 8, &  ! number of stresses/GP
                          nnode, & ! number of nodes per element
                          nreqs, & ! number of GP for hist. output
                          narch    ! number of output unit
      LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                 !  defined or not updated
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: quad            ! .TRUE. -> quadratic approach for membrane
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      INTEGER :: locax           ! local x definition option
      REAL (kind=8) ::  angdf    ! angle between X_1 and orthotropic dir 1
      TYPE (ele03), POINTER    :: head, tail !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      REAL (kind=8) ::  stabq    ! stabilization factor for shear
      TYPE (ele03_set), POINTER :: next
    END TYPE ele03_set
    TYPE (ele03_set), POINTER :: head
    TYPE (ele03_set), POINTER :: tail

  CONTAINS
    SUBROUTINE ini_ele03 (head, tail)
      !initialize a list of SHELT sets

      TYPE (ele03_set), POINTER :: head, tail

      NULLIFY (head, tail)

    END SUBROUTINE ini_ele03

    SUBROUTINE add_ele03 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele03_set), POINTER :: new, head, tail

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
    END SUBROUTINE add_ele03

    SUBROUTINE srch_ele03 (head, anter, posic, name, found)
      !This subroutine searches for a set named "name"
      !Dummy arguments
      LOGICAL :: found
      CHARACTER (len=*) :: name ! set name
      TYPE (ele03_set), POINTER :: head, anter, posic

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
    END SUBROUTINE srch_ele03

    SUBROUTINE del_ele03 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele03_set), POINTER :: head, anter, posic

      TYPE (ele03), POINTER :: ea,ep
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
        CALL del_ele03e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO

     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele03


   SUBROUTINE ini_ele03e (head, tail)
     !initialize a list of ELE03 elements

     TYPE (ele03), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele03e

   SUBROUTINE add_ele03e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele03), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele03e

   SUBROUTINE srch_ele03e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele03), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele03e

   SUBROUTINE del_ele03e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele03), POINTER :: head, tail, anter, posic
     TYPE (ele03), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     IF(ASSOCIATED(posic%nangl))DEALLOCATE (posic%nangl)
     IF(ASSOCIATED(posic%jas  ))DEALLOCATE (posic%jas)
     IF(ASSOCIATED(posic%strsg))DEALLOCATE (posic%strsg)
     IF(ASSOCIATED(posic%ehist))DEALLOCATE (posic%ehist,posic%strap)
     IF(ASSOCIATED(posic%stres))DEALLOCATE (posic%stres)
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele03e

   SUBROUTINE cut_ele03e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele03), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele03e

   SUBROUTINE new_ele03e(elm)
   !Create a new element of ELE07 sets

     TYPE(ele03),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvolu = 0d0      !     "     initial area
     elm%cartd = 0d0      !     "     cartesyan derivatives
     elm%jacin = 0d0      !     "     inverse jacobians
     elm%stra0 = 0d0      !     "     initial strains
     elm%gamm0 = 0d0      !     "     initial distorsions
     elm%qbar  = 0d0      !     "     equivalent shears
     elm%ambda = 0d0      !     "     thickness ratio
     NULLIFY(elm%ehist,elm%strap,elm%stres,elm%strsg,elm%nangl,elm%jas,elm%lnods)
     NULLIFY(elm%next)
   END SUBROUTINE new_ele03e


   INCLUDE 'actu03.fi'
   INCLUDE 'acvd03.fi'
   INCLUDE 'axep03.fi'
   INCLUDE 'bmat03.fi'
   INCLUDE 'bshe03.fi'
   INCLUDE 'btdbpr2.fi'
   INCLUDE 'comm03.fi'
   INCLUDE 'comp_ang03.fi'
   INCLUDE 'dmat03.fi'
   INCLUDE 'dump03.fi'
   INCLUDE 'elmd03.fi'
   !INCLUDE 'expo03.fi'
   INCLUDE 'gaus03.fi'
   INCLUDE 'impo03.fi'
   INCLUDE 'intr03.fi'
   INCLUDE 'kgeo03.fi'
   INCLUDE 'kgsh03.fi'
   INCLUDE 'load03.fi'
   INCLUDE 'mase03.fi'
   INCLUDE 'masm03.fi'
   !INCLUDE 'nods03.fi'
   INCLUDE 'outp03.fi'
   INCLUDE 'rest03.fi'
   INCLUDE 'resv03.fi'
   !INCLUDE 'secd03.fi'
   INCLUDE 'setg03.fi'
   INCLUDE 'stif03.fi'
   INCLUDE 'stra03.fi'
   INCLUDE 'tanm03.fi'
   INCLUDE 'toar03.fi'
  END MODULE ele03_db
