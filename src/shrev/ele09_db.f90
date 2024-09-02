  MODULE ele09_db
    USE param_db,ONLY: mnam,midn,mlin
    USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,elastiff
    USE c_input
    IMPLICIT NONE

    SAVE

   TYPE ele09
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4), POINTER :: lnods(:)  ! Conectivities
     REAL (kind=8), POINTER :: phil(:), & ! angles of local system
                               jac(:),  & ! Gauss associated length
                               r0(:)   ,& !
                               stra0(:,:),  & ! initial strains
                               stres(:,:),  & ! present forces and moments
                               sede0(:),    & ! (ngaus)initial area ratio
                               sedef(:),    & ! (ngaus)area ratio
                               !plastictiy variables
                               !ehist = (/ epstr-actual,eqstr,epstr-lastconverged,consist-param(2) /)
                               ehist(:,:),   & !(5,ngaus)       first set
                               strap(:,:,:), & !(nstre,ngaus,2) second set
                               lstre(:,:,:), & !(3,nlayr,ngaus) third set
                               !auxiliar variables
                               auxil(:)        !
     INTEGER(kind=4) :: prpro                  ! 0  No print
                                               ! 1-> print stress
                                               ! 2-> print displacement profile
                                               ! 3-> print stress and displacement profile
     TYPE (ele09), POINTER :: next              !pointer to next element
   END TYPE ele09

    ! Derived type for a set of SHREV elements
    TYPE ele09_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) ::   &
                          ndofe, & ! number of DOFs per node
                          nnode, & ! number of nodes per element
                          ngaus, & ! number of GP
                          nstre, & ! number of stresses/GP
                          nelem, & ! number of elements
                          axesc, & ! flag for local system
                          nreqs, & ! number of GP for hist. output
                          narch    ! number of output unit
      REAL (kind=8), POINTER :: posgp(:),   & !Gauss point position
                                shape(:,:), & !shape functions
                                deriv(:,:), & !shape functions derivatives
                                weigh(:)      !gauss points weigths
      REAL (kind=8), POINTER :: estr(:,:)     !Element strains
      INTEGER (kind=4), POINTER :: esecs(:)   !Element secs
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      LOGICAL :: zigzpp = .FALSE.  ! .TRUE. use superimposed zig-z++ function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      TYPE (ele09), POINTER    :: head, tail !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      TYPE (ele09_set), POINTER :: next
    END TYPE ele09_set
    TYPE (ele09_set), POINTER :: head
    TYPE (ele09_set), POINTER :: tail

  CONTAINS
    SUBROUTINE ini_ele09 (head, tail)
      !initialize a list of SHELQ sets

      TYPE (ele09_set), POINTER :: head, tail

      NULLIFY (head, tail)

    END SUBROUTINE ini_ele09

    SUBROUTINE add_ele09 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele09_set), POINTER :: new, head, tail

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
    END SUBROUTINE add_ele09

    SUBROUTINE srch_ele09 (head, anter, posic, name, found)
      !This subroutine searches for a set named "name"
      !Dummy arguments
      LOGICAL :: found
      CHARACTER (len=*) :: name ! set name
      TYPE (ele09_set), POINTER :: head, anter, posic

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
    END SUBROUTINE srch_ele09

    SUBROUTINE del_ele09 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele09_set), POINTER :: head, anter, posic

      TYPE (ele09), POINTER :: ea,ep
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
        CALL del_ele09e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO

     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele09


   SUBROUTINE new_ele09(elset)
   !Create a new element of ELE06 sets

     TYPE(ele09_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%ndofe = 3        !     "     number of dof per node
     elset%nnode = 2        !     "     number of nodes per element
     elset%ngaus = 1        !     "     number of integration points
     elset%nstre = 5        !     "     number of stress measures per gauss point
     elset%nelem = 0        !     "     number of elements
     elset%axesc = 0        !     "     flag for local system (nodal)
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%zigzag= .FALSE.  !     "     no additional in plane dofs
     elset%zigzpp= .FALSE.  !     "     no additional in plane dofs
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     NULLIFY(elset%head,elset%tail,elset%ngrqs)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele09

   SUBROUTINE ini_ele09e (head, tail)
     !initialize a list of ELE09 elements

     TYPE (ele09), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele09e

   SUBROUTINE add_ele09e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele09), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele09e

   SUBROUTINE srch_ele09e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele09), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele09e

   SUBROUTINE del_ele09e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele09), POINTER :: head, tail, anter, posic
     TYPE (ele09), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE( posic%phil,posic%jac,posic%r0,posic%stra0,posic%stres, &
                 posic%sede0,posic%sedef)
     IF(ASSOCIATED(posic%ehist))DEALLOCATE (posic%ehist,posic%strap)
     IF(ASSOCIATED(posic%lstre))DEALLOCATE (posic%lstre)
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele09e

   SUBROUTINE cut_ele09e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele09), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele09e

   INCLUDE 'actua9.fi'
   INCLUDE 'acvdf9.fi'
   INCLUDE 'bmatx9.fi'
   INCLUDE 'btstr9.fi'
   INCLUDE 'commv9.fi'
   INCLUDE 'dump09.fi'
   INCLUDE 'elmda9.fi'
   INCLUDE 'gauss9.fi'
   INCLUDE 'gauss9c.fi'
   INCLUDE 'intrf9.fi'
   INCLUDE 'intrf9z.fi'
   INCLUDE 'istgp9.fi'
   INCLUDE 'istgp9z.fi'
   INCLUDE 'kgeom9.fi'
   INCLUDE 'loadp9.fi'
   INCLUDE 'locla9.fi'
   INCLUDE 'masel9.fi'
   INCLUDE 'masmt9.fi'
   INCLUDE 'nodxy9.fi'
   INCLUDE 'outpu9.fi'
   INCLUDE 'rest09.fi'
   INCLUDE 'resvp9.fi'
   INCLUDE 'resvp9c.fi'
   INCLUDE 'rztmat09n.fi'
   INCLUDE 'rztmat09m.fi'
   INCLUDE 'rztmat09o.fi'
   INCLUDE 'setga9.fi'
   INCLUDE 'stiff9.fi'
   INCLUDE 'stran9.fi'
   INCLUDE 'tanma9.fi'
   INCLUDE 'tanma9c.fi'
   INCLUDE 'tanma9z.fi'
   INCLUDE 'zigzag_pro1.fi'
   INCLUDE 'zigzag_pro2.fi'

  END MODULE ele09_db
