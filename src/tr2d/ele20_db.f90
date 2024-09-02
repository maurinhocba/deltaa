MODULE ele20_db
  USE param_db,ONLY: mnam,midn
  USE c_input
  USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,snn
  IMPLICIT NONE

  ! Derived type for an ELE20 element

  ! Total Lagrangial (Log strain) 2-D triangular element

  !  Reference:  F.Flores "A two-dimensional linear assumed strain
  !              triangular element for finite deformation analysis",
  !              ASME Journal of Applied Mechanics

  INTEGER, PARAMETER :: nnode = 3, & !number of nodes per element
                        ngaus = 1    !number of integration points

  INTEGER(kind=4), PARAMETER  :: kk(2,3) = RESHAPE((/ 3,2, 1,3, 2,1 /), &
                                                    (/2,3/) )

  TYPE ele20
    INTEGER (kind=4) :: numel  ! label of element
    INTEGER (kind=4) :: matno  ! Material number
    INTEGER (kind=4) :: lnods(6)  ! Conectivities
    REAL (kind=8) :: Area1,       & ! Initial area
                     angle,       & ! angle between dir 1 and orthotropic dir 1
                     cd(6,3),     & ! cartesian derivatives
                     stint(4)       ! Actual stresses
    REAL (kind=8), POINTER :: gausv(:)       !Gauss-point internal variables
    TYPE (ele20), POINTER :: next            !pointer to next element
  END TYPE ele20

  ! Derived type for a set of ELE20 elements
  TYPE ele20_set
    CHARACTER (len=mnam) :: sname ! set name
    INTEGER (kind=4) :: nelem, &  ! number of elements
                        nreqs, &  ! number of GP for hist. output
                        narch     ! number of output unit

    LOGICAL :: gauss, &        ! .FALSE. -> Initial constants not defined or not updated
               lside, &        ! .FALSE. -> patch node not computed yet
               eulrf, &        ! .TRUE. for Eulerian Formulation (Garino J2 algorithm)
               swapc           ! .FALSE. not swap the corner alone elements (default is TRUE)

    INTEGER :: plstr           ! compute Plastic Strain Flag
        ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
    REAL (kind=8) :: angdf     ! angle between X_1 and orthotropic dir 1
    TYPE (ele20), POINTER    :: head, tail !pointer to first and last elm.
    INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
    TYPE (ele20_set), POINTER :: next      !pointer to next set
  END TYPE ele20_set

  TYPE(ele20_set),POINTER,SAVE:: head=>NULL(), tail=>NULL()  !first and last elements sets

 CONTAINS

  SUBROUTINE new_ele20(elms)
  !Create a new ELE20 sets of the list

    !--- Dummy variables
    TYPE(ele20_set),POINTER:: elms

    ALLOCATE(elms)   !Create ele20_set element
    elms%sname = ''        !Initialise set name
    elms%nelem = 0         !     "     number of elements
    elms%nreqs = 0         !     "     number of GP for hist. output
    elms%narch = 0         !     "     number of output unit
    elms%gauss = .FALSE.   !     "     flag for initial constants defined or updated
    elms%lside = .FALSE.   !     "     flag to compute Gauss constants
    elms%plstr = 0         !     "     compute Plastic Strain Flag
    elms%angdf = 0d0       !     "     angle between X_1 and orthotropic dir 1
    NULLIFY(elms%head, elms%tail, elms%ngrqs)  !Initialises pointer
    NULLIFY(elms%next)   !Initialises pointer to next element of the list

  RETURN
  END SUBROUTINE new_ele20

  SUBROUTINE add_ele20(new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ele20_set), POINTER :: new, head, tail

    !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add an element to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    ENDIF
  END SUBROUTINE add_ele20

  SUBROUTINE srch_ele20(head, anter, posic, name, found)
    !This subroutine searches for a set named "name"
    !Dummy arguments
    LOGICAL :: found
    CHARACTER (len=*) :: name ! set name
    TYPE (ele20_set), POINTER :: head, anter, posic

    found = .FALSE.                     !initializes flag
    NULLIFY (posic,anter)               !initializes pointers
    !Check if a list is empty
    IF (ASSOCIATED (head)) THEN         !if there are sets
      posic => head                     !point to first set
      DO
        IF (TRIM(posic%sname) == TRIM(name)) THEN    !check name
          found = .TRUE.                !set flag to .TRUE.
          EXIT                          !O.K. exit search
        END IF
        IF (ASSOCIATED(posic%next) ) THEN   !there is a next set
          anter => posic                    !point anter to present
          posic => posic%next               !point present to next
        ELSE
          EXIT                              !list exhausted, Exit search
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)   !set not found, null pointers
  END SUBROUTINE srch_ele20

  SUBROUTINE del_ele20(head, anter, posic)
    !This subroutine deletes a set pointed with posic

    TYPE (ele20_set), POINTER :: head, anter, posic

    TYPE (ele20), POINTER :: ea,ep
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
      CALL del_ele20e (posic%head,posic%tail, ea, ep )  !deletes element
    END DO

    NULLIFY (posic,anter)          !point to nothing
  END SUBROUTINE del_ele20

  ! ******* functions for a list of elements in a set ********

  SUBROUTINE new_ele20e(elm)
  !Create a new ELE20 sets of the list

    !--- Dummy variables
    TYPE(ele20),POINTER:: elm

    ALLOCATE(elm)   !Create ele20 element
    elm%numel = 0     !Initialise label of element
    elm%matno = 0     !     "     material number
    elm%lnods = 0     !     "     conectivities
    elm%Area1 = 0d0   !     "     initial area
    elm%angle = 0d0   !     "     angle between dir 1 and orthotropic dir 1
    elm%cd = 0d0      !     "     cartesian derivatives
    elm%stint = 0d0   !     "     actual stresses
    NULLIFY(elm%gausv)  !Initialises pointer
    NULLIFY(elm%next)   !Initialises pointer to next element of the list

  RETURN
  END SUBROUTINE new_ele20e

  SUBROUTINE add_ele20e(new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ele20), POINTER :: new, head, tail

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
  END SUBROUTINE add_ele20e

  SUBROUTINE srch_ele20e(head, anter, posic, kelem, found)
    !This subroutine searches for an element labeled "kelem"
    !Dummy arguments
    LOGICAL :: found
    INTEGER (kind=4) :: kelem
    TYPE (ele20), POINTER :: head, anter, posic

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
  END SUBROUTINE srch_ele20e

  SUBROUTINE del_ele20e(head, tail, anter, posic)
    !This subroutine deletes element pointed with posic

    TYPE (ele20), POINTER :: head, tail, anter, posic
    TYPE (ele20), POINTER :: e

    IF (.NOT.ASSOCIATED (anter)) THEN    !
      head => posic%next
    ELSE
      anter%next => posic%next
    END IF
    e => posic%next                       !keep pointer to next element
    IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
    IF( ASSOCIATED( posic%gausv ) )DEALLOCATE (posic%gausv)                    !deallocate internal variables
    !IF( ASSOCIATED( posic%naxis ) )DEALLOCATE( posic%naxis )
    DEALLOCATE (posic)                    !deallocate fixed space
    posic => e                            !point to next element
    ! NULLIFY (posic,anter)
    RETURN
  END SUBROUTINE del_ele20e

  SUBROUTINE cut_ele20e(head, anter, posic)
    !This subroutine deletes a set pointed with posic
    ! without nullifying anter    ???? what for ????
    TYPE (ele20), POINTER :: head, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    NULLIFY (posic)
  END SUBROUTINE cut_ele20e

  SUBROUTINE delete_ele20e(head, tail)
  !This subroutine deletes element pointed with posic
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(ele20),POINTER:: head, tail
    !--- Local variables
    TYPE(ele20),POINTER:: elm

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      elm => head%next                  !keep pointer to next element
      DEALLOCATE(head%gausv)            !deallocate internal variables
      DEALLOCATE(head)                  !deallocate element
      head => elm                       !point to next element
    END DO
  RETURN
  END SUBROUTINE delete_ele20e

   INCLUDE 'Actu20.fi'
   INCLUDE 'Acvd20.fi'
   INCLUDE 'Axep20.fi'
   INCLUDE 'Bmat20.fi'
   INCLUDE 'Cdac20.fi'
   INCLUDE 'Comm20.fi'
   INCLUDE 'Dump20.fi'
   INCLUDE 'Elmd20.fi'
   INCLUDE 'Expo20.fi'
   INCLUDE 'Gaus20.fi'
   INCLUDE 'impo20.fi'
   INCLUDE 'Kgmm20.fi'
   INCLUDE 'Load20.fi'
   INCLUDE 'Mase20.fi'
   INCLUDE 'Masm20.fi'
   INCLUDE 'nods20.fi'
   INCLUDE 'Outp20.fi'
   INCLUDE 'Rest20.fi'
   INCLUDE 'Resv20.fi'
   INCLUDE 'secd20.fi'
   INCLUDE 'Stif20.fi'
   INCLUDE 'Stra20.fi'
   INCLUDE 'Toar20.fi'

END MODULE ele20_db

