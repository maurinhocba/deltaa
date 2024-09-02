 MODULE ele04_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,mater,postv
   USE c_input
   IMPLICIT NONE

   ! Derived type for an ELE04 element

   ! Total Lagrangial (Log strain) 3-D tetrahedral element

   !  Reference:  F.Flores  (unpublished)

   INTEGER, PARAMETER :: nnode = 4, & !number of nodes per element
                         ngaus = 1, & !number of integration points
                         nvare = 14   !number of variables per integration point

   INTEGER(kind=4), PARAMETER  :: kk(3,4) = RESHAPE((/ 2,3,4, 1,4,3, &
                     1,2,4,   2,1,3 /),  (/3,4/) )

 REAL (kind=8), PARAMETER :: nd0(4,3) =      & !natural derivatives of 4-node element
                (/ -1d0, 1d0, 0d0, 0d0,      &
                   -1d0, 0d0, 1d0, 0d0,      &
                   -1d0, 0d0, 0d0, 1d0/)
!                (/ 1d0, 0d0, 0d0,-1d0,      &
!                   0d0, 1d0, 0d0,-1d0,      &
!                   0d0, 0d0, 1d0,-1d0/)


   TYPE ele04
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(8)  ! Conectivities
     REAL (kind=8) :: dvol,        & ! Initial volume
                      angle(3),    & ! angle between dir 1 and orthotropic dir 1
                      cd(4,3),     & ! cartesian derivatives
                      facs(0:4),   & ! Original metric tensor
                      stint(6)       ! Actual stresses
     REAL (kind=8), POINTER :: gausv(:)       !Gauss-point internal variables
     TYPE (ele04), POINTER :: next            !pointer to next element
   END TYPE ele04

   ! Derived type for a set of ELE04 elements
   TYPE ele04_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for hist. output
                         narch     ! number of output unit
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: linear          ! .TRUE. -> use linear approach
     LOGICAL :: lside           ! .FALSE. -> patch node not computed yet
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf(3)  ! angle between X_1 and orthotropic dir 1
     TYPE (ele04), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele04_set), POINTER :: next      !pointer to next set
   END TYPE ele04_set
   TYPE (ele04_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS
   SUBROUTINE ini_ele04 (head, tail)
     !initialize a list of ELE04 sets

     TYPE (ele04_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele04

   SUBROUTINE add_ele04 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele04_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele04

   SUBROUTINE srch_ele04 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=mnam) :: name ! set name
     TYPE (ele04_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF(posic%sname == name) THEN    !check name
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
   END SUBROUTINE srch_ele04

   SUBROUTINE del_ele04 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele04_set), POINTER :: head, anter, posic

     TYPE (ele04), POINTER :: ea,ep
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
       CALL del_ele04e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele04

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele04e (head, tail)
     !initialize a list of ELE04 elements

     TYPE (ele04), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele04e

   SUBROUTINE add_ele04e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele04), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele04e

   SUBROUTINE srch_ele04e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele04), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele04e

   SUBROUTINE del_ele04e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele04), POINTER :: head, tail, anter, posic
     TYPE (ele04), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%gausv)                    !deallocate internal variables
     !IF( ASSOCIATED( posic%naxis ) )DEALLOCATE( posic%naxis )
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele04e

   SUBROUTINE cut_ele04e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele04), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele04e

   !INCLUDE 'Actu04.fi'
   INCLUDE 'Acvd04.fi'
   INCLUDE 'Axep04.fi'
   INCLUDE 'Bmat04.fi'
   INCLUDE 'Comm04.fi'
   !INCLUDE 'Dump04.fi'
   INCLUDE 'Elmd04.fi'
   INCLUDE 'Gaus04.fi'
   !INCLUDE 'Kgmm04.fi'
   INCLUDE 'Load04.fi'
   INCLUDE 'Mase04.fi'
   INCLUDE 'Masm04.fi'
   INCLUDE 'Outp04.fi'
   !INCLUDE 'Rest04.fi'
   INCLUDE 'resv04.fi'
   INCLUDE 'Stif04.fi'
   !INCLUDE 'Stra04.fi'
   INCLUDE 'Toar04.fi'

 END MODULE ele04_db
