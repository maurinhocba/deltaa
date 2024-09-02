 MODULE ele30_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,pmater,snn
   USE c_input
   IMPLICIT NONE

   ! Derived type for an ELE30 element
   ! Butterfly Bezier Basic thin Shell Triangle (BBST)

   ! Reference: To be developed

   INTEGER (kind=4), PARAMETER :: hh(2,3) = RESHAPE( &
                     [ 3,2, 1,3, 2,1 ], [ 2,3 ] ), &    !nodes of each side
                                  kk(3,3) = RESHAPE( &
                     [ 4,3,2, 5,1,3, 6,2,1 ], [3,3] ) !side element connectivities
   REAL (kind=8) :: ar(12,6),  areg(36,18)              !A matrix at regular elements
   REAL (kind=8) :: mb(12,6,6),abou(36,18,6)            !A matrix at boundaries
   LOGICAL :: modif(3,6)=RESHAPE((/.TRUE., .TRUE., .TRUE., &
                                   .TRUE., .TRUE., .TRUE., &
                                   .TRUE., .TRUE.,.FALSE., &
                                   .TRUE., .TRUE.,.FALSE., &
                                   .TRUE., .TRUE.,.FALSE., &
                                   .TRUE.,.FALSE.,.FALSE. /),(/3,6/))

   TYPE ele30
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(12) ! Conectivities
     INTEGER (kind=4) :: lside(3)  ! Neighbour elements, used to compute
                                   ! transvere shear forces deriving moments
     INTEGER (kind=4) :: bcode     ! boundary case
     REAL (kind=8) :: Area1,     & ! Initial area
                      lb,        & ! Thickness ratio
                      angle,     & ! angle between dir 1 and orthotropic dir 1
                      a(1:3),    & ! side projections in dir 1
                      b(1:3),    & ! side projections in dir 2
                      cdn(6,3),  & ! side normal derivatives
                      stra0(3),  & ! Original curvature tensor
                      stra1(6)     ! Actual mid-surface metric tensor and curvatures
     REAL (kind=8), POINTER :: nab(:,:)       !normal at boundary
     REAL (kind=8), POINTER :: gausv(:,:)     !layer values
     TYPE (ele30), POINTER :: neig1,neig2,neig3 !pointer to neighbour elements
     TYPE (ele30), POINTER :: next            !pointer to next element
   END TYPE ele30

   ! Derived type for a set of ELE30 elements
   TYPE ele30_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem  ! number of elements
     INTEGER (kind=4) :: nreqs  ! number of GP for hist. output
     INTEGER (kind=4) :: narch  ! number of output unit
     LOGICAL :: logst           ! use logarithmic strain
     LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                !  defined or not updated
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) ::  angdf     ! angle between X_1 and orthotropic dir 1
     REAL (kind=8), POINTER :: stint(:,:)   !moments, forces and shears
     TYPE (ele30), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele30_set), POINTER :: next      !pointer to next set
   END TYPE ele30_set
   TYPE (ele30_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS
   SUBROUTINE new_ele30(elset)
   !Create a new element of ELE30 sets

     TYPE(ele30_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%logst = .FALSE.  !     "     use logarithmic strain
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele30

   SUBROUTINE add_ele30 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele30_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele30

   SUBROUTINE srch_ele30 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele30_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF(TRIM(posic%sname) == TRIM(name)) THEN    !check name
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
   END SUBROUTINE srch_ele30

   SUBROUTINE del_ele30 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele30_set), POINTER :: head, anter, posic

     TYPE (ele30), POINTER :: ea,ep
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
       CALL del_ele30e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele30

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE new_ele30e(elm)
   !Create a new element of ELE30 sets

     TYPE(ele30),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:12) = 0  !     "     conectivities
     elm%lside(1:3) = 0   !     "     neighbour elements
     elm%bcode = 0        !     "     to inner all-regular nodes
     elm%area1 = 0d0      !     "     initial area
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%lb = 0d0         !     "     initial area
     elm%a = 0d0          !     "     side projections in dir 1
     elm%b = 0d0          !     "       "       "       "  "  2
     elm%cdn = 0d0        !     "     side normal derivatives
     elm%stra0(1:3) = 0d0 !     "     original curvature tensor
     elm%stra1(1:6) = (/ 1d0,1d0,0d0,0d0,0d0,0d0 /) !     "     actual mid-surface metric tensor
     NULLIFY(elm%gausv)
     NULLIFY(elm%next,elm%neig1,elm%neig2,elm%neig3,elm%nab)

   RETURN
   END SUBROUTINE new_ele30e

   SUBROUTINE add_ele30e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele30), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele30e

   SUBROUTINE ins_ele30e(new, top, tail)
   !This subroutine is simmilar to add_ele, but insert a member list between two pointers

     !--- Dummy variables
     TYPE(ele30),POINTER:: new, top, tail

     IF (ASSOCIATED(top%next)) THEN
       new%next => top%next
       top%next => new
     ELSE
       NULLIFY(new%next)
       top%next => new
       !It should update the tail of the list if necessary
       tail => new
     ENDIF

   RETURN
   END SUBROUTINE ins_ele30e


   SUBROUTINE srch_ele30e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele30), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele30e

   SUBROUTINE del_ele30e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele30), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     IF( .NOT.ASSOCIATED(posic%next) )tail => anter !last element in list
     IF(ASSOCIATED(posic%gausv)) DEALLOCATE (posic%gausv) !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     IF( ASSOCIATED( anter) )posic => anter%next                   !point to next element
     RETURN
   END SUBROUTINE del_ele30e


!  INCLUDE 'Actu30.fi'
   INCLUDE 'Acvd30.fi'
   INCLUDE 'Axep30.fi'
   INCLUDE 'Bfle30.fi'
   INCLUDE 'bmem30.fi'
   INCLUDE 'Comm30.fi'
!  INCLUDE 'Dump30.fi'
   INCLUDE 'Elmd30.fi'
!  INCLUDE 'expo30.fi'
   INCLUDE 'Gaus30.fi'
!  INCLUDE 'impo30.fi'
   INCLUDE 'kgmm30.fi'
   INCLUDE 'load30.fi'
   INCLUDE 'mase30.fi'
   INCLUDE 'Masm30.fi'
   INCLUDE 'modvec30.fi'
   INCLUDE 'modmat30.fi'
!  INCLUDE 'nods30.fi'
   INCLUDE 'Outp30.fi'
!  INCLUDE 'Rest30.fi'
   INCLUDE 'Resv30.fi'
!  INCLUDE 'secd30.fi'
   INCLUDE 'Shap30.fi'
   INCLUDE 'Stif30.fi'
   INCLUDE 'stra30.fi'
!  INCLUDE 'Stst30.fi'
   INCLUDE 'Top_array30.fi'
 END MODULE ele30_db
