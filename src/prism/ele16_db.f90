 MODULE ele16_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input
  IMPLICIT NONE


  REAL(kind=8), PARAMETER :: &
   gpa3(2,3) =  RESHAPE( (/  0.5d0, 0.0d0, 0.5d0, 0.5d0, 0.0d0, 0.5d0 /), (/2,3/)), &
   gpa8(2,8) =  RESHAPE( (/ 0.2113248654052, 0.0000000000000, 0.7886751345948, 0.0000000000000, 0.7886751345948, 0.2113248654052, &
                            0.2113248654052, 0.7886751345948, 0.0000000000000, 0.7886751345948, 0.0000000000000, 0.2113248654052, &
                            0.3333333333333, 0.3333333333333, 0.3333333333333, 0.3333333333333 /),(/2,8/)),             &
   amaq(8,8) = RESHAPE( (/                                                                                              &
 1.3660254037844,-1.7320508075689,-3.7320508075689, 1.7320508075689, 0.0000000000000, 1.3660254037844, 0.0000000000000, &
-2.3660254037844,-0.3660254037844, 1.7320508075689,-0.2679491924311,-1.7320508075689, 0.0000000000000,-0.3660254037844, &
 0.0000000000000,-0.6339745962156, 0.0000000000000, 0.0000000000000, 1.4142135623731,-3.3460652149512, 0.0000000000000, &
-1.4142135623731, 0.0000000000000, 0.8965754721681, 0.0000000000000, 0.0000000000000, 1.4142135623731,-0.8965754721681, &
 0.0000000000000,-1.4142135623731, 0.0000000000000, 3.3460652149512, 0.0000000000000, 0.0000000000000,-0.3660254037844, &
-0.6339745962156,-0.3660254037844,-0.2679491924311, 1.7320508075689,-1.7320508075689, 0.0000000000000, 0.0000000000000, &
 1.3660254037844,-2.3660254037844, 1.3660254037844,-3.7320508075689,-1.7320508075689, 1.7320508075689, 0.0000000000000, &
 0.0000000000000, 6.0000000000000,-3.0000000000000, 0.0000000000000,-3.0000000000000, 0.0000000000000, 6.0000000000000, &
 0.0000000000000, 0.0000000000000,-3.0000000000000, 6.0000000000000, 0.0000000000000, 6.0000000000000, 0.0000000000000, &
-3.0000000000000 /),  (/8,8/) ),                                                                                        &
   amal(6,6) = RESHAPE( (/                                                                                              &
 1.3660254037844,-1.7320508075689,-1.3660254037844, 0.0000000000000,-0.3660254037844, 0.0000000000000,                  &
-0.3660254037844, 1.7320508075689, 0.3660254037844, 0.0000000000000, 1.3660254037844, 0.0000000000000,                  &
 0.0000000000000, 0.0000000000000, 0.5176380902050, 0.0000000000000, 1.9318516525781, 0.0000000000000,                  &
 0.0000000000000, 0.0000000000000,-1.9318516525781, 0.0000000000000,-0.5176380902050, 0.0000000000000,                  &
 0.0000000000000, 0.0000000000000, 1.3660254037844,-0.3660254037844, 0.3660254037844, 1.7320508075689,                  &
 0.0000000000000, 0.0000000000000,-0.3660254037844, 1.3660254037844,-1.3660254037844,-1.7320508075689 /),  (/6,6/) )


   ! Derived type for an ELE16 element
   ! Total Lagrangial (Log strain) 3-D prism element
   INTEGER, PARAMETER :: nvare = 14   !number of internal variables per integration point
   TYPE ele16
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4), POINTER :: lnods(:)  ! Conectivities 6-12-15
     REAL (kind=8) :: angle(3)                  ! Euler angle between XYZ and orthotropic
     REAL (kind=8), POINTER :: dvol(:),       & ! Initial volume
                               cartd(:,:,:),  & !cartesian derivatives of Shape Functions
                               cdq(:,:,:,:),  & !cartesian derivatives for QUAD approach
                               stint(:,:),    & !Actual Kirchhoff stresses
                               gausv(:,:),    & !Gauss-point internal variables
                               nfdas(:,:,:),  & ! Nodal Function Derivatives at the Assumed Strain points
                               jacin(:,:,:),  & ! in-plane inverse jacobian
                               stra0(:,:)
     TYPE (ele16), POINTER :: next              !pointer to next element
   END TYPE ele16

   ! Derived type for a set of ELE16 elements
   TYPE ele16_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nnode, &  ! number of nodes per element 6,12,15,18
                         nreqs, &  ! number of GP for History output
                         narch, &  ! number of output unit
                         nassp, &  ! number of output unit
                         ngaus     !number of integration points
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: small           ! .TRUE. -> use Green strain tensor
                                ! .FALSE. -> Use log strains (Default)
     LOGICAL :: quad            ! .TRUE. -> use QUADratic approximation
                                ! .FALSE.   based on neighbour elements
     LOGICAL :: lface           ! .TRUE. -> if extended connectivities have been computed
                                ! .FALSE.
     LOGICAL :: bezier          ! .TRUE. -> if bezier polynomials will be used
                                ! .FALSE.   Lagrange polynomials wil be used
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf(3)    ! Default Euler angles between
                                  ! Global X-Y-Z and orthotropic system
     LOGICAL :: bbar            ! .TRUE. -> use averge volumetric deformation
                                ! .FALSE. -> use standard formulation
     LOGICAL :: shell           ! .TRUE. -> modify transverse shear
                                ! .FALSE. -> use standard formulation
     INTEGER :: locax           ! local x definition option
     REAL(kind=8), ALLOCATABLE :: gpa(:,:),    & !2,nassp
                                  nfnda(:,:,:),& !nnode,nassp,2
                                  amat(:,:),   & !nassp,nassp
                                  pag(:,:,:)     !2,nassp,ngaus/2
     TYPE (ele16), POINTER    :: head, tail  !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)   !gauss points for output
     TYPE (ele16_set), POINTER :: next       !pointer to next set
   END TYPE ele16_set
   TYPE (ele16_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE ini_ele16 (head, tail)
     !initialize a list of ELE16 sets

     !Dummy arguments
     TYPE (ele16_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele16

   SUBROUTINE add_ele16 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele16_set), POINTER :: new, head, tail

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

     END IF
   END SUBROUTINE add_ele16

   SUBROUTINE srch_ele16 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele16_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele16

   SUBROUTINE del_ele16 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele16_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele16), POINTER :: ea,ep
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
       CALL del_ele16e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele16

   !----------- Element management routines

   SUBROUTINE ini_ele16e (head, tail)
     !initialize a list of ELE16 elements

     ! dummy arguments
     TYPE (ele16), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele16e

   SUBROUTINE new_ele16e(elm)
   !Create a new element of ELE18 sets

     TYPE(ele16),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     NULLIFY(elm%lnods)   !     "     conectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     NULLIFY(elm%dvol,elm%cartd,elm%cdq,elm%stint,elm%gausv,elm%nfdas,elm%jacin)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele16e

   SUBROUTINE add_ele16e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele16), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele16e

   SUBROUTINE srch_ele16e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele16), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele16e

   SUBROUTINE del_ele16e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele16), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele16), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%dvol)               !deallocate vol arrays
     DEALLOCATE (posic%stint)              !deallocate stress array
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic%lnods)              !deallocate fixed space
     DEALLOCATE (posic%cartd)              !deallocate fixed space
     IF(ASSOCIATED(posic%cdq)) DEALLOCATE (posic%cdq)  !deallocate fixed space
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele16e

   SUBROUTINE cut_ele16e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele16), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele16e

   INCLUDE 'actu16.fi'
   INCLUDE 'acvd16.fi'
   INCLUDE 'bbar16.fi'
   INCLUDE 'bmat16.fi'
   INCLUDE 'bmat16q.fi'
   INCLUDE 'bsma16.fi'
   INCLUDE 'comm16.fi'
   INCLUDE 'dump16.fi'
   INCLUDE 'elmd16.fi'
   INCLUDE 'expo16.fi'
   INCLUDE 'gaus16.fi'
   INCLUDE 'impo16.fi'
   INCLUDE 'jacob16.fi'
   INCLUDE 'kgmm16.fi'
   INCLUDE 'kgmm16q.fi'
   INCLUDE 'kgmm16t.fi'
   INCLUDE 'lcsy16.fi'
   INCLUDE 'load16.fi'
   INCLUDE 'mase16.fi'
   INCLUDE 'masm16.fi'
   INCLUDE 'nodx16.fi'
   INCLUDE 'outp16.fi'
   INCLUDE 'rest16.fi'
   INCLUDE 'resv16.fi'
   INCLUDE 'stif16.fi'
   INCLUDE 'stra16.fi'
   INCLUDE 'toar16.fi'

 END MODULE ele16_db
