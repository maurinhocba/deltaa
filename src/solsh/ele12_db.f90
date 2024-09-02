 MODULE ele12_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn,elastiff,rotor3d
   USE c_input  !All?
  IMPLICIT NONE

   ! Derived type for an ELE12 element
   ! Total Lagrangial  3-D prism Laminated Solid-Shell  element
   INTEGER, PARAMETER :: nstre =  6, &! number of stress measures
                         ngaud =  1, &! number of Gauss point to define volume
                         ngaus =  2   !number of integration points
   TYPE ele12
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER, POINTER  :: lnods(:)  ! Conectivities
     REAL (kind=8) :: angle,            &  ! Euler angle between local t1-t2 and orthotropic
                      dvol,             &  ! Initial volume
                      jacin(2,2,2),     &  ! in-plane inverse jacobian for shear ANS
                      stres(6,2),       &  !Actual 2ndPK stresses at faces
                      alpha,ka             ! EAS variable, stiffness and integrated stresses
     REAL (kind=8), POINTER :: cartd(:),      & ! cartesian der (y3) of Shape Functions at center (EAS)
                               nfdas(:,:,:),  & ! Nodal Function Derivatives at the Assumed Strain points
                               cdq(:,:,:,:),  & ! cartesian derivatives for QUAD approach (AS)
                               area(:,:),     & ! jacobians at mid side points and its sum
                               h(:),          & ! submatrix
                               nangl(:,:),    & ! nodal angles for zigzag functions
                               asazz(:,:),    & ! Assumed strain side derivatives for zigzag functions
                               stint(:),      & !integrated stresses
                               se(:)            !to compute geometric stiffness

     TYPE (ele12), POINTER :: next              !pointer to next element
   END TYPE ele12

   ! Derived type for a set of ELE12 elements
   TYPE ele12_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nnb,   &  ! number of basic nodes per element
                         nnode, &  ! number of nodes per element
                         nstra, &  ! number of shell strain components
                         nreqs, &  ! number of GP for History output
                         narch     ! number of output unit
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: check           ! .TRUE. -> check connectivities orientation
                                ! .FALSE. -> does not check element orientation (Default)
     LOGICAL :: lface           ! .TRUE. -> if extended connectivities have been computed
                                ! .FALSE.
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf     ! Default Euler angles between
                                ! Global X-Y-Z and orthotropic system
     LOGICAL :: quad            ! .TRUE. -> use in-plane quadratic approach
                                ! .FALSE. -> use standard formulation
     INTEGER :: locax           ! local x definition option
     LOGICAL :: zigzag          ! use zigzag approach
     REAL (kind=8) :: beta(2)   ! stabilization factors

     TYPE (ele12), POINTER     :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele12_set), POINTER :: next      !pointer to next set
   END TYPE ele12_set
   TYPE (ele12_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE new_ele12(elset)
   !Create a new element of ELE12 sets

     TYPE(ele12_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nnb   = 6        !     "     number of nodes per elements
     elset%nnode = elset%nnb!     "     number of nodes per elements
     elset%nstra = 8        !     "     number of shell strain components
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%lface = .FALSE.  !     "     flag to generate extended conns
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     elset%quad  = .FALSE.  !     "     to consider quadratic approach
     elset%locax = 3        !     "     local x definition
     elset%zigzag= .FALSE.  !     "     hierarchical in-plane DOFs
     NULLIFY(elset%head,elset%tail,elset%ngrqs)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele12

   SUBROUTINE add_ele12 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele12_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele12

   SUBROUTINE srch_ele12 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele12_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele12

   SUBROUTINE del_ele12 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele12_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele12), POINTER :: ea,ep
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
       CALL del_ele12e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele12

   !----------- Element management routines

   SUBROUTINE ini_ele12e (head, tail)
     !initialize a list of ELE12 elements

     ! dummy arguments
     TYPE (ele12), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele12e

   SUBROUTINE new_ele12e(elm)
   !Create a new element of ELE12 sets

     TYPE(ele12),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvol  = 0d0      !     "     Element volume
     elm%jacin = 0d0      !     "     In-plane Inverse jacobian for shear ANS
     elm%stres = 0d0      !     "     stresses at faces (for output only)
     elm%alpha = 0d0      !     "     EAS parameter
     elm%ka    = 1d0      !     "     EAS parameter stiffness

     NULLIFY(elm%cartd,elm%nfdas,elm%lnods,elm%area,elm%cdq,elm%h,elm%nangl, &
             elm%asazz,elm%stint,elm%se)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele12e

   SUBROUTINE add_ele12e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele12), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele12e

   SUBROUTINE srch_ele12e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele12), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele12e

   SUBROUTINE del_ele12e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele12), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele12), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%lnods)              !deallocate connectivities array
     DEALLOCATE (posic%cdq)                !deallocate in-plane cartesyan derivatives array
     DEALLOCATE (posic%h)                  !deallocate EAS H array
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele12e

   SUBROUTINE cut_ele12e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele12), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele12e

   INCLUDE 'acvd12.fi'
   INCLUDE 'bmat12.fi'
   INCLUDE 'bmat12p.fi'
   INCLUDE 'bmat12q.fi'
   INCLUDE 'bmat12s.fi'
   INCLUDE 'bphi12.fi'
   INCLUDE 'bphi12s.fi'
   INCLUDE 'bsma12.fi'
   INCLUDE 'bsma12s.fi'
   INCLUDE 'comm12.fi'
   INCLUDE 'dump12.fi'
   INCLUDE 'elmd12.fi'
!   INCLUDE 'expo12.fi'
   INCLUDE 'gaus12.fi'
!   INCLUDE 'impo12.fi'
   INCLUDE 'jacob12p.fi'
   INCLUDE 'jacob12s.fi'
   INCLUDE 'kgmm12.fi'
   INCLUDE 'kgmm12q.fi'
   INCLUDE 'kgmm12s.fi'
   INCLUDE 'kgms12p.fi'
!   INCLUDE 'kgms12s.fi'
!   INCLUDE 'kgmt12p.fi'
!   INCLUDE 'kgmt12s.fi'
   INCLUDE 'lcsy12.fi'
   INCLUDE 'load12.fi'
   INCLUDE 'mase12.fi'
   INCLUDE 'masm12.fi'
   INCLUDE 'outp12.fi'
   INCLUDE 'rest12.fi'
   INCLUDE 'resv12p.fi'
   INCLUDE 'resv12s.fi'
   INCLUDE 'resv12si.fi'
   INCLUDE 'stif12p.fi'
   INCLUDE 'stif12s.fi'
   INCLUDE 'stif12si.fi'
   INCLUDE 'toar12.fi'
   INCLUDE 'zigzag_pro.fi'

 END MODULE ele12_db
