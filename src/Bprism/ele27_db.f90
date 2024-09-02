! que cosas puedo probar para tratar de darme cuenta que pasa:
! usar 4 puntos de integración en el plano (esto no me parece)
! usar los puntos de muestreo membranales de Bathe et al
! evaluar las deformaciones membranales en el interior del elemento y no en las superficies externas
! Son el problema los puntos de control ?
! ahora todo esto se tendría que poder evaluar en PRISM y luego ver que pasa aquí
!
! OPCIONES DEL ELEMENTO
!
! Posición de los puntos de integración: GPINT
!      para 3 puntos de integración
!   A- a la mitad de cada lado  DEFAULT
!   B- puntos internos
!      gpcoo(2,ngaup),weigp(ngaup)
!
! funciones de interpolación: BEZIER/STANDAR
!   A- polinomios de Bezier DEFAULT
!   B- polinomios cuadráticos estándar
!      shape(nnod1,ngaup),deriv(nnod1,3,ngaup)
!      notar que en la dirección 3, depende de la posición de cada superficie de muestreo
!      para obtener las funciones hay que multiplicar por L1 y L2
!      para obtener las derivadas en el plano hay que multiplicar por L1 y L2
!      para obtener las derivadas en dicha dirección hay que multiplicar por -1/2 y +1/2
!
! parte membranal ANSMM: Se evalúa en dos superficies y se interpola en el interior
!   estas superficies pueden ser:  MMSURF  gpm(2)   valuado en INPD27
!      -las superficies externas zeta = +/- 1         DEFAULT
!      -las superficies internas zeta = +/- 1/raiz(3)
!   en dichas superficies:
!   0 :standard en desplazamientos, evaluar el tensor métrico en el plano directamente
!      para ello usando shape y deriv hay que calcular y guardar IPCDM(nnode,2,nface,ngaup) Evaluada en JACOB27
!   1 :usar los puntos de Kin y Bathe para MITC
!   2 :usar los puntos a la mitad de cada lado de cada subtriángulo para ANS
!      MPG(2,nasmm)  incluye los puntos de evaluación del gradiente  valuado en INPD27
!      NTAN2(nnode,nasmm,nface)  incluye las derivadas naturales para evaluar el gradiente en las direcciones correspondientes
!      requieren invertir A y usar las interpolaciones a los puntos de Gauss, hecho en GAUSS 27
!      PA2(3,nasmm,ngaup)    Para las componentes de C
!      PA2B(3,nasmm,ngaup)   Para las componentes de delta C (variaciones)
!      JACIM(2,2,3,2) se requiere para ANSMM /= 0 en los superficies de evaluación
!
! parte corte Tra ANSSH: Se evalúa en dos superficies y se interpola en el interior
!   estas superficies pueden ser: SHSURF   GPSH(2)  valuada en INPD27
!      -las superficies externas zeta = +/- 1
!      -las superficies internas zeta = +/- 1/raiz(3)
!   en dichas superficies:
!   0 :standard en desplazamientos, evaluar las componentes cartesianas directamente
!      para ello usando shape y deriv hay que calcular y guardar
!      CARTD(nnode,3,ngaup,nface)  valuada en JACOB27
!      IPCDM(nnode,2,nface,ngaup) Evaluada en JACOB27 hay que distinguir con ANSMM = 0!!!
!   1 :usar los puntos de Kin y Bathe para MITC
!   2 :usar los puntos de Oñate et al para ANS
!     ambas opciones en GPS(3,nnas2) valuadas en INPD27
!                       NTAN(nnode,nnas2,nface) valuadas en GAUS27
!     NFDAS(nnode,nnass,nface)  valuadas en JACOB27 a partir de DERS2
!     JACIS(2,2,3,2) se requiere para ANSSH /= 0 en los superficies de evaluación,
!
! deformación axial transversal EASTS:
!    : evaluar en dos superficies que pueden ser
!   0: -las superficies externas zeta = +/- 1
!   3: -las superficies internas zeta = +/- 1/raiz(3)
!   1: evaluar sólo en la superficie media
!   2: evaluar en la superficie media + EAS
!      NF3(nnode,faces,ngaup) contiene las derivadas para evalur f3 en las dos superficies o en la media
!

 MODULE ele27_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input  !All?
  IMPLICIT NONE

   ! Total Lagrangial (Log strain) 3-D Bezier prism Solid-Shell element
   INTEGER, PARAMETER :: nvare = 14, &!number of internal variables per integration point
                         nstre =  6, &! number of stress measures
                         nnode = 12, &! number of nodes per element
                         ndofe = 36, &! number of DOFs per element
                         nasmm =  9, &! number of Assumed Strain for membrane
                         nnas1 =  6, &! number of shear assumed strain points (lineal OZTS)
                         nnas2 =  8, &! number of shear assumed strain points (lineal + Quad)
                         nface =  2, &! number of faces
                         ngaup =  3, &! number of in-plane Gauss points
                         ngaud =  3   ! number of Gauss point to define volume

!  gauss points positions for jacobian computation
   REAL(kind=8), PARAMETER :: zface = 0.7071067811865  ! 5.77350269189626D-01  ! 1d0
   REAL(kind=8), PARAMETER :: gpzv(ngaud) = (/  -zface , +zface , 0d0 /)
!  In-plane Gauss points (xita,eta) and weights
  REAL(kind=8) :: gpcoo(2,ngaup),wp(ngaup)
!  Position of Surfaces for membrane, shear and transvers evaluation
   REAL(kind=8) :: gpzm(nface),  & ! Gauss Point z for membrane components
                   gpzs(nface),  & ! Gauss Point z for shear components
                   gpzt(nface)     ! Gauss Point z for transverse components

! PARAMETERS FOR MEMBRANE FORMULATION
!  Gradient computation at sampling points
   REAL (kind=8) :: ansmp(2,nasmm)              !sampling points positions
   REAL (kind=8) :: mdtan(nnode,nasmm,nface)    !natural derivative (tangent) at sampling points

! interpolates the assumed strain values to Gauss points
   REAL (kind=8) :: pagm(3,nasmm,ngaup)          ! PA matrix from sampling to integration points (C-components)
   REAL (kind=8) :: pagm2(3,nasmm,ngaup)         ! PA matrix from sampling to integration points (delta C-components)

! PARAMETERS FOR SHEAR FORMULATION (mixed assumed shear)

!  Gauss points position of Assumed Shear Strain Sampling Points
 REAL (kind=8) :: anssp(2,nnas2)
 REAL (kind=8), ALLOCATABLE :: ams(:,:)       ! A-matrix (Lineal)
!  tangent derivatives at each assumed strain point
 REAL (kind=8)  :: sdtan(nnode,nnas2,nface)

!  PA Matrix Product at in-plane Gauss Points (incomplete quadratic approach)
  REAL (kind=8) :: pags(2,nnas2,ngaup)


   ! Derived type for an ELE27 element

   TYPE ele27
     INTEGER (kind=4) :: numel                     ! element label
     INTEGER (kind=4) :: matno                     ! Material number
     INTEGER (kind=4) :: lnods(nnode)              ! Conectivities
     REAL (kind=8) :: angle,                &      ! Euler angle between local t1-t2 and orthotropic
                      dvol(ngaud,ngaup),    &      ! Initial volume
                      sem(3,ngaup,2),       &      ! integrated stresses for membrane geometric stiffness
                      ses(nnas2,2)                 ! integrated stresses for shear geometric stiffness
     REAL (kind=8) :: cqi(3,2,ngaup),           &  ! initial metric tensor at each face at each in-plane gauss point
                      ccasi(2,2,ngaup)             ! Initial shear stain at each face at each in-plane gauss point
     REAL (kind=8), POINTER :: ipcdm(:,:,:,:),  &  ! InPlaneCartesianDerivativeMembrane (nnod1,2,nface,ngaup)
                               cartd(:,:,:,:),  &  ! cartesian ders for transverse shear at gauss points
                               nfdas(:,:,:),    &  ! Nodal Function Derivatives at the Assumed Strain points (Shear)(nnode,nassp,nface)
                               nf3(:,:,:),      &  ! Nodal Function Derivatives at the Gauss points (transverse axialShear)(nnode,ngaup,nface)
                               jacim(:,:,:,:),  &  ! in-plane inverse jacobian for membrane ANS (2,2,nface,ngaud)
                               jacis(:,:,:,:),  &  ! in-plane inverse jacobian for shear ANS (2,2,nface,ngaud)
                               alpha(:),ka(:),  &  ! EAS variable & stiffness at each Gauss point (ngaup)
                               h(:,:),          &  ! h vector at each Gauss point (ndofe,ngaup)
                               stint(:,:,:),    &  !Actual Kirchhoff stresses (6, ngaus, ngaup)
                               c33i(:,:), set(:,:),   &  ! initial strains and integrated stresses for transvers stress (ngaup, 1 or 2
                               gausv(:,:,:)     !Gauss-point internal variables  (7, ngaus, ngaup)
     TYPE (ele27), POINTER :: next              !pointer to next element
   END TYPE ele27

   ! Derived type for a set of ELE27 elements
   TYPE ele27_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for History output
                         narch, &  ! number of output unit
                         ngaus, &  ! number of integration points
                         nassp, &  ! number of Assumed Strain for Shear
                         anssh, &  ! shear option 0: displac.  1: linear 2:imcomplete quadratic
                         easts, &  ! EAS transvers strain option  0: displacemente  1:1 DOF EAS
                         ansmm     ! membrane option
                                   ! 0:standard
                                   ! 1:ANS at corner nodes (not used)
                                   ! 2:ANS at mid-side points of each subtriangle
                                   ! 3:ANS 6 values at sides and 3 at element center (not implemented)
     LOGICAL :: bezier     ! .TRUE. -> use Bezier Functions
     LOGICAL :: gauss      ! .FALSE. -> Initial constants not defined or not updated
     LOGICAL :: small      ! .TRUE. -> use Green strain tensor
                           ! .FALSE. -> Use log strains (Default)
     LOGICAL :: quad       ! output .False. 6 node prism, .TRUE. 15-node prism
     INTEGER :: plstr      ! compute Plastic Strain Flag
                           ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax            ! local x definition option
     REAL (kind=8) :: angdf      ! Default Euler angles between
                                 ! Global X-Y-Z and orthotropic system
     REAL (kind=8) :: fpsg(2,2)   ! factors for output at external faces
     INTEGER(kind=4) :: isg(2,2) ! inner gauss points for output

     TYPE (ele27), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele27_set), POINTER :: next      !pointer to next set
   END TYPE ele27_set
   TYPE (ele27_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE ini_ele27 (head, tail)
     !initialize a list of ELE27 sets

     !Dummy arguments
     TYPE (ele27_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele27

   SUBROUTINE add_ele27 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele27_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele27

   SUBROUTINE srch_ele27 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele27_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele27

   SUBROUTINE del_ele27 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele27_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele27), POINTER :: ea,ep
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
       CALL del_ele27e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele27

   !----------- Element management routines

   SUBROUTINE ini_ele27e (head, tail)
     !initialize a list of ELE27 elements

     ! dummy arguments
     TYPE (ele27), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele27e

   SUBROUTINE new_ele27e(elm)
   !Create a new element of ELE27 sets

     TYPE(ele27),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods = 0d0      !Initializes connectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvol  = 0d0      !     "     elemental volumes
     elm%sem = 0d0        !Initializes equivalent stresses for membrane
     elm%ses = 0d0        !Initializes equivalent stresses for shear
     NULLIFY(elm%cartd,elm%ipcdm,elm%nfdas,elm%jacim,elm%jacis,elm%nf3,elm%stint,elm%gausv) !initializes pointers to allocatable arrays
     NULLIFY(elm%c33i,elm%set) !initializes pointers to allocatable arrays
     NULLIFY( elm%alpha,elm%ka,elm%h)  !nullifies EAS variables
     NULLIFY(elm%next)    !Initializes pointer to next element in the list
   RETURN
   END SUBROUTINE new_ele27e

   SUBROUTINE add_ele27e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele27), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele27e

   SUBROUTINE srch_ele27e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele27), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele27e

   SUBROUTINE del_ele27e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele27), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele27), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%stint)              !deallocate stress array
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele27e

   SUBROUTINE cut_ele27e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele27), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele27e

    INCLUDE 'actu27.fi'
    INCLUDE 'acvd27.fi'
   INCLUDE 'bmat27.fi'
   INCLUDE 'bmai27.fi'
   INCLUDE 'bmma27.fi'
   INCLUDE 'btma27.fi'
   INCLUDE 'bsma27.fi'
   INCLUDE 'comm27.fi'
!!   INCLUDE 'dump27.fi'
   INCLUDE 'elmd27.fi'
!!   INCLUDE 'expo27.fi'
   INCLUDE 'gaus27.fi'
!!   INCLUDE 'impo27.fi'
   INCLUDE 'jacob27.fi'
   INCLUDE 'kgmm27.fi'
   INCLUDE 'kgms27.fi'
   INCLUDE 'kgmt27.fi'
   INCLUDE 'lcsy27.fi'
   INCLUDE 'load27.fi'
   INCLUDE 'mase27.fi'
   INCLUDE 'masm27.fi'
   INCLUDE 'nodx27.fi'
   INCLUDE 'outp27.fi'
!!   INCLUDE 'rest27.fi'
   INCLUDE 'resv27.fi'
   INCLUDE 'stif27.fi'
   INCLUDE 'stra27.fi'

 END MODULE ele27_db
