 LOGICAL :: ttti,    &! TRUE if Trougth The Thickness Integration is necessary
            pflag     ! TRUE if plastic flow in the step

 LOGICAL :: newmt,   &! TRUE if material constant computation necessary
            found,   &! TRUE if section exists
            natst,   &! TRUE for large strain analysis
            elast,   &! TRUE if material is strictly elastic
            shell,   &! TRUE if bending included
            plast     ! TRUE if plasticity is to be considered

 INTEGER (kind=4) ielem,    & !element number
                  i,j,k,n,l,jn,kn !different indexes

 INTEGER (kind=4) isec,  & !associated material
                  nlayr, & !number of layers
                  mtype, & !associated material type
                  secty, & !section type
                  !numpt, & !number of points in curve
                  osec,  & !associated material of previous element
                  nvar     !number of internal variables per layer

 REAL (kind=8) stran(3),  & !C=U^2  also Log strains
               stres(3),  & !layer stresses or t.t.t integrated stresses (forces)
               stra1(6),  & !first and second fundamental forms
               lambd(7),  & !consistency parameters
               u2(3),     & !C^-1
               r1,r2,     & !eigevenctor components in local system
               lb(3),     & !eigenvalues
               thnew,     & !present thickness
               zk,zk2,    & !distance to mid surface
               aux,       & !auxiliar value
               voli,area1,& !volume and area element
               !! t0,t1,j0,  & !thermical dilatation coeff
               efpst

 REAL (kind=8) thick,     & !thickness (original)
               alpha,     & !thermical dilatation coeff
               propi(13), & !Plastic Material properties
               chi(12),   & !Hill 48 coefficients
               c(4),      & !Elastic constitutive matrix for plane stress
               dm(21),    & !Elastic integrated constitutive matrix
               dmatx(6,6), &
               daux(21),   &
               cm(3,3),    &
               d(3,3),     &
               !deatht,    & !end time for plasticity
               minstr,     & !minimum strain to integrate trougth the thickness
               dummy

 ! Gauss points throught the thickness
 REAL (kind=8) :: thf(7),wei(7)

 !REAL (kind=8), ALLOCATABLE :: cm(:,:), prop(:,:), volfr(:), rr(:,:)
 !REAL (kind=8), POINTER :: val(:,:)
 REAL (kind=8), ALLOCATABLE :: varin(:)                              !internal variables

 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (mater), POINTER :: mat    !pointer to a section data

 INTERFACE
   INCLUDE 'rubberps.h'
 END INTERFACE

