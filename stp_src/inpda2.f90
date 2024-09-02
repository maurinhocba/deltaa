 SUBROUTINE inpda2( nel, sname )
 !
 !  read data for a TRUSS set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,m,i,j,k,ngaus
 TYPE (truss), POINTER :: eset

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 IF( truss_sets == 0 )THEN   !for the first TRUSS set
   truss_head => eset
   truss_nvarg = 0           !initializes number of variables at Gauss points
   IF( truss_force > 1 ) truss_nvarg = truss_nvarg + 1
   IF( truss_stres > 1 ) truss_nvarg = truss_nvarg + 1
   IF( truss_eqpst > 1 ) truss_nvarg = truss_nvarg + 1
   truss_nvarn = 0           !initializes number of variables at Nodal points
   IF( MOD(truss_force,2) == 1 ) truss_nvarn = truss_nvarn + 1
   IF( MOD(truss_stres,2) == 1 ) truss_nvarn = truss_nvarn + 1
   IF( MOD(truss_eqpst,2) == 1 ) truss_nvarn = truss_nvarn + 1
 ELSE                        !for subsequent sets
   truss_tail%next => eset
 END IF
 truss_tail => eset          !last set position

 truss_sets = truss_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set
 ngaus = 1                     !number of gauss points per element

 READ(17) eset%nnode,eset%nstre !number of nodes per element and
                                !number of variables to read per Gauss point (3)
 eset%ngaus = 1                 !initializes

 ALLOCATE ( eset%lnods(eset%nnode,nel) )   !connectivities
 IF(truss_nvarg > 0 ) ALLOCATE (eset%elvar(truss_nvarg,ngaus,nel)) !Gauss point variables

 ! read connectivities
 DO ielem=1,nel
   READ(17) eset%matno,(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO

 RETURN
 END SUBROUTINE inpda2
