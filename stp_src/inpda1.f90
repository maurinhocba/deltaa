 SUBROUTINE inpda1( nel, sname )
 !
 !  read data for a SPOT set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,m,i,j,k,ngaus
 TYPE (spot), POINTER :: eset

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 IF( spot_sets == 0 )THEN   !for the first SPOT set
   spot_head => eset
   spot_nvarg = 0           !initializes number of variables at Gauss points
   IF( spot_force > 1 ) spot_nvarg = spot_nvarg + 1
   IF( spot_shear > 1 ) spot_nvarg = spot_nvarg + 1
   IF( spot_eqpst > 1 ) spot_nvarg = spot_nvarg + 1
 ELSE                        !for subsequent sets
   spot_tail%next => eset
 END IF
 eset%first_l = last_label         !auxiliar labels

 spot_tail => eset           !last set position

 spot_sets = spot_sets + 1   !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set
 ngaus = 1                     !number of gauss points per element

 READ(17) eset%nnode,eset%nstre !number of nodes per element and
                                !number of variables to read per Gauss point (3)
 eset%ngaus = 1                 !initializes

 ALLOCATE ( eset%lnods(eset%nnode,nel) )   !connectivities
 IF(spot_nvarg > 0 ) ALLOCATE (eset%elvar(spot_nvarg,ngaus,nel)) !Gauss point variables

 ! read connectivities
 IF( eset%nnode == 6 ) THEN
   ALLOCATE (eset%a_x(ndime,2*eset%nelem,2), eset%lc(3,2,eset%nelem))
   eset%a_x = 0d0
   i = 0
   DO ielem=1,nel
     READ(17) eset%matno,(eset%lnods(n,ielem),n=1,eset%nnode)
     READ(17) eset%lc(1:2,:,ielem)
     j = 0
     DO m=1,2
       eset%lc(3,m,ielem) = 1d0 - eset%lc(1,m,ielem) - eset%lc(2,m,ielem)
       i = i+1
       DO n=1,3
         j = j + 1
         eset%a_x(:,i,1) = eset%a_x(:,i,1) + coord(:,eset%lnods(j,ielem))*eset%lc(n,m,ielem) !original coordinates
         eset%a_x(:,i,2) = eset%a_x(:,i,2) + coors(:,eset%lnods(j,ielem))*eset%lc(n,m,ielem) !stage coordinates
       END DO
     END DO
   END DO
   last_label = eset%first_l + 2*eset%nelem
 ELSE
   m = 0                              !counter of auxiliar nodes
   DO ielem=1,nel                     !loop over elements in the set
     READ(17) eset%matno,(eset%lnods(n,ielem),n=1,eset%nnode)   !read connectivities
     IF( eset%lnods(2,ielem) == 0 )THEN   !one node spot
       m = m+1                              !increase number of auxiliar nodes
       eset%lnods(2,ielem) = -m             !store as a flag
     END IF
   END DO
   IF( m > 0 )THEN                        !if auxiliar nodes defined
     ALLOCATE (eset%a_x(ndime,m,1))       !reserve memory for nodes
     last_label = last_label + m                        !increase first label for next set
     m = 0                                  !initializes counter
     DO ielem=1,nel                       !loop again
       IF( eset%lnods(2,ielem) < 0 )THEN    !if an auxiliar node
         m = m+1                            !correct counter
         READ(17)(eset%a_x(i,m,1),i=1,ndime)   !read original coordinates
         !eset%a_x(:,m,1) = coord(:,eset%lnods(1,ielem))   !keep original coordinates
       END IF
     END DO
   END IF
 END IF
 spot_auxiliar_nodes = last_label - label(npoin)

 RETURN
 END SUBROUTINE inpda1
