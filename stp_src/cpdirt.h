 SUBROUTINE  cpdirt(ndime,nnode,nelem,lnods,dirt,coord)

 ! computes thicknes direction for each element

 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ndime, & !=2 problem dimension
                                 nnode, & !3 or 4 number of nodes per element
                                 nelem, & !number of elements in the mesh
                                 lnods(:,:)     !connectivities
 REAL (kind=8), INTENT(IN) :: coord(:,:)        !mesh coordinates
 REAL (kind=8), INTENT(OUT) :: dirt(:,:)        !thickness direction
 END SUBROUTINE cpdirt
