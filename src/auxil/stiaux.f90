 SUBROUTINE stiaux(nnode,lnods,nvarl,ndofe,stiff,force,gstif )

! auxiliar  Routine to assemble stiffness and rhs vector

 USE kinc_db, ONLY : nvelr,velor,nn
 USE npo_db, ONLY : ifpre
 USE ctrl_db, ONLY : neulr
 INTEGER (kind=4), INTENT(IN) :: nnode,        & !number of element nodes
                                 lnods(nnode), & !element connectivities
                                 nvarl,        & !number of element DOFs
                                 ndofe           !number of node DOFs
 REAL(kind=8), INTENT(IN) :: stiff(1)            !nvarl*(nvarl+1)/2 element matrix
 REAL(kind=8), INTENT(IN OUT) :: force(1),     & !global rhs vector
                                 gstif(1)        !global matrix

 LOGICAL :: presd,addof
 INTEGER (kind=4) :: i,j,k,l,n,ij,ii
 INTEGER (kind=4) :: lm(nvarl)
 REAL (kind=8) :: predp,auxil(nvarl)

 INTEGER (kind=4) poesti
 poesti(i,j,n) = (2*n-i)*(i-1)/2+j    ! position i,j in stiff(nxn)


 addof = ndofe == 7 .AND. neulr == 9  !additional DOFs for 3D shells
 !to consider different number of DOFs per node
 presd = .FALSE.              !initializes to no prescribed displacments
 l = 0                        !initializes DOF counter
 DO n = 1,nnode               !for each node
   j = lnods(n)                 !node number
   IF( j > 0)THEN               !if node exist
     DO i = 1,ndofe               !for each possible DOF
       ii = i
       IF( addof .AND. i > 5)ii=i+1  !only for zigzag in 3D problems
       k = ifpre(ii,j)               !DOF number
       l = l+1    !update
       IF(k >= -nn) THEN
         lm(l) = k
       ELSE
         !               prescribed displacement
         lm(l) = 0
         predp = velor(-k-nn,nvelr+1)
         IF (predp /= 0d0)THEN
           IF(.NOT.presd) auxil = 0d0
           DO ij = 1,l
             auxil(ij)= auxil(ij)+stiff(poesti(ij,l,nvarl))*predp
           END DO
           DO ij = l+1,nvarl
             auxil(ij)= auxil(ij)+stiff(poesti(l,ij,nvarl))*predp
           END DO
           presd = .TRUE.
         END IF
       END IF
     END DO
   ELSE
     lm(l+1:l+ndofe) = 0
     l = l+ndofe
   END IF
 END DO
 CALL ensmat(nvarl,lm(1),stiff(1),gstif(1))
 IF(presd)CALL ensvec(nvarl,lm(1),auxil(1),force(1))

 RETURN
 END SUBROUTINE stiaux
