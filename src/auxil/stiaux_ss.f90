 SUBROUTINE stiaux_ss(nnode,lnods,nvarl,zigzag,stiff,force,gstif )

! auxiliar  Routine to assemble stiffness and rhs vector for SOLSH + ZIGZAG
 !DOFs are ordered as follows
 ! node 1-8 1-3             1-24
 ! node 1-4 7-8  if zigzag 25-32

 USE kinc_db, ONLY : nvelr,velor,nn
 USE npo_db, ONLY : ifpre
 INTEGER (kind=4), INTENT(IN) :: nnode,        & !number of element nodes
                                 lnods(nnode), & !element connectivities
                                 nvarl           !number of element DOFs
 LOGICAL, INTENT(IN) ::          zigzag          !if additional DOFs
 REAL(kind=8), INTENT(IN) :: stiff(1)            !nvarl*(nvarl+1)/2 element matrix
 REAL(kind=8), INTENT(IN OUT) :: force(1),     & !global rhs vector
                                 gstif(1)        !global matrix

 INTEGER (kind=4) :: i,j,k,l,n,ij,ia,le
 INTEGER (kind=4) :: lm(nvarl)
 REAL (kind=8) :: auxil(nvarl),predp
 LOGICAL :: presd
 INTEGER (kind=4) poesti
 poesti(i,j,n) = (2*n-i)*(i-1)/2+j    ! position i,j in stiff(nxn)

 ! first Gather DOFs

 l = 1                   !initializes position in array LM
 DO n=1,nnode            !for each basic node
   j = lnods(n)                !node number
   lm(l:l+2) = ifpre(1:3,j)    !pass DOFs
   l = l+3                     !update counter
 END DO
 IF( zigzag ) THEN
   ia = 25                 !first position of additional DOFs
   DO n=1,nnode/2
     j = lnods(n)                !node number
     lm(ia:ia+1) = ifpre(7:8,j)  !additiional DOFs if exist
     ia = ia+2
   END DO
 END IF
 ! check prescribed displacements and generate AUXIL
 presd = .FALSE.              !initializes to no prescribed displacments
 DO l=1,nvarl
   IF( lm(l) >= -nn )CYCLE
   predp = velor(-lm(l)-nn,nvelr+1)
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
 END DO

 CALL ensmat(nvarl,lm(1),stiff(1),gstif(1))           !assemble stiffnes matrix
 IF(presd)CALL ensvec(nvarl,lm(1),auxil(1),force(1))  !assemble rhs if exists

 RETURN
 END SUBROUTINE stiaux_ss
