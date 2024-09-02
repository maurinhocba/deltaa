 SUBROUTINE stiaux_ss_NS(nnode,lnods,nvarl,zigzag,estif,force,lst,ust )

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
 REAL(kind=8), INTENT(IN) :: estif(nvarl,nvarl)  !non-symmetric element matrix
 REAL(kind=8), INTENT(IN OUT) :: force(1),lst(1),ust(1)   !global rhs & lhs vector

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
 i = 0                        !initializes DOF counter and column/file
 DO i = 1,nvarl               !for each node
   gdof = lm(i)               !global DOF
   IF(gdof < -nn) THEN     !if DOF exist
     !               prescribed displacement
     predp = velor(-gdof-nn,nvelr+1)  !prescribed value (velocity)
     IF (predp /= 0d0)THEN               !if non-zero
       IF(.NOT.presd) auxil = 0d0        !initializes if necessary
       auxil = auxil + predp*estif(:,i)  !column i x prescribed value
       presd = .TRUE.                    !assemble to RHS
     END IF
   END IF
 END DO

 IF(presd)CALL ensvec(nvarl,lm(1),auxil(1),force(1))

 CALL ensmatN(nvarl,lm(1),estif(1,1),lst(1),ust(1)) !assemble stiffnes matrix

 RETURN
 END SUBROUTINE stiaux_ss_NS
