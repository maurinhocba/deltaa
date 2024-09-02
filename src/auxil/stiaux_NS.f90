 SUBROUTINE stiaux_NS(nnode,lnods,nvarl,ndofe,estif,force,lst,ust )

! auxiliar  Routine to assemble stiffness and rhs vector

 USE kinc_db, ONLY : nvelr,velor,nn
 USE npo_db, ONLY : ifpre
 INTEGER (kind=4), INTENT(IN) :: nnode,        & !number of element nodes
                                 lnods(nnode), & !element connectivities
                                 nvarl,        & !number of element DOFs
                                 ndofe           !number of node DOFs
 REAL(kind=8), INTENT(IN) :: estif(nvarl,nvarl)  !element matrix
 REAL(kind=8), INTENT(IN OUT) :: force(1),lst(1),ust(1)   !global rhs & lhs vector

 LOGICAL :: presd
 INTEGER (kind=4) :: i,j,k,nnu,n,idof,gdof
 INTEGER (kind=4) :: lm(nvarl)
 REAL (kind=8) :: predp,auxil(nvarl)

 ! rewrite as two arrays
 presd = .FALSE.              !initializes to no prescribed displacments
 i = 0                        !initializes DOF counter and column/file
 DO n = 1,nnode               !for each node
   nnu = lnods(n)                 !node number
   DO idof = 1,ndofe               !for each possible DOF
     gdof = ifpre(idof,nnu)              !DOF number
     i = i+1    !update DOF counter
     IF(gdof >= -nn) THEN     !if DOF exist
       lm(i) = gdof             !pass equation number to array
     ELSE
       !               prescribed displacement
       lm(i) = 0                !no associated equation
       predp = velor(-gdof-nn,nvelr+1)  !prescribed value (velocity)
       IF (predp /= 0d0)THEN               !if non-zero
         IF(.NOT.presd) auxil = 0d0        !initializes if necessary
         auxil = auxil + predp*estif(:,i)  !column i x prescribed value
         presd = .TRUE.                    !assemble to RHS
       END IF
     END IF
   END DO
 END DO

 IF(presd)CALL ensvec(nvarl,lm(1),auxil(1),force(1))

 CALL ensmatN(nvarl,lm(1),estif(1,1),lst(1),ust(1))
 RETURN
 END SUBROUTINE stiaux_NS
