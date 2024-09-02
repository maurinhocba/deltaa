 SUBROUTINE addlhs(ien,ndime,nen,force,stiff,ustif,nsymm,lstif,astif,nee)

 !.... add element stiffness into global stiffness
 !  Symetric part      = (Upper + Lower )/2
 !  Anti-Symetric part = (Lower - Upper )/2
 !  Upper part = Symetric part - Antisymmetric part
 !  Lower part = Symetric part + Antisymmetric part

 USE kinc_db, ONLY: nn,nv1,velor
 USE npo_db, ONLY : ifpre
 IMPLICIT NONE
 !       arguments
 INTEGER (kind=4) nen,       &  !number of nodes included
                  ien(nen),  &  !connectivities
                  nee,       &  !number of DOFs (size of local stiffness) = nen*ndime
                  ndime,     &  !number of DOFs per node, Problem dimension
                  nsymm         !=0 symmetric matrix only  =1 unsymmetric matrix availabe
 REAL (kind=8)    force(*),  &  !force vector due prescribed displacements
                  stiff(*),  &  !global stiffness (symmetric part)
                  ustif(*),  &  !global stiffness (asymmetric part)
                  lstif(*),  &  !local stiffness (symmetric part)
                  astif(*)      !local stiffness (asymmetric part)
 !       local variables
 LOGICAL presd !flag
 INTEGER (kind=4) k,n,j,i,ij,iposn,lm(nee),kk
 REAL (kind=8) erest(nee),predp

 presd = .FALSE.  ! initializes prescribed displacement flag
 k = 0            ! initializes DOF value (column in local stiffness)
 DO n=1,nen       ! for each node
   j = ien(n)     ! associated node
   DO i=1,ndime      !for each DOF in the node
     k = k+1         !update DOF value of local stiffness
     iposn = ifpre(i,j)      !associated DOF
     IF(iposn > -nn) THEN    !free or null DOF
       lm(k) = iposn           !Equation number of 0
     ELSE                    !prescribed displacement
       IF(.NOT.presd) erest = 0d0       !initializes array for first time
       predp = velor(-iposn-nn,nv1)     !recover prescribed velocity
       kk = k                           !initializes KK to first column, row K
       ! upper triangle values (transpose row)
       DO ij=1,k-1                      !for each Previous DOF
         !                       Upper = Symmetric - Antisymmetric
         erest(ij) = erest(ij) + (lstif(kk)-astif(kk-ij))*predp
         kk = kk + nee - ij             !update KK to next column, row K
       END DO
       erest(k) = erest(k) + lstif(kk)*predp  !diagonal value
       ! lower triangle values (direct column)
       DO ij=k+1,nee                    !for each remaining DOF
         kk = kk + 1                    !next position in column
         !                Lower = Symmetric + Antisymmetric
         erest(ij) = erest(ij) + (lstif(kk)+astif(kk-k))*predp
       END DO
       lm(k) = 0                        !no associated DOF
       presd = .TRUE.                   !prescribed displacement exists
     END IF
   END DO
 END DO
 !     assembles symmetric component
 CALL ensmat(nee,lm,lstif,stiff)
 !     assembles asymmetric component
 !      kk = nee*(nee-1)/2
 ! IF(nsymm == 1 .AND. ANY(astif(1:kk) /= 0d0))                      &
 !    CALL ensmau(nee,lm,astif,ustif)
 !     assembles prescribed displacement loads
 IF(presd) CALL ensvec(nee,lm,erest,force)
 RETURN

 END SUBROUTINE addlhs
