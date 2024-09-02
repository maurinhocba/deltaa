 SUBROUTINE gauss1( iset )
 !
 !  read gauss point information
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN OUT) :: iset

 REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
 TYPE (spot), POINTER, SAVE :: eset
 INTEGER :: ielem,g,i,iv,j,n,m,nnode
 INTEGER, ALLOCATABLE :: lnods(:)
 REAL (kind=8), ALLOCATABLE :: elvar(:,:)

 iset = iset + 1   !update set number

 IF( iset == 1)THEN
   eset => spot_head  !for first set, point the head
 ELSE                  !point to next set
   eset => eset%next
 END IF

 nnode = eset%nnode

 ALLOCATE(elvar(eset%nstre,eset%ngaus))   !get memory for auxiliar arrays

 DO ielem=1,eset%nelem               !process all elements
   DO g=1,eset%ngaus                 !for each Gauss point
     READ(16,END=100) (elvar(i,g),i=1,eset%nstre)   !read variables from disk
     !                  check values are in tolerances
     DO i=1,eset%nstre
       IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
       IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) elvar(i,g) = valor2
     END DO

     iv = 0       !initializes pointer to Gauss Values
     !                              Process Axial Forces
     IF( spot_force /= 0 )THEN
       IF( spot_force > 1)THEN     !if Gauss Force desired
         iv = iv+1                  !update pointer to Gaussian values
         eset%elvar(iv,g,ielem) = elvar(1,g)  !assign
       END IF
     END IF
     !                              Process Section stress
     IF( spot_shear /= 0 )THEN
       IF( spot_shear > 1 )THEN    !if Gauss Shear desired
         iv = iv+1                  !update pointer to Gaussian values
         eset%elvar(iv,g,ielem) = elvar(2,g)   !assign
       END IF
     END IF
     !                              Process Section Equivalent plastic strain
     IF( spot_eqpst /= 0 )THEN
       IF( spot_eqpst > 1 )THEN   !if Gauss point EqPlSt desired
         iv = iv+1                 !update pointer to Gaussian values
         eset%elvar(iv,g,ielem) = elvar(3,g)   !assign
       END IF
     END IF
   END DO

 END DO
 RETURN
 100 fin = .TRUE.         !abnormal end of file detected
 RETURN

 END SUBROUTINE gauss1
