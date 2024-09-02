 SUBROUTINE gauss9c( iset )
 !
 ! Process Gauss information and compute nodal values for
 ! 2-D classic beam elements
 !
 USE data_db
 IMPLICIT NONE
 ! dummy arguments
 INTEGER, INTENT(IN OUT) :: iset   !set number

 ! local variables
 REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
 TYPE (shrev), POINTER, SAVE :: eset
 INTEGER :: ielem,g,i,nelem,nnode,nstre,ngaus,n
 REAL (kind=8) :: facta,pos,d(2,2),x0(2,2),disp(2)
 INTEGER :: lnods(4)
 REAL (kind=8), ALLOCATABLE :: elvar(:)
 REAL (kind=8), POINTER :: vargs(:,:)

 INTEGER, POINTER :: nodes(:,:)
 vargs => shrev_vargs         !pointers to use short names

 iset = iset + 1   !update set number

 IF( iset == 1)THEN
   eset => shrev_head  !for first set point the head and reserve memory
 ELSE                  !point to next set
   eset => eset%next
 END IF

 nelem = eset%nelem      !number of elements in the set
 nnode = eset%nnode      !number of nodes per element
 nstre = eset%nstre      !number of variables at each Gauss point
 ngaus = eset%ngaus      !number of integration points

                ! get memory for auxiliar arrays
 ALLOCATE( elvar(nstre))


 DO ielem=1,nelem                    !for each element
   lnods = eset%lnods(:,ielem)       !element connectivities
   x0 = coord(1:2,lnods(1:2))    !original coordinates
   d(:,1) = x0(:,2) - x0(:,1)
   facta = SQRT(DOT_PRODUCT(d(:,1),d(:,1)))
   d(:,1) = d(:,1)/facta
   d(1,2) = -d(2,1)
   d(2,2) =  d(1,1)
   n = lnods(3)
   DO g=1,ngaus
     READ(16,END=100) (elvar(i),i=1,nstre)     !read variables from disk
     !  modify values read considering tolerances
     DO i=1,nstre
       IF (ABS(elvar(i)) > valor1) elvar(i) = SIGN(valor1,elvar(i))
       IF (ABS(elvar(i)) < valor2 .AND. elvar(i) /= 0.0) &
                                     elvar(i) = valor2
     END DO
     eset%elvar(1:nstre,g,ielem) = elvar(1:nstre)  !assign
     disp = MATMUL(d,elvar(4:5))
     dispa(:,n) = disp
     n = n + 1
   END DO

 END DO
 DEALLOCATE (elvar)
 RETURN
 100 fin = .TRUE.      !abnormal end of file detected
 RETURN
 END SUBROUTINE gauss9c
