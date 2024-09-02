SUBROUTINE gauss2( iset )
!
!  read gauss point information and computes smoothing arrays
!
USE data_db
IMPLICIT NONE
INTEGER, INTENT(IN OUT) :: iset

REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
TYPE (truss), POINTER, SAVE :: eset
INTEGER :: ielem,g,i,iv,iw,j,n,nnode
INTEGER, ALLOCATABLE :: lnods(:)
REAL (kind=8) :: leng,dist(ndime)
REAL (kind=8), ALLOCATABLE :: elvar(:,:)
INTEGER, POINTER :: nodes(:,:)
REAL (kind=8), POINTER :: vargs(:,:),accpn(:)

iset = iset + 1   !update set number
nodes => truss_nodes    !point to short name arrays
vargs => truss_vargs
accpn => truss_accpn

IF( iset == 1)THEN
  eset => truss_head  !for first set, point the head
  IF( truss_nvarn > 0 )THEN  !initializes smoothing arrays
    accpn = 0.0
    vargs = 0.0
  END IF
ELSE                  !point to next set
  eset => eset%next
END IF

ALLOCATE(elvar(eset%nstre,eset%ngaus))   !get memory for auxiliar arrays

nnode = eset%nnode          !Number of nodes per element
ALLOCATE (lnods(nnode))     !memory for 1-element connectivities

DO ielem=1,eset%nelem               !process all elements
  lnods = eset%lnods(:,ielem)       !element connectivities
  dist = coord(:,lnods(1)) - coord(:,lnods(2))  !element orientation
  leng = SQRT(DOT_PRODUCT(dist,dist))           !element length
  IF( leng == 0d0 )leng = 1d0                   !for dashpots
  leng = 1d0/leng                               !smoothing factor
  DO g=1,eset%ngaus                 !for each Gauss point
    READ(16,END=100) (elvar(i,g),i=1,eset%nstre)   !read variables from disk
    !                  check values are in tolerances
    DO i=1,eset%nstre
      IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
      IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) &
                                    elvar(i,g) = valor2
    END DO

    iv = 0       !initializes pointer to Gauss Values
    iw = 0       !initializes pointer to Nodal Values
    !                              Process Section Forces
    IF( truss_force /= 0 )THEN
      IF( truss_force > 1)THEN     !if Gauss point Force desired
        iv = iv+1                  !update pointer to Gaussian values
        eset%elvar(iv,g,ielem) = elvar(1,g)  !assign
      END IF
      IF( MOD(truss_force,2) == 1)THEN    !if nodal section force desired
        iw = iw+1                         !update pointer to vargs
        DO j=1,nnode                      !for each node
          n = nodes(lnods(j),1)             !node order in vargs
          vargs(iw,n) = vargs(iw,n) + elvar(1,g)*leng  !sum on smoothing array
          IF(iw == 1) accpn(n) = accpn(n) + leng         !sum on nodal weigths
        END DO
      END IF
    END IF
    !                              Process Section stress
    IF( truss_stres /= 0 )THEN
      IF( truss_stres > 1 )THEN    !if Gauss point stress desired
        iv = iv+1                  !update pointer to Gaussian values
        eset%elvar(iv,g,ielem) = elvar(2,g)   !assign
      END IF
      IF( MOD(truss_stres,2) == 1)THEN    !if nodal section stress desired
        iw = iw+1                         !update pointer to vargs
        DO j=1,nnode                      !for each node
          n = nodes(lnods(j),1)             !node order in vargs
          vargs(iw,n) = vargs(iw,n) + elvar(2,g)*leng !sum on smoothing array
          IF(iw == 1) accpn(n) = accpn(n) + leng         !sum on nodal weigths
        END DO
      END IF
    END IF
    !                              Process Section Equivalent plastic strain
    IF( truss_eqpst /= 0 )THEN
      IF( truss_eqpst > 1 )THEN   !if Gauss point EqPlSt desired
        iv = iv+1                 !update pointer to Gaussian values
        eset%elvar(iv,g,ielem) = elvar(3,g)   !assign
      END IF
      IF( MOD(truss_eqpst,2) == 1)THEN    !if nodal Equiv. Pl Str. desired
        iw = iw+1                         !update pointer to vargs
        DO j=1,nnode                      !for each node
          n = nodes(lnods(j),1)             !node order in vargs
          vargs(iw,n) = vargs(iw,n) + elvar(3,g)*leng !sum on smoothing array
          IF(iw == 1) accpn(n) = accpn(n) + leng         !sum on nodal weigths
        END DO
      END IF
    END IF
  END DO

END DO
RETURN
100 fin = .TRUE.         !abnormal end of file detected
RETURN

END SUBROUTINE gauss2
