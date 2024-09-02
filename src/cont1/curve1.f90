SUBROUTINE curve1(x,surf,np)
!
! Compute nodal normals and segment Z-increments of a surface
!
USE cont1_db
IMPLICIT NONE
INTEGER (KIND=4), INTENT (IN) :: np
REAL (KIND=8), INTENT(IN) :: x(:,:)
TYPE (surf1_db) :: surf

!local variables
INTEGER (KIND=4), PARAMETER :: nxn(2) = (/ 2,1 /)
INTEGER (KIND=4) :: ncnod,nsegm,i,j,k,n,nn(2)
INTEGER (KIND=4), ALLOCATABLE :: nr(:)
REAL (KIND=8) :: a(2),t(2),long,tt(2,2)

LOGICAL :: master,slave
REAL (KIND=8), POINTER :: tn(:,:),ta(:,:)


ncnod = surf%ncnod            !number of nodes in the surface
nsegm = surf%nsegm            !number of segments defining surface

slave =  surf%iscod == 1      !If acts as slave
master = surf%imcod == 1      !If acts as master

IF(master) ALLOCATE(surf%cu(nsegm),ta(2,nsegm))
ALLOCATE( tn(2,ncnod),nr(np) )
tn = 0d0                      !average normals initialization
nr = 0                        !inverse relation initialization

DO i=1,ncnod              !generate inverse relation
  ! it stores at the global position the local position in LCNOD
  nr(surf%lcnod(i)) = i
END DO
! Compute average normals
DO n=1,nsegm                  !for each segment
  nn(:) = surf%lcseg(:,n)     !connectivities of the segment
  t(1) = -x(2,nn(2)) + x(2,nn(1))  !normal*l
  t(2) =  x(1,nn(2)) - x(1,nn(1))  !normal*l

  DO i=1,2                    !sums on average normal weighted by the area
    k = nr(nn(i))             !local numeration of the node
    tn(:,k) = tn(:,k) + t     !sums
  END DO
  IF(master)ta(:,n) = t  !store normals Weigthed by the length
END DO
! compute unitary normals
DO n=1,ncnod
  t = tn(:,n)                       !weighted normal at nodes
  long = SQRT(DOT_PRODUCT(t,t))     !length
  tn(:,n) = t/long                  !unit vector
END DO

!If surface acts as MASTER compute normal increments at the side
IF( master)THEN

  DO n=1,nsegm                      !for each segment
    t = ta(:,n)                !restore normal
    long = SQRT(DOT_PRODUCT(t,t))   !compute length (twice the area)
    t = t/long                      !unit vector
    nn(:) = surf%lcseg(:,n)         !connectivities of the segment
    DO i=1,2
      k = nr(nn(i))                 !local numeration of the node
      tt(:,i) = tn(:,k)             !nodal averaged normal
      !COS of angle between nodal average normal and segment normal
      a(i) = DOT_PRODUCT(tt(:,i),t)
    END DO

    !         Compute Z increments  along sides
    surf%cu(n) = DOT_PRODUCT(tt(:,2)-tt(:,1),(/ ta(2,n), -ta(1,n) /) )

    k = 0
    DO i=1,2
      !modify Z-incrementes for plane segments
      IF(ABS(a(i) - 1d0) < 1D-5)THEN      !same normal for node and segment
        j = nxn(i)                        !side nodes
        !if a side element exist
        IF(surf%nhseg(i,n) > 0 .OR. surf%nhseg(j,n) > 0) surf%cu(n) = 0d0
        EXIT
      END IF
    END DO
  END DO
  DEALLOCATE (ta)
END IF

IF( slave )THEN
  surf%tn => tn     !keep normals if the surface acts as SLAVE
ELSE
  DEALLOCATE( tn, surf%lcnod ) !release memory
END IF

DEALLOCATE (nr)     !release temporary memory

RETURN
END SUBROUTINE curve1

