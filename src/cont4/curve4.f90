SUBROUTINE curve4(x,surf,np)
!
! Compute nodal normals and segment Z-increments of a surface
!
USE cont4_db
IMPLICIT NONE
INTEGER (KIND=4), INTENT (IN) :: np
REAL (KIND=8), INTENT(IN) :: x(:,:)
TYPE (surf4_db) :: surf

!local variables
INTEGER (KIND=4), PARAMETER :: nxn(3) = (/ 2,3,1 /)
INTEGER (KIND=4) :: ncnod,nsegm,i,j,jj,k,n,nn(3)
INTEGER (KIND=4), ALLOCATABLE :: nr(:)
REAL (KIND=8) :: ls(3,3),a(3),t(3),long,tt(3,3)

LOGICAL :: master,slave
REAL (KIND=8), POINTER :: tn(:,:)


ncnod = surf%ncnod            !number of nodes in the surface
nsegm = surf%nsegm            !number of segments defining surface

slave =  surf%iscod == 1      !If acts as slave
master = surf%imcod == 1      !If acts as master

IF(master) ALLOCATE(surf%cu(3,nsegm))
ALLOCATE(tn(3,ncnod),nr(np))
tn = 0d0                      !average normals initialization
nr = 0                        !inverse relation initialization

DO i=1,ncnod              !generate inverse relation
  ! it stores at the global position the local position in LCNOD
  nr(surf%lcnod(i)) = i
END DO
! Compute average normals
DO n=1,nsegm                  !for each segment
  nn(:) = surf%lcseg(:,n)     !connectivities of the segment
  DO i=1,2                    !first two side elements
    k = nxn(i)                              !next node
    ls(:,i) = x(:,nn(nxn(k))) - x(:,nn(k))  !opposite side
  END DO
  t(1) = ls(2,1)*ls(3,2) - ls(2,2)*ls(3,1)   !normal*a2
  t(2) = ls(3,1)*ls(1,2) - ls(3,2)*ls(1,1)
  t(3) = ls(1,1)*ls(2,2) - ls(1,2)*ls(2,1)
  DO i=1,3                    !sums on average normal weighted by the area
    k = nr(nn(i))             !local numeration of the node
    tn(:,k) = tn(:,k) + t     !sums
  END DO
  IF(master)surf%cu(:,n) = t  !store normals (temporary)
END DO
! compute unitary normals
DO n=1,ncnod
  t = tn(:,n)                       !weighted normal
  long = SQRT(DOT_PRODUCT(t,t))     !length
  tn(:,n) = t/long                  !unit vector
END DO

!If surface acts as MASTER compute normal increments at the side
IF( master)THEN

  DO n=1,nsegm                      !for each segment
    t = surf%cu(:,n)                !restore normal
    long = SQRT(DOT_PRODUCT(t,t))   !compute length
    t = t/long                      !unit vector
    nn(:) = surf%lcseg(:,n)         !connectivities of the segment
    DO i=1,3
      k = nr(nn(i))                 !local numeration of the node
      tt(:,i) = tn(:,k)             !nodal average normal
      !COS of angle between nodal average normal and segment normal
      a(i) = DOT_PRODUCT(tt(:,i),t)
      j = nxn(i)                    !next node
      ls(:,i) = x(:,nn(nxn(j))) - x(:,nn(j))  !opposite side
    END DO

    !         Compute Z increments  along sides
    DO i=1,3
      k = nxn(i)                    !next node
      surf%cu(i,n) = DOT_PRODUCT(tt(:,nxn(k))-tt(:,k),ls(:,i))
    END DO

    k = 0
    DO i=1,3
      !modify Z-incrementes for plane segments
      IF(ABS(a(i) - 1d0) < 1D-8)THEN      !same normal for node and segment
        j = nxn(i)                        !side nodes
        jj = nxn(j)
        !if a side element exist
        IF(surf%nhseg(j,n) /= 0 .OR. surf%nhseg(jj,n) /=0) surf%cu(i,n) = 0d0
        k = k+1
      END IF
    END DO
    IF( k > 1) surf%cu(1:3,n) = 0d0
    !write(55,"(3e15.3)")surf%cu(:,n)
  END DO

END IF

IF( slave )THEN
  surf%tn => tn     !keep normals if the surface acts as SLAVE
ELSE
  DEALLOCATE( tn, surf%lcnod ) !release memory
END IF

DEALLOCATE (nr)     !release temporary memory

RETURN
END SUBROUTINE curve4
