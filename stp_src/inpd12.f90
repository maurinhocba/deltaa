SUBROUTINE inpd12( nel, sname, etype )
!
!  read data for a 3-D SHELL (BST) set
!
USE data_db
USE flc_db,ONLY: nflc
IMPLICIT NONE

  !Dummy variables
  INTEGER(kind=4),INTENT(IN)::  nel,  & !number of elements
                                etype   !element type
  CHARACTER(len=30),INTENT(IN):: sname
  !Local variables
  INTEGER(kind=4):: ielem, olbflc,nvar,nvarn,i,matno
  REAL(kind=8):: angle
  LOGICAL:: found
  TYPE(bst),POINTER:: eset
  TYPE (udmat), POINTER :: postd

  ALLOCATE (eset)               !get memory for this set
  NULLIFY (eset%next)           !nullify pointer to next set

  IF( bst_sets == 0 )THEN   !for the first 3-D SHELL set
    bst_head => eset        !initializes head pointer
    bst_nvarg = 0           !initializes number of variables at Gauss Points
    IF( bst_force > 1 ) bst_nvarg = bst_nvarg + 3
    IF( bst_momen > 1 ) bst_nvarg = bst_nvarg + 3
    IF( bst_shear > 1 ) bst_nvarg = bst_nvarg + 2
    IF( bst_logst > 1 ) bst_nvarg = bst_nvarg + 3
    IF( bst_curva > 1 ) bst_nvarg = bst_nvarg + 3
    IF( bst_thrat > 1 ) bst_nvarg = bst_nvarg + 1
    IF( bst_eqpst > 1 ) bst_nvarg = bst_nvarg + 2
    IF( bst_vmise > 1 ) bst_nvarg = bst_nvarg + 2
    IF(nflc == 0)THEN  !check if FLC is available
      bst_fldma = 0
      bst_wfldFZ = .FALSE.
      bst_wfldSZ = .FALSE.
    END IF
    IF( bst_fldma > 1 ) bst_nvarg = bst_nvarg + 1
    IF( bst_wfldFZ .OR. bst_wfldSZ ) bst_nvarg = bst_nvarg + 2
    bst_nvarn = 0           !initializes number of variables at Nodal Points
    IF( MOD(bst_force,2) == 1 ) bst_nvarn = bst_nvarn + 3
    IF( MOD(bst_momen,2) == 1 ) bst_nvarn = bst_nvarn + 3
    IF( MOD(bst_shear,2) == 1 ) bst_nvarn = bst_nvarn + 2
    IF( MOD(bst_logst,2) == 1 ) bst_nvarn = bst_nvarn + 3
    IF( MOD(bst_curva,2) == 1 ) bst_nvarn = bst_nvarn + 3
    IF( MOD(bst_thrat,2) == 1 ) bst_nvarn = bst_nvarn + 1
    IF( MOD(bst_eqpst,2) == 1 ) bst_nvarn = bst_nvarn + 2
    IF( MOD(bst_vmise,2) == 1 ) bst_nvarn = bst_nvarn + 2
    IF( MOD(bst_fldma,2) == 1 ) bst_nvarn = bst_nvarn + 1
    IF( bst_wfldFZ .OR. bst_wfldSZ ) bst_nvarn = bst_nvarn + 2
  ELSE                        !for subsequent sets
    bst_tail%next => eset
  END IF
  bst_tail => eset          !last set position

  bst_sets = bst_sets + 1 !increase number of sets
  eset%set = nsets        !set position (possibly unnecessary)

  eset%sname = sname            !set name
  eset%nelem = nel              !number of elements in the set

  READ(17) eset%nnode,   &   !number of nodes per element  = 3 or 4
           eset%ngaus,   &   !number of gauss points       = 1
           eset%nstre,   &   !number of stress variables to read = 13
           eset%locax        !local system definition

  ! read user defined data
  IF( eset%nstre > 13 )THEN
    nvarn = eset%nstre - 13
    READ (17) nvar,matno
    CALL new_post_data(postd,nvar)
    postd%matno = matno
    postd%nvarv = nvarn
    DO i=1,nvar
      READ (17) postd%type(i),postd%dim(i),postd%name(1:nvarn+1,i)
      IF( postd%type(i) == 2)   & !TENSOR
      postd%dim(i) = postd%dim(i) * (postd%dim(i)+1) /2
    END DO
    bst_nvarg = bst_nvarg + nvarn
  END IF

  IF(bst_nvarg > 0 )ALLOCATE(eset%elvar(bst_nvarg,eset%ngaus,nel))   !Gauss point variables

  ! read connectivities
  IF( eset%nnode == 3 )THEN
    ALLOCATE ( eset%lnods(9,nel) ) !patch Connectivities + elements
  ELSE !IF ( eset%nnode = 4 )THEN
    ALLOCATE ( eset%lnods(4,nel) ) !element Connectivities
  END IF
  ALLOCATE(eset%matno(nel))  !associated material
  DO ielem=1,nel
    IF( eset%nnode == 3 )THEN
      READ(17) eset%matno(ielem),eset%lnods(1:6,ielem)
    ELSE
      READ(17) eset%matno(ielem),eset%lnods(1:4,ielem)
    END IF
  END DO
  IF( eset%nstre > 13 )postd%matno = eset%matno(1)

  ALLOCATE( eset%angle(nel) )
  DO ielem=1,nel
    READ(17) angle
    eset%angle(ielem) = angle
  END DO

RETURN
END SUBROUTINE inpd12
