 SUBROUTINE inpda6( nel, sname )
 !
 !  read data for a 3-D SHELL (SIMO Theory) set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,i,nvarn,nvar,matno
 TYPE (shl3d), POINTER :: eset
 REAL (kind=8) :: angle
 TYPE (udmat), POINTER :: postd

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 IF( shl3d_sets == 0 )THEN   !for the first 3-D SHELL set
   shl3d_head => eset        !initializes head pointer
   shl3d_nvarg = 0           !initializes number of variables
   IF( shl3d_force > 1 ) shl3d_nvarg = shl3d_nvarg + 3
   IF( shl3d_momen > 1 ) shl3d_nvarg = shl3d_nvarg + 3
   IF( shl3d_shear > 1 ) shl3d_nvarg = shl3d_nvarg + 2
   IF( shl3d_logst > 1 ) shl3d_nvarg = shl3d_nvarg + 3
   IF( shl3d_eqpst > 1 ) shl3d_nvarg = shl3d_nvarg + 1
   IF( shl3d_vmise > 1 ) shl3d_nvarg = shl3d_nvarg + 1
   IF( shl3d_thrat > 1 ) shl3d_nvarg = shl3d_nvarg + 1
   shl3d_nvarn = 0           !initializes number of variables
   IF( MOD(shl3d_force,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 3
   IF( MOD(shl3d_momen,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 3
   IF( MOD(shl3d_shear,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 2
   IF( MOD(shl3d_logst,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 3
   IF( MOD(shl3d_eqpst,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 1
   IF( MOD(shl3d_vmise,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 1
   IF( MOD(shl3d_thrat,2) == 1 ) shl3d_nvarn = shl3d_nvarn + 1
 ELSE                        !for subsequent sets
   shl3d_tail%next => eset
 END IF
 shl3d_tail => eset          !last set position

 shl3d_sets = shl3d_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set

 READ(17) eset%nnode,   &   !number of nodes per element
          eset%ngaus,   &   !number of gauss points
          eset%nstre        !number of stress variables to read

 ALLOCATE ( eset%lnods(eset%nnode,nel),eset%matno(nel) )    !Connectivities

  ! read user defined data
  IF( eset%nstre > 11 )THEN
    nvarn = eset%nstre - 11
    READ (17) nvar,matno
    CALL new_post_data(postd,nvar)
    postd%matno = matno
    postd%nvarv = nvarn
    DO i=1,nvar
      READ (17) postd%type(i),postd%dim(i),postd%name(1:nvarn+1,i)
      IF( postd%type(i) == 2)   & !TENSOR
      postd%dim(i) = postd%dim(i) * postd%dim(i)
    END DO
    shl3d_nvarg = shl3d_nvarg + nvarn
  END IF

 IF(shl3d_nvarg > 0) ALLOCATE(eset%elvar(shl3d_nvarg,ABS(eset%ngaus),nel))   !Gauss point variables

 ! read connectivities (first element material only)
 DO ielem=1,nel
   READ(17) eset%matno(ielem),(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO
 DO ielem=1,nel
   READ(17) angle !not used yet
 END DO

 RETURN
 END SUBROUTINE inpda6
