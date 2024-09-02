 SUBROUTINE inpda9( nel, sname )
 !
 !  read data for a 2-D SHELL/BEAM set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,mtype
 TYPE (shrev), POINTER :: eset

 ALLOCATE (eset)             !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 READ(17) eset%nnode,   &   !number of nodes per element
          eset%nstre,   &   !number of Gauss variables to read
          eset%ngaus,   &   !number of gauss points
          mtype             !problem type

 IF( shrev_sets == 0 )THEN   !for the first 2-D SHELL/BEAM set
   shrev_head => eset        !initializes head pointer
   shrev_nvarg = 0           !initializes number of variables
   IF( shrev_force > 1 ) THEN
     IF( mtype == 1 .OR. mtype == 5)THEN
       shrev_nvarg = shrev_nvarg + 1
     ELSE
       shrev_nvarg = shrev_nvarg + 2
     END IF
   END IF
   IF( shrev_momen > 1 ) THEN
     IF( mtype == 1 .OR. mtype ==5 )THEN
       shrev_nvarg = shrev_nvarg + 1
     ELSE
       shrev_nvarg = shrev_nvarg + 2
     END IF
     IF( (mtype == 1 .AND. eset%nstre == 5 ) .OR. eset%nstre == 7) shrev_nvarg = shrev_nvarg + 1
   END IF
   IF( shrev_shear > 1 ) THEN
     shrev_nvarg = shrev_nvarg + 1
     IF( (mtype == 1 .AND. eset%nstre == 5 ) .OR. eset%nstre == 7) shrev_nvarg = shrev_nvarg + 1
   END IF
   IF( shrev_eqpst > 1 ) shrev_nvarg = shrev_nvarg + 1
   IF( shrev_vmise > 1 ) shrev_nvarg = shrev_nvarg + 1
   IF( shrev_thrat > 1 ) shrev_nvarg = shrev_nvarg + 1
   shrev_nvarn = 0           !initializes number of variables
   IF( MOD(shrev_force,2) == 1 ) THEN
     IF( mtype == 1 .OR. mtype == 5 )THEN
       shrev_nvarn = shrev_nvarn + 1
     ELSE
       shrev_nvarn = shrev_nvarn + 2
     END IF
   END IF
   IF( MOD(shrev_shear,2) == 1 ) shrev_nvarn = shrev_nvarn + 1
   IF( MOD(shrev_momen,2) == 1 ) THEN
     IF( mtype == 1 )THEN
       shrev_nvarn = shrev_nvarn + 1
     ELSE
       shrev_nvarn = shrev_nvarn + 2
     END IF
   END IF
   IF( MOD(shrev_eqpst,2) == 1 ) shrev_nvarn = shrev_nvarn + 1
   IF( MOD(shrev_vmise,2) == 1 ) shrev_nvarn = shrev_nvarn + 1
   IF( MOD(shrev_thrat,2) == 1 ) shrev_nvarn = shrev_nvarn + 1
 ELSE                        !for subsequent sets
   shrev_tail%next => eset
 END IF
 shrev_tail => eset          !last set position

 shrev_sets = shrev_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set

 IF( mtype == 4 )ngp = ngp+eset%nelem*eset%ngaus

 IF( ntype == 0 )THEN
   ntype = mtype
   IF( mtype == 4 )THEN  ! modify flags
     shrev_force = 3  !only in Gauss points
     shrev_shear = 3
     shrev_momen = 3
     shrev_nvarg = 5
     shrev_nvarn = 5
     shrev_eqpst = 0
     shrev_vmise = 0
     shrev_thrat = 0
   END IF
 ELSE IF( ntype /= mtype )THEN
   WRITE (*,*) ' mixed problem types for 2-D ????? '
 END IF

 IF( mtype == 4 )THEN  ! modify flags
   ALLOCATE ( eset%lnods(4,nel) )   !Connectivities
 ELSE
   ALLOCATE ( eset%lnods(eset%nnode,nel) )   !Connectivities
 END IF
 IF(shrev_nvarg > 0)ALLOCATE(eset%elvar(shrev_nvarg,eset%ngaus,nel)) !Gauss point variables

 ! read connectivities (first element material only)
 DO ielem=1,nel
   READ(17) eset%matno,(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO
 RETURN
 END SUBROUTINE inpda9
