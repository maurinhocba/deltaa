 SUBROUTINE inpda8( nel, sname )
 !
 !  read data for a 3-D BEAME set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n
 TYPE (beame), POINTER :: eset

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 IF( beame_sets == 0 )THEN   !for the first 3-D BEAME set
   beame_head => eset        !initializes head pointer
   beame_nvarg = 0           !initializes number of variables
   IF( beame_force > 1 ) beame_nvarg = beame_nvarg + 3
   IF( beame_momen > 1 ) beame_nvarg = beame_nvarg + 3
   IF( beame_eqpst > 1 ) beame_nvarg = beame_nvarg + 1
   IF( beame_vmise > 1 ) beame_nvarg = beame_nvarg + 1
   IF( beame_arrat > 1 ) beame_nvarg = beame_nvarg + 1
   beame_nvarn = 0           !initializes number of variables
   IF( MOD(beame_force,2) == 1 ) beame_nvarn = beame_nvarn + 3
   IF( MOD(beame_momen,2) == 1 ) beame_nvarn = beame_nvarn + 3
   IF( MOD(beame_eqpst,2) == 1 ) beame_nvarn = beame_nvarn + 1
   IF( MOD(beame_vmise,2) == 1 ) beame_nvarn = beame_nvarn + 1
   IF( MOD(beame_arrat,2) == 1 ) beame_nvarn = beame_nvarn + 1
 ELSE                        !for subsequent sets
   beame_tail%next => eset
 END IF
 beame_tail => eset          !last set position

 beame_sets = beame_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set

 READ(17) eset%nnode,   &   !number of nodes per element
          eset%nstre,   &   !number of stress variables to read
          eset%ngaus        !number of gauss points

 ALLOCATE ( eset%lnods(eset%nnode,nel) ) !Connectivities
 IF(beame_nvarg > 0) ALLOCATE(eset%elvar(beame_nvarg,eset%ngaus,nel))   !Gauss point variables

 ! read connectivities (first element material only)
 DO ielem=1,nel
   READ(17) eset%matno,(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO
 RETURN
 END SUBROUTINE inpda8
