 SUBROUTINE inpda3( nel, sname, etype )
 !
 !  read data for a 2-D SOLID set
 !
 USE data_db
 USE flc_db,ONLY: flc_tp, hpflc, tpflc, new_flc, add_flc, srch_flc
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel,  & !number of elements
                         etype   !element type
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,i,mgaus,mtype,ngaus
 REAL (kind=8) :: angle
 TYPE (sol2d), POINTER :: eset
 TYPE(flc_tp),POINTER:: flc
 LOGICAL :: found

! INTEGER(kind=4), PARAMETER :: kk(2,3) = RESHAPE((/ 3,2, 1,3, 2,1 /), &
!                                                    (/2,3/) )

 INTERFACE
   INCLUDE 'cpdirt.h'
 END INTERFACE

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set

 IF( sol2d_sets == 0 )THEN   !for the first 2-D SOLID set
   sol2d_head => eset
   sol2d_nvarg = 0           !initializes
   IF( sol2d_stres > 1 ) sol2d_nvarg = sol2d_nvarg + 4
   IF( sol2d_logst > 1 ) sol2d_nvarg = sol2d_nvarg + 3
   IF( sol2d_shtst > 1 ) sol2d_nvarg = sol2d_nvarg + 3
   IF( sol2d_thrat > 1 ) sol2d_nvarg = sol2d_nvarg + 1
   IF( sol2d_eqpst > 1 ) sol2d_nvarg = sol2d_nvarg + 1
   IF( sol2d_vmise > 1 ) sol2d_nvarg = sol2d_nvarg + 1
   IF( sol2d_fldma > 1 ) sol2d_nvarg = sol2d_nvarg + 1
   IF( sol2d_wfldFZ .OR. sol2d_wfldSZ ) sol2d_nvarg = sol2d_nvarg + 2
   sol2d_nvarn = 0           !initializes
   IF( MOD(sol2d_stres,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 4
   IF( MOD(sol2d_logst,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 3
   IF( MOD(sol2d_shtst,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 3
   IF( MOD(sol2d_thrat,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 1
   IF( MOD(sol2d_eqpst,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 1
   IF( MOD(sol2d_vmise,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 1
   IF( MOD(sol2d_fldma,2) == 1 ) sol2d_nvarn = sol2d_nvarn + 1
   IF( sol2d_wfldFZ .OR. sol2d_wfldSZ ) sol2d_nvarn = sol2d_nvarn + 2
 ELSE                        !for subsequent sets
   sol2d_tail%next => eset
 END IF
 sol2d_tail => eset          !last set position

 sol2d_sets = sol2d_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set

 READ(17) mtype,        &   !problem type
          mgaus,        &   !number of gauss point in each space direction
          eset%nnode,   &   !number of nodes per element and
          eset%nstre,   &   !number of variables per Gauss point
          eset%ver          !GP position

 IF( ntype == 0 )THEN
   ntype = mtype
 ELSE IF( ntype /= mtype )THEN
   WRITE (*,*) ' mixed problem types for 2-D ????? '
 END IF

 IF( etype == 17 )THEN
   eset%ngaus = mgaus*mgaus
 ELSE
   eset%ngaus = mgaus
 END IF

 ALLOCATE ( eset%lnods(eset%nnode,nel), eset%matno(nel) ) !Connectivities
 IF(sol2d_nvarg > 0) ALLOCATE(eset%elvar(sol2d_nvarg,eset%ngaus,nel))   !Gauss point variables

 ! read connectivities (first element material only)
 DO ielem=1,nel
   READ(17) eset%matno(ielem),(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO

 DO ielem=1,nel
   READ(17) angle
 END DO

 IF(sol2d_thrat > 0 .OR. sol2d_shtst) THEN
   ALLOCATE(eset%dirt(ndime,nel))   !thickness direction
   CALL cpdirt(ndime,eset%nnode,nel,eset%lnods,eset%dirt,coord)
 END IF

 IF( etype == 19 )THEN
   IF( .NOT.ASSOCIATED (sol2d_cpx) )THEN !allocate and initializes
     ALLOCATE(sol2d_cpx(2,npoin))
     sol2d_cpx = 0
   END IF

   DO ielem=1,nel
     DO i=1,3
       n = eset%lnods(i+3,ielem)
       IF(sol2d_cpx(1,n) /= 0)CYCLE
       sol2d_cpx(1,n) = eset%lnods(i,ielem)
       sol2d_cpx(2,n) = eset%lnods(MOD(i,3)+1,ielem)
     END DO
   END DO
 END IF

 RETURN
 END SUBROUTINE inpda3
