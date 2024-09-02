 SUBROUTINE pcte10 (flag,iload)
 !
 !  Write mesh information for TECPLOT for RIGID elements
 !
 USE data_db
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: flag     !.FALSE. mesh only, .TRUE. results included
 INTEGER(kind=4), INTENT(IN) :: iload

 TYPE( rigid ), POINTER :: e
 INTEGER :: iset,iel,ngv,i,l,nnode,n
 CHARACTER (len=15) ::   etype
 CHARACTER (len=32) :: sname
 LOGICAL :: one                       !flag
 REAL(kind=4), POINTER :: gv(:,:)     !global variables

 INTERFACE
   INCLUDE 'tecgetgv.h'
 END INTERFACE

 CALL tecgetgv(flag,gv,rigid_nodes,rigid_nps,ngv)  !get global variables
 IF( flag) CALL tecgenvarshls(8)   !Generate Variables Share List and Passive Variables list
 e => rigid_head                 !point to first element set

 one = .TRUE.                    !true for first set of beame elemenents
 DO iset=1,rigid_sets            !for each element set
   nnode = rigid_head%nnode
   SELECT CASE(e%ntype)
   CASE (1)
     ! do not print POINTS
     e => e%next
     CYCLE
   CASE (2)
     SELECT CASE (nnode)
     CASE (2)
       etype = 'FELINESEG      '
     CASE (3)
       etype = 'FETRIANGLE     '
     CASE (4)
       etype = 'FEQUADRILATERAL'
     END SELECT
   CASE (3)
     SELECT CASE (nnode)
     CASE (2)
       etype = 'FELINESEG      '
     CASE (4)
       etype = 'FETETRAHEDRON  ' !'FETETRAHEDRAL  '
     CASE (6,8)
       etype = 'FEBRICK        '
     END SELECT
   END SELECT

   l = LEN_TRIM(e%sname)
   sname(1:l+2) = '"'//e%sname(1:l)//'"'     !set name with added apostrophes
   WRITE(11,"('ZONE T=',a,' Nodes=',i6,' Elements=',i6,' ZONETYPE=',a15,    &
              ' DATAPACKING=BLOCK  SOLUTIONTIME=',e12.4 )") &   !,' STRANDID=',i4
              sname(1:l+2),rigid_nps, e%nelem, etype, ttime     !, iload

   IF( flag )THEN !SOLUTION
     IF( t_npv > 0)WRITE(11,"(a)")TRIM(ps_var)  !print list of passive variables
     IF( one )THEN                  !only print for first set
       DO i=1,ngv    !print global variables
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,rigid_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO

     ELSE
       WRITE(11,"(a)")TRIM(sh_var)  !print lis of common variables
     END IF

   ELSE           !GRID
     IF( one )THEN               !only print for first set
       DO i=1,ngv
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,rigid_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
     ELSE
       WRITE(11,*)TRIM(sh_var)
     END IF
     DO iel = 1,e%nelem            !for each element in the set
        WRITE(11,"(8i8)") rigid_nodes(e%lnods(:,iel),1)  !print local numeration
     END DO
   END IF
   one = .FALSE.
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcte10
