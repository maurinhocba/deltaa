 SUBROUTINE pcte12 (flag,iload)
 !
 !  Write mesh information for TECPLOT for BST elements
 !
 USE data_db
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: flag  !.FALSE. mesh only, .TRUE. results included
 INTEGER(kind=4), INTENT(IN) :: iload

 TYPE( bst ), POINTER :: e
 INTEGER :: ngv,i,n,iset,l,iel,nlgv
 CHARACTER (len=15) ::   etype
 CHARACTER (len=32) :: sname

 LOGICAL :: one                       !flag
 REAL(kind=4), POINTER :: gv(:,:)     !global variables
 REAL(kind=4), ALLOCATABLE :: lgv(:)  !local gauss variables

 INTERFACE
   INCLUDE 'tecgetgv.h'
 END INTERFACE

 ! tasks to do
 ! 1-PRINT mesh (flag=.FALSE.)
 !   a-print nodal coordinates (coord or coors if stage displements) of Truss sets only
 !   b-print connectivities of each truss set  using VARSHARELIST=([1,2,3]) (ndime)
 ! 2-PRINT results (flag=.TRUE.)
 !   variables are of three 'Groups'
 !         I-Global nodal variables: displacements, velocities, accelerations, pressure and temperatures (1:ngvar)
 !        II-Nodal Smoothed variables (truss_nvarn) for this type of elements (n1:n1-1+truss_nvarn), n1 must be determined
 !       III-Gaussian variables (truss_nvarg) for this type of elements (n2:n2-1+truss_nvarg), n2 must be determined
 !   a-print nodal displacements, velocities, etc. and smoothed variables at the first zone
 !   b-print Gaussian variables at each zone using VARLOCATION VARLOCATION=([?,?]=CELLCENTERED)
 !   c-it must be determined what to include  VARSHARELIST    (type I and II)
 !   d-it must be determined what to include  PASSIVEVARLIST  (thouse no included in type I, II and III)

 CALL tecgetgv(flag,gv,bst_nodes,bst_nps,ngv)  !get global variables
 IF( flag) CALL tecgenvarshls(9)   !Generate Variables Share List and Passive Variables list
 nlgv = t_ev(2,9)                  !number of Local Gauss Variables

 e => bst_head                 !point to first element set
 one = .TRUE.                    !true for first set of bst elemenents
 DO iset=1,bst_sets            !for each element set
   IF( e%nnode == 4)THEN
     etype = 'FEQUADRILATERAL'
   ELSE                          !
     etype = 'FETRIANGLE     '
   END IF

   l = LEN_TRIM(e%sname)
   sname(1:l+2) = '"'//e%sname(1:l)//'"'     !set name with added apostrophes
   WRITE(11,"('ZONE T=',a,' Nodes=',i6,' Elements=',i6,' ZONETYPE=',a15,    &
              ' DATAPACKING=BLOCK  SOLUTIONTIME=',e12.4  )") &   !,' STRANDID=',i4
               sname(1:l+2),  bst_nps, e%nelem, etype, ttime     !, iload

   IF( flag )THEN !SOLUTION
     IF( t_npv > 0)WRITE(11,"(a)")TRIM(ps_var)  !print list of passive variables
     IF( nlgv > 0) WRITE(11,"(a)")TRIM(lc_var)  !print variables location
     IF( one )THEN                  !only print for first set
       DO i=1,ngv    !print global variables
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,bst_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DO i=1,bst_nvarn !print smoothed variables
         WRITE(11,"(10e15.6)")(bst_vargs(i,n),n=1,bst_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO

     ELSE
       WRITE(11,"(a)")TRIM(sh_var)  !print lis of common variables
     END IF
     ! print Gaussian Variables
     IF( nlgv > 0 )THEN
       ALLOCATE( lgv(e%nelem) )   !get memory for average element values
       DO i=1,nlgv                ! for each variabls
         DO iel=1,e%nelem         ! compute cell-center values
           lgv(iel) = e%elvar(i,1,iel)
         END DO
         WRITE(11,"(10e15.6)") (lgv(iel),iel=1,e%nelem)  !print gauss values
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DEALLOCATE( lgv )
     END IF

   ELSE           !GRID
     IF( one )THEN               !only print for first set
       DO i=1,ngv
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,bst_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
     ELSE
       WRITE(11,*)TRIM(sh_var)
     END IF
     DO iel = 1,e%nelem            !for each element in the set
       WRITE(11,"(4i8)") bst_nodes(e%lnods(1:e%nnode,iel),1)  !print local numeration
     END DO
   END IF
   one = .FALSE.
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcte12
