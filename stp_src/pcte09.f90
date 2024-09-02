 SUBROUTINE pcte09 (flag,iload)
 !
 !  Write mesh information for TECPLOT for 2_D Shell elements
 !
 USE data_db
 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: flag  !.FALSE. mesh only, .TRUE. results included
 INTEGER(kind=4), INTENT(IN) :: iload

 TYPE( shrev ), POINTER :: e
 INTEGER :: ngv,i,n,iset,l,iel,nlgv,nelem
 CHARACTER (len=15) :: etype = 'FELINESEG      '
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

 CALL tecgetgv(flag,gv,shrev_nodes,shrev_nps,ngv)  !get global variables
 IF( flag) CALL tecgenvarshls(7)   !Generate Variables Share List and Passive Variables list
 nlgv = t_ev(2,7)                  !number of Local Gauss Variables

 e => shrev_head                 !point to first element set
 one = .TRUE.                    !true for first set of shrev elemenents
 DO iset=1,shrev_sets            !for each element set
   IF( e%nnode == 2)THEN
     nelem = e%nelem
   ELSE                          !
     nelem = 2*nelem
   END IF

   l = LEN_TRIM(e%sname)
   sname(1:l+2) = '"'//e%sname(1:l)//'"'     !set name with added apostrophes
   WRITE(11,"('ZONE T=',a,' Nodes=',i6,' Elements=',i6,' ZONETYPE=',a15,    &
              ' DATAPACKING=BLOCK  SOLUTIONTIME=',e12.4 )") &  !,' STRANDID=',i4
               sname(1:l+2),shrev_nps, nelem, etype, ttime     !, iload

   IF( flag )THEN !SOLUTION
     IF( t_npv > 0) WRITE(11,"(a)")TRIM(ps_var)  !print list of passive variables
     IF( nlgv > 0) WRITE(11,"(a)")TRIM(lc_var)  !print variables location
     IF( one )THEN                  !only print for first set
       DO i=1,ngv    !print global variables
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,shrev_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DO i=1,shrev_nvarn !print smoothed variables
         WRITE(11,"(10e15.6)")(shrev_vargs(i,n),n=1,shrev_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO

     ELSE
       WRITE(11,"(a)")TRIM(sh_var)  !print lis of common variables
     END IF
     ! print Gaussian Variables
     IF( nlgv > 0 )THEN
       ALLOCATE( lgv(nelem) )   !get memory for average element values
       DO i=1,nlgv                ! for each variabls
         IF( e%nnode == 2 )THEN
           DO iel=1,e%nelem         ! compute cell-center values
             lgv(iel) = SUM(e%elvar(i,:,iel))/e%ngaus
           END DO
         ELSE
           l = 0
           DO iel=1,e%nelem         ! compute cell-center values
             IF( e%ngaus == 2 )THEN
               lgv(l+1) = e%elvar(i,1,iel)
               lgv(l+2) = e%elvar(i,2,iel)
             ELSE
               lgv(l+1) = 0.65*e%elvar(i,1,iel) + 0.35*e%elvar(i,2,iel)
               lgv(l+2) = 0.35*e%elvar(i,2,iel) + 0.65*e%elvar(i,3,iel)
             END IF
             l = l+1
           END DO
         END IF
         WRITE(11,"(10e15.6)") (lgv(iel),iel=1,nelem)  !print gauss values
         WRITE(11,"()")            !leave a line between data blocks
       END DO
       DEALLOCATE( lgv )
     END IF

   ELSE           !GRID
     IF( one )THEN               !only print for first set
       DO i=1,ngv
         WRITE(11,"(10e15.6)")(gv(i,n),n=1,shrev_nps)
         WRITE(11,"()")            !leave a line between data blocks
       END DO
     ELSE
       WRITE(11,*)TRIM(sh_var)
     END IF
     DO iel = 1,nelem            !for each element in the set
       IF( e%nnode == 2 )THEN
         WRITE(11,"(4i8)") shrev_nodes(e%lnods(:,iel),1)  !print local numeration
       ELSE
         WRITE(11,"(4i8)") shrev_nodes(e%lnods(1:2,iel),1)  !print local numeration
         WRITE(11,"(4i8)") shrev_nodes(e%lnods(2:3,iel),1)  !print local numeration
       END IF
     END DO
   END IF
   one = .FALSE.
   e => e%next
 END DO

 RETURN

 END SUBROUTINE pcte09
