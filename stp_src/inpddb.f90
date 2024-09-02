 SUBROUTINE inpddb()
 !
 !  read data for a Drawbead sheets
 !
 USE data_db   !give access to all the data
 IMPLICIT NONE

 INTEGER (kind=4) ielem,n,ndb             !counters
 TYPE (drawb), POINTER :: eset,oset       !pointers to drawbeads
 ! for search
 TYPE (bst), POINTER :: bst_e          !pointer to a rotation-free shell
 TYPE (shl3d), POINTER :: shl3d_e      !pointer to a Simot type shell
 LOGICAL :: found                      !search flag

 ALLOCATE (drawb_shnel(drawb_sets))    !number of elements in each sheet
 drawb_shnel = 0                       !initializes
 DO ndb=1,drawb_sets           !loop over all drawbeads

   ALLOCATE (eset)             !get memory for this set
   NULLIFY (eset%next)         !nullify pointer to next set

   IF( ndb == 1 )THEN   !for the first drawbead sheet set
     drawb_head => eset        !associate head
   ELSE                        !for subsequent sets
     drawb_tail%next => eset   !generate link
   END IF

   drawb_tail => eset          !last set position
   ! read drawbead general information
   READ (43) eset%pname,eset%sname,eset%nnode,eset%nelem,eset%icode
   ! get memory for connectivities and affected zone variable
   ALLOCATE(eset%lnods(eset%nnode,eset%nelem),eset%daz(eset%nelem))
   ! read drawbead connectivities
   READ(43) ((eset%lnods(n,ielem),n=1,eset%nnode),ielem=1,eset%nelem)
   eset%daz = .FALSE.  !initializes
   eset%eall = .FALSE. !flag indicates if all the elements in the sheet are included
   ! search if an associated shell type set exist
   IF( eset%icode == 1 ) THEN !if an associated surface set exists
     ! search in rotation-free data base
     found = .FALSE.  !initializes
     eset%icode = 0   !initializes
     IF( bst_sets > 0 ) THEN   !if rotation-free shell exist
       bst_e => bst_head           !point to first set
       DO n=1,bst_sets                 !loop over all the sets
         IF( TRIM(bst_e%sname) == TRIM(eset%sname))THEN  !compare names
           found = .TRUE.                           !sheet found
           eset%icode = n                           !position in list
           IF( bst_e%nelem == eset%nelem ) THEN     !compare number of elements
             eset%eall = .TRUE.            !all the elements included
           ELSE
             eset%eall = .FALSE.           !partially included
             ALLOCATE( eset%ass(eset%nelem))    !get memory for auxiliar array
             ! relate elements included with elements in sheet => ass
             CALL assign_element(eset%lnods,bst_e%lnods,eset%nelem,bst_e%nelem,eset%nnode,3,eset%ass,found)
           END IF
           drawb_shnel(ndb) = bst_e%nelem  !keep number of elements in the sheet
           EXIT               !exit loop if sheet found
         END IF
         bst_e => bst_e%next  !point to next set
       END DO
     END IF
     IF( .NOT.found ) THEN    !if sheet is not a rotation-free set
       ! search in Simo shells data base
       IF( shl3d_sets > 0 ) THEN     !if this type of shell exists
         shl3d_e => shl3d_head          !point to first set
         DO n=1,shl3d_sets                 !loop over all the sets
           IF( TRIM(shl3d_e%sname) == TRIM(eset%sname))THEN  !compare names
             found = .TRUE.                             !sheet found
             eset%icode = n+100                         !position in list
             IF( shl3d_e%nelem == eset%nelem ) THEN     !compare number of elements
               eset%eall = .TRUE.              !all the elements included
             ELSE
               eset%eall = .FALSE.             !partially included
               ALLOCATE( eset%ass(eset%nelem))     !get memory for auxiliar array
               ! relate elements included with elements in sheet => ass
               CALL assign_element(eset%lnods,shl3d_e%lnods,eset%nelem,shl3d_e%nelem,eset%nnode,shl3d_e%nnode,eset%ass,found)
             END IF
             drawb_shnel(ndb) = shl3d_e%nelem  !keep number of elements in the sheet
             EXIT                          !exit loop if found
           END IF
           shl3d_e => shl3d_e%next       !point to next sheet
         END DO
       END IF
     END IF
   END IF

 END DO
 ! compute number of associated and non-associated sheets
 ! only results in associated sheets are printed presently
 drawb_shs = 0      !initializes number of associated sheets
 drawb_shn = 0      !initializes number of non-associated sheets
 eset => drawb_head !point to first drawbead
 DO ndb=1,drawb_sets      !loop
   found = .FALSE.
   IF( eset%icode > 0 )THEN    !if surface is associated
     IF(drawb_shs == 0) THEN   !first associated sheet
       drawb_shs = 1
     ELSE                      !second ..
       oset => drawb_head        !point to head
       DO n=1,ndb-1                 !check previos drawbeads
         IF( oset%icode == eset%icode )THEN    !if sheet already used
           found = .TRUE.                      !exit loop
           EXIT
         ELSE
           oset => oset%next                   !point to next db and continue search
         END IF
       END DO
       IF(.NOT.found )THEN  !if sheet no found among previous db
         drawb_shs = drawb_shs + 1             !increase number of associated sheets
       END IF
     END IF
   ELSE                        !if surface is not associated
     IF(drawb_shn == 0) THEN   ! for the first non-associated sheet
       drawb_shn = 1                !initializes
       eset%icode = -drawb_shn      !use a negative value
     ELSE                      ! for second ..
       oset => drawb_head           !point to head
       DO n=1,ndb-1                   !loop over previous DB
         IF( TRIM(oset%sname) == TRIM(eset%sname) )THEN   !compare sheet names
           found = .TRUE.             !sheet already counted
           eset%icode = oset%icode    !pass the same code
           EXIT                       !exit loop
         ELSE
           oset => oset%next        !point to next and continue search
         END IF
       END DO
       IF(.NOT.found )THEN    ! if sheet was not used by previous db
         drawb_shn = drawb_shn + 1   !increment number of non-associated sheets
         eset%icode = -drawb_shn     !pass the negative order
       END IF
     END IF
   END IF
   eset => eset%next
 END DO
 RETURN
 END SUBROUTINE inpddb
 !---------------------------------------------------------
 SUBROUTINE assign_element(ln1,ln2,nel1,nel2,nn1,nn2,rel,ok)
 !---------------------------------------------------------
 ! relates included elements in a DB with the hole SHEET
 !---------------------------------------------------------
 IMPLICIT NONE
 ! dummy arguments
 INTEGER(kind=4), INTENT(IN) :: nel1,  & !number of included elements in a DB
                                nel2,  & !number of elements in the sheet > NEL1
                                nn1,   & !number of nodes in each db element
                                nn2      !number of nodes in each sheet element
 INTEGER(kind=4), INTENT(IN) :: ln1(nn1,nel1), & ! DB included conns
                                ln2(nn2,nel2)    ! sheet connectivities
 INTEGER(kind=4), POINTER :: rel(:)      !(nel1) relates ln1 with ln2
 LOGICAL, INTENT(OUT) :: ok              !if array REL constructed correctly
 ! local variables
 INTEGER(kind=4) :: ie1,ie2,i
 LOGICAL, ALLOCATABLE :: used(:)
 LOGICAL :: found

 ALLOCATE(used(nel2))    !get auxiliar memory
 used = .FALSE.          !initializes
 rel = 0          !initializes
 DO ie1=1,nel1     !for each DB element
   found = .FALSE.   !initializes flag
   !  search SHEET element with the same nodes
   DO ie2=1,nel2  !for each element
     IF( used(ie2) )CYCLE !element already paired, do not check
     DO i=1,nn1   !for each DB element node
       IF( .NOT.ANY(ln2(:,ie2) == ln1(i,ie1) ))EXIT  !compare with all the nodes
       IF( i == nn1 )found = .TRUE.   !if all nodes have been checked
     END DO
     IF(found) THEN       !if found, pair elements
        rel(ie1) = ie2    !associate element
        used(ie2) = .TRUE.   !do not check this element any longer
        EXIT              !exit loop
     END IF
   END DO
 END DO
 ok = ALL(rel /= 0)   !verify that all elements have been paired
 DEALLOCATE(used)     !release memory
 RETURN

 END SUBROUTINE assign_element
