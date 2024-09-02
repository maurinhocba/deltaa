 MODULE ift_db
   ! this MODULE contains a data base to manage Temperature boundary conditions

   USE ctrl_db, ONLY : ndoft  !number of temperature DOFs per node
   USE param_db,ONLY: mnam    !lenght of labels

   IMPLICIT NONE
   INTEGER(kind=4), PARAMETER :: nmf=3  !(>= NDOFT)
   SAVE

   !******* for prescribed values (changing in time with a known curve)

  ! node in a list
  TYPE rpt_nod
    INTEGER (kind=4) :: node        !node label
    REAL (kind=8) :: v(nmf)         !fixed temperatures
    TYPE (rpt_nod), POINTER :: next !next position
  END TYPE rpt_nod

  ! set (list of nodes)
  TYPE rpt_set
    CHARACTER (len=mnam) :: lc                !associated curpt label
    INTEGER (kind=4) :: nrv                   !number of nodes in the set
    REAL (kind=8) :: factor                   !curpt factor
    TYPE (rpt_nod), POINTER :: head,tail      !pointers to first and last
    TYPE (rpt_set), POINTER :: next           !pointer to next set
  END TYPE rpt_set

  ! global pointers
  TYPE (rpt_set), POINTER :: headv,tailv   !pointers to first & last sets

  INTEGER (kind=4) :: &
    nprev,            &  !number of prescribed values
    npret                !number of constrained temperature sets

  INTEGER (kind=4),POINTER :: &
    lctmp(:)      ! associated curves to prescribed temp. sets

  REAL (kind=8),POINTER :: &
    prtmp(:,:)     ! prescribed temperatures in time


 CONTAINS


     ! functions for rpt nodes

     SUBROUTINE ini_rptn (head, tail)
       !initialize a list of nodes

       !Dummy arguments
       TYPE (rpt_nod), POINTER :: head, tail

       NULLIFY (head, tail)

     END SUBROUTINE ini_rptn

     SUBROUTINE add_rptn (new, head, tail)
       !This subroutine adds data to the end of the list of nodes
       !Dummy arguments
       TYPE (rpt_nod), POINTER :: new, head, tail

       !Check if a list is empty
       IF (.NOT. ASSOCIATED (head)) THEN
         !list is empty, start it
         head => new
         tail => new
         NULLIFY (tail%next)

       ELSE
         !add a node to the list
         tail%next => new
         NULLIFY (new%next)
         tail => new

       END IF
       RETURN
     END SUBROUTINE add_rptn

     ! functions for rpt sets

     SUBROUTINE ini_rpts (head, tail)
       !initialize a list

       !Dummy arguments
       TYPE (rpt_set), POINTER :: head, tail

       NULLIFY (head, tail)

     END SUBROUTINE ini_rpts

     SUBROUTINE add_rpts (new, head, tail)
       !This subroutine adds a set to the end of the list of sets
       !Dummy arguments
       TYPE (rpt_set), POINTER :: new, head, tail

       !Check if a list is empty
       IF (.NOT. ASSOCIATED (head)) THEN
         !list is empty, start it
         head => new
         tail => new
         NULLIFY (tail%next)

       ELSE
         !add a set to the list
         tail%next => new
         NULLIFY (new%next)
         tail => new

       END IF
       RETURN
     END SUBROUTINE add_rpts

     SUBROUTINE srch_rpts (head, anter, posic, lc, found)
       !This subroutine searches for a set associated to curpt LC
       !Dummy arguments
       LOGICAL :: found            ! Answer
       CHARACTER (len=mnam) :: lc  ! curpt label
       TYPE (rpt_set), POINTER :: head, anter, posic

       found = .FALSE.                   !initializes
       NULLIFY (posic,anter)
       !Check if a list is empty
       IF (ASSOCIATED (head)) THEN
         posic => head                   !point to head of the list
         DO
           IF(posic%lc == lc) THEN       !compare with associated curpt
             found = .TRUE.              !found
             EXIT                        !leave search
           END IF
           IF (ASSOCIATED(posic%next) ) THEN     !more sets to check
             anter => posic                      !keep previous
             posic => posic%next                 !point to next
           ELSE
             EXIT                        !no more sets => exit
           END IF
         END DO
       ENDIF
       IF (.NOT.found) NULLIFY (posic,anter)   !not found => nullify pointers
       RETURN
     END SUBROUTINE srch_rpts

     SUBROUTINE del_rpts (head, tail, anter, posic)
       !deletes a set pointed with posic
       TYPE (rpt_set), POINTER :: head, tail, anter, posic

       IF (.NOT.ASSOCIATED (anter)) THEN  !if first set
         head => posic%next               !new head
       ELSE
         anter%next => posic%next         !skip posic
       END IF
       ! if posic == tail    (last set)
       IF (.NOT.ASSOCIATED (posic%next) ) tail => anter  !new last set
       CALL dalloc_rpts (posic)           !release memory
       NULLIFY (anter)                    !both anter & posic are null now
     END SUBROUTINE del_rpts

     SUBROUTINE dalloc_rpts (rpts)
       ! deallocates a rpt set (release memory)
       TYPE (rpt_set), POINTER :: rpts
       TYPE (rpt_nod), POINTER :: rptn, rptnaux

       rptn => rpts%head    !point to first node
       DO
         IF (.NOT.ASSOCIATED (rptn) ) EXIT
         rptnaux => rptn%next    !keep next pointer
         DEALLOCATE (rptn)       !release memory of the node
         rptn => rptnaux         !point to next
       END DO

       DEALLOCATE (rpts)         !release rest of the vars. lc nrv & factor
       RETURN
     END SUBROUTINE dalloc_rpts

     SUBROUTINE store_rpt (head,vel,nods,ndoft)
       !This subroutine stores the data

       !Dummy argument
       TYPE (rpt_nod), POINTER :: head
       INTEGER (kind=4) :: ndoft, nods(:)
       REAL (kind=8) :: vel(:,:)

       !Local variables
       TYPE (rpt_nod), POINTER :: ptr
       INTEGER :: n

       IF (ASSOCIATED(head))THEN
         ptr => head

         n = 0
         DO
           n = n + 1
           nods(n) = ptr%node
           vel(1:ndoft,n) = ptr%v(1:ndoft)
           ptr => ptr%next
           IF (.NOT.ASSOCIATED(ptr)) EXIT
         END DO
       ENDIF
     END SUBROUTINE store_rpt

 END MODULE ift_db
