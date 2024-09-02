MODULE flc_db
IMPLICIT NONE

  TYPE flc_tp
    INTEGER(kind=4):: lbl,   & !FLC label
                      sttyp, & !strain type
                      npt      !Number of points that define FLC curve
    REAL(kind=8),POINTER:: cv(:,:)   !FLC curve defined by points
    REAL(kind=8):: LmM,      & !Marginal tolerance
                   LmPS,     & !Plane strain tolerance
                   LmLS,     & !Low strain tolerance
                   LmT         !Excessive thinning
    TYPE(flc_tp),POINTER:: next
  END TYPE flc_tp

  !--- General variables
  INTEGER(kind=4):: nflc=0
  INTEGER, ALLOCATABLE :: lblflc(:)
  TYPE(flc_tp),POINTER:: hpflc=>NULL(), tpflc=>NULL()

CONTAINS

  !------------------------------------------------
  SUBROUTINE new_flc(new)
  !Create a new 'flc_tp' pointer
  IMPLICIT NONE
  TYPE(flc_tp),POINTER:: new

    ALLOCATE(new)
    new%lbl = 0
    new%sttyp = 1
    new%npt = 0
    new%LmM  = 0d0
    new%LmPS = 0d0
    new%LmLS = 0d0
    new%LmT  = 0d0
    NULLIFY(new%cv,new%next)

  RETURN
  END SUBROUTINE new_flc

  !------------------------------------------------
  SUBROUTINE add_flc(new,head,tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE
  TYPE(flc_tp),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
    ELSE
      !add a segment to the list
      tail%next => new
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_flc

  !------------------------------------------------
  SUBROUTINE del_flc(elm,head,tail)
  !Delete a list
  TYPE(flc_tp),POINTER:: elm, head, tail
  !Local variables
  TYPE(flc_tp),POINTER:: ant, aux

    NULLIFY(ant)
    aux => head
    DO
      IF (.NOT.ASSOCIATED(aux)) EXIT
      IF (aux%lbl == elm%lbl) EXIT
      ant => aux
      aux => aux%next
    END DO

    IF (ASSOCIATED(aux)) THEN
      IF (.NOT.ASSOCIATED(ant)) THEN
        head => elm%next
      ELSE IF (.NOT.ASSOCIATED(elm%next)) THEN
        tail => ant
        NULLIFY(ant%next)
      ELSE
        ant%next => elm%next
      END IF
      NULLIFY(aux)
    END IF

    IF (ASSOCIATED(elm%cv)) DEALLOCATE(elm%cv)
    DEALLOCATE(elm)

  RETURN
  END SUBROUTINE del_flc

!!  !------------------------------------------------
!!  SUBROUTINE dalloc_flc(head,tail)
!!  !Delete a list
!!  TYPE(flc_tp),POINTER:: head, tail
!!  !Local variables
!!  TYPE(flc_tp),POINTER:: aux
!!
!!    NULLIFY(tail)
!!    DO
!!      IF (.NOT.ASSOCIATED(head)) EXIT
!!      aux => head%next
!!      IF (ASSOCIATED(head%cv)) DEALLOCATE(head%cv)
!!      DEALLOCATE(head)
!!      head => aux
!!    END DO
!!
!!  RETURN
!!  END SUBROUTINE dalloc_flc

  !------------------------------------------------
  SUBROUTINE srch_flc(head,lbl,found,flc)
  !This subroutine searches for a FLC curve labed "lbl"
  LOGICAL,INTENT(OUT):: found
  INTEGER(kind=4),INTENT(IN):: lbl
  TYPE(flc_tp),POINTER:: head, flc

    found = .FALSE.
    !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      flc => head
      DO
        IF (.NOT.ASSOCIATED(flc)) EXIT
        IF (flc%lbl == lbl) THEN
          found = .TRUE.
          EXIT
        END IF
        flc => flc%next
      END DO
    ENDIF
    IF (.NOT.found) NULLIFY(flc)

  RETURN
  END SUBROUTINE srch_flc

END MODULE flc_db
