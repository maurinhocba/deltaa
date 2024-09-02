MODULE bonds_db

  ! scratch database to store connected spheres - contact bonds

  IMPLICIT NONE
  PRIVATE

  ! element (derived type storing 2 connected nodes & contact forces)
  TYPE elem
    PRIVATE
    INTEGER (kind=4) :: nodes(2)
    REAL    (kind=8) :: f(2)
    TYPE (elem), POINTER :: next
  END TYPE elem

  ! list of elements
  TYPE list
    PRIVATE
    INTEGER (kind=4)  :: length
    TYPE (elem), POINTER :: head,tail
    TYPE (elem), POINTER :: posic,anter
  END TYPE list

  PUBLIC :: list,          & !
            add_at_end,    & ! insert an elem at the end of the list
            delete_elem,   & ! deletes the elem pointed with posic
            end_ls,        & ! test End of List for posic
            get_nodes,     & ! gets nodes associated to posic
            get_f,         & ! gets nodes associated to posic
            get_length,    & ! gets number of elements stored in the list
            ls_head,       & ! sets the pointer posic on the head of the list
            ls_ini,        & ! Initializes a list of elem
            ls_next,       & ! moves the pointer posic on the next element of the list
            dalloc_list      ! deallocates the list

CONTAINS

  ! ********   PUBLIC procedures for list of elem  ************


  SUBROUTINE add_at_end (ls, nodes, f)
    ! insert a elem at the end of the list

    !Dummy arguments
    TYPE (list), INTENT(INOUT) :: ls
    INTEGER (kind=4), INTENT(INOUT) :: nodes(2)
    REAL    (kind=8), INTENT(INOUT) :: f(2)

    ! local
    TYPE (elem), POINTER :: new

    new => n_alloc (nodes, f)            !get memory and assign
    IF ( ls%length == 0 ) THEN         !if list is empty
      ls%head => new                   !initializes list
    ELSE
      ls%anter => ls%tail              !keep pointer to previous last elem
      ls%tail%next => new              !new last elem
    END IF
    ls%tail => new                     !updates pointer to last elem
    ls%length = ls%length + 1          !updates number of elem in the list
    ls%posic => new                    !updates pointer to current elem

  END SUBROUTINE add_at_end


  SUBROUTINE delete_elem (ls, iostat)
    ! deletes the elem pointed with posic
    ! posic is redirected on the succeeding elem, or nullified if end_ls

    !Dummy arguments
    TYPE (list), INTENT(IN OUT) :: ls
    INTEGER (kind=4),  INTENT(OUT), OPTIONAL :: iostat

    ! local
    TYPE (elem), POINTER :: addr, addr_next

    IF ( end_ls(ls) ) THEN       !If at the end of the list
      iostat = -1                ! error
    ELSE
      addr => ls%posic           !point to current position
      addr_next => addr%next     !point to succeeding elem
      IF ( .NOT.ASSOCIATED (ls%anter) ) THEN       !If no previous
        ls%head => addr_next                       ! head  of the list
      ELSE
        ls%anter%next => addr_next                 !skip posic
      END IF
      ls%posic => addr_next      ! posic redirected on the succeeding elem
      IF ( .NOT.ASSOCIATED (addr_next) )  &  !IF last elem in the list deleted
        ls%tail => ls%anter                  !point tail to anter
      ls%length = ls%length - 1              !update number of elem in list
      CALL n_free(addr)                      !release memory
      IF (PRESENT(iostat)) iostat = 0        !O.K.
    END IF

  END SUBROUTINE delete_elem


  LOGICAL FUNCTION end_ls (ls)
    ! test End of List for posic

    !Dummy arguments
    TYPE (list), INTENT(IN) :: ls

    end_ls = .NOT.ASSOCIATED (ls%posic)

  END FUNCTION end_ls


  FUNCTION get_nodes (ls, iostat)
    ! gets nodes associated to posic

    INTEGER (kind=4), DIMENSION(2) :: get_nodes
    !Dummy arguments
    TYPE (list), INTENT(IN) :: ls
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_ls (ls) ) THEN              !if at the end of the list
      IF (PRESENT(iostat)) iostat = -1  !error
    ELSE
      get_nodes = ls%posic%nodes        !nodes
      IF (PRESENT(iostat)) iostat = 0   !O.K.
    END IF

  END FUNCTION get_nodes


  FUNCTION get_f (ls, iostat)
    ! gets contact forces associated to posic

    REAL (kind=8), DIMENSION(2) :: get_f
    !Dummy arguments
    TYPE (list), INTENT(IN) :: ls
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_ls (ls) ) THEN              !if at the end of the list
      IF (PRESENT(iostat)) iostat = -1  !error
    ELSE
      get_f = ls%posic%f                !contact forces
      IF (PRESENT(iostat)) iostat = 0   !O.K.
    END IF

  END FUNCTION get_f


  FUNCTION get_length (ls)
    ! gets number of elements stored in the list

    INTEGER (kind=4) :: get_length
    !Dummy arguments
    TYPE (list), INTENT(IN) :: ls

    get_length = ls%length

  END FUNCTION get_length


  SUBROUTINE ls_head (ls)
    ! sets the pointer posic on the head of the list

    !Dummy arguments
    TYPE (list), INTENT(IN OUT) :: ls

    NULLIFY (ls%anter)   ! head has no preceding element
    ls%posic => ls%head  ! point to first

  END SUBROUTINE ls_head


  SUBROUTINE ls_ini (ls)
    ! Initializes a list of elem
    !Dummy arguments
    TYPE (list), POINTER :: ls

    ALLOCATE (ls)              !get memory
    NULLIFY (ls%head, ls%tail, ls%anter, ls%posic) !initializes pointers
    ls%length = 0              !no elem yet

  END SUBROUTINE ls_ini


  SUBROUTINE ls_next (ls, iostat)
    ! moves the pointer posic on the next element of the list

    !Dummy arguments
    TYPE (list), INTENT(IN OUT) :: ls
    INTEGER (kind=4), INTENT(OUT), OPTIONAL :: iostat

    IF (end_ls (ls) ) THEN        !if at the end of the list
      iostat = -1                 !error
    ELSE
      ls%anter => ls%posic        !keep previous position
      ls%posic => ls%posic%next   !point to next
      IF (PRESENT(iostat)) iostat = 0  !O.K.
    END IF

  END SUBROUTINE ls_next


  SUBROUTINE dalloc_list (ls)
    ! deallocates all the list - deleting element from end to head

    !Dummy arguments
    TYPE (list), INTENT(INOUT) :: ls

    ! local
    INTEGER (kind=4) :: i, iostat, length

    length = ls%length
    ls%posic => ls%head  
    DO i = 1,length
      CALL delete_elem (ls, iostat)
      IF ( iostat == -1 ) EXIT
    END DO

  END SUBROUTINE dalloc_list




  ! *********** PRIVATE Procedures for list of elem  ************

  FUNCTION n_alloc (nodes, f)
    ! allocates a elem and assign value

    !Dummy arguments
    TYPE (elem),  POINTER        :: n_alloc
    INTEGER (kind=4), INTENT(INOUT) :: nodes(2)
    REAL    (kind=8), INTENT(INOUT) :: f(2)

    ALLOCATE (n_alloc)              ! gets memory
    n_alloc%nodes = nodes           ! assign value
    n_alloc%f     = f               ! assign value
    NULLIFY (n_alloc%next)

  END FUNCTION n_alloc


  SUBROUTINE n_free (adr)
    ! frees memory used with one elem

    !Dummy arguments
    TYPE (elem),  POINTER :: adr

    DEALLOCATE (adr)

  END SUBROUTINE n_free


  LOGICAL FUNCTION empty (ls)
    ! test if List is empty

    !Dummy arguments
    TYPE (list), INTENT(IN) :: ls

    empty = ls%length == 0       ! no elem ==> TRUE

  END FUNCTION empty


END MODULE bonds_db
