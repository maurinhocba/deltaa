MODULE cont4_db
  IMPLICIT NONE
  SAVE

  ! Global Variables for contact 4 Algorithm

  INTEGER (kind=4), PARAMETER :: &
     nisdb = 4, &! number of integer values for each slave node
     nrsdb = 9   ! number of real values for each slave node
  INTEGER (kind=4) :: &
     nsurf, &! Number of contact surfaces
     npair   ! Number of contact surfaces pairs
  LOGICAL ::  wear     !.TRUE. if friction work is computed
  REAL (kind=8) :: &
    oldis, & ! maximum displacement expected
    disma, & ! maximum displacement expected
    ffdis, & ! factor to update displacement
    ctime    ! time increment for force computation
  REAL (kind=8), ALLOCATABLE :: &
    surtf(:,:) !(ndime,npair) Total forces on each pair interface
  REAL (kind=8), POINTER :: &
          wwear(:)   !wwear(npoin) friction work

!** Derived type for contact 4 database **!

  TYPE pair4_db
    ! For each contact pair the following variables are allocated

    CHARACTER (len=6) :: &
      pname,  &! pair name
      master, &! master surface name
      slave    ! slave surface name
    INTEGER (kind=4)  :: &
      imast,  &! Master surface internal number
      islav,  &! Slave surface internal number
      indcon, &! contact type 0 1 2 forces on both, slave only, master only
      ncnod,  &! number of nodes in slave surface
      mtsur,  &! bottom, reversed, central or top surface for master
      slsur,  &! bottom, reversed, central or top surface for slave
      freq     ! Frequency of global search
    REAL (kind=8) :: &
      npenal, &! Normal penalty coeff
      tpenal, &! Tangential penalty coeff
      static, &! Static Friction coefficient
            kinet , &! Kinetic Friction coefficient
      cutoff, &! Maximum Addmisible gap
            gapinc, &! Maximum incremental gap
      bhforc, &! Blank Holder Force
      start,  &! Activation time
      end      ! Deactivation time
    LOGICAL :: prev,  & !.TRUE. if pair active in previous step
               press, & !.TRUE. if nodal forces are stored
               wrink    !.TRUE. if wrinkles control is desired

    ! Integer values for slave Data_base
    INTEGER (kind=4), POINTER :: issdb(:,:)  ! ISSDB(nisdb,ncnod)
               !(1) = nears    segment with projection
               !(2) = nearo    previous segment with proy
               !(3) < 0  number of steps without contact
               !    1 penetration and stuck   (Static friccion)
               !    2 penetration and sliding (Kinetic friccion)
               !(4) >= 0         iterative penetration
    REAL (kind=8), POINTER :: presn(:) !presn(ncnod) normal nodal force
    REAL (kind=8), POINTER :: mingp(:) !mingp(ncnod) minimum gap

    ! Real values for slave Data_base
    REAL (kind=8), POINTER :: rssdb(:,:) !RSSDB(nrsdb,ncnod)
                      !(1) = gap
                      !(2) = vnx
                      !(3) = vny
                      !(4) = vnz
                      !(5) = rp       onset local coord
                      !(6) = sp
                      !(7) = rpo      onset local coord (last converged)
                      !(8) = spo
                      !(9) = gap1     converged penetration

    TYPE (pair4_db), POINTER :: next  ! Pointer to next pair

  END TYPE pair4_db

  TYPE (pair4_db), POINTER :: &
      headp,  &! Pointer to first pair
      tailp    ! Pointer to last pair

  ! Derived type for the surface database
  TYPE surf4_db
    CHARACTER (len=6) :: sname    ! surface name
    LOGICAL :: cxc,    & !.TRUE. compute triangle center coordinates
               bottom, & !.TRUE. bottom surface used by some pair
               confor, & !.TRUE. conforming surface
               iwrit,  & ! .T. : Save surface for Post Process
               press,  & ! .T. : binder pressure computed for some pair
               imcod,  & ! Code for master surface
               iscod     ! Code for slave surface
               ! .F., Does not act as a master/slave surf. for any pair
               ! .T., Act as a master/slave surf. for some pair
    LOGICAL :: curved !treat surface as faceted or curved

    INTEGER (kind=4)  :: &
       ncnod,   &! Number of nodes defining a surface
       nsegm     ! Number of segments in the surf
    REAL (kind=8) :: density !surface density, to compute mass

    INTEGER (kind=4), POINTER :: &
       lcnod(:),   &!(ncnod) list of nodes in the surface
       lcseg(:,:), &!(3,nsegm) surface connectivities
       nhseg(:,:), &!(3,nsegm) connected segments to each segment
       lcseb(:,:), &!(3,nsegm) inverted surface connectivities (bottom)
       nhseb(:,:)   !(3,nsegm) inverted connection (bottom)

    REAL (kind=8), POINTER :: &
       xc(:,:),  &!(3,nsegm)  coordinates of the segment center
       cu(:,:),  &!(3,nsegm)  surface curvatures
       tn(:,:)    !(3,ncnod)  normal (outward) at nodes

    TYPE (surf4_db), POINTER :: next        ! pointer to next surface
  END TYPE surf4_db

  TYPE (surf4_db), POINTER :: &
      shead,  &! pointer to first surface
      stail    ! pointer to last surface

CONTAINS

  !************    pair managment routines ************

  SUBROUTINE ini_cont4 (head, tail)
    !Initialize the contact 4 PAIRS database
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair4_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini_cont4

  SUBROUTINE add_pair4 (new, head, tail)
    !This subroutine adds a pair dbase to the end of the list
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair4_db), POINTER :: new, head, tail

       !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE   !add a pair to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF
    RETURN
  END SUBROUTINE add_pair4

  SUBROUTINE srch_pair4 (head, anter, posic, name, found)
    !Searches for a pair named "name"
    IMPLICIT NONE
       !Dummy arguments
    LOGICAL :: found
    CHARACTER (len=6) :: name ! set name
    TYPE (pair4_db), POINTER :: head, anter, posic

    found = .FALSE.
    NULLIFY (posic,anter)
       !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head
      DO
        IF(posic%pname == name) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)
    RETURN
  END SUBROUTINE srch_pair4

  SUBROUTINE del_pair4 (head, tail, anter, posic)
    !Deletes a pair pointed with posic
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair4_db), POINTER :: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    END IF
       !if posic == tail
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
    CALL dalloc_pair4 (posic)
    NULLIFY (anter)
    RETURN
  END SUBROUTINE del_pair4

  SUBROUTINE dalloc_pair4 (pair)
    !Deallocates a pair
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair4_db), POINTER :: pair

    DEALLOCATE ( pair%issdb, pair%rssdb )
    IF( pair%press ) DEALLOCATE( pair%presn )
    IF( pair%wrink ) DEALLOCATE( pair%mingp )
    DEALLOCATE (pair)
    RETURN
  END SUBROUTINE dalloc_pair4

  SUBROUTINE new_pair4 (pair)
    !Allocates a pair
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair4_db), POINTER :: pair

    ALLOCATE (pair)
    NULLIFY ( pair%issdb, pair%rssdb )
    RETURN
  END SUBROUTINE new_pair4

  !************    surface managment routines ************

  SUBROUTINE ini_srf4 (head, tail)
    !Initialize the contact 4 SURFACES database
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini_srf4

  SUBROUTINE ini4_srf (head, tail)
    !Initialize a list of surfaces
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini4_srf

  SUBROUTINE add4_srf (new, head, tail)
    !Adds a surface to the end of the list
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: new, head, tail

       !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE   !add a surface to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF
    RETURN
  END SUBROUTINE add4_srf

  SUBROUTINE srch4_srf (head, anter, posic, name, found)
    !This subroutine searches for a surface named "name"
    IMPLICIT NONE
       !Dummy arguments
    LOGICAL, INTENT(OUT) :: found
    CHARACTER (len=6), INTENT(IN) :: name ! set name
    TYPE (surf4_db), POINTER :: head, anter, posic
    ! INTENT(IN) :: head  begining of the data base
    ! INTENT(OUT) :: posic  pointer to searched surface
    ! INTENT(OUT) :: anter  pointer to previous surface

    found = .FALSE.    !initializes
    NULLIFY (posic,anter)
       !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head
      DO
        IF(posic%sname == name) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)
    RETURN
  END SUBROUTINE srch4_srf

  SUBROUTINE del4_srf (head, tail, anter, posic)
    !Deletes a surface pointed with posic
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: head, tail, anter, posic
    ! INTENT(IN OUT) :: head  begining of the data base
    ! INTENT(IN OUT) :: tail  end of the data base
    ! INTENT(IN) :: posic  pointer to surface to delete
    ! INTENT(IN OUT) :: anter  pointer to previous surface

    IF (.NOT.ASSOCIATED (anter)) THEN   !IF deleled surface is the first
      head => posic%next                !head ==> 2nd surface
    ELSE
      anter%next => posic%next          !previous points to next
    END IF
       !if posic == tail                !If deleted surface is the last
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter  !
    CALL dallo4_srf (posic)   !release memory
    NULLIFY (anter)           !what for ?
    RETURN
  END SUBROUTINE del4_srf

  SUBROUTINE dallo4_srf (surface)
    !Deallocates a surface
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: surface

    IF( ASSOCIATED ( surface%lcnod )) DEALLOCATE ( surface%lcnod )
    IF( ASSOCIATED ( surface%lcseg )) DEALLOCATE ( surface%lcseg )
    IF( ASSOCIATED ( surface%nhseg )) DEALLOCATE ( surface%nhseg )
    IF( ASSOCIATED ( surface%lcseb )) DEALLOCATE ( surface%lcseb )
    IF( ASSOCIATED ( surface%nhseb )) DEALLOCATE ( surface%nhseb )
    IF( ASSOCIATED ( surface%xc    )) DEALLOCATE ( surface%xc    )
    IF( ASSOCIATED ( surface%cu    )) DEALLOCATE ( surface%cu    )
    IF( ASSOCIATED ( surface%tn    )) DEALLOCATE ( surface%tn    )

    DEALLOCATE (surface)
    RETURN
  END SUBROUTINE dallo4_srf

  SUBROUTINE new_surf4 (surf)
    !Allocates a surface
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf4_db), POINTER :: surf

    ALLOCATE (surf)
    surf%bottom = .FALSE.
    surf%press  = .FALSE.
    NULLIFY ( surf%lcnod, surf%lcseg, surf%nhseg, surf%lcseb, &
              surf%nhseb, surf%xc, surf%cu, surf%tn )
    RETURN
  END SUBROUTINE new_surf4

END MODULE cont4_db
