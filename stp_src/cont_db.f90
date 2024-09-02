MODULE cont_db
  IMPLICIT NONE

  ! Global Variables for contact

  INTEGER :: &
     nsurf  ! Number of contact surfaces
  LOGICAL :: wwear,  & !flag to print friction work
             wwrink, & !flag to print wrinkles (gaps) information
             wpress, & !flag to print binder press
             press,  & !flag to read binder press information
             wear,   & !flag to read tools wearing information
             wrink     !flag to read contact gap information

  CHARACTER(len=15) :: wear_l(3),  wear_u,  & !label for friction work
                       wrink_l(4), wrink_u, & !label for wrinkles (gaps) information
                       press_l(3), press_u    !label for binder press

  REAL (kind=8), ALLOCATABLE :: wwork(:),   & !(npoin) friction work
                                areas(:),   & !(npoin) nodal areas
                                presn(:,:) ,& !(2,npoin) nodal press 1:     2:
                                wrinkn(:,:)   !(2,npoin) nodal gap   1:     2:

  REAL (kind=8)              :: wear_f,  &  !factor for friction work
                                wrink_f, &  !factor for wrinkles (gaps) information
                                press_f     !factor for binder press


!** Derived type for contact  database **!


  ! Derived type for the surface database
  TYPE surf_db
    CHARACTER (len=30) :: sname  !surface name

    LOGICAL :: press, & !read contact press for this surface
               cpress,& !not a binder press for this surface
               wrink    !read wrinkles (gap) info for this surface

    INTEGER :: &
       ncnod,   &! Number of nodes defining a surface
       nsegm     ! Number of segments in the surf

    INTEGER, POINTER :: &
       lcnod(:),   & !(ncnod) list of nodes in the surface
       lcseg(:,:)    !(3,nsegm) surface connectivities

    TYPE (surf_db), POINTER :: next        ! pointer to next surface
  END TYPE surf_db

  TYPE (surf_db), POINTER :: &
      shead,  & ! pointer to first surface
      stail     ! pointer to last surface

CONTAINS


  !************    surface managment routines ************

  SUBROUTINE ini_srf (head, tail)
    !Initialize the contact SURFACES database
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini_srf

  SUBROUTINE add_srf (new, head, tail)
    !Adds a surface to the end of the list
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf_db), POINTER :: new, head, tail

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
  END SUBROUTINE add_srf

  SUBROUTINE new_surf (surf)
    !Allocates a surface
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf_db), POINTER :: surf

    ALLOCATE (surf)
    NULLIFY ( surf%lcnod, surf%lcseg )
    RETURN
  END SUBROUTINE new_surf

  SUBROUTINE readn (nsurf)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nsurf

    INTEGER :: i
    TYPE (surf_db), POINTER :: surf

    surf => shead
    DO i=1,nsurf
      IF( surf%press )READ(44)
      IF( surf%wrink )READ(44)
      IF( surf%wrink )READ(44)
      surf => surf%next
    END DO
  RETURN

  END SUBROUTINE readn

END MODULE cont_db
