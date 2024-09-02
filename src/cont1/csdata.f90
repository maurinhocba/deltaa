      SUBROUTINE csdata(maxve,iprint,label,npoin,coord,surf) !,esets,eset)

!.... read contact surface data and generate database

      USE cont1_db

      IMPLICIT NONE
!     arguments
      INTEGER (kind=4), INTENT(IN) :: iprint,npoin,label(:) !,esets,eset(esets)
      INTEGER (kind=4), INTENT(IN) :: maxve
      REAL (kind=8), INTENT(IN) :: coord(:,:)
      TYPE (surf1_db), POINTER :: surf  !INTENT(OUT)
!     local variables
      INTEGER (kind=4), PARAMETER ::  nnseg = 2
      INTEGER (kind=4) ncnod,nsegm,ncnxx,nsexx,maxcc
      CHARACTER (len=6) :: sname
      INTEGER (kind=4), ALLOCATABLE :: lcnod(:),lcseg(:,:)

      INTERFACE
        INCLUDE 'csdat1.h'
        INCLUDE 'mastdb.h'
        INCLUDE 'curve1.h'
      END INTERFACE

!.... set control surface DATA

      maxcc  = maxve               !maximum work space available

      ncnod = surf%ncnod           !flag to read nodes
      nsegm = surf%nsegm           !flag to read segments

!.... READ node numbers and segment connectivities

      ! set maximum values for nodes(ncnxx) and segments (nsexx) to be read
      IF( nsegm > 0 )THEN        !if segments are to be read
        ncnxx = MAX( ncnod, maxcc/2, nsegm/2+100 )
        nsexx = MAX( ncnxx,nsegm)
      ELSE                       !only nodes will be read
        ncnxx = MAX( ncnod, maxcc )
        nsexx = 1
      END IF
      ALLOCATE( lcnod(ncnxx), lcseg(nnseg,nsexx) )  !temporary space

      ! read into lcnod <== nodes   lcseg <== segments
      sname = surf%sname
      CALL csdat1(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iprint, &
                  label,npoin,coord,sname) !surf%sname,esets,eset)

!.... set memory pointers

      IF(surf%iscod == 1 .OR. surf%curved )THEN  !for slave surface keep nodes
        ALLOCATE (surf%lcnod(ncnod))  !reserve space to list of nodes
        surf%lcnod = lcnod(1:ncnod)   !assign list of nodes
      END IF

      IF( nsegm > 0 )THEN          !if segments read keep them
        ALLOCATE ( surf%lcseg(nnseg,nsegm) ) !reserve space for connectivities
        surf%lcseg = lcseg(1:nnseg,1:nsegm)  !assign connectivities
      END IF
      DEALLOCATE( lcnod, lcseg )

!.... search for neighbours segments if the surface will act as a
!.... master surface for some pair ==>  NHSEG

      IF(surf%imcod == 1) THEN
        ALLOCATE( surf%nhseg(nnseg,nsegm) )
        !   form the master surface connection
        CALL mastdb (nsegm,nnseg,surf%lcseg,surf%nhseg,surf%confor,coord)
        ALLOCATE( surf%xc(2,nsegm) )
        surf%cxc = .FALSE.
      END IF

      surf%ncnod = ncnod     !keep number of nodes
      surf%nsegm = nsegm     !keep number of segments

      IF( surf%curved )CALL curve1(coord,surf,npoin)

      RETURN
      END SUBROUTINE csdata
