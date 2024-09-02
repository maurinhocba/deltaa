SUBROUTINE elmdat( input,leng )
!       Read element data
USE data_db
IMPLICIT NONE

  !Dummy variables
  INTEGER(kind=4),INTENT(IN):: leng
  CHARACTER(len=*):: input

  !Local variables
  INTEGER(kind=4)::  etype, nel, iset, iel, i, j, n, ist, ng
  REAL (kind=8) :: x0(ndime,2),f,df
  CHARACTER(len=30):: sname
  TYPE(spot ),POINTER:: espot
  TYPE(truss),POINTER:: etruss
  TYPE(sol2d),POINTER:: esol2d
  TYPE(sol3d),POINTER:: esol3d
  TYPE(shl3d),POINTER:: eshl3d
  TYPE(beame),POINTER:: ebeame
  TYPE(shrev),POINTER:: eshrev
  TYPE(rigid),POINTER:: erigid
  TYPE(bst  ),POINTER:: ebst

  ! initializes sets
  spot_sets  = 0
  truss_sets = 0
  sol2d_sets = 0
  sol3d_sets = 0
  shl3d_sets = 0
  beame_sets = 0
  shrev_sets = 0
  bst_sets   = 0
  rigid_sets = 0
  setyp = 0
  nsets = 0

  ntype = 0      !initializes 2-D problem type
  ngp = 0

  DO
    READ(17) etype, nel, sname !read element type and number of elements in the set
    IF(etype == 0)EXIT         !etype = 0 is a flag to exit loop of element sets

    IF(nel == 0)CYCLE         !this is unnecessary
    nsets = nsets + 1         !increase number of sets
    setyp(nsets) = etype      !keep element type

    SELECT CASE ( etype )
    CASE (1)      ! read initial data for SPOT elements
      CALL inpda1 (nel,sname)

    CASE (2)      ! read initial data for TRUSS elements
      CALL inpda2 (nel,sname)

    CASE (17,19,20)  ! read initial data for 2-D SOLID elements
      CALL inpda3 (nel,sname,etype)

    CASE (4,5,12,16,18,27,29)  ! read initial data for 3-D SOLID elements
      CALL inpda5 (nel,sname,etype)

    CASE (6,7)    ! read initial data for 3-D SHELL (Simo Theory) elements
      CALL inpda6 (nel,sname)

    CASE (8)      ! read initial data for 3-D BEAM (Simo Theory) element
      CALL inpda8 (nel,sname)

    CASE (9,11)   ! read initial data for 2-D SHELL/BEAM (Simo Theory) element
      CALL inpda9 (nel,sname)
      IF( ntype == 4 ) ip=2 !GID-ASCII for linear 2-D frames
    CASE (10)     ! read initial data for RIGID element
      CALL inpd10 (nel,sname,.FALSE.)

    CASE (13:15,25)  ! read initial data for Rotation-free SHELL element
      CALL inpd12 (nel,sname,etype)

    END SELECT
  END DO

  IF( rigid_show )THEN
    OPEN(18,FILE=input(1:leng)//'.f18',STATUS='OLD',IOSTAT=ist, &
         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    IF (ist == 0) THEN
      DO
        READ(18,IOSTAT=ist) nel,sname  !read number of elements in the set, and Surf_name
        IF (ist /= 0) EXIT
        CALL inpd10 (nel,sname, .TRUE. )
      END DO
      CLOSE(18)
    END IF
  END IF

  IF( drawb_wafz )THEN
    OPEN(43,FILE=input(1:leng)//'.dbl',STATUS='OLD',IOSTAT=ist, &
         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    IF (ist == 0) THEN
      READ(43,IOSTAT=ist) drawb_sets !read number of drawbead sets
      CALL inpddb ()              !read drawbead sheet data
    END IF
  END IF

  IF( spot_sets > 0 )THEN
    ALLOCATE( spot_nodes(npoin,2) )  !get memory for auxiliar array
    spot_nodes = 0                 !initializes array
    espot => spot_head
    DO iset=1,spot_sets            !for each element set
      IF( espot%nnode == 2 )THEN     !only for two nodes spots
        DO iel = 1,espot%nelem            !for each element in the set
          !                           add nodes to the list
          DO i=1,espot%nnode          !for each node in the element
            n = espot%lnods(i,iel)    !global node
            IF( n < 0 )CYCLE             !for one node spot
            IF( spot_nodes(n,1) == 0 )THEN   !if node not considered yet
              spot_nodes(n,1) = 1            !set 1 as a flag
              spot_nodes(n,2) = iset         !keep set
            END IF
          END DO
        END DO
      END IF
      espot => espot%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(spot_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the spot mesh
      spot_nodes(n,1) = i               !local number
    END DO
    spot_nps = i                      !transfer number of nodes in the mesh

  END IF

  IF( truss_sets > 0 )THEN
    ALLOCATE( truss_nodes(npoin,2) )  !get memory for auxiliar array
    truss_nodes = 0                 !initializes array
    etruss => truss_head
    DO iset=1,truss_sets            !for each element set
      DO iel = 1,etruss%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,etruss%nnode          !for each node in the element
          n = etruss%lnods(i,iel)    !global node
          IF( truss_nodes(n,1) == 0 )THEN   !if node not considered yet
            truss_nodes(n,1) = 1            !set 1 as a flag
            truss_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      etruss => etruss%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(truss_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the truss mesh
      truss_nodes(n,1) = i               !local number
    END DO
    truss_nps = i                      !transfer number of nodes in the mesh
    IF(truss_nps > 0 .AND. truss_nvarn > 0 )  &
         ALLOCATE( truss_accpn(truss_nps), truss_vargs(truss_nvarn,truss_nps))
  END IF

  IF( sol2d_sets > 0 )THEN
    ALLOCATE( sol2d_nodes(npoin,2) )  !get memory for auxiliar array
    sol2d_nodes = 0                 !initializes array
    esol2d => sol2d_head
    DO iset=1,sol2d_sets            !for each element set
      DO iel = 1,esol2d%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,esol2d%nnode          !for each node in the element
          n = esol2d%lnods(i,iel)    !global node
          IF( sol2d_nodes(n,1) == 0 )THEN   !if node not considered yet
            sol2d_nodes(n,1) = 1            !set 1 as a flag
            sol2d_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      esol2d => esol2d%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(sol2d_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      sol2d_nodes(n,1) = i               !local number
    END DO
    sol2d_nps = i                      !transfer number of nodes in the mesh
    IF(sol2d_nps > 0 .AND. sol2d_nvarn > 0 )  &
          ALLOCATE( sol2d_accpn(sol2d_nps), sol2d_vargs(sol2d_nvarn,sol2d_nps))
    ! compute mid-side node coordinates from control points
    IF( ASSOCIATED( sol2d_cpx) )THEN
      DO n=1,npoin                       !for each node in the global mesh
        IF(sol2d_cpx(1,n) == 0 )CYCLE    !It is not a mid-side node
        coord(:,n) =  coord(:,n)/2d0 + coord(:,sol2d_cpx(1,n))/4d0 + coord(:,sol2d_cpx(2,n))/4d0
        coors(:,n) =  coors(:,n)/2d0 + coors(:,sol2d_cpx(1,n))/4d0 + coors(:,sol2d_cpx(2,n))/4d0
      END DO
    END IF
  END IF

  IF( sol3d_sets > 0 )THEN
    ALLOCATE( sol3d_nodes(npoin,2) )  !get memory for auxiliar array
    sol3d_nodes = 0                 !initializes array
    esol3d => sol3d_head
    DO iset=1,sol3d_sets            !for each element set
      DO iel = 1,esol3d%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,esol3d%nnode          !for each node in the element
          n = esol3d%lnods(i,iel)    !global node
          IF( n > npoin ) CYCLE
          IF( sol3d_nodes(n,1) == 0 )THEN   !if node not considered yet
            sol3d_nodes(n,1) = 1            !set 1 as a flag
            sol3d_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      esol3d => esol3d%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(sol3d_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      sol3d_nodes(n,1) = i               !local number
    END DO
    sol3d_nps = i                      !transfer number of nodes in the mesh
    IF(sol3d_nps > 0 .AND. sol3d_nvarn > 0 )  &
          ALLOCATE( sol3d_accpn(sol3d_nps), sol3d_vargs(sol3d_nvarn,sol3d_nps))
  END IF

  IF( shl3d_sets > 0 )THEN
    ALLOCATE( shl3d_nodes(npoin,2) )  !get memory for auxiliar array
    shl3d_nodes = 0                 !initializes array
    eshl3d => shl3d_head
    DO iset=1,shl3d_sets            !for each element set
      DO iel = 1,eshl3d%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,eshl3d%nnode          !for each node in the element
          n = eshl3d%lnods(i,iel)    !global node
          IF( shl3d_nodes(n,1) == 0 )THEN   !if node not considered yet
            shl3d_nodes(n,1) = 1            !set 1 as a flag
            shl3d_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      eshl3d => eshl3d%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(shl3d_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      shl3d_nodes(n,1) = i               !local number
    END DO
    shl3d_nps = i                      !transfer number of nodes in the mesh
    IF(shl3d_nps > 0 .AND. shl3d_nvarn > 0 )  &
          ALLOCATE( shl3d_accpn(shl3d_nps), shl3d_vargs(shl3d_nvarn,shl3d_nps))
  END IF

  IF( beame_sets > 0 )THEN
    ALLOCATE( beame_nodes(npoin,2) )  !get memory for auxiliar array
    beame_nodes = 0                 !initializes array
    ebeame => beame_head
    DO iset=1,beame_sets            !for each element set
      DO iel = 1,ebeame%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,ebeame%nnode          !for each node in the element
          n = ebeame%lnods(i,iel)    !global node
          IF( beame_nodes(n,1) == 0 )THEN   !if node not considered yet
            beame_nodes(n,1) = 1            !set 1 as a flag
            beame_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      ebeame => ebeame%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(beame_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      beame_nodes(n,1) = i               !local number
    END DO
    beame_nps = i                      !transfer number of nodes in the mesh
    IF(beame_nps > 0 .AND. beame_nvarn > 0 )  &
          ALLOCATE( beame_accpn(beame_nps), beame_vargs(beame_nvarn,beame_nps))
  END IF

  IF( shrev_sets > 0 )THEN
    ALLOCATE( shrev_nodes(npoin,2) )  !get memory for auxiliar array
    shrev_nodes = 0                 !initializes array
    eshrev => shrev_head
    DO iset=1,shrev_sets            !for each element set
      DO iel = 1,eshrev%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,eshrev%nnode          !for each node in the element
          n = eshrev%lnods(i,iel)    !global node
          IF( shrev_nodes(n,1) == 0 )THEN   !if node not considered yet
            shrev_nodes(n,1) = 1            !set 1 as a flag
            shrev_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      eshrev => eshrev%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(shrev_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      shrev_nodes(n,1) = i               !local number
    END DO
    shrev_nps = i                      !transfer number of nodes in the mesh
    IF(shrev_nps > 0 .AND. shrev_nvarn > 0 )  &
          ALLOCATE( shrev_accpn(shrev_nps), shrev_vargs(shrev_nvarn,shrev_nps))
    IF(ntype == 4 )THEN
      IF( ngp > 0 ) ALLOCATE (coora(ndime,ngp), labea(ngp), dispa(ndime,ngp))
      n = MAXVAL(label)
      eshrev => shrev_head
      i = 0
      DO iset=1,shrev_sets            !for each element set
        ng = eshrev%ngaus
        df = 1d0/(ng-1)
        DO iel = 1,eshrev%nelem            !for each element in the set
          x0 = coord(:,eshrev%lnods(1:2,iel))
          x0(:,2) = x0(:,2) - x0(:,1)        !axial vector
          eshrev%lnods(3,iel) = i+1          !additional nodes
          eshrev%lnods(4,iel) = i+ng
          f = 0d0
          DO j=1,ng
            i=i+1
            coora(:,i) = x0(:,1) + f*x0(:,2)
            f = f+df
            labea(i) = n+j
          END DO
          n = n+ng
        END DO
        eshrev => eshrev%next
      END DO
    END IF
  END IF

  IF( rigid_sets > 0 )THEN
    ALLOCATE( rigid_nodes(npoin,2) )  !get memory for auxiliar array
    rigid_nodes = 0                 !initializes array
    erigid => rigid_head
    DO iset=1,rigid_sets            !for each element set
      DO iel = 1,erigid%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,erigid%nnode          !for each node in the element
          n = erigid%lnods(i,iel)    !global node
          IF( rigid_nodes(n,1) == 0 )THEN   !if node not considered yet
            rigid_nodes(n,1) = 1            !set 1 as a flag
            rigid_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      n = erigid%nmast           !global node
      IF( n > 0 )THEN
        IF( rigid_nodes(n,1) == 0 )THEN   !if node not considered yet
          rigid_nodes(n,1) = 1            !set 1 as a flag
          rigid_nodes(n,2) = iset+10000   !keep set
        END IF
      END IF
      erigid => erigid%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(rigid_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
      rigid_nodes(n,1) = i               !local number
    END DO
    rigid_nps = i                      !transfer number of nodes in the mesh
  END IF

  IF( bst_sets > 0 )THEN
    ALLOCATE(   bst_nodes(npoin,2) )  !get memory for auxiliar array
    bst_nodes = 0                 !initializes array
    ebst => bst_head
    DO iset=1,bst_sets            !for each element set
      DO iel = 1,ebst%nelem            !for each element in the set
        !                           add nodes to the list
        DO i=1,ebst%nnode          !for each node in the element
          n = ebst%lnods(i,iel)    !global node
          IF(   bst_nodes(n,1) == 0 )THEN   !if node not considered yet
              bst_nodes(n,1) = 1            !set 1 as a flag
              bst_nodes(n,2) = iset         !keep set
          END IF
        END DO
      END DO
      ebst => ebst%next
    END DO
    ! compute local numeration for this element type
    i = 0                              !initializes
    DO n=1,npoin                       !for each node in the global mesh
      IF(  bst_nodes(n,1) == 0 )CYCLE    !does not belong to this mesh
      i = i + 1                        !increase nodes in the 2-D solid mesh
        bst_nodes(n,1) = i               !local number
    END DO
    bst_nps = i                      !transfer number of nodes in the mesh
    IF(bst_nps > 0 .AND. bst_nvarn > 0 )  &
          ALLOCATE( bst_accpn(bst_nps), bst_vargs(bst_nvarn,bst_nps))
  END IF

RETURN
END SUBROUTINE elmdat
