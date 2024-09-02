SUBROUTINE pcggib05( )
!  Write mesh information for GiD for 3-D-Solids
USE data_db,ONLY: sol3d, sol3d_head, sol3d_sets, sol3d_nvarg,given
IMPLICIT NONE

  TYPE(sol3d),POINTER:: e
  INTEGER(kind=4):: etype,iset,ipos,ngaus,g
  REAL(kind=8) :: a,b
  LOGICAL:: gauss
  CHARACTER(len=32):: sname
  CHARACTER(len=34):: gpname
  REAL(kind=8), ALLOCATABLE :: p(:),w(:)

  a = 1d0/SQRT(3D0)
  b = 1d0/3D0

  gauss =  sol3d_nvarg > 0
  IF (.NOT.gauss) RETURN

  e => sol3d_head                 !point to first element set
  DO iset=1,sol3d_sets            !for each element set
    sname = TRIM(e%sname)
    gpname = 'GP'//TRIM(e%sname)

    !-----------------------------------------------------------------------------
    ipos = 1             !internal
    ngaus = e%ngaus      !
    SELECT CASE (e%nnode)
    CASE(4,10)          !Tetrahedra
      etype = 5

    CASE(8,20)          !Hexahedra
      etype = 6
      IF( e%ngaus == 4 .OR.  e%ngaus == 2) THEN
        IF( given ) THEN
          ipos = 0  !given
        ELSE
          ngaus = 8
        END IF
      END IF

    CASE(6,15)            !Prism
      etype = 7
      IF( e%ngaus == 2)THEN
        IF(given) THEN
          ipos = 0
        ELSE
          ngaus = 6
        END IF
      ELSE
        ipos = 0
        ALLOCATE( p(ngaus),w(ngaus))
        CALL gaussq(ngaus,p,w)
      END IF
    END SELECT

    CALL GID_BEGINGAUSSPOINT(TRIM(gpname),etype,TRIM(sname),ngaus,0,ipos)

    IF( ipos == 0 )THEN    !non-standard positions
      IF( etype == 6  .AND.  ngaus == 4 )THEN     !Hexahedra with 4 GP
        CALL GiD_WriteGaussPoint3D( a ,-a ,-a)
        CALL GiD_WriteGaussPoint3D(-a , a ,-a)
        CALL GiD_WriteGaussPoint3D(-a ,-a , a)
        CALL GiD_WriteGaussPoint3D( a , a , a)
      ELSE IF( ngaus == 2 )THEN                   !2 GP along element axis
        IF ( etype == 6 )THEN      !nnode = 8       !Hexahedra
          CALL GiD_WriteGaussPoint3D( 0 , 0 ,-1d0)    !bottom Face
          CALL GiD_WriteGaussPoint3D( 0 , 0 , 1d0)    !top Face
        ELSE !IF(etype == 7 )THEN  !nnode = 6       !Prisma
          IF( e%etype == 16 )THEN       !PRISM
            CALL GiD_WriteGaussPoint3D( b , b ,-a)
            CALL GiD_WriteGaussPoint3D( b , b , a)
          ELSE ! IF( e%etype == 5 )THEN       !SPRISM
            CALL GiD_WriteGaussPoint3D( b , b ,-1d0)
            CALL GiD_WriteGaussPoint3D( b , b , 1d0)
          END IF
        END IF
      END IF
    END IF
    CALL GID_ENDGAUSSPOINT()

    e => e%next
  END DO

RETURN
END SUBROUTINE pcggib05
