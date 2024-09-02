SUBROUTINE csdat1(ncnod,nsegm,ncnxx,nsexx,nnseg,lcnod,lcseg,iwrit, &
                  label,npoin,coord,sname)

!.... READ node numbers and segment connectivities

USE c_input
USE surf_db
IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: nnseg,iwrit,npoin,                 &
                  label(:),ncnxx,nsexx
INTEGER (kind=4), INTENT(IN OUT) :: ncnod,nsegm
INTEGER (kind=4), INTENT(OUT) :: lcnod(:),lcseg(:,:)
REAL(kind=8), INTENT(IN) :: coord(:,:)
CHARACTER (len=6), INTENT(IN) :: sname
!     Local variables
INTEGER (kind=4) i,ii,j,k,l,chnode,nodes(4),noden(4),isegm
REAL(kind=8) :: c1(3), c2(3), l1, l2

LOGICAL :: istop,found,flag
INTERFACE
  INCLUDE 'elemnt.h'
  INCLUDE 'getnod.h'
END INTERFACE

!.... READ and WRITE node numbers
IF(ncnod > 0)THEN
  CALL listen('CSDAT1')
  IF (.NOT.exists('NODEDE'))                                       &
    CALL runend('CSDAT1:NODE_DESCRIPTION CARD EXPECT')

  CALL rdfrin('CSDAT1',lcnod,ii,ncnxx)
  IF(ii > ncnxx)CALL RUNEND('CSDAT1: Increase CONTACT MEMORY   !')
  ncnod = ii
  CALL listen('CSDAT1')
  IF (.NOT.exists('ENDNOD'))                                       &
      CALL runend('CSDAT1:END_NODE_DESCRIPTION CARD   ')
  IF(iwrit /= 0) THEN
    WRITE(lures,"('      Nodes on contact surface ',a7,/)")sname
    DO i=1,ncnod,8
      WRITE(lures,"(8x,8i9)") lcnod(i:MIN(i+7,ncnod))
    END DO
  END IF
  flag = .TRUE.
END IF
!.... READ and WRITE nodal segment connectivities
IF (nsegm > 0)THEN
  CALL listen('CSDAT1')
  IF (.NOT.exists('ELEMEN'))THEN
     backs = .TRUE.
    ! surface definition based on an element set
    flag = .FALSE.
    found = .FALSE.
    ALLOCATE (surfa)
    surfa%sname = sname
    CALL elemnt ('SURFAC',iwrit=iwrit,flag2=found)
    IF (.NOT.found) CALL runend('RDSURF:ELEMENT SET NOT FOUND       ')
    CALL store_segs (surfa%head,lcseg,nnseg,nsegm)
    CALL dalloc_srf (surfa)

  ELSE
    nsegm = 0
    flag = .TRUE.
    outer : DO
      CALL listen('CSDAT1')
      IF (exists('ENDELE'))EXIT outer
      nsegm = nsegm+1
      IF(nsegm > nsexx)CALL RUNEND('CSDAT1: Increase CONTACT MEMORY   !')
      DO i=1,nnpar
        nodes(i) = INT(param(i))
        noden(i) = chnode(nodes(i))
      END DO
      IF((nnpar < nnseg) .OR. (nnpar > nnseg))THEN
        IF( nnpar == 4 .AND. nnseg == 3)THEN !divide QUAD in 2 triangles
          !first check if no repeated nodes
          DO i=1,nnseg
            DO j=i+1,nnpar
              IF(nodes(i) == nodes(j))THEN
                l = 0
                DO k=1,nnpar
                  IF(k /= j) THEN
                    l = l + 1
                    lcseg(l,nsegm) = nodes(k)
                  END IF
                END DO
                CYCLE outer
              END IF
            END DO
          END DO
          ! choose the smallest diagonal to split
          c1 = coord(:,noden(3)) - coord(:,noden(1))
          c2 = coord(:,noden(4)) - coord(:,noden(2))
          l1 = DOT_PRODUCT(c1,c1)
          l2 = DOT_PRODUCT(c2,c2)
          nsegm = nsegm+1
          IF(nsegm > nsexx)CALL RUNEND('CSDAT1: Increase CONTACT MEMORY   !')
          IF( l1 > l2 )THEN
            lcseg(1:3,nsegm-1) = (/ nodes(1),nodes(2),nodes(4) /)
            lcseg(1:3,nsegm) = nodes(2:4)
          ELSE
            lcseg(1:3,nsegm-1) = nodes(1:3)
            lcseg(1:3,nsegm) = (/ nodes(1), nodes(3), nodes(4) /)
          END IF
        ELSE
          WRITE(lures,"(' nnpar',i5,' nnseg',i5)")nnpar,nnseg
          CALL runend('CSDAT1:erroneous number of nodes   ')
        END IF
      ELSE
        lcseg(1:nnseg,nsegm) = nodes(1:nnseg)
      END IF
    END DO outer
  END IF
  istop = .FALSE.
  !.... WRITE NODAL SEGMENT CONNECTIVITIES
  if (lcseg(1,1) == 0)CALL runend ('CSDAT1: Contact segment missing  ')
  DO isegm=1,nsegm
    IF(iwrit > 0 .AND. flag)THEN
      IF(MOD(isegm,50) == 1) THEN
        WRITE(lures,"(/,'  Nodal Segment  Connectivities ',        &
             &   ' of contact surface ',a7)") sname
        WRITE(lures,"(/,8x,'Segment No.  ',4('Node',i2,5x,:)//)")  &
             &  (i,i=1,nnseg)
      END IF
      WRITE(lures,"(i15,1x,4i11)")isegm,(lcseg(i,isegm),i=1,nnseg)
    END IF
    outer1 : DO i=1,nnseg-1
      k = lcseg(i,isegm)
      DO j=i+1,nnseg
        IF(k == lcseg(j,isegm))THEN
          WRITE(lures,"('repeated nodes in segment',4i10)")lcseg(1:nnseg,isegm)
          istop = .TRUE.
          EXIT outer1
        END IF
      END DO
    END DO outer1
  END DO
  IF(istop) CALL runend('CSDAT1: repeated nodes in segment  ')

  ii = ncnod
  CALL getnod(ncnod,ncnxx,nsegm,nnseg,lcnod,lcseg)
  IF(iwrit /= 0 .AND. ncnod > ii .AND. flag) THEN
    WRITE(lures,"('  Nodes added to surface ',a6,/)")sname
    DO i=ii+1,ncnod,8
      WRITE(lures,"(8x,8i9)") lcnod(i:MIN(i+7,ncnod))
    END DO
  END IF

  IF( flag )THEN
    DO isegm=1,nsegm
      DO i=1,nnseg
        lcseg(i,isegm) = chnode(lcseg(i,isegm))
      END DO
    END DO
  END IF
END IF
CALL listen('CSDAT1')
IF (.NOT.exists('ENDSUR'))CALL runend('CSDAT1:END_SURFACE CARD EXPECTED   ')
IF(flag)THEN
  DO ii=1,ncnod
    lcnod(ii) = chnode(lcnod(ii))
  END DO
END IF

RETURN
END SUBROUTINE csdat1
