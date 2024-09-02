 SUBROUTINE prgidb ( )
 !
 !  Writes drawbead affected zones
 !  only for associated sheets presently
 !  for GiD tool in both ASCII and Binary types
 !
  ! give access to all data
  USE data_db
  IMPLICIT NONE
  ! local variable
  TYPE (drawb), POINTER :: eset        !pointer to a drawbead sheet
  INTEGER(kind=4) :: iset,ielem,iel,ic,nelem,is,ik   !different indexes
  REAL(kind=8) :: value                ! 0 not affected  1 affected
  CHARACTER (len=20) :: saux           !variable name
  CHARACTER (len=6) :: caux='      '   !component name
  CHARACTER (len=34) :: gpname         !Gauss Point Set name
  LOGICAL, ALLOCATABLE :: sdaz(:)      !variable at each SHEET element
  INTEGER, ALLOCATABLE :: icode(:)     !associated surface to each db
  LOGICAL :: ssheet                    !same sheet ?

  INTERFACE
    INCLUDE 'cab_gid_bin.h'
    INCLUDE 'cab_gid.h'
  END INTERFACE

  ic = 0                      !initializes sheet identifier
  IF( drawb_sets > 0 )THEN      !if drawbead exist
    ! generate auxiliar array with sheet codes
    ALLOCATE(icode(drawb_sets) )  !get memory
    eset => drawb_head            !point to first
    DO iset=1,drawb_sets             !loop over each drawbead
      icode(iset) = eset%icode       !pass sheet code
      eset => eset%next              !point to next
    END DO
    ! loop over associated sheet
    DO is=1,drawb_shs
      eset => drawb_head            !point to first DB
      DO iset=1,drawb_sets             !loop over each drawbead
        IF( icode(iset) > 0 )THEN          !if an associated element set
          IF( ic == 0 )THEN                  !check
            ic = icode(iset)                   !associate code
            nelem =drawb_shnel(iset)           !number of elements in the sheet
            ALLOCATE( sdaz(nelem) )            !get memory
            sdaz = .FALSE.                     !initializes
            ssheet = .TRUE.                    !this sheet
            gpname = 'GP'//eset%sname          !Gauss Point Name
          ELSE
            ssheet = ic == icode(iset)       !same sheet
          END IF
          IF( ssheet )THEN             !for the associated sheet
            saux = 'Drawbead Mark-'//eset%pname(1:6)
            IF( ip == 2 )THEN  !ASCII type header
              CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,rrt_db%name)
            ELSE               !Binary type header
              CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,rrt_db%name)
            END IF
            DO ielem=1,eset%nelem
              value = 0d0
              IF( eset%daz(ielem) ) value = 1d0
              IF( eset%eall )THEN
                iel = ielem
              ELSE
                iel = eset%ass(ielem)
              END IF
              IF( ip == 2 )THEN   !ASCII type result
                WRITE(13,"(i8,(e15.5))")iel,value
              ELSE                !Binary type result
                CALL GID_WRITESCALAR(iel,value)
              END IF
              sdaz(iel) = sdaz(iel) .OR. eset%daz(ielem)
            END DO
            IF( ip == 2 )THEN
              WRITE(13,"('End Values')")
            ELSE
              CALL GID_ENDRESULT()
            END IF
          END IF
          icode(iset) = 0
        END IF
        eset => eset%next
      END DO
      ! prints values for all the sheet
      IF( ic /= 0 )THEN
        saux = 'Drawbead Mark-'//'ALL_DB'
        IF( ip == 2 )THEN   !ASCII type header
          CALL cab_gid(1,2,saux,caux,1,loadstep,ttime,gpname,rrt_db%name)
        ELSE                !Binary type header
          CALL cab_gid_bin(1,2,saux,caux,1,loadstep,ttime,gpname,rrt_db%name)
        END IF
        DO iel=1,nelem
          value = 0d0
          IF( sdaz(iel) ) value = 1d0
          IF( ip == 2 )THEN !ASCII type result
            WRITE(13,"(i8,(e15.5))")iel,value
          ELSE              !Binary type result
            CALL GID_WRITESCALAR(iel,value)
          END IF
        END DO
        IF( ip == 2 )THEN   !close list
          WRITE(13,"('End Values')")
        ELSE
          CALL GID_ENDRESULT()
        END IF
        ic = 0
        DEALLOCATE(sdaz)    !release memory
      END IF
   END DO
   DEALLOCATE(icode)
 END IF
 RETURN
 END SUBROUTINE
