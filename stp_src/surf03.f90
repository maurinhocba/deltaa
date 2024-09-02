      SUBROUTINE surf03 ( lnods, nnode, nelem, bound, ne)
      !******************************************************************
      !
      !*** get the boundary definition from the element set
      !    segments are oriented (outward normal)
      !
      !******************************************************************
      IMPLICIT NONE
      !dummy arguments
      INTEGER (kind=4), INTENT(IN) :: nelem      !number of elements
      INTEGER (kind=4), INTENT(IN) :: nnode      !number of nodes/element
      INTEGER (kind=4), INTENT(IN) :: lnods(nnode,nelem) !connectivities
      INTEGER (kind=4), POINTER :: bound(:,:)
      INTEGER (kind=4), INTENT(OUT) :: ne

      ! local variables
      INTEGER (kind=4) :: in,ie,n1,n2,je,jn,m1,m2
      INTEGER (kind=4), ALLOCATABLE :: lside(:,:),order(:,:)

      ALLOCATE (lside(nnode,nelem))  !this arrays stores -1 if inner side
                                     !  > 1 if part of the boundary
      lside = 0                      ! initializes to unknown

      !   SEARCH FOR THE 'TWIN' SIDES
      ne = 0                     !initializes number of boundary sides
      DO ie=1,nelem              !for each element
        jsearch : DO in=1,nnode        !for each side in the element
          IF( lside(in,ie) /=0 )CYCLE  !if processed side then cycle
          n1=lnods(in,ie)          !first node of the side
          n2=lnods(MOD(in,nnode)+1,ie)      !second node of the side
          IF( n1 == n2 ) THEN      !if both node are the same ==> not a side
            lside(in,ie) = -1      !set as if where an inner side
            CYCLE                  !next side
          END IF
          !   search for the 'twin' side of n1-n2
          DO je=ie+1,nelem  !for each remaining element
            DO jn=1,nnode                 !for each side
              IF( lside(jn,je) /= 0 )CYCLE   !if processed side then cycle
              m1=lnods(jn,je)                  !first node of the side
              m2=lnods(MOD(jn,nnode)+1,je)     !second node of the side
              IF( m1 == n2 .AND. m2 == n1) THEN  !if this is the TWIN side
                lside(in,ie) = -1         !set to Inner side
                lside(jn,je) = -1         !set to Inner side
                CYCLE jsearch              !TWIN side found, exit search
              END IF
            END DO
          END DO
          ne = ne + 1         !increase number of boundary sides
          lside(in,ie) = ne   !assign number to boundary side
        END DO jsearch
      END DO

      ! SORT sides according to connectivities
      ALLOCATE ( order(3,ne) )  !(1)


      order = 0   !initializes
      DO ie=1,nelem              !for each element
        DO in=1,nnode                !for each side in the element
          IF( lside(in,ie) < 0 )CYCLE  !if an inner side then cycle
          n1=lside(in,ie)          !segment order
          n2=lnods(in,ie)          !second node of the side
          !   search for the 'next' side of n1-n2
          nsearch : DO je=1,nelem  !for each remaining element (including IE)
            DO jn=1,nnode                 !for each side
              IF( lside(jn,je) <  0 )CYCLE   !if an inner side then cycle
              m1=lnods(MOD(jn,nnode)+1,je)   !first node of the side
              IF( n2 == m1) THEN  !if this is the next side
                m2=lside(jn,je)            !segment order
                order(1,n1) = in           !node   (this side)(local order)
                order(2,n1) = ie           !element(this side)
                order(3,n1) = m2           !node   (next side in order)
                EXIT nsearch               !next segment found, exit search
              END IF
            END DO
          END DO nsearch
        END DO
      END DO

      !GENERATE SURFACE data_base
      !Initialize empty list
      ALLOCATE( bound(2,ne) )
      in = 1
      DO ie=1,ne     !for each side
        IF(order(1,in) == 0)THEN  !find a non-stored side
          in = 1         !initializes
          DO
            IF(order(1,in) /= 0)EXIT  !non-stored side
            in = in+1                 !next side
          END DO
        END IF
        jn = order(1,in)     !side order in element connectivities
        je = order(2,in)     !element
        ! store segment in list
        bound(2,ie) = lnods(jn,je)                !first node
        bound(1,ie) = lnods(MOD(jn,nnode)+1,je)   !second node
        order(1,in) = 0      !set as stored or non-existent
        in = order(3,in)     !next side
      END DO
      DEALLOCATE (lside, order)               !release memory
      RETURN
      END SUBROUTINE surf03
