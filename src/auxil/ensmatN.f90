SUBROUTINE ensmatN(nvarl,lm,locst,gl,gu)
!*************************************************************************
!
!     assembles a local non-symmetric matrix into the global symmetric matrices
!     local matrix is stored as a square matrix
!     both global matrices (triangles) are stored as arrays
!     THIS ROUTINE MUST BE CHECHED FOR SLAVE/MASTER DOFS
!
!*************************************************************************
USE kinc_db, ONLY : maxa,maxav ,nn,npsdf,nesdf,ftsdf,ftsd0
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl)
REAL (kind=8),INTENT(IN) :: locst(nvarl,nvarl)
REAL (kind=8),INTENT(IN OUT) :: gl(1:*),gu(1:*)


INTEGER (kind=4) :: i,j,k,neqi,neqj,l,m,ib,ie,jb,je,neql,neqm
REAL    (kind=8) :: stk,fl,fr,stkm


DO i = 1,nvarl                                   !for each row
  neqi = lm(i)                                   !assoc. equation
  DO j = 1,nvarl                               !for each column
    neqj = lm(j)                               !assoc. equation
    stk = locst(i,j)                           !value to assemble
    IF(neqi > 0) THEN                            !if active DOF
      IF(neqj > 0) THEN                          !   BOTH Active DOFs
        IF(neqj <= neqi) THEN                    !if neq(j) <= neq(i) (lower triangle including Diag)
          k = maxav(neqi) + neqi - neqj          !post. in lower global matrix
          gl(k)=gl(k)+stk                        !sums on lower global matrix
        ELSE                                     !if neq(j) >  neq(i) (upper triangle)
          k = maxav(neqj) + neqj - neqi          !post. in upper global matrix
          gu(k)=gu(k)+stk                        !sums on upper global matrix
        END IF
      ELSE IF(neqj < 0 .AND. neqj > -nn) THEN    !if J is a slave DOF              !
        jb = npsdf(-neqj)                        !first position in array          !
        je = npsdf(-neqj+1)-1                    !last position in array           !
        DO l=jb,je                               !for each master DOF              !
          neql=nesdf(l)                          !assoc. equation                  !
          IF(neql > 0) THEN                      !if DOF active                    !
            stkm = stk*ftsdf(l)                                                    !
            IF(neql <= neqi) THEN                !if neq(l) <= neq(i)              !
              k = maxav(neqi) + neqi - neql      !post. in lower global matrix     !
              gl(k)=gl(k)+stkm                   !sums on lower global matrix      !
            ELSE                                 !if neq(l) >  neq(i)              !
              k = maxav(neql) + neql - neqi      !post. in upper global matrix     !
              gu(k)=gu(k)+stkm                   !sums on upper global matrix      !
            END IF                                                                 !
          END IF                                                                   !
        END DO                                   !l=jb,je                          !
      END IF
    ELSE IF(neqi < 0 .AND. neqi > -nn) THEN        !if slave DOF                     !
      ib = npsdf(-neqi)                            !first position in array          !
      ie = npsdf(-neqi+1)-1                        !last position in array           !
      DO m=ib,ie                                                                   !
        fl = ftsd0(m)                            !assoc. factor                    !
        neqm=nesdf(m)                            !assoc. equation                  !
        IF(neqm > 0) THEN                        !if active DOF                    !
          stkm = stk*fl                          !modified value                   !
          IF(neqj > 0) THEN                      !if active DOF                    !
            IF(neqj <= neqm) THEN                !if neq(j) <= neq(m)              !
              k = maxav(neqm) + neqm - neqj      !post. in lower global matrix     !
              gl(k)=gl(k)+stkm                   !sums on lower global matrix      !
            ELSE                                 !if neq(j) >  neq(m)              !
              k = maxav(neqj) + neqj - neqm      !post. in upper global matrix     !
              gu(k)=gu(k)+stkm                   !sums on upper global matrix      !
            END IF                                                                 !
          ELSE IF(neqj < 0 .AND. neqj > -nn) THEN!if slave DOF too                 !
            jb = npsdf(-neqj)                    !first equation in array          !
            je = npsdf(-neqj+1)-1                !last equation in array               !
            DO l=jb,je                           !for each master DOF              !
              fr = ftsdf(l)                      !assoc. factor                    !
              neql=nesdf(l)                      !assoc. equation                  !
              IF(neql > 0) THEN                  !if DOF inactive                  !
                stkm = stk*fr*fl                 !modified value                   !
                IF(neql <= neqm) THEN            !if neq(l) <= neq(m)              !
                  k= maxav(neqm)+neqm-neql       !post. in lower global matrix     !
                  gl(k)=gl(k)+stkm               !sums on lower global matrix      !
                ELSE                             !if neq(l) >  neq(m)              !
                  k= maxav(neql)+neql-neqm       !post. in upper global matrix     !
                  gu(k)=gu(k)+stkm               !sums on upper global matrix      !
                END IF                                                             !
              END IF                                                               !
            END DO                               !l=jb,je                          !
          END IF                                 !neqj ?                           !
        END IF                                   !neqm > 0                         !
      END DO                                     !m=ib,ie                          !
    END IF
  END DO                                         !j=1,nvarl                        !
END DO                                           !i=1,nvarl

RETURN
END SUBROUTINE ensmatN
