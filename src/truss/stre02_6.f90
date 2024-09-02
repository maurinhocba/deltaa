 SUBROUTINE stre02_6 (mat,l,stres,newm)

 ! computes the internal nodal forces  1D (truss elements)
 ! for rubbers (material type 6)
 USE mat_dba, ONLY : mater
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: l  !lambda
 REAL (kind=8), INTENT(IN OUT) :: stres(:)
 TYPE (mater), POINTER :: mat
 LOGICAL, INTENT(IN OUT) :: newm
 !  local variables
 REAL (kind=8) le,beta,db !dl,yistr,signo,efpst
 REAL (kind=8), SAVE :: pr(12)  !,uniax,hards,k
 LOGICAL, SAVE :: plast

 !     evaluates incremental and total stress

 IF( newm )THEN
   pr    = mat%prope(7:18)
   plast = .FALSE.
   !plast =  mat%matdef(3) > 1
   !IF( plast )THEN
   !  uniax = mat%propp(1)
   !  hards = mat%propp(2)
   !END IF
   newm = .FALSE.
 END IF

 !IF ( plast ) THEN          ! plastic
 !  ! implicit analysis (last converged values)
 !  stres(4:6) = stres(8:10)     !plastic strain, back stress, eq.pl.st
 !  stres(7) = 0d0               !dl
 !  le = l/exp(stres(4))         !elastic strain
 !  efpst = stres(6)             !plastic strain
 !  yistr = uniax + hards*efpst  !linear hardening
 !ELSE
   le = l
 !END IF

 CALL rubber1d(pr(1),le,plast,mat%matdef(8),beta,db)

 !IF( plast )THEN
 !  dl = ABS(beta) - yistr       !yield function
 !  IF( dl > 0d0 )THEN
 !    signo = 1d0                                  !strain sign
 !    IF(beta < 0d0) signo = -1d0
     !DO
     !  dl = dl / (db+hards)     !consistency param.
     !  le = le/EXP(dl*signo)      !elastic lambda
     !  CALL  rubber1d(pr(1),le,.TRUE.,mat%matdef(8),beta,db)
     !  efpst = efpst + dl
     !  yistr = uniax + hards*efpst    !linear hardening
     !  dl = ABS(beta) - yistr        !yield function
     !  IF( ABS(dl) < 1d-4 )EXIT
     !END DO
     !stres(4) = LOG(l/le)                        !plastic strain
     !stres(5) = stres(5) + kinet*stres(4)        !back stress
     !stres(6) = ABS(stres(4))                    !effective plastic strain
     !stres(7) = dl                               !implicit analysis
 !  END IF
 !END IF

 stres(2) = beta        !converged stress
 !IF( k > 0d0 )THEN
 !  j = EXP(beta/2d0/k)    !volumetric
 !ELSE
 !  j = 1d0
 !END IF

 RETURN
 END SUBROUTINE stre02_6
