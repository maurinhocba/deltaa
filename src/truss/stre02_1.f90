 SUBROUTINE stre02_1 (mat,e,stres,hyper,newm,j)

 USE mat_dba, ONLY : mater,inte_cr

 ! computes the internal nodal forces  1D (truss elements)

 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: e               !strain
 REAL (kind=8), INTENT(IN OUT) :: stres(:)    !(1:7)stre0,stres,?,ep,p,eps,dl,(ep,p,eps)(LC)
 REAL (kind=8), INTENT(OUT) :: j              !jacobian
 TYPE (mater), POINTER :: mat
 LOGICAL, INTENT(IN) :: hyper
 LOGICAL, INTENT(IN OUT) :: newm
 !  local variables
 REAL (kind=8) stran,dsttr,dlamb,yistr,signo,efpst
 INTEGER (kind=4), SAVE :: is,ik,i
 REAL (kind=8), SAVE :: young,c0,c1,c2,c3,kinmod,kinsat,&
                        aprim,kinet,poiss
 LOGICAL, SAVE :: plast

 !     evaluates incremental and total stress

 IF( newm )THEN
   newm = .FALSE.
   young = mat%prope(1)
   poiss = mat%prope(2)
   poiss = 1d0-2d0*poiss
   plast =  mat%matdef(3) > 1
   IF( plast )THEN
     is    =  mat%matdef(4)
     c0   = mat%propp(1)        !C0 constant or Initial Yield
     c1   = mat%propp(2)        !Efref or Hardening constant
     c2   = mat%propp(3)        !exponent
     c3   = mat%propp(4)        !residual flow stress

     ik    =  mat%matdef(5)
     SELECT CASE (ik)
     CASE (1)
       kinet = 0d0
     CASE (2)
       kinmod = mat%propp(6)
       kinet = kinmod
     CASE (3)
       kinmod = mat%propp(6)
       kinsat = mat%propp(7)
     END SELECT
   END IF
 END IF

   ! implicit analysis (last converged values)
 IF ( plast ) THEN
   stres(4:6) = stres(8:10)     !plastic strain, back stress, eq.pl.st
   stres(7) = 0d0               !dlamb
 END IF

 IF( hyper )THEN  !for hyper-elastic models
   IF( plast )THEN
     stran = e - stres(4)                         !trial elastic strain
     stres(2) = young*stran + stres(1)
   ELSE
     stres(2) = young*e + stres(1)
   END IF
 ELSE  !for hipo-elastic models
   stres(2) = young*e + stres(2)
 END IF

 IF ( plast ) THEN        ! plastic

   !  compute present yield stress
   efpst = stres(6)

   IF( is == 5 )THEN
     i = 1   !begin at first interval
     c0 = inte_cr (mat%chead%val,mat%chead%np,efpst,i)
     c1 = mat%chead%val(3,i)
     c0 = c0 - c1 * efpst
     CALL isoha14(2,yistr,aprim,efpst,c0,c1,c2,c3)
   ELSE
     CALL isoha14(is,yistr,aprim,efpst,c0,c1,c2,c3)
   END IF

   signo = 1d0                                  !strain sign
   IF(stran < 0d0) signo = -1d0
   SELECT CASE ( ik ) !kinematic hardening model
   CASE (1) ! no hardening
     dsttr = ABS(stres(2)) - yistr       !yield function
   CASE (2) ! linear kinematic hardening
     dsttr = ABS(stres(2)-stres(5)) - yistr       !yield function
   CASE (3) ! saturation law for kinematic hardening
     dsttr = ABS(stres(2)-stres(5)) - yistr       !yield function
     kinet = kinmod       !??
   END SELECT

   IF(dsttr > 0) THEN                          !if greater than zero
     dlamb = dsttr / (young+aprim+kinet)     !consistency param.
     stres(4) = stres(4) + dlamb*signo              !plastic strain
     stres(5) = stres(5) + dlamb*kinet*signo        !back stress
     stres(6) = stres(6) + dlamb                    !effective plastic strain
     stres(2) = stres(2) - young*dlamb*signo        !stress
     stres(7) = dlamb                               !implicit analysis
   END IF
   stran = stres(2)/young
 ELSE
   stran = e
 END IF
 j = 1d0 + stran*poiss
 RETURN
 END SUBROUTINE stre02_1
