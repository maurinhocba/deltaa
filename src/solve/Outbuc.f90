 SUBROUTINE outbuc(npoin,ndofn,istep,ifpre,x,lambd,buckl,flag)

 !     print eigevector into output file

 USE param_db, ONLY: mich
 USE ctrl_db, ONLY : ndime
 USE outp_db, ONLY : iwrit
 USE kinc_db, ONLY : neq,npsdf,nesdf,ftsdf,nn
 USE lispa0
 USE npo_db, ONLY : label,cpx
 USE c_input,ONLY : openfi
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN)::npoin,ndofn,istep,ifpre(:,:),flag
 REAL (kind=8),INTENT(IN) :: lambd,buckl
 REAL (kind=8),INTENT(IN OUT) :: x(:)
 CHARACTER(len=mich):: inttoch

 INTEGER (kind=4) :: i,j,k,l,m,n
 REAL (kind=8)    :: y(ndofn),pi2,omega,period


 SELECT CASE (flag)
 CASE (:0) !linear buckling analysis
   IF( iwrit == 1 ) &
     WRITE(lures,"(//' linear buckling analysis from a non-linear path '/     &
        &   ' step ',i3,'  actual load', e15.6,'  critical load ',e15.6,//, &
        &   ' eigenvector components ',//,'   node       ',                 &
        &   'DOF-1       DOF-2       DOF-3  ')") istep,lambd,buckl
   WRITE(10,"('GiD Post Results File 1.0',/)")

   WRITE(10,"(//'Result ""Buck_mode_',a2,'"" ""Load Anal.""   0.1000E+01 Vector OnNodes ',/,  &
        &       'ComponentNames ""X_Disp"",""Y_Disp"",""Z_Disp"" ',/,                      &
        &       'Values')")TRIM(inttoch(-flag,2))
 CASE (1)    !first derivative of asyntotic analysis
   WRITE(lures,"(/,' First derivatives of displacements ',/           &
        &      '   node     DOF-1       DOF-2       DOF-3  ')")
 CASE (2)    !second derivative of asyntotic analysis
   WRITE(lures,"(/,' Second derivatives of displacements ',           &
        &      '   node     DOF-1       DOF-2       DOF-3  ')")
 CASE (3)    !natural mode of dynamic analysis
   pi2 = ATAN(1d0)*8d0      !2*Pi
   omega = SQRT(buckl)      !natural frequency
   period = pi2/omega       !period
   WRITE(lures,"(//' linear dynamic eigen-mode analysis  '/     &
        &   'Order ',i3,'  natural freq ',e15.6,' Period ', e15.6)") istep,omega,period
   IF( iwrit == 1 ) &
     WRITE(lures,"(// ' eigenvector components ',//,'   node       ',                 &
        &   'DOF-1       DOF-2       DOF-3  ')")
   WRITE(10,"(//'Result ""Eigen_mode"" ""Load Anal.""',e15.5,' Vector OnNodes ',/,  &
        &       'ComponentNames ""X_Disp"",""Y_Disp"",""Z_Disp"" ',/,               &
        &       'Values')")REAL(istep)
 END SELECT

 DO i =1,npoin       !for each node
   DO j=1,ndofn        !for each DOF
     k = ifpre(j,i)      !associated global DOF
     SELECT CASE (k)        !according to value
     CASE (1:)                !active DOF
       y(j) = x(k)
     CASE (-nn+1:-1)          !slave DOF
       l = npsdf(-k)
       m = npsdf(-k+1) - 1
       y(j) = 0d0
       DO n=l,m
         IF( nesdf(n) > 0 ) &
           y(j) = y(j)+ftsdf(n)*x(nesdf(n))
       END DO
     CASE (:-nn,0)            !prescribed DOF or inactive
       y(j) = 0d0
     END SELECT
   END DO
   IF( ASSOCIATED (cpx) )THEN
     IF( cpx(1,i) /= 0 )THEN
       DO k=1,ndime
         m = ifpre(k,cpx(2,i))
         n = ifpre(k,cpx(3,i))
         y(k) = y(k)/2d0
         IF( m > 0 ) y(k) = y(k) + x(m) /4d0
         IF( n > 0 ) y(k) = y(k) + x(n) /4d0
       END DO
     END IF
   END IF
   IF(ANY(ifpre(1:ndime,i) /= 0))THEN
     IF( iwrit == 1 )WRITE(lures,"(i7,3x,8E15.5)")label(i),y
     IF( flag <= 0 .OR. flag == 3)WRITE(10,"(i7,3x,8E15.5)")label(i),y(1:ndime)
   END IF
 END DO
 IF( flag <= 0 .OR. flag == 3 )THEN
  WRITE(10,"('End Values')")
  !CLOSE(10)
 END IF
 RETURN

 END SUBROUTINE outbuc
