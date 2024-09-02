 SUBROUTINE colsol(l,maxa,neq,task,iw,ish,nsymm,v,u)
 !.....................................................................
 !.                                                                   .
 !.   P R O G R A M                                                   .
 !.     TO SOLVE EQUILIBRIUM EQUATIONS FROM FINITE ELEMENTS           .
 !.     IN FAST MEMORY USING COMPACT STORAGE WITH SKYLINE SCHEME      .
 !.                                                                   .
 !.  --INPUT VARIABLE                                                 .
 !.     L(:)     = STIFFNES MATRIX (SYMMETRIC PART) IN COMPACT FORM   .
 !.     V(NEQ)   = LOADING VECTOR                                     .
 !.     MAXA(NEQ+1)= VECTOR CONTAINING THE DIAGONAL POSITIONS IN A    .
 !.     NEQ      = NUMBER OF EQUATION                                 .
 !.     TASK     = FLAG                                               .
 !.          EQ. 1 MATRIX FACTORIZATION                               .
 !.          EQ. 2 BACK SUBSTITUTION OF VECTOR V                      .
 !.     IW       = NUMBER ASSOCIATED TO OUTPUT DEVICE FOR MESSAGES    .
 !.     ISH      = FLAG TO INDICATE IF NON-POSSITIVE MATRICES ALLOWED .
 !.          EQ. 0   ONLY POSSITIVE MATRICES ALLOWED                  .
 !.          EQ. 1   ANY NON-SINGULAR MATRIZ ALLOWED                  .
 !.     NSYMM    = INDICATES IF MATRIZ IS SYMMETRIC                   .
 !.          EQ. 0   SYMMETRIC MATRIX                                 .
 !.          EQ. 1   NON-SYMMETRIC MATRIX                             .
 !.     IN THE LAST CASE                                              .
 !.        U(:)  = NON-SYMMETRIC PART OF THE STIFFNESS MATRIX         .
 !.                                                                   .
 !.-- OUTPUT                                                          .
 !.     TASK = 1                                                      .
 !.        L(:)  = D & L - FACTOR OF THE STIFFNES MATRIX              .
 !.        U(:)  = UPPER TRIANGULAR FACTOR FOR NON-SYMMETRIC MATRIX   .
 !.     TASK = 2                                                      .
 !.        V(NEQ) = EQUATION SOLUTION                                 .
 !.                                                                   .
 !.....................................................................
 IMPLICIT NONE
 !.....................................................................
 INTEGER (kind=4),INTENT(IN) :: neq,task,iw,ish,maxa(:),nsymm
 REAL (kind=8),INTENT(IN OUT) :: l(:)
 REAL (kind=8),INTENT(IN OUT), OPTIONAL :: u(:),v(:)
 !Local Variables
 INTEGER (kind=4) :: n,kn,kl,ku,k,ic,klt,ki,nd,m,kh,j,kk
 REAL (kind=8) :: b,c,d !,e(neq,neq)
 !
 IF(task == 1) THEN

   ! L*D*U Factorization of Stiffness Matrix

   IF(nsymm == 0 .OR. nsymm == -1) THEN   ! U = L(T)
      !e = 0d0 !initializes
      !DO n=1,neq                           !para cada fila/columna i==N
      !  kn = maxa(n)                       !pos. diagonal (i,i) en L y U
      !  kl = kn + 1                        !pos. de K(i,i-1) y U(i-1,i), primer punto del bucle
      !  ku = maxa(n+1) - 1                 !pos.de K(Jmin,i) y U(i,Jmin), ultimo punto del bucle
      !  j=n                                !initializes to diagonal value
      !  e(j,j) = l(kn)                     !u(kn) diagonal value
      !  DO kn=kl,ku                        !for each point in the column
      !    j=j-1                            !previos value
      !    e(j,n) = l(kn)                   !along row
      !    e(n,j) = l(kn)                   !along column
      !  END DO
      !END DO
      !DO n=1,neq                           !para cada fila/columna i==N
      !  !WRITE(58,"(24E15.6)")e(n,:)
      !  WRITE(58,"(35(E19.12,','),E19.12,';')")e(n,:)
      !END DO
     DO n=1,neq                           !para cada fila i==N
       kn = maxa(n)                       !pos. de K(i,i)
       kl = kn+1                          !pos. de K(i,i-1)
       ku = maxa(n+1)-1                   !pos.del 1ro de la fila K(i,Jmin)
       kh = ku-kl                         !long de la fila sin extremos
       IF(kh > 0) THEN                    !si hay val.de g(i,j) a calcular
         k = n - kh                       !ec. asoc. al 1ro de la fila Jmin
         ic = 0                           !contador
         klt = ku                         !pos del el (i,j) = (i,Jmin)
         DO j=1,kh                        !calculo de los g(i,j)
           ic = ic + 1                    !num de terminos de la sumatoria
           klt = klt - 1                  !pos de g(i,j)
           ki = maxa(k)                   !pos de d(j,j)
           nd = maxa(k+1) - ki - 1        !long. fila j - diagonal
           IF(nd > 0) THEN                !num de terminos de la sumatoria
             kk = MIN(ic,nd)              !numero definitivo de terminos
             c = 0.d0                     !pone a 0 el acumulador
             DO m=1,kk                    ! KI+M==(j,j-m)  KLT+M==(i,j-m)
               c = c + l(ki+m)*l(klt+m)   !+l(j,k)*g(i,k)
             END DO
             l(klt) = l(klt) - c          !g(i,j) --> A
           END IF
           k = k + 1                      !pasa a la siguiente col. k==j
         END DO
       END IF
       IF(kh >= 0) THEN                   !si hay el. fuera de la diagonal
         k = n                            !K==j   N==i  --> j=i
         b = 0.d0                         !pone a 0 el acumulador diag.
         DO kk=kl,ku                      !para c/ter. no diag.
           k = k - 1                      !decrementa j (hasta Jmin)
           ki = maxa(k)                   !pos. de d(j,j)
           c = l(kk)/l(ki)                !l(i,j)=g(i,j)/d(j,j)
           b = b + c*l(kk)                !+l(i,j)*g(i,j)
           l(kk) = c                      !l(j,i) --> A
         END DO
         l(kn) = l(kn) - b                !d(n,n) --> A
       END IF
       IF(l(kn) <= 0) THEN
         IF(ish == 0) THEN
           WRITE(iw,8000)
           WRITE(iw,8020)n,l(kn)
           STOP
         ELSE
           WRITE(iw,8020) n,l(kn)
           IF(l(kn) == 0) l(kn) = 1.d-16
         END IF
       END IF
     END DO
   ELSE                                   ! U /= L(T) Non-symmetric matrix
     !! write the whole matrix
     !e = 0d0 !initializes
     !DO n=1,neq                           !para cada fila/columna i==N
     !  kn = maxa(n)                       !pos. diagonal (i,i) en L y U
     !  kl = kn + 1                        !pos. de K(i,i-1) y U(i-1,i), primer punto del bucle
     !  ku = maxa(n+1) - 1                 !pos.de K(Jmin,i) y U(i,Jmin), ultimo punto del bucle
     !  j=n                                !initializes to diagonal value
     !  e(j,j) = l(kn)                     !u(kn) diagonal value
     !  DO kn=kl,ku                        !for each point in the column
     !    j=j-1                            !previos value
     !    e(j,n) = l(kn)                   !along row
     !    e(n,j) = u(kn)                   !along column
     !  END DO
     !END DO
     !DO n=1,neq                           !para cada fila/columna i==N
     !  WRITE(58,"(24E18.9)")e(n,:)
     !END DO
     DO n=1,neq                           !para cada fila/columna i==N
       kn = maxa(n)                       !pos. diagonal (i,i) en L y U
       kl = kn + 1                        !pos. de K(i,i-1) y U(i-1,i), primer punto del bucle
       ku = maxa(n+1) - 1                 !pos.de K(Jmin,i) y U(i,Jmin), ultimo punto del bucle
       kh = ku - kl                       !long. de la fila/columna sin extremos
       IF(kh > 0) THEN                    !si hay val.de g(i,j) a calcular ¿?
         k = n - kh                       !ec. asoc. al 1ro ¿o 2do? de la fila/columna Jmin
         ic = 0                           !contador
         klt = ku                         !pos del el (i,j) = (i,Jmin)
         DO j=1,kh                        !calculo de los g(i,j) y f(j,i)
           ic = ic + 1                    !num de terminos de la sumatoria
           klt = klt - 1                  !pos de g(i,j) y f(j.i)
           ki = maxa(k)                   !pos de d(j,j)
           nd = maxa(k+1) - ki - 1        !long. fila j - diagonal
           IF(nd > 0) THEN                !num de terminos de la sumatoria
             kk = MIN(ic,nd)              !numero definitivo de terminos
             c = 0.d0                     !pone a 0 el acumulador de g(i.j)
             d = 0.d0                     !pone a 0 el acumulador de f(j.i)
             DO m=1,kk                    ! KI+M==(j,j-m)  KLT+M==(i,j-m)
               c = c + u(ki+m)*l(klt+m)   !+u(k,j)*g(i,k)
               d = d + l(ki+m)*u(klt+m)   !+l(j,k)*f(k,i)
             END DO
             l(klt) = l(klt) - c          !g(i,j) --> A
             u(klt) = u(klt) - d          !f(j,i) --> U
           END IF
           k = k + 1                      !pasa a la siguiente col. k==j
         END DO
       END IF
       IF(kh >= 0) THEN                   !si hay el. fuera de la diagonal
         k = n                            !K==j   N==i  --> j=i
         b = 0.d0                         !pone a 0 el acumulador diag.
         DO kk=kl,ku                      !para c/ter. no diag.
           k = k - 1                      !decrementa j (hasta Jmin)
           ki = maxa(k)                   !pos. de d(j,j)
           c = l(kk)/l(ki)                !l(i,j)=g(i,j)/d(j,j)
           d = u(kk)/l(ki)                !u(j,i)=f(j,i)/d(j,j)
           b = b + c*u(kk)                !+l(i,j)*f(j,i)
           l(kk) = c                      !l(i,j) --> A
           u(kk) = d                      !u(j,i) --> U
         END DO
         l(kn) = l(kn) - b                !d(i,i) --> A
       END IF
       IF(l(kn) <= 0) THEN
         IF(ish == 0) THEN
           WRITE(iw,8000)
           WRITE(iw,8020) n,l(kn)
           STOP
         ELSE
           WRITE(iw,8020) n,l(kn)
           IF(l(kn) == 0) l(kn) = 1.d-16
         END IF
       END IF
     END DO
   END IF
   !       d = MINVAL(l(maxa(1:neq)))
   !       WRITE(IW,"(' minimum pivot =',e15.4)")d
 ELSE

   !   Reduce Independent vector

   DO n=1,neq
     kl = maxa(n)+1
     ku = maxa(n+1)-1
     IF(ku-kl >= 0) THEN
       k = n
       c = 0.d0
       DO kk=kl,ku
         k = k-1
         c = c + l(kk)*v(k)
       END DO
       v(n) = v(n) - c
     END IF
   END DO

   !       Back-Substitution

   DO n=1,neq
     k = maxa(n)
     v(n) = v(n)/l(k)
   END DO
   n = neq
   IF(nsymm == 0 .OR. nsymm == -1) THEN   ! Symmetric Matrix
     DO m=2,neq
       kl = maxa(n)+1
       ku = maxa(n+1)-1
       k  = n
       DO kk=kl,ku
         k = k-1
         v(k) = v(k) - l(kk)*v(n)
       END DO
       n = n-1
     END DO
   ELSE                                   !Non-symmetric matrix
     DO m=2,neq
       kl = maxa(n)+1
       ku = maxa(n+1)-1
       k = n
       DO kk=kl,ku
         k = k - 1
         v(k) = v(k) - u(kk)*v(n)
       END DO
       n = n-1
     END DO
   END IF

 END IF
 8000 FORMAT(//,'   STOP - Non-Positive Stiffness Matrix ',//)
 8020 FORMAT('   Non-positive pivot for equation',I5,'  PIVOT=',E20.12)

 RETURN
 END SUBROUTINE colsol
