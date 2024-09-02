      SUBROUTINE SSPACE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,BLO,&
     &BUPC,NN,NNM,NWK,NWM,NROOT,RTOL,NC,NNC,NITEM,IFSS,IFPR,NSTIF,IOUT,IFAC)
!***********************************************************************
!                                                                      *
!  P R O G R A M A                                                     *
!       PARA CALCULAR LOS AUTOVALORES MAS PEQUENOS Y LOS CORRESPONDIEN-*
!       TES AUTOVECTORES EN UN PROBLEMA GENERALIZADO USANDO EL METODO  *
!       DE ITERACION DE SUBESPACIOS                                    *
!                                                                      *
!  -  -  VARIABLES DE ENTRADA  -  -                                    *
!        A(NWK)    = MATRIZ DE RIGIDEZ EN FORMA COMPACTA ( ASUMIDA     *
!                    POSITIVA DEFINIDA )                               *
!        B(NWN)    = MATRIZ DE MASA EN FORMA COMPACTA                  *
!        MAXA(NNM) = VECTOR CONTENIENDO LAS DIRECCIONES DE LOS ELEMEN- *
!                    TOS DE LA DIAGONAL EN LA MATRIZ DE RIGIDEZ A      *
!        R(NN,NC)  = AUTOVECTORES EN LA SALIDA DEL PROBLEMA            *
!        EIGV(NC)  = AUTOVALORES EN LA SALIDA DEL PROBLEMA             *
!        TT(NN)    = VECTOR DE TRABAJO                                 *
!        W(NN)     = VECTOR DE TRABAJO                                 *
!        AR(NNC)   = MATRIZ DE TRABAJO QUE GUARDA LA PROJECCION DE K   *
!        BR(NNC)   = MATRIZ DE TRABAJO QUE GUARDA LA PROJECCION DE M   *
!        VEC(NC,NC)= VECTOR DE TRABAJO                                 *
!        D(NC)     = VECTOR DE TRABAJO                                 *
!        RTOLV(NC) = VECTOR DE TRABAJO                                 *
!        BUP(NC)   = VECTOR DE TRABAJO                                 *
!        BLO(NC)   = VECTOR DE TRABAJO                                 *
!        BUPC(NC)  = VECTOR DE TRABAJO                                 *
!        NN        = ORDEN DE LA MATRICES DE MASA Y RIGIDEZ            *
!        NNM       = NN + 1                                            *
!        NWK       = NUMERO DE ELEMENTOS BAJO EL SKYLINE DE LA         *
!                    MATRIZ DE RIGIDEZ                                 *
!        NWM       = NUMERO DE ELEMENTOS BAJO EL SKYLINE DE LA         *
!                    MATRIZ DE MASA                                    *
!                      I. E. NWM=NWK PARA MATRIZ DE MASA CONSISTENTE   *
!                            NWM=NN PARA MATRIZ DE MASA REDUCIDA       *
!        NROOT     = NUMERO DE AUTOVALORES Y AUTOVECTORES PEDIDOS      *
!        RTOL      = TOLERANCIA DE CONVERGENCIA EN LOS AUTOVALORES     *
!                    (1.E-06 OR SMALLER )                              *
!        NC        = NUMERO DE VECTORES DE ITERACION USADOS            *
!                    (GENERALMENTE FIJADO EN MIN(2*NROOT, NROOT+8),PERO*
!                    NC NO PUEDE SER MAYOR QUE EL NUMERO DE GRADOS DE  *
!                    LIBERTAD MASICOS)                                 *
!        NNC       = NC*(NC+1)/2 DIMENSION DE LOS ARREGLOS DE ALMACENA-*
!                    MIENTO AR, BR                                     *
!        NITEM     = NUMERO MAXIMO DE ITERACIONES DE SUBESPACIOS PERMI-*
!                    TIDOS (USUALMENTE FIJADO EN 16)                   *
!                    LOS PARAMETROS NC Y/O NITEM DEBEN SER AUMENTADOS  *
!                    SI UNA SOLUCION NO CONVERGE                       *
!        IFSS      = BANDERA PARA VERIFICACION DE LA SECUENCIA DE STURN*
!                    EQ. 0 NO  VERIFICAR                               *
!                    EQ. 1 VERIFICAR                                   *
!        IFPR      = BANDERA PARA IMPRIMIR DURANTE LA ITERACION        *
!                    EQ. 0 NO  IMPRIME                                 *
!                    EQ. 1 IMPRIME                                     *
!        NSTIF     = ARCHIVO SECUENCIAL PARA GUARDAR LA MATRIZ DE RIG. *
!        IOUT      = ARCHIVO DE IMPRESION DE SALIDA                    *
!        IFAC      = INDICA (1) QUE LA MATRIZ YA ESTA FACTORIZADA      *
!                                                                      *
!   -  -  SALIDA  - -                                                  *
!        EIGV(NROOT)= AUTOVALORES                                      *
!        R(NN,NROOT)=AUTOVECTORES                                      *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!     argumentos
      INTEGER(kind=4), INTENT(IN) :: NN,NNM,NWK,NWM,NROOT,NC,NNC,NITEM, &
     &                               IFPR,NSTIF,IOUT,IFAC
      REAL(kind=8), INTENT(IN OUT) ::  RTOL,                            &
     &          A(NWK),B(NWM),R(0:NN,NC),TT(NN),W(NN),EIGV(NC),         &
     &          D(NC),VEC(NC,NC),AR(NNC),BR(NNC),RTOLV(NC),BUP(NC),     &
     &          BLO(NC),BUPC(NC)
      INTEGER(kind=4), INTENT(IN) :: MAXA(NNM)
      INTEGER(kind=4), INTENT(IN OUT) ::IFSS
!     Variables locales
      INTEGER (kind=4) :: iconv,max,nsch,nsmax,n1,nc1,i,nd,ii,j,l,ij,   &
     &                    ish,nite,k,ind,itemp,is,nei,nmis
      REAL (kind=8) :: tolj,rt,art,brt,eigvt,bt,dif,vnorm,wnorm,shift
      INTERFACE
        INCLUDE 'colsol.h'
      END INTERFACE
!***********************************************************************
!                                                                      *
!        ESTABLECE TOLERANCIA PARA LA ITERACION DE JACOBI              *
      TOLJ=0.000000000001
!
!        INICIALIZACION
!
      max = maxa(nnm)         !tamaño de la matriz B
      ICONV=0                 !inicializa convergencia
      NSCH=0
      NSMAX=12
      n1 = nc+1               !No de vectores de iteracion + 1
      nc1= nc-1               !No de vectores de iteracion - 1
      IF( ifac == 0 )THEN     !Matriz sin factorizar
        REWIND NSTIF
        WRITE(NSTIF) A        !guarda en disco
      END IF
      d   =0d0                !inicializa el vector de trabajo
!
!     ESTABLECE LOS VECTORES DE ITERACION INICIALES
!
      nd=nn/nc            !relacion entre el No de GdL y el No de vectores de iteracion
      IF(nwm > nn) THEN   !Si la matriz B no es diagonal
        DO I=1,NN
          II=MAXA(I)
          R(I,1)= B(II)
          W(I)=B(II)/A(II)
        END DO
      ELSE                ! B diagonal
        J=0            !inicializa numero de valores diagonales no-nulos en B
        DO I =1 ,NN    !para cada ecuación
          II=MAXA(I)     !posicion diagonal en A
          R(I,1)=B(I)    !valor diagonal en B
          IF(B(I).GT.0) J=J+1 !aumenta el contador de valores no-nulos
          W(I)=B(I)/A(II)     !relación entre valores diagonales
        END DO
        IF(NC > J)THEN       !verifica que haya + pivotes positivos que autovalores a calcular
          WRITE(IOUT,1007)
          STOP                  ! se detiene si no hay suficientes pivotes no nulos
        END IF
      END IF
      DO  J=2,NC
        DO I=1,NN
          R(I,J)=0
        END DO
      END DO
!
      L=NN-ND
      DO J=2,NC
        RT=0.
        DO I=1,L
          IF(W(I).LT.RT)CYCLE
          RT=W(I)
          IJ=I
        END DO
        DO I=L,NN
          IF(W(I).LE.RT)CYCLE
          RT=W(I)
          IJ=I
        END DO
        TT(J)= FLOAT(IJ)
        W(IJ)=0.
        L=L-ND
        R(IJ,J)=1.
      END DO
!
!     WRITE(IOUT,1008)
!     WRITE(IOUT,1002)(TT(J),J=2,NC)
!
!     FACTORIZA LA MATRIZ A EN (L)*(D)*(L(T))
!
      ISH=0    !SOLO MATRIZ DE RIGIDEZ POSITIVA DEFINIDA
      !ISH=1    !PERMITE MATRIZ DE RIGIDEZ NO-DEFINIDA POSITIVA
      !CALL COLSOL(A,W,MAXA,NN,NWK,NNM,1,IOUT,ISH)
      IF(IFAC == 0) CALL colsol(a,maxa,nn,1,iout,ish,0)  !NUEVA VERSION DE COLSOL
!
!  --- C O M I E N Z A   E L   C I C L O   D E   I T E R A C I O N ----
!
      NITE=0
      DO
        NITE=NITE+1
        write(*,'(''+'',T5,'' ITER ='',i5,''   AUTOV ='',e20.14)')        &
     &            nite,EIGV(NROOT)
        IF(IFPR/=0)WRITE(IOUT,1010)NITE
!
!       CALCULA LAS PROYECCIONES DE A Y B
!
  90    IJ=0
        DO J=1,NC
          DO K=1,NN
            TT(K)=R(K,J)
          END DO
          CALL REDBAK(A,TT,MAXA,NN,max)
          DO I=J,NC
            ART =0.
            DO K=1,NN
              ART=ART+R(K,I)*TT(K)
            END DO
            IJ=IJ+1
            AR(IJ)=ART
          END DO
          DO K=1,NN
            R(K,J)=TT(K)
          END DO
        END DO
        IJ=0
        DO J=1,NC
          CALL MULT (TT,B,R(1,J),MAXA,NN,NWM)
          DO I=J,NC
            BRT=0.
            DO K=1,NN
              BRT = BRT +R(K,I)*TT(K)
            END DO
            IJ=IJ+1
            BR(IJ)=BRT
          END DO
          IF(ICONV.GT.0) EXIT
          DO K=1,NN
            R(K,J)=TT(K)
          END DO
        END DO
!
!       RESUELVE EL AUTOSISTEMA DE OPERADORES DEL SUBESPACIO
!
        DO IND=1,2
          IF(IFPR/=0)THEN
            WRITE(IOUT,1020)
            II=1
            DO I=1,NC
              ITEMP =II+NC-I
              WRITE(IOUT,1005)(AR(J),J=II,ITEMP)
              II=II+N1-I
            END DO
            WRITE(IOUT,1030)
            II=1
            DO I=1,NC
              ITEMP= II+NC-I
              WRITE(IOUT,1005)(BR(J),J=II,ITEMP)
              II=II+N1-I
            END DO
          END IF
          IF(IND.EQ.2)EXIT

          CALL JACOBI (AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX,IFPR,IOUT)

          WRITE(IOUT,1040)

        END DO
!
!   ACOMODA LOS AUTOVALORES EN ORDEN ASCENDENTE
!
 350    IS=0
        II=1
        DO 360 I=1,NC1
          ITEMP=II+N1-I
          IF(EIGV(I+1).GE.EIGV(I))GO TO 360
          IS=IS+1
          EIGVT=EIGV(I+1)
          EIGV(I+1)=EIGV(I)
          EIGV(I)=EIGVT
          BT=BR(ITEMP)
          BR(ITEMP)=BR(II)
          BR(II)=BT
          DO K=1,NC
            RT=VEC(K,I+1)
            VEC(K,I+1)=VEC(K,I)
            VEC(K,I)=RT
          END DO
 360      II=ITEMP
        IF(IS.GT.0)GO TO 350
        IF(IFPR /= 0)THEN
          WRITE(IOUT,1035)
          WRITE(IOUT,1006)(EIGV(I),I=1,NC)
        END IF
!
!   CALCULA B VECES LOS AUTOVECTORES APROXIMADOS (ICONV.EQ.0)
!   O APROXIMACIONES A LOS AUTOVECTORES FINALES (ICONV.GT.0)
!
        DO 420 I=1,NN
        DO 422 J=1,NC
  422   TT(J)=R(I,J)
        DO 424 K=1,NC
        RT=0.
        DO 430 L=1,NC
  430   RT=RT+TT(L)*VEC(L,K)
  424   R(I,K)= RT
  420   CONTINUE
        IF(ICONV.GT.0) GO TO 500
!
!       VERIFICA LA CONVERGENCIA DE LOS AUTOVALORES
!
        DO 380 I=1,NC
        DIF=DABS(EIGV(I)-D(I))
  380   RTOLV(I)=DIF/EIGV(I)
        IF(IFPR.EQ.0)GO TO 385
        WRITE(IOUT,1050)
        WRITE(IOUT,1005)(RTOLV(I),I=1,NC)
!
 385    DO 390 I=1,NROOT
        IF(RTOLV(I).GT.RTOL)GO TO 400
  390   CONTINUE
        WRITE(IOUT,1060) RTOL
        ICONV=1
        CYCLE
  400   IF(NITE.LT.NITEM)GO TO 410
        WRITE(IOUT,1070)
        ICONV=2
        IFSS=0
        CYCLE
!
  410   DO 440 I=1,NC
  440   D(I)=EIGV(I)
      END DO
!
!  - - -   F I N   D E L   C I C L O   D  E   I T E R A C I O N - - -
!
  500 continue
!      WRITE(IOUT,1100)
!     WRITE(IOUT,1006)(EIGV(I),I=1,NROOT)
!     WRITE(IOUT,1110)
!     DO 530 J=1,NROOT
!530  WRITE(IOUT,1005) (R(K,J),K=1,NN)
!
!    CALCULA E IMPRIME LA NORMA DE LOS ERRORES
!
      REWIND NSTIF
      READ(NSTIF) A
!
      DO 580 L=1,NROOT
      RT=EIGV(L)
      CALL MULT(TT,A,R(1,L),MAXA,NN,NWK)
      VNORM=0
      DO 590 I=1,NN
  590 VNORM=VNORM+TT(I)*TT(I)
      CALL MULT(W,B,R(1,L),MAXA,NN,NWM)
      WNORM=0.
      DO 600 I=1,NN
      TT(I)=TT(I)-RT*W(I)
  600 WNORM=WNORM+TT(I)*TT(I)
      VNORM=DSQRT(VNORM)
      WNORM=DSQRT(WNORM)
      D(L)=WNORM/VNORM
  580 CONTINUE
      WRITE(IOUT,1115)
      WRITE(IOUT,1006)(D(I),I=1,NROOT)
!
!     APLICA LA VERIFICACION DE LA SECUENCIA DE STURN
!
      IF(IFSS.EQ.0)GO TO 700
      CALL SCHECK(EIGV,RTOLV,BUP,BLO,BUPC,D,NC,NEI,RTOL,SHIFT)
!
      WRITE(IOUT,1120)SHIFT
!
!     DESPLAZA LA MATRIZ A
!
      REWIND NSTIF
      READ(NSTIF)A
      IF(NWM.GT.NN)GO TO 645
      DO 640 I=1,NN
      II=MAXA(I)
  640 A(II)=A(II)-B(I)*SHIFT
      GO TO 660
 645  DO 650 I=1,NWK
 650  A(I)=A(I)-B(I)*SHIFT
!
!     FACTORIZA LA MATRIZ DESPLAZADA
!
  660 ISH=1
      !CALL COLSOL(A,W,MAXA,NN,NWK,NNM,1,IOUT,ISH)
      CALL colsol(a,maxa,nn,1,iout,ish,0)
!
!     CUENTA EL NUMERO DE ELEMENTOS NEGATIVOS EN LA DIAGONAL
!
      NSCH=0
      DO 664 I=1,NN
      II=MAXA(I)
      IF(A(II).LT.0)NSCH=NSCH+1
  664 CONTINUE
      IF(NSCH.EQ.NEI)GO TO 670
      NMIS=NSCH-NEI
      WRITE(IOUT,1130)NMIS
      GO TO 700
  670 WRITE(IOUT,1140)NSCH
  700 RETURN
!
 1002 FORMAT(1H0,10F10.0)
 1005 FORMAT(1H ,11E11.4)
 1006 FORMAT(1H0,6E21.13)
 1007 FORMAT(///'  STOP ,  NC ES MAYOR QUE EL NUMERO DE FILAS INDEPEN', &
     &'DIENTES DE LA MATRIZ  B ')
 1008 FORMAT(///,'  GRADOS DE LIBERTAD EXCITADOS POR VECTORES UNITARIO',&
     &'S INICIALES ')
 1010 FORMAT('1',' N U M E R O   D E   I T E R A C I O N ' , I4)
 1020 FORMAT(' PROYECCION DE A (MATRIZ AR)')
 1030 FORMAT(' PROYECCION DE B (MATRIZ BR)')
 1035 FORMAT(' AUTOVALORES DE AR-LAMBDA *BR')
 1040 FORMAT('0 AR Y BR LUEGO DE LA DIAGONALIZACION')
 1050 FORMAT('0 TOLERANCIA RELATIVA ALCANZADA EN LOS AUTOVALORES')
 1060 FORMAT(///,'CONVERGENCIA ALCANZADA PARA RTOL  ',E10.4)
 1070 FORMAT('1','***NO CONVERGE EN EL NUMERO MAXIMO DE ITERACIONES ',  &
     &'PERMITIDO'/,' ACEPTAMOS LOS VALORES DE ITERACION ACTUALES'/,     &
     &' NO SE REALIZA LA VERIFICACION DE LA SECUENCIA DE STURN ')
 1100 FORMAT(///' LOS AUTOVALORES CALCULADOS SON')
 1115 FORMAT(//,' IMPRIME LAS NORMAS DE LOS ERRORES DE LOS AUTOVALORES')
 1110 FORMAT(//,'LOS AUTOVECTORES CALCULADOS SON'//)
 1120 FORMAT(///,' VERIF. APLICADA AL DESPLAZ.  ',E22.14)
 1130 FORMAT(//,' HAY ',I4,' AUTOVALORES FALTANTES ')
 1140 FORMAT(// ' ENCONTRAMOS LOS ', I4,'  PRIMEROS AUTOVALORES')
 1150 FORMAT(/6E12.4/6E12.4)
!
      END
      SUBROUTINE REDBAK (A,V,MAXA,NN,max)
!**********************************************************************
!                                                                     *
!     P R O G R A M A                                                 *
!         QUE REDUCE Y SUSTITUYE HACIA ATRAS LOS VECTORES DE ITERACION*
!                                                                     *
!**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXA(NN+1),V(NN),A(max)
!
      DO N=1,NN
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU < KL)CYCLE
        K=N
        C=0.
        DO KK=KL,KU
          K=K-1
          C=C+A(KK)*V(K)
        END DO
        V(N)=V(N)-C
      END DO
!
      DO N=1,NN
        K=MAXA(N)
        V(N)=V(N)/A(K)
      END DO
      DO N=NN,2,-1
        KL=MAXA(N)+1
        KU=MAXA(N+1)-1
        IF(KU < KL)CYCLE
        K=N
        DO KK=KL,KU
          K=K-1
          V(K)=V(K)-A(KK)*V(N)
        END DO
      END DO
!
      RETURN
      END
      SUBROUTINE MULT (TT,B,RR,MAXA,NN,NWM)
!**********************************************************************
!                                                                     *
!     P R O G R A M A                                                 *
!         QUE EVALUA EL PRODUCTO DE RR VECES B Y LO GUARDA EN TT      *
!                                                                     *
!**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TT(1),B(1),MAXA(1),RR(1)
!
      IF(NWM.GT.NN)GO TO 20
      DO 10 I=1,NN
   10 TT(I)=B(I)*RR(I)
      RETURN

   20 DO 40 I =1,NN
   40 TT(I)=0.
      DO 100 I=1,NN
      KL=MAXA(I)
      KU=MAXA(I+1)-1
      II=I+1
      CC=RR(I)
      DO 100 KK=KL,KU
      II=II-1
  100 TT(II)=TT(II)+B(KK)*CC
      IF(NN.EQ.1)RETURN
      DO 200 I=2,NN
      KL=MAXA(I)+1
      KU=MAXA(I+1)-1
      IF(KU-KL)200,210,210
  210 II=I
      AA=0.
      DO 220 KK=KL,KU
      II=II-1
  220 AA=AA+B(KK)*RR(II)
      TT(I)=TT(I)+AA
  200 CONTINUE
!
      RETURN
      END
      SUBROUTINE SCHECK(EIGV,RTOLV,BUP,BLO,BUPC,NEIV,NC,NEI,RTOL,SHIFT)
!**********************************************************************
!                                                                     *
!     P R O G R A M A                                                 *
!         QUE EVALUA EL DESPLAZAMIENTO PARA VERIF.DE LA SEC.DE STURN  *
!                                                                     *
!**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EIGV(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),NEIV(NC)
!
      FTOL=0.01
!
      DO 100 I=1,NC
      BUP(I)=EIGV(I)*(1.+FTOL)
  100 BLO(I)=EIGV(I)*(1.-FTOL)
      NROOT=0
      DO 120 I=1,NC
  120 IF(RTOLV(I).LT.RTOL)NROOT=NROOT+1
      IF(NROOT.GE.1)GO TO 200
      WRITE(IOUT,1010)
      WRITE(6,1010)
      STOP
!
!  ENCUENTRA LIMITES SUPERIORES................
!
  200 DO 240 I=1,NROOT
  240 NEIV(I)=1
      IF(NROOT.NE.1)GO TO 260
      BUPC(1)=BUP(1)
      LM=1
      L=1
      I=2
      GO TO 295
  260 L=1
      I=2
  270 IF(BUP(I-1).LE.BLO(I))GO TO 280
      NEIV(L)=NEIV(L)+1
      I=I+1
      IF(I.LE.NROOT)GO TO 270
  280 BUPC(L)=BUP(I-1)
      IF(I.GT.NROOT)GO TO 290
      L=L+1
      I=I+1
      IF(I.LE.NROOT)GO TO 270
      BUPC(L)=BUP(I-1)
  290 LM=L
      IF(NROOT.EQ.NC)GO TO 300
  295 IF(BUP(I-1).LE.BLO(I))GO TO 300
      IF(RTOLV(I).GT.RTOL)GO TO 300
      BUPC(L)=BUP(I)
      NEIV(L)=NEIV(L)+1
      NROOT=NROOT+1
      IF(NROOT.EQ.NC)GO TO 300
      I=I+1
      GO TO 295
!
!    ENCUENTRA EL DESPLAZAMIENTO
!
  300 WRITE(IOUT,1020)
      WRITE(IOUT,1005)(BUPC(I),I=1,LM)
      WRITE(IOUT,1030)
      WRITE(IOUT,1006) (NEIV(I),I=1,LM)
      LL=LM-1
      IF(LM.EQ.1)GO TO 310
  330 DO 320 I=1,LL
  320 NEIV(L)=NEIV(L)+NEIV(I)
      L=L-1
      LL=LL-1
      IF(L.NE.1)GO TO 330
  310 WRITE(IOUT,1040)
      WRITE(IOUT,1006)(NEIV(I),I=1,LM)
      L=0
      DO 340 I=1,LM
      L=L+1
      IF(NEIV(I).GE.NROOT)GO TO 350
  340 CONTINUE
  350 SHIFT=BUPC(L)
      NEI=NEIV(L)
!
      RETURN
!
 1005 FORMAT('0',6E22.14)
 1006 FORMAT('0',6I22)
 1010 FORMAT('0 **** ERROR SOLUCION DETENIDA EN SCHECK',/12X,           &
     &'NO DE AUTOVALORES ENCONTRADOS.',/1X)
 1020 FORMAT(///,'  LIMITES SUPERIORES DE LOS AUTOVALORES ')
 1030 FORMAT(' NO DE AUTOVALORES EN CADA INTERVALO')
 1040 FORMAT(' NO DE AUTOVALORES MENORES QUE LOS LIMITES SUPERIORES')
      END
      SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)
!***********************************************************************
!                                                                      *
!    P R O G R A M A                                                   *
!       PARA RESOLVER UN AUTOPROBLEMA GENERALIZADO USANDO LA           *
!       ITERACION DE JACOBI GENERALIZADA                               *
!***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)
!     ABS(X)=DABS(X)
!     SQRT(X)=DSQRT(X)
!
!INICIALIZA LOS AUTOVALORES Y  MATRICES DE AUTOVECTORES
!
      N1=N+1
      II=1
      DO 10 I=1,N
      IF(A(II).GT.0..AND.B(II).GT.0)GO TO 4
      WRITE(IOUT,2020)II,A(II),B(II)
      STOP
    4 D(I)=A(II)/B(II)
      EIGV(I)=D(I)
   10 II=II+N1-I
      DO 30 I=1,N
      DO 20 J=1,N
   20 X(I,J)=0.
   30 X(I,I)=1.
      IF(N.EQ.1)RETURN
!
!INCIALIZA EL CONTADOR DE BARRIDAS Y EMPIEZA LA ITERACION
!
      NSWEEP=0
      NR=N-1
   40 NSWEEP=NSWEEP+1
      IF(IFPR.EQ.1)WRITE(IOUT,2000)NSWEEP
!
!       VER. SI ESTE ELEMENTO NO DIAGONAL ES SUFICIENTEMENTE GRANDE
!     COMO PARA REQUERIR QUE SE LO ANULE
      EPS=(.01**NSWEEP)**2
      DO 210 J=1,NR
      JP1=J+1
      JM1=J-1
      LJK=JM1*N-JM1*J/2
      JJ=LJK+J
      DO 210 K=JP1,N
      KP1=K+1
      KM1=K-1
      JK=LJK+K
      KK=KM1*N-KM1*K/2+K
      EPTOLA=(A(JK)*A(JK))/(A(JJ)*A(KK))
      EPTOLB=(B(JK)*B(JK))/(B(JJ)*B(KK))
      IF((EPTOLA.LT.EPS).AND.(EPTOLB.LT.EPS))GO TO 210
!
!  SI SE REQUIERE ANULAR ,CALCULE LOS ELEM.DE LAS MAT.DE ROT. CA Y CG
!
      AKK=A(KK)*B(JK)-B(KK)*A(JK)
      AJJ=A(JJ)*B(JK)-B(JJ)*A(JK)
      AB =A(JJ)*B(KK)-A(KK)*B(JJ)
      CHECK=(AB*AB+4.*AKK*AJJ)/4.
      IF(CHECK)50,60,60
   50 WRITE(IOUT,2020)
      STOP
   60 SQCH=DSQRT(CHECK)
      D1=AB/2.+SQCH
      D2=AB/2.-SQCH
      DEN=D1
      IF(DABS(D2).GT.DABS(D1))DEN=D2
      IF(DEN)80,70,80
   70 CA=0.
      CG=-A(JK)/A(KK)
      GO TO 90
   80 CA=AKK/DEN
      CG=-AJJ/DEN
!
! REALIZA LA ROTACION GENERALIZADA PARA ANULAR EL PRESENTE ELEMENTO
!
   90 IF(N-2)100,190,100
  100 IF(JM1-1)130,110,110
  110 DO 120 I=1,JM1
      IM1=I-1
      IJ=IM1*N-IM1*I/2+J
      IK=IM1*N-IM1*I/2+K
      AJ=A(IJ)
      BJ=B(IJ)
      AK=A(IK)
      BK=B(IK)
      A(IJ)=AJ+CG*AK
      B(IJ)=BJ+CG*BK
      A(IK)=AK+CA*AJ
  120 B(IK)=BK+CA*BJ
  130 IF(KP1-N)140,140,160
  140 LJI=JM1*N-JM1*J/2
      LKI=KM1*N-KM1*K/2
      DO 150 I=KP1,N
      JI=LJI+I
      KI=LKI+I
      AJ=A(JI)
      BJ=B(JI)
      AK=A(KI)
      BK=B(KI)
      A(JI)=AJ+CG*AK
      B(JI)=BJ+CG*BK
      A(KI)=AK+CA*AJ
  150 B(KI)=BK+CA*BJ
  160 IF(JP1-KM1)170,170,190
  170 LJI=JM1*N-JM1*J/2
      DO 180 I=JP1,KM1
      JI=LJI+I
      IM1=I-1
      IK=IM1*N-IM1*I/2+K
      AJ=A(JI)
      BJ=B(JI)
      AK=A(IK)
      BK=B(IK)
      A(JI)=AJ+CG*AK
      B(JI)=BJ+CG*BK
      A(IK)=AK+CA*AJ
  180 B(IK)=BK+CA*BJ
  190 AK=A(KK)
      BK=B(KK)
      A(KK)=AK+2.*CA*A(JK)+CA*CA*A(JJ)
      B(KK)=BK+2.*CA*B(JK)+CA*CA*B(JJ)
      A(JJ)=A(JJ)+2.*CG*A(JK)+CG*CG*AK
      B(JJ)=B(JJ)+2.*CG*B(JK)+CG*CG*BK
      A(JK)=0.
      B(JK)=0.
!
!     ACTUALIZA LA MATRIZ DE AUTOVECTORES LUEGO DE CADA ROTACION
!
      DO 200 I=1,N
      XJ=X(I,J)
      XK=X(I,K)
      X(I,J)=XJ+CG*XK
  200 X(I,K)=XK+CA*XJ
  210 CONTINUE
!
!   ACTUALIZA LOS AUTOVECTORES LUEGO DE CADA SALIDA
!
      II=1
      DO 220 I=1,N
      IF(A(II).GT.0..AND.B(II).GT.0.) GO TO 215
      WRITE(IOUT,2020)II,A(II),B(II)
      STOP
  215 EIGV(I)=A(II)/B(II)
  220 II=II+N1-I
      IF(IFPR.EQ.0)GO TO 230
      WRITE(IOUT,2030)
      WRITE(IOUT,2010)(EIGV(I),I=1,N)
!
!    VERIFICA CONVERGENCIA
!
  230 DO 240 I=1,N
      TOL=RTOL*D(I)
      DIF=DABS(EIGV(I)-D(I))
      IF(DIF.GT.TOL)GO TO 280
  240 CONTINUE
!
!  VER LOS OTROS ELEM. NO-DIAG. Y VER SI OTRA BARRIDA ES NECESARIA
!
      EPS=RTOL**2
      DO 250 J=1,NR
      JM1=J-1
      JP1=J+1
      LJK=JM1*N-JM1*J/2
      JJ=LJK+J
      DO 250 K=JP1,N
      KM1=K-1
      JK=LJK+K
      KK=KM1*N-KM1*K/2+K
      EPSA=(A(JK)*A(JK))/(A(JJ)*A(KK))
      EPSB=(B(JK)*B(JK))/(B(JJ)*B(KK))
      IF((EPSA.LT.EPS).AND.(EPSB.LT.EPS))GO TO 250
      GO TO 280
  250 CONTINUE
!
!   ANULA EL TRIANGULO INFERIOR DE LAS MATRICES RESULTANTES Y ESCALE
!     LOS AUTOVECTORES
  255 II=1
      DO 275 I=1,N
      BB=DSQRT(B(II))
      DO 270 K=1,N
  270 X(K,I)=X(K,I)/BB
  275 II=II+N1-I
      RETURN
!
!  ACTUALIZE LA MATRIZ D Y EMPIEZE UN NUEVO BARRIDO , SI ESTA PERMITIDO
!
  280 DO 290 I=1,N
  290 D(I)=EIGV(I)
      IF(NSWEEP.LT.NSMAX)GO TO 40
      GO TO 255
 2000 FORMAT('0 NUMERO DE BARRIDA EN *JACOBI* =',I4)
 2010 FORMAT('1',6E20.12)
 2020 FORMAT('0 *** ERROR EN LA SOLUCION STOP'/                         &
     &' MATRICES NO DEFINIDA-POSITIVA'/                                 &
     &'  II=',I4,' A(II)=',E20.12,' B(II)=',E20.12)
 2030 FORMAT(' LOS AUTOVALORES ACTUALES EN *JACOBI*  SON '/)
       END
