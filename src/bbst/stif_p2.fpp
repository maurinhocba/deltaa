   IF( ASSOCIATED(e%gausv) )THEN
     ttti = .TRUE.
   ELSE IF( plast )THEN
     aux = ABS(stra1(1)-1d0)+ABS(stra1(2)-1d0)+ABS(stra1(3))+ &
          (ABS(stra1(4))+ABS(stra1(5))+ABS(stra1(6)))*thick/2d0
     ttti = aux > minstr
   ELSE
     ttti = mtype == 6  ! it is neccessary for rubbers
   END IF

   IF( secty == 12 )THEN              !for standard solid section
     IF( ttti )THEN
     !IF( .FALSE. )THEN
       ! Integrate throught the thickness Elasto-plastic consistent matrix

       thnew = thick*e%lb             !thickness at new configuration

       IF( mtype == 1 )THEN
         IF( ASSOCIATED(e%gausv) )THEN
           lambd(1:nlayr) = e%gausv(8,1:nlayr)           !consistency parameter
           pflag = ANY( lambd /= 0d0 )
         ELSE
           pflag = .FALSE.
         END IF
       END IF

       dmatx = 0d0                         ! Initializes D matrix

       DO l=1,nlayr                        ! for each layer
         zk = thf(l)*thnew                 ! Z coordinate
         zk2 = zk*zk                       ! Z coord squared
         IF( mtype == 6 )THEN     !hyperelastic materials
           stran = stra1(1:3)+stra1(4:6)*zk  ! layer U^2
           CALL lgst14(stran,r1,r2,lb(1),'STIF14',found) !compute eigenvalues
           lb(3) = 1d0/lb(1)/lb(2)
           CALL rubberps(chi,lb,mat%matdef(8),stran,mat=cm,r1=r1,r2=r2)
         ELSE
           IF( logst )THEN
             stran = stra1(1:3)+stra1(4:6)*zk  ! layer U^2
             ! inverse of U^2
             aux  = (stran(1)*stran(2)-stran(3)**2)
             u2(1) = stran(2)/aux
             u2(2) = stran(1)/aux
             u2(3) = -stran(3)/aux
           END IF

           ! compute layer constitutive matrix (tangent algorithmic)
           IF( mtype == 1 )THEN   ! for an isotropic homogeneous material
             IF( pflag )THEN
               efpst = e%gausv(4,l) + lambd(l)    !Total Eff. plastic strain
               stres(1:3) = e%gausv(9:11,l)       !stresses (Hencky or 2PK)
             ELSE
               efpst = 0d0
             END IF
             CALL dmat14(dummy,dummy,dummy,stres(1),lambd(l),efpst,d(1,1),.FALSE.)
           ELSE IF( mtype == 5 )THEN     ! Orthotropic material
             !d(1,1) = c(1); d(1,2) = c(2); d(2,2) = c(3); d(3,3) = c(4)
           END IF

           IF( logst )THEN           !large strain (modify using the metric tensor)
             cm(1,1) = d(1,1)*u2(1)*u2(1) + 2d0*d(1,3)*u2(1)*u2(3) &
                     + d(3,3)*u2(3)*u2(3)

             cm(1,2) = d(1,3)*u2(1)*u2(3) + d(1,2)*u2(1)*u2(2)     &
                     + d(3,3)*u2(3)*u2(3) + d(2,3)*u2(3)*u2(2)

             cm(1,3) = ( d(1,1)*u2(1)*u2(3) + d(1,2)*u2(1)*u2(3)   &
                     + d(1,3)*( u2(1)*u2(1) + u2(1)*u2(2) + u2(3)*u2(3) )   &
                     + d(2,3)*u2(3)*u2(3)                                   &
                     + d(3,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/2d0

             cm(2,2) = d(2,2)*u2(2)*u2(2) + 2d0*d(2,3)*u2(2)*u2(3) &
                     + d(3,3)*u2(3)*u2(3)

             cm(2,3) = ( d(1,2)*u2(2)*u2(3) + d(1,3)*u2(3)*u2(3)   &
                     + d(2,2)*u2(2)*u2(3)                                   &
                     + d(2,3)*( u2(1)*u2(2) + u2(2)*u2(2) +u2(3)*u2(3) )    &
                     + d(3,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/2d0

             cm(3,3) = (d(3,3)*(u2(1)*u2(1)                        &
                     + 2d0*u2(1)*u2(2)+u2(2)*u2(2))                         &
                     + d(1,1)*u2(3)*u2(3) + 2d0*d(1,2)*u2(3)*u2(3)          &
                     + 2d0*d(1,3)*( u2(1)*u2(3) + u2(2)*u2(3) )             &
                     + d(2,2)*u2(3)*u2(3)                                   &
                     + 2d0*d(2,3)*( u2(1)*u2(3) + u2(2)*u2(3) ) )/4d0
           ELSE      !small strain
             cm      = d
           END IF
         END IF

         cm = cm*wei(l)
         DO i=1,3
           DO j=i,3
             dmatx(i,j)     = dmatx(i,j)     + cm(i,j)
             dmatx(i,j+3)   = dmatx(i,j+3)   + cm(i,j)*zk
             dmatx(i+3,j+3) = dmatx(i+3,j+3) + cm(i,j)*zk2
           END DO
         END DO
       END DO

       dmatx(2,4) = dmatx(1,5)
       dmatx(3,4) = dmatx(1,6)
       dmatx(3,5) = dmatx(2,6)

       voli = area1*thick   ! Initial Vol
       k = 0
       DO i=1,6
         DO j=i,6
           k = k+1
           daux(k) = voli*dmatx(i,j)
         END DO
       END DO

     ELSE     !use elastic integrated matrix
       daux = dm*area1
     END IF
   ELSE ! secty == 13
     daux = dm*area1
   END IF
