 SUBROUTINE gauss5( iset )
 !
 ! Process Gauss information and compute nodal values for
 ! 3-D solid elements
 !
 USE data_db
 IMPLICIT NONE
 ! dummy arguments
 INTEGER, INTENT(IN OUT) :: iset   !set number

 ! local variables
 REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
 TYPE (sol3d), POINTER, SAVE :: eset
 INTEGER :: ielem,g,i,j,n,iv,nelem,nnode,nstre,ngaus,iw,nvarn,nps,nn,naux,etype
 REAL (kind=8) :: f(3,3),stran(6),lb(3),fld,facta,st(3), &
                  r,r1,r2,re(3,3)
 INTEGER, ALLOCATABLE :: lnods(:)
 REAL (kind=8), ALLOCATABLE :: elvar(:,:),posgp(:,:),shape(:,:),deriv(:,:,:), &
                               x0(:,:),x(:,:),cartd(:,:),dvol(:),strsg(:,:),  &
                               strea(:,:),weigh(:),shapr(:,:)
 INTEGER, POINTER :: nodes(:,:)
 REAL (kind=8), POINTER :: vargs(:,:),accpn(:)
 LOGICAL :: locsy


 iset = iset + 1   !update set number
 vargs => sol3d_vargs         !pointers to use short names
 accpn => sol3d_accpn
 nodes => sol3d_nodes
 nps   = sol3d_nps
 nvarn = sol3d_nvarn

 IF( iset == 1)THEN
   eset => sol3d_head  !for first set point the head and reserve memory
   IF( nvarn > 0 )THEN !initializes smoothing arrays
     accpn = 0.0
     vargs = 0.0
   END IF
 ELSE                  !point to next set
   eset => eset%next
 END IF

 nelem = eset%nelem      !number of elements in the set
 nnode = eset%nnode      !number of nodes per element
 nstre = eset%nstre      !number of variables at each Gauss point
 ngaus = eset%ngaus      !number of integration points
 naux = nstre - 8        !number of additional internal variables
 nn = nnode              !number of vertex nodes
 IF( sol3d_prsmo .AND. nnode >  8 )THEN   !special smoothing for 15-node prism
   ALLOCATE( shapr(nnode,ngaus) )
   CALL shapri(shapr,ngaus,nnode)   !compute GP shape functions for 15-node prism
 ELSE   !standard smoothing
   IF( nn == 20 ) nn = 8
   IF( nn == 15 ) nn = 6
   IF( nn == 12 ) nn = 6
 END IF

                ! get memory for auxiliar arrays
 ALLOCATE( elvar(nstre,ngaus), posgp(3,ngaus), shape(nn,ngaus), &
           deriv(nn,3,ngaus), cartd(nn,3), dvol(ngaus),      &
           x0(3,nn), x(3,nn), lnods(nn), weigh(ngaus))

 IF( nvarn > 0 )ALLOCATE( strsg(nvarn,ngaus), strea(nvarn,nn) ) !get memory

 CALL shapes(nn,ngaus,shape,deriv,posgp,weigh,3,eset%etype)   !compute shape function derivatives

 ! Modify shape functions for smoothing purposes
 ! transform then to extrapolation functions
 SELECT CASE (nn)
 CASE (8) !Hexahedra
   SELECT CASE (ngaus)
   CASE (1)
     shape(1:8,1) = 1d0
   CASE (2)
     shape(1:4,1) = 1d0
     shape(5:8,1) = 0d0
     shape(1:4,2) = 0d0
     shape(5:8,2) = 1d0
   CASE (4)
     shape(1,1:4)   =(/ 0.683012702,    0.683012702,    0.683012702,   -1.049038106 /)
     shape(2,1:4)   =(/ 1.549038106,   -0.183012702,   -0.183012702,   -0.183012702 /)
     shape(3,1:4)   =(/-0.183012702,    1.549038106,   -0.183012702,   -0.183012702 /)
     shape(4,1:4)   =(/ 0.683012702,    0.683012702,   -1.049038106,    0.683012702 /)
     shape(5,1:4)   =(/-0.183012702,   -0.183012702,    1.549038106,   -0.183012702 /)
     shape(6,1:4)   =(/ 0.683012702,   -1.049038106,    0.683012702,    0.683012702 /)
     shape(7,1:4)   =(/-1.049038106,    0.683012702,    0.683012702,    0.683012702 /)
     shape(8,1:4)   =(/-0.183012702,   -0.183012702,   -0.183012702,    1.549038106 /)
   CASE (8)
     shape(1,1:8)   = (/ 2.549038,  -0.683013,  -0.683013,    0.183013, -0.683013,    0.183013,   0.183013, -0.049038 /)
     shape(2,1:8)   = (/-0.683013,   2.549038,   0.183013,   -0.683013,  0.183013,   -0.683013,  -0.049038,  0.183013 /)
     shape(3,1:8)   = (/ 0.183013,  -0.683013,  -0.683013,    2.549038, -0.049038,    0.183013,   0.183013, -0.683013 /)
     shape(4,1:8)   = (/-0.683013,   0.183013,   2.549038,   -0.683013,  0.183013,   -0.049038,  -0.683013,  0.183013 /)
     shape(4,1:8)   = (/-0.683013,   0.183013,   2.549038,   -0.683013,  0.183013,   -0.049038,  -0.683013,  0.183013 /)
     shape(5,1:8)   = (/-0.683013,   0.183013,   0.183013,   -0.049038,  2.549038,   -0.683013,  -0.683013,  0.183013 /)
     shape(6,1:8)   = (/ 0.183013,  -0.683013,  -0.049038,    0.183013, -0.683013,    2.549038,   0.183013, -0.683013 /)
     shape(7,1:8)   = (/-0.049038,   0.183013,   0.183013,   -0.683013,  0.183013,   -0.683013,  -0.683013,  2.549038 /)
     shape(8,1:8)   = (/ 0.183013,  -0.049038,  -0.683013,    0.183013, -0.683013,    0.183013,   2.549038, -0.683013 /)
   CASE (27)
   END SELECT

 CASE (6)    !Prisma
   SELECT CASE (ngaus)
   CASE (1)
     shape(1:6,1) = 1d0
   CASE (2)
     IF( etype == 16 )THEN
       shape(1:3,1) =  1.366025404
       shape(4:6,1) = -0.366025404
       shape(1:3,2) = -0.366025404
       shape(4:6,2) =  1.366025404
     ELSE
       shape(1:3,1) = 1d0
       shape(4:6,1) = 0d0
       shape(1:3,2) = 0d0
       shape(4:6,2) = 1d0
     END IF
   CASE (6)
     IF( type27gp )THEN  
       !      using mid-side Gauss points
       shape(1,1:6)   = (/  1.000000,  0.000000, -1.000000,  0.000000,  1.000000,  0.000000 /)    !
       shape(2,1:6)   = (/  1.000000,  0.000000,  1.000000,  0.000000, -1.000000,  0.000000 /)    !
       shape(3,1:6)   = (/ -1.000000,  0.000000,  1.000000,  0.000000,  1.000000,  0.000000 /)    !
       shape(4,1:6)   = (/  0.000000,  1.000000,  0.000000, -1.000000,  0.000000,  1.000000 /)    !
       shape(5,1:6)   = (/  0.000000,  1.000000,  0.000000,  1.000000,  0.000000, -1.000000 /)    !
       shape(6,1:6)   = (/  0.000000, -1.000000,  0.000000,  1.000000,  0.000000,  1.000000 /)    !
     ELSE
       !    using internal Gauss points
       shape(1,1:6)   = (/ -0.455342,   -0.455342,   2.276709,    0.122008,    0.122008,  -0.610042  /)
       shape(2,1:6)   = (/  2.276709,   -0.455342,  -0.455342,   -0.610042,    0.122008,   0.122008  /)
       shape(3,1:6)   = (/ -0.455342,    2.276709,  -0.455342,    0.122008,   -0.610042,   0.122008  /)
       shape(4,1:6)   = (/  0.122008,    0.122008,  -0.610042,   -0.455342,   -0.455342,   2.276709  /)
       shape(5,1:6)   = (/ -0.610042,    0.122008,   0.122008,    2.276709,   -0.455342,  -0.455342  /)
       shape(6,1:6)   = (/  0.122008,   -0.610042,   0.122008,   -0.455342,    2.276709,  -0.455342  /)
     END IF
   CASE (7)
   CASE (8)
   CASE (9)
   END SELECT

 CASE (4)     !Tetrahedra
   SELECT CASE (ngaus)
   CASE (1)
     shape(1:4,1) = 1d0
   END SELECT

 END SELECT

 locsy =  ( sol3d_thrat /= 0 .OR. sol3d_logst /= 0 )
 etype = eset%etype

 DO ielem=1,nelem                    !for each element
   lnods = eset%lnods(1:nn,ielem)       !element connectivities
   x0 = coord(1:3,lnods)    !original coordinates
   x  = coorf(1:3,lnods)    !present coordinates
   IF( locsy )THEN
     re  = eset%dirt(:,:,ielem) !rotation matrix
   END IF
   DO g=1,ngaus
     READ(16,END=100) (elvar(i,g),i=1,nstre)  !read variables from disk
     !  modify values read considering tolerances
     DO i=1,nstre
       IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
       IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) &
                                     elvar(i,g) = valor2
     END DO
     iw = 0       !initializes pointer to nodal variables
     iv = 0       !initializes pointer to Gaussian variables
     !  Calculates transformation iso-p jacobian.
     CALL jacob5(cartd,deriv(1,1,g),dvol(g),x0,nn)
     dvol(g) = dvol(g) * weigh(g)
     !                            Process stresses
     IF( sol3d_stres /= 0 )THEN
       r = elvar(6,g) ; elvar(6,g) = elvar(5,g); elvar(5,g) = r !swap S_xz & S_yz
       IF( sol3d_stres > 1 )THEN   !Gauss points stresses desired?
         iv = iv+6                 !update pointer
         eset%elvar(iv-5:iv,g,ielem) = elvar(1:6,g)  !assign
       END IF
       IF( MOD(sol3d_stres,2) == 1)THEN   !nodal stresses desired
         iw = iw + 6                      !update pointer
         strsg(iw-5:iw,g)= elvar(1:6,g)   !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process log strains
     IF( sol3d_logst /= 0 .OR. sol3d_thrat /= 0 )THEN
       !compute strain
       cartd = MATMUL(cartd,re)
       f = MATMUL(x,cartd)                      !Deformation gradient
       stran(1) = DOT_PRODUCT(f(:,1),f(:,1))    !U^2(1,1)
       stran(2) = DOT_PRODUCT(f(:,2),f(:,2))    !U^2(2,2)
       stran(3) = DOT_PRODUCT(f(:,3),f(:,3))    !U^2(3,3)
       stran(4) = DOT_PRODUCT(f(:,1),f(:,2))    !U^2(1,2)
       stran(5) = DOT_PRODUCT(f(:,1),f(:,3))    !U^2(1,3)
       stran(6) = DOT_PRODUCT(f(:,2),f(:,3))    !U^2(2,3)
       IF( sol3d_thrat > 0)THEN            !thickness ratio desired
         !If thickness direction exists, compute log strain on sheet plane
         st = (/ stran(1), stran(2), stran(4) /)
         CALL lgst2d(st,r1,r2,lb(1:2))   !Log strains STRAN = three components
         r = SQRT( stran(3) )
         lb(3) = LOG(r)
       ELSE
         CALL lgst3d(stran,lb)   !Log strains STRAN = six components
       END IF
       IF( sol3d_logst > 1)THEN            !Gauss points strains desired?
         iv = iv+3                         !update pointer
         eset%elvar(iv-2:iv,g,ielem) = lb    !assign
       END IF
       IF( MOD(sol3d_logst,2) == 1)THEN    !nodal strains desired
         iw = iw + 3                       !update pointer
         strsg(iw-2:iw,g)= lb              !assign to auxiliar array for smoothing
       END IF
       IF( sol3d_thrat > 0)THEN            !thickness ratio desired
         IF( sol3d_thrat > 1)THEN
           iv = iv+1                       !update pointer
           eset%elvar(iv,g,ielem) = r  !assign
         END IF
         IF( MOD(sol3d_thrat,2) == 1)THEN  !nodal thickness ratio desired
           iw = iw + 1                     !update pointer
           strsg(iw,g)= r                  !assign to auxiliar array for smoothing
         END IF
       END IF
     END IF
     !                            Process Equivalent Plastic Strain
     IF( sol3d_eqpst /= 0 )THEN
       IF( sol3d_eqpst > 1 )THEN          !Gauss points Eq. Pl strain desired?
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = elvar(7,g)   !assign
       END IF
       IF( MOD(sol3d_eqpst,2) == 1)THEN   !nodal Eq. Pl strain desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= elvar(7,g)          !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process Equivalent Von Mises Stress
     IF( sol3d_vmise /= 0 )THEN
       IF( sol3d_vmise > 1 )THEN          !Gauss points Von Mises stress desired?
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = elvar(8,g)   !assign
       END IF
       IF( MOD(sol3d_vmise,2) == 1)THEN   !nodal Eq. Mises stress desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= elvar(8,g)          !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process Forming Limit Diagram
     IF( sol3d_fldma /= 0 )THEN
       ! compute fld_map ????
       fld = 0d0          !initialization only (FF)
       IF( sol3d_fldma > 1 )THEN          !Gauss points Fld Map desired?
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = fld     !assign
       END IF
       IF( MOD(sol3d_fldma,2) == 1)THEN   !nodal Fld Maps desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= fld                 !assign to auxiliar array for smoothing
       END IF
     END IF

     IF( naux > 0 )THEN                  !if additionao variables
       n  = iv+1                         !auxiliar to first value
       iv = iv+naux                      !update pointer (last value)
       eset%elvar(n:iv,g,ielem) = elvar(9:nstre,g)    !assign
     END IF

   END DO


   IF( nvarn == 0 )CYCLE                  !no variables to smooth, cycle

   !     smoothing computations


   !IF( nn > 8 )THEN
   !  CALL smotpr(ngaus,nvarn,dvol,strsg,strea,shapr,facta,nn)
   !ELSE IF( nnode == 6 )THEN   ! standard smoothing with lumped diagonal inverse
   !  CALL smotpr6(ngaus,nvarn,dvol,strsg,strea,shape,facta,nn)
   !ELSE   ! standard smoothing with lumped diagonal inverse
   !  CALL smot05(ngaus,nvarn,dvol,strsg,strea,shape,facta,nn)
   !END IF
     facta = 1d0/SUM(dvol(1:ngaus))
     strea = 0d0
     DO n=1,nn
       DO g=1,ngaus
         strea(1:nvarn,n) = strea(1:nvarn,n) + shape(n,g)*strsg(1:nvarn,g)
       END DO
     END DO

   !     assemble global smoothing matrix & rhs for stresses

   DO i=1,nn                           !for each nodal vertex
     n = nodes(lnods(i),1)                  !associated node
     accpn(n) = accpn(n) + facta  !add nodal factor (area)
     DO j=1,nvarn                         !for each smoothed variable
       vargs(j,n) = vargs(j,n) + strea(j,i)*facta  !sum
     END DO
   END DO

 END DO

 !           release auxiliar arrays
 DEALLOCATE( elvar, posgp, shape, deriv, cartd, dvol, x0, x, lnods, weigh)
 IF( nvarn > 0 ) DEALLOCATE( strsg , strea )

 RETURN
 100 fin = .TRUE.      !abnormal end of file detected
 RETURN
 END SUBROUTINE gauss5
