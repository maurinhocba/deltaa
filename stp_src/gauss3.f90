 SUBROUTINE gauss3( iset )
 !
 ! Process Gauss information and compute nodal values for
 ! 2-D solid elements
 !
 USE data_db
 USE flc_db,ONLY: flc_tp, hpflc, srch_flc, lblflc
 IMPLICIT NONE
 ! dummy arguments
 INTEGER, INTENT(IN OUT) :: iset   !set number

 ! local variables
 REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
 TYPE (sol2d), POINTER, SAVE :: eset
 INTEGER :: ielem,g,i,j,n,iv,nelem,nnode,nstre,ngaus,iw,nvarn,nvert,nps,ver
 REAL (kind=8) :: f(2,2), stran(3), r1,r2,lb(3),r0,r,fld,facta,dir(2) !,p
 INTEGER, ALLOCATABLE :: lnods(:)
 REAL (kind=8), ALLOCATABLE :: elvar(:,:),posgp(:,:),shape(:,:),deriv(:,:,:), &
                               x0(:,:),x(:,:),cartd(:,:),dvol(:),strsg(:,:),  &
                               strea(:,:),weigh(:),shap1(:,:)
 INTEGER, POINTER :: nodes(:,:)
 REAL (kind=8), POINTER :: vargs(:,:),accpn(:)
 LOGICAL ::  needst,found
 TYPE(flc_tp),POINTER:: flc
 INTEGER :: isec,osec,lbl,etype


 iset = iset + 1   !update set number
 vargs => sol2d_vargs         !pointers to use short names
 accpn => sol2d_accpn
 nodes => sol2d_nodes
 nps   = sol2d_nps
 nvarn = sol2d_nvarn


 IF( iset == 1)THEN
   eset => sol2d_head  !for first set point the head and reserve memory
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
 ver   = eset%ver        !Gauss Point Positions

 SELECT CASE (nnode)     !compute number of vertex
 CASE(3,4)               !for linear elements
   nvert = nnode
 CASE(6,8,9)             !for quadratic elements
   nvert = nnode/2
 END SELECT
                ! get memory for auxiliar arrays
 ALLOCATE( elvar(nstre,ngaus), posgp(2,ngaus), shape(nnode,ngaus), &
           deriv(nnode,2,ngaus), cartd(nnode,2), dvol(ngaus),      &
           x0(2,nnode), x(2,nnode), lnods(nnode), weigh(ngaus),    &
           shap1(nvert,ngaus) )

 IF( nvarn > 0 )ALLOCATE( strsg(nvarn,ngaus), strea(nvarn,nvert) ) !get memory

 etype = 3
 IF( ver == 2 )etype = -3
 CALL shapes(nvert,ngaus,shap1,deriv,posgp,weigh,2,etype)   !compute shape functions
 CALL shapes(nnode,ngaus,shape,deriv,posgp,weigh,2,etype)   !compute shape functions

 lbl = -1
 DO ielem=1,nelem                    !for each element
   IF( lbl /= lblflc(eset%matno(ielem)))THEN
     lbl = lblflc(eset%matno(ielem))
     CALL srch_flc(hpflc,lbl,found,flc)
   END IF
   lnods = eset%lnods(:,ielem)       !element connectivities
   x0 = coord(1:2,lnods)    !original coordinates
   x  = coorf(1:2,lnods)    !present coordinates
   IF( sol2d_thrat /= 0 )dir = eset%dirt(:,ielem)  !normal direction (TLF)
   DO g=1,ngaus
     READ(16,END=100) (elvar(i,g),i=1,nstre)  !read variables from disk
     !  modify values read considering tolerances
     DO i=1,nstre
       IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
       IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) &
                                     elvar(i,g) = valor2
     END DO
     !  Calculates transformation iso-p jacobian.
     CALL jacob3(cartd,deriv(1,1,g),dvol(g),x0,nnode)
     dvol(g) = dvol(g) * weigh(g)
     iw = 0       !initializes pointer to nodal variables
     iv = 0       !initializes pointer to Gaussian variables
     !                            Process stresses
     IF( sol2d_stres /= 0 )THEN
       !! next line pp deviator stresses in the plane and pressure
       !!p = (elvar(1,g)+elvar(2,g)+elvar(4,g))/3d0
       !!elvar(1,g) = elvar(1,g)-p
       !!elvar(2,g) = elvar(2,g)-p
       !!elvar(4,g) = p
       IF( sol2d_stres > 1 )THEN   !Gauss points stresses desired?
         iv = iv+4                 !update pointer
         eset%elvar(iv-3:iv,g,ielem) = elvar(1:4,g)  !assign
       END IF
       IF( MOD(sol2d_stres,2) == 1)THEN   !nodal stresses desired
         iw = iw + 4                      !update pointer
         strsg(iw-3:iw,g)= elvar(1:4,g)   !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process log strains
     needst = (sol2d_fldma /= 0) .OR. sol2d_wfldFZ .OR. sol2d_wfldSZ
     IF( sol2d_logst /= 0 .OR. sol2d_shtst /= 0 .OR. sol2d_thrat /= 0 .OR. needst)THEN
       !compute strain
       f = MATMUL(x,cartd)                          !Deformation gradient
       stran(1) = f(1,1)*f(1,1) + f(2,1)*f(2,1)     !U^2
       stran(2) = f(1,2)*f(1,2) + f(2,2)*f(2,2)
       stran(3) = f(1,1)*f(1,2) + f(2,1)*f(2,2)
       CALL lgst2d(stran,r1,r2,lb(1:2))   !Log strains STRAN = three components
                                          !            lb = principal values
       SELECT CASE (ntype)   !compute lb(3)  according to problem type
       CASE (1)    !for plane stress problems
         lb(3) = -lb(1) - lb(2)           !isochoric approximation
       CASE (2)    !for Plane strain problems
         lb(3) = 0d0
       CASE (3)    !for axilsymmetric problems
         r0 = DOT_PRODUCT(x0(1,:),shape(:,g))     !original X_1 coordinates
         r  = DOT_PRODUCT(x(1,:),shape(:,g))      !present X_1 coordinates
         lb(3) = LOG( r/r0 )                      !log strain (dir Z)
       END SELECT

       IF( sol2d_logst > 1)THEN            !Gauss points strains desired?
         iv = iv+3                         !update pointer
         eset%elvar(iv-2:iv,g,ielem) = lb    !assign
       END IF
       IF( MOD(sol2d_logst,2) == 1)THEN    !nodal strains desired
         iw = iw + 3                       !update pointer
         strsg(iw-2:iw,g)= lb              !assign to auxiliar array for smoothing
       END IF
       IF( sol2d_shtst > 0)THEN            !Gauss points strains desired?
         !    U component along dir  ( dir . U . dir )
         r0 = r1*dir(1) + r2*dir(2)        ! r1 . dir
         r  = r1*dir(2) - r2*dir(1)        ! r2 . dir
         stran(2) = lb(1)*r**2 + lb(2)*r0**2 ! x strain
         stran(3) = lb(1)*r0**2 + lb(2)*r**2 ! thickness strain
         IF( lb(3) >= stran(2) )THEN
           stran(1) = lb(3)                  ! circunferential strain
         ELSE
           stran(1) = stran(2)
           stran(2) = lb(3)
         END IF
         IF( sol2d_shtst > 1)THEN            !Gauss points strains desired?
           iv = iv+3                         !update pointer
           eset%elvar(iv-2:iv,g,ielem) = stran    !assign
         END IF
         IF( MOD(sol2d_logst,2) == 1)THEN    !nodal strains desired
           iw = iw + 3                       !update pointer
           strsg(iw-2:iw,g)= stran              !assign to auxiliar array for smoothing
         END IF
       END IF
       IF( sol2d_thrat > 0)THEN            !thickness ratio desired
         !    U component along dir  ( dir . U . dir )
         r0 = r1*dir(1) + r2*dir(2)        ! r1 . dir
         r  = r1*dir(2) - r2*dir(1)        ! r2 . dir
         r  = EXP(lb(1))*r0**2 + EXP(lb(2))*r**2  ! dir . (ri Li ri) . dir
         IF( sol2d_thrat > 1)THEN
           iv = iv+1                       !update pointer
           eset%elvar(iv,g,ielem) = r  !assign
         END IF
         IF( MOD(sol2d_thrat,2) == 1)THEN  !nodal thickness ratio desired
           iw = iw + 1                     !update pointer
           strsg(iw,g)= r                  !assign to auxiliar array for smoothing
         END IF
       END IF
     END IF
     !                            Process Equivalent Plastic Strain
     IF( sol2d_eqpst /= 0 )THEN
       IF( sol2d_eqpst > 1 )THEN          !Gauss points Eq. Pl strain desired?
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = elvar(5,g)   !assign
       END IF
       IF( MOD(sol2d_eqpst,2) == 1)THEN   !nodal Eq. Pl strain desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= elvar(5,g)          !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process Equivalent Von Mises Stress
     IF( sol2d_vmise /= 0 )THEN
       IF( sol2d_vmise > 1 )THEN          !Gauss points Von Mises stress desired?
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = elvar(6,g)   !assign
       END IF
       IF( MOD(sol2d_vmise,2) == 1)THEN   !nodal Eq. Mises stress desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= elvar(6,g)          !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process Forming Limit Diagram
     IF( sol2d_fldma /= 0 )THEN
       ! Principal strains to pass must be adequately choosed
       !CALL CalDisFLD(flc%npt,flc%cv,lb(2),lb(1),fld)
       fld = 0d0
       IF( sol2d_fldma > 1 )THEN            !Gauss points Fld Map desired
         iv = iv+1                        !update pointer
         eset%elvar(iv,g,ielem) = fld     !assign
       END IF
       IF( MOD(sol2d_fldma,2) == 1)THEN   !nodal Fld Maps desired
         iw = iw + 1                      !update pointer
         strsg(iw,g)= fld                   !assign to auxiliar array for smoothing
       END IF
     END IF
     !                            Process Forming Zone & Safety Zone
     IF( sol2d_wfldFZ .OR. sol2d_wfldSZ )THEN
       iv = iv + 2                       !update pointer
       eset%elvar(iv-1:iv,g,ielem) = fld !assign
       iw = iw + 2                       !update pointer
       strsg(iw-1:iw,g)= lb(1:2)           !assign to auxiliar array for smoothing
     END IF
   END DO

   IF( nvarn == 0 )CYCLE                  !no variables to smooth, cycle

   !     smoothing computations

   CALL smot06(ngaus,nvarn,dvol,strsg,strea,shap1,facta,nvert)

   !     assemble global smoothing matrix & rhs for stresses

   DO i=1,nvert                           !for each nodal vertex
     n = nodes(lnods(i),1)                  !associated node
     accpn(n) = accpn(n) + facta  !add nodal factor (area)
     DO j=1,nvarn                         !for each smoothed variable
       vargs(j,n) = vargs(j,n) + strea(j,i)*facta  !sum
     END DO
   END DO

 END DO

 !           release auxiliar arrays
 DEALLOCATE( elvar, posgp, shape, deriv, cartd, dvol, x0, x, lnods, weigh, shap1 )
 IF( nvarn > 0 ) DEALLOCATE( strsg , strea )

 RETURN
 100 fin = .TRUE.      !abnormal end of file detected
 RETURN
 END SUBROUTINE gauss3
