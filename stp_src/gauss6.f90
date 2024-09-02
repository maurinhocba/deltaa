SUBROUTINE gauss6( iset )
!
! Process Gauss information and compute nodal values for
! 3-D shell (Simo) elements
!
USE data_db
IMPLICIT NONE
! dummy arguments
INTEGER, INTENT(IN OUT) :: iset   !set number

! local variables
REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
TYPE (shl3d), POINTER, SAVE :: eset
INTEGER :: ielem,g,i,j,n,iv,nelem,nnode,nstre,ngaus,iw,nvarn,nvert,nps,naux
REAL (kind=8) :: facta,t(3),stran(3),r1,r2,lb(3)
INTEGER, ALLOCATABLE :: lnods(:)
REAL (kind=8), ALLOCATABLE :: elvar(:,:),posgp(:,:),shape(:,:),deriv(:,:,:), &
                              x0(:,:),d0(:,:),dvol(:),strsg(:,:),strea(:,:),  &
                              weigh(:),shap1(:,:),x(:,:),d(:,:)
INTEGER, POINTER :: nodes(:,:)
REAL (kind=8), POINTER :: vargs(:,:),accpn(:)
INTEGER(kind=8) :: ety


iset = iset + 1   !update set number
vargs => shl3d_vargs         !pointers to use short names
accpn => shl3d_accpn
nodes => shl3d_nodes
nps   = shl3d_nps
nvarn = shl3d_nvarn

IF( iset == 1)THEN
  eset => shl3d_head  !for first set point the head and reserve memory
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
ngaus = ABS(eset%ngaus)      !number of integration points
nvert = nnode           !number of vertex nodes
IF(nnode == 6) nvert = 3 !quadratic triangles
naux = nstre - 11       !number of additional internal variables

               ! get memory for auxiliar arrays
ALLOCATE( elvar(nstre,ngaus), posgp(2,ngaus), shape(nnode,ngaus),      &
          deriv(nnode,2,ngaus), d0(ndime,2), dvol(ngaus),weigh(ngaus), &
          x0(3,nnode), lnods(nnode),shap1(nvert,ngaus), x(3,nnode),    &
          d(ndime,2))

IF( nvarn > 0 )ALLOCATE( strsg(nvarn,ngaus), strea(nvarn,nvert) ) !get memory

ety = 6
IF(eset%ngaus < 0 ) ety = 3
CALL shapes(nvert,ngaus,shap1,deriv,posgp,weigh,2,ety)   !compute shape for lower order nodes
CALL shapes(nnode,ngaus,shape,deriv,posgp,weigh,2,ety)   !compute shape functions

DO ielem=1,nelem                    !for each element
  lnods = eset%lnods(:,ielem)       !element connectivities
  x0 = coord(1:3,lnods)    !original coordinates
  x  = coorf(1:3,lnods)    !present coordinates
  DO g=1,ngaus
    READ(16,END=100) (elvar(i,g),i=1,nstre)     !read variables from disk
    !  modify values read considering tolerances
    DO i=1,nstre
      IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
      IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) &
                                    elvar(i,g) = valor2
    END DO
    iw = 0       !initializes pointer to nodal variables
    iv = 0       !initializes pointer to Gaussian variables
    !  compute Gauss point volume
    d0 = MATMUL(x0,deriv(:,:,g))                !tangent vectors
    CALL vecpro(d0(1,1),d0(1,2),t(1))           !normal vector
    CALL vecuni(3,t(1),r1)                      !unit normal vector
    dvol(g) = r1*weigh(g)  !gauss point length
    !                            Process membrane forces
    IF( shl3d_force /= 0 )THEN
      IF( shl3d_force > 1 )THEN   !Gauss points forces desired?
        iv = iv+3                 !update pointer
        eset%elvar(iv-2:iv,g,ielem) = elvar(1:3,g)  !assign
      END IF
      IF( MOD(shl3d_force,2) == 1)THEN   !nodal forces desired
        iw = iw + 3                      !update pointer
        strsg(iw-2:iw,g)= elvar(1:3,g)   !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process moments
    IF( shl3d_momen /= 0 )THEN
      IF( shl3d_momen > 1 )THEN   !Gauss points moments desired?
        iv = iv+3                 !update pointer
        eset%elvar(iv-2:iv,g,ielem) = elvar(4:6,g)  !assign
      END IF
      IF( MOD(shl3d_momen,2) == 1)THEN   !nodal momens desired
        iw = iw + 3                      !update pointer
        strsg(iw-2:iw,g)= elvar(4:6,g)   !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process shear forces
    IF( shl3d_shear /= 0 )THEN
      IF( shl3d_shear > 1 )THEN   !Gauss points shears desired?
        iv = iv+2                 !update pointer
        eset%elvar(iv-1:iv,g,ielem) = elvar(7:8,g)  !assign
      END IF
      IF( MOD(shl3d_shear,2) == 1)THEN   !nodal shears desired
        iw = iw + 2                      !update pointer
        strsg(iw-1:iw,g)= elvar(7:8,g)   !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Compute log strains
    IF( shl3d_logst /= 0 )THEN
      d = MATMUL(x,deriv(:,:,g))                !tangent vectors
      CALL stran6(d0,d,t,r1,stran) !compute metric tensor F^T F

      CALL lgst2d(stran(1:3),r1,r2,lb(1:2))   !Log strains STRAN = three components
      lb(3) = -lb(1) - lb(2)             !isochoric approximation

      IF( shl3d_logst > 1)THEN            !Gauss points strains desired?
        iv = iv+3                         !update pointer
        eset%elvar(iv-2:iv,g,ielem) = lb    !assign
      END IF
      IF( MOD(shl3d_logst,2) == 1)THEN      !nodal strains desired
        iw = iw + 3                       !update pointer
        strsg(iw-2:iw,g)= lb              !assign to auxiliar array for smoothing
      END IF
    END IF

    !                            Process Equivalent Plastic Strain
    IF( shl3d_eqpst /= 0 )THEN
      IF( shl3d_eqpst > 1 )THEN          !Gauss points Eq. Pl strain desired?
        iv = iv+1                        !update pointer
        eset%elvar(iv,g,ielem) = elvar(9,g)   !assign
      END IF
      IF( MOD(shl3d_eqpst,2) == 1)THEN   !nodal Eq. Pl strain desired
        iw = iw + 1                      !update pointer
        strsg(iw,g)= elvar(9,g)      !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process Equivalent Von Mises Stress
    IF( shl3d_vmise /= 0 )THEN
      IF( shl3d_vmise > 1 )THEN          !Gauss points Von Mises stress desired?
        iv = iv+1                        !update pointer
        eset%elvar(iv,g,ielem) = elvar(10,g)   !assign
      END IF
      IF( MOD(shl3d_vmise,2) == 1)THEN   !nodal Eq. Mises stress desired
        iw = iw + 1                      !update pointer
        strsg(iw,g)= elvar(10,g)         !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process thickness ratio
    IF( shl3d_thrat /= 0 )THEN
      IF( shl3d_thrat > 1 )THEN          !Gauss points thickness ratio desired?
        iv = iv+1                        !update pointer
        eset%elvar(iv,g,ielem) = elvar(11,g)   !assign
      END IF
      IF( MOD(shl3d_thrat,2) == 1)THEN   !nodal thickness ratio desired
        iw = iw + 1                      !update pointer
        strsg(iw,g)= elvar(11,g)      !assign to auxiliar array for smoothing
      END IF
    END IF

    IF( naux > 0 )THEN                  !if additional variables
      n  = iv+1                         !auxiliar to first value
      iv = iv+naux                      !update pointer (last value)
      eset%elvar(n:iv,g,ielem) = elvar(12:nstre,g)  !assign
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
DEALLOCATE( elvar, posgp, shape, deriv, d, dvol, x0, lnods, weigh, shap1 )
IF( nvarn > 0 ) DEALLOCATE( strsg , strea )

RETURN
100 fin = .TRUE.      !abnormal end of file detected
RETURN
END SUBROUTINE gauss6
