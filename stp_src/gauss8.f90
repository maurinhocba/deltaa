SUBROUTINE gauss8( iset )
!
! Process Gauss information and compute nodal values for
! 3-D /beam elements
!
USE data_db
IMPLICIT NONE
! dummy arguments
INTEGER, INTENT(IN OUT) :: iset   !set number

! local variables
REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
TYPE (beame), POINTER, SAVE :: eset
INTEGER :: ielem,g,i,j,n,iv,nelem,nnode,nstre,ngaus,iw,nvarn,nps
REAL (kind=8) :: facta
INTEGER, ALLOCATABLE :: lnods(:)
REAL (kind=8), ALLOCATABLE :: elvar(:,:),posgp(:,:),shape(:,:),deriv(:,:,:), &
                              x0(:,:),d(:),dvol(:),strsg(:,:),strea(:,:),    &
                              weigh(:)
INTEGER, POINTER :: nodes(:,:)
REAL (kind=8), POINTER :: vargs(:,:),accpn(:)


iset = iset + 1   !update set number
vargs => beame_vargs         !pointers to use short names
accpn => beame_accpn
nodes => beame_nodes
nps   = beame_nps
nvarn = beame_nvarn

IF( iset == 1)THEN
  eset => beame_head  !for first set point the head and reserve memory
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

               ! get memory for auxiliar arrays
ALLOCATE( elvar(nstre,ngaus), posgp(1,ngaus), shape(nnode,ngaus),     &
          deriv(nnode,1,ngaus), d(ndime), dvol(ngaus),weigh(ngaus),   &
          x0(ndime,nnode), lnods(nnode) )

IF( nvarn > 0 )ALLOCATE( strsg(nvarn,ngaus), strea(nvarn,nnode) ) !get memory

CALL shapes(nnode,ngaus,shape,deriv,posgp,weigh,1,8)   !compute shape functions

!                  compute a special type of weigthing functions
IF(ngaus == 1)THEN
  shape(1:nnode,1) = 1d0
ELSE IF( ngaus == 2 )THEN
  IF( nnode == 2 )THEN
    shape(1:2,1:2) = RESHAPE((/ 1.732, -0.732, &
                               -0.732,  1.732  /),(/2,2/) )
  ELSE     !nnode = 3
    shape(1:3,1:2) = RESHAPE((/ 1.732, 0.5, -0.732, &
                               -0.732, 0.5,  1.732  /),(/3,2/))
  END IF
ELSE ! ngaus = 3   nnode = 3
  shape(1:3,1:3) = RESHAPE((/ 1d0, 0d0, 0d0, &     !just an approximation
                              0d0, 1d0, 0d0, &
                              0d0, 0d0, 1d0 /),(/3,3/))
END IF

DO ielem=1,nelem                    !for each element
  lnods = eset%lnods(:,ielem)       !element connectivities
  x0 = coord(1:ndime,lnods)    !original coordinates
  !x  = coorf(1:ndime,lnods)    !present coordinates
  DO g=1,ngaus
    !  modify values read considering tolerances
    DO i=1,nstre
      IF (ABS(elvar(i,g)) > valor1) elvar(i,g) = SIGN(valor1,elvar(i,g))
      IF (ABS(elvar(i,g)) < valor2 .AND. elvar(i,g) /= 0.0) &
                                    elvar(i,g) = valor2
    END DO
    iw = 0       !initializes pointer to nodal variables
    iv = 0       !initializes pointer to Gaussian variables
    ! compute Gauss point volume
    d  = MATMUL(x0,deriv(:,1,g))                !tangent vector
    dvol(g) = SQRT(DOT_PRODUCT(d,d))*weigh(g)  !gauss point length
    READ(16,END=100) (elvar(i,g),i=1,nstre)     !read variables from disk
    !                            Process membrane forces
    IF( beame_force /= 0 )THEN
      IF( beame_force > 1 )THEN   !Gauss points forces desired?
        iv = iv+3                 !update pointer
        eset%elvar(iv-2:iv,g,ielem) = elvar(1:3,g)  !assign
      END IF
      IF( MOD(beame_force,2) == 1)THEN   !nodal forces desired
        iw = iw + 3                      !update pointer
        strsg(iw-2:iw,g)= elvar(1:3,g)   !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process moments
    IF( beame_momen /= 0 )THEN
      IF( beame_momen > 1 )THEN   !Gauss points moments desired?
        iv = iv+3                 !update pointer
        eset%elvar(iv-2:iv,g,ielem) = elvar(4:6,g)  !assign
      END IF
      IF( MOD(beame_momen,2) == 1)THEN   !nodal momens desired
        iw = iw + 3                      !update pointer
        strsg(iw-2:iw,g)= elvar(4:6,g)   !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process Equivalent Plastic Strain
    IF( beame_eqpst /= 0 )THEN
      IF( beame_eqpst > 1 )THEN          !Gauss points Eq. Pl strain desired?
        iv = iv+1                        !update pointer
        eset%elvar(iv,g,ielem) = elvar(7,g)   !assign
      END IF
      IF( MOD(beame_eqpst,2) == 1)THEN   !nodal Eq. Pl strain desired
        iw = iw + 1                      !update pointer
        strsg(iw,g)= elvar(7,g)      !assign to auxiliar array for smoothing
      END IF
    END IF
    !                            Process Equivalent Von Mises Stress
    !IF( beame_vmise /= 0 )THEN
    !  IF( beame_vmise > 1 )THEN          !Gauss points Von Mises stress desired?
    !    iv = iv+1                        !update pointer
    !    eset%elvar(iv,g,ielem) = elvar(8,g)   !assign
    !  END IF
    !  IF( MOD(beame_vmise,2) == 1)THEN   !nodal Eq. Mises stress desired
    !    iw = iw + 1                      !update pointer
    !    strsg(iw,g)= elvar(8,g)    !assign to auxiliar array for smoothing
    !  END IF
    !END IF
    !                            Process area ratio
    IF( beame_arrat /= 0 )THEN
      IF( beame_arrat > 1 )THEN          !Gauss points area ratio desired?
        iv = iv+1                        !update pointer
        eset%elvar(iv,g,ielem) = elvar(nstre,g)   !assign
      END IF
      IF( MOD(beame_arrat,2) == 1)THEN   !nodal area ratio desired
        iw = iw + 1                      !update pointer
        strsg(iw,g)= elvar(nstre,g)      !assign to auxiliar array for smoothing
      END IF
    END IF
  END DO

  IF( nvarn == 0 )CYCLE                  !no variables to smooth, cycle

  !     smoothing computations

  CALL smot08(ngaus,nvarn,nnode,dvol,strsg,strea,shape,facta)

  !     assemble global smoothing matrix & rhs for stresses

  DO i=1,nnode                           !for each nodal vertex
    n = nodes(lnods(i),1)                  !associated node
    accpn(n) = accpn(n) + facta  !add nodal factor (area)
    DO j=1,nvarn                         !for each smoothed variable
      vargs(j,n) = vargs(j,n) + strea(j,i)*facta  !sum
    END DO
  END DO
END DO

!           release auxiliar arrays
DEALLOCATE( elvar, posgp, shape, deriv, d, dvol, x0, lnods, weigh )
IF( nvarn > 0 ) DEALLOCATE( strsg , strea )

RETURN
100 fin = .TRUE.      !abnormal end of file detected
RETURN
END SUBROUTINE gauss8
