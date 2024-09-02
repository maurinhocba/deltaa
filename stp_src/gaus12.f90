 SUBROUTINE gaus12( iset, etype )
 !
 ! Process Gauss information and compute nodal values for
 ! 3-D BST Shell elements
 !
 USE data_db
 USE flc_db,ONLY: flc_tp, hpflc, srch_flc, lblflc
 IMPLICIT NONE
 ! dummy arguments
 INTEGER, INTENT(IN OUT) :: iset   !set number
 INTEGER, INTENT(IN) :: etype      !element type

 ! local variables
 REAL (kind=8), PARAMETER :: valor1 = 0.9999D99, valor2 = 0.999D-99
 TYPE (bst), POINTER, SAVE :: eset
 INTEGER :: ielem,i,j,n,l,iv,nelem,nnode,nstre,ngaus,iw,nvarn,nps, &
            lnods(6),naux,nno,olbl,locax
 REAL (kind=8) :: stran(6), r1,r2,lb(3),fld,s(3),str(3),            &
                  ls(3,2),gs(3,2),t(3),l11,l12,l21,l22,long,        &
                  x0(3,6),x(3,6),a(3,4),b(3,4),tt(3,3),             &
                  stra0(6),area2,angle
 REAL (kind=8), ALLOCATABLE :: elvar(:),strsg(:)
 INTEGER, POINTER :: nodes(:,:)
 REAL (kind=8), POINTER :: vargs(:,:),accpn(:)
 LOGICAL :: sides(3), found, needst
 TYPE(flc_tp),POINTER:: flc


 iset = iset + 1   !update set number
 vargs => bst_vargs         !pointers to use short names
 accpn => bst_accpn
 nodes => bst_nodes
 nps   = bst_nps
 nvarn = bst_nvarn

 IF( iset == 1)THEN
   eset => bst_head  !for first set point the head and reserve memory
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
 locax = eset%locax      !local system option 
 naux = nstre - 13       !number of additional internal variables

 IF( nnode == 3 )THEN !triangles
   nno = 6  !patch nodes
 ELSE !IF ( nnode = 4 )THEN
   nno = 4  !element nodes only
 END IF
                ! get memory for auxiliar arrays
 ALLOCATE( elvar(nstre) )

 IF( nvarn > 0 )ALLOCATE( strsg(nvarn) ) !get memory

 olbl = -1
 DO ielem=1,nelem                    !for each element

   IF( lblflc(eset%matno(ielem)) /= olbl)THEN
     olbl = lblflc(eset%matno(ielem))
     CALL srch_flc(hpflc,olbl,found,flc)
   END IF
   lnods(1:nno) = eset%lnods(:,ielem)       !element connectivities
   DO n=1,nno     !for each node
     l =lnods(n)           !node
     IF( l > 0 )THEN       !if node exists
       DO i=1,ndime             !for each direction
         x0(i,n) = coord(i,l)   !original coord.
         x(i,n) = coorf(i,l)    !actual coord
       END DO
       IF(n > 3) sides(n-3) = .TRUE.
     ELSE
       sides(n-3) = .FALSE.
     END IF
   END DO

   READ(16,END=100) (elvar(i),i=1,nstre)  !read variables from disk
   !  modify values read considering tolerances
   DO i=1,nstre
     IF (ABS(elvar(i)) > valor1) elvar(i) = SIGN(valor1,elvar(i))
     IF (ABS(elvar(i)) < valor2 .AND. elvar(i) /= 0.0) elvar(i) = valor2
   END DO
   iw = 0       !initializes pointer to nodal variables
   iv = 0       !initializes pointer to Gaussian variables
   !         compute local systems
   IF( etype == 12 )THEN  !updated Lagrangean elements
     ls(1:3,1) = x(1:3,2) - x(1:3,1) !first side vector
     ls(1:3,2) = x(1:3,3) - x(1:3,1) !second side vector
     CALL vecpro(ls(1,1),ls(1,2),t(1))       !element normal
     CALL vecuni(3,ls(1,1),long)             !Unit local X vector
     CALL vecuni(3,t(1),long)                !Unit normal vector
     CALL vecpro(t(1),ls(1,1),ls(1,2))       !Unit local Y vector

     !         SELECT local x=t1 in the global xy plane
     gs(1:3,1) = (/ -t(2), t(1) , 0d0 /)     !intersection between planes
     !         of course local y = t2 = t3 x t1
     gs(1:3,2) = (/ -t(1)*t(3), -t(2)*t(3), (t(1)*t(1)+t(2)*t(2))/)
     IF(ABS( gs(3,2) ) < 1.0d-5) THEN      !If planes are almost paralell
       gs(1:3,1) = (/  1d0, 0d0, 0d0 /)      !Local X = Global X
       gs(1:3,2) = (/  0d0, 1d0, 0d0 /)      !Local Y = Global Y
     ELSE                                  !If they are different
       CALL vecuni(ndime,gs(1,1),long)       !normalizes t1
       CALL vecuni(ndime,gs(1,2),long)       !normalizes t2
     END IF

     !compute Change-of-Base matriz (Lambda)
     l11 = DOT_PRODUCT(gs(:,1),ls(:,1))     ! Gx . Lx
     l12 = DOT_PRODUCT(gs(:,1),ls(:,2))     ! Gx . Ly
     l21 = DOT_PRODUCT(gs(:,2),ls(:,1))     ! Gy . Lx
     l22 = DOT_PRODUCT(gs(:,2),ls(:,2))     ! Gy . Ly
     ! compute a and b factors, local system and 2A
     CALL axes12(x0,a,b,sides,tt,area2)
   ELSE !( etype == 14 .OR. etype == 25)THEN !Total Lagrangean elements
     !compute Change-of-Base matriz (Lambda)
     angle = eset%angle(ielem)
     l11 = COS(angle)
     l12 =-SIN(angle)
     l21 =-l12
     l22 = l11
     IF ( etype == 14 )THEN !Total Lagrangean elements
       ! compute a and b factors, local system and 2A
       CALL axes14(x0,a,b,sides,tt,angle,area2,locax)
     ELSE !( etype == 25 )THEN !Total Lagrangean quadrilateral
       ! compute derivatives at the center
       CALL axes25(x0,x,a,b,tt,angle,area2)
     END IF
   END IF
   !                            Process forces
   IF( bst_force /= 0 )THEN
     s   = elvar(1:3)                 ! Membrane stresses
     str(1) = l11**2*s(1)+2d0*l12*l11*s(3)+l12**2*s(2)
     str(2) = l21**2*s(1)+2d0*l21*l11*s(3)+l22**2*s(2)
     str(3) = l11*l21*s(1)+(l22*l11+l21*l12)*s(3)+l22*l12*s(2)
     IF( bst_force > 1 )THEN   !Gauss points forces desired?
       iv = iv+3               !update pointer
       eset%elvar(iv-2:iv,1,ielem) = str  !assign
     END IF
     IF( MOD(bst_force,2) == 1)THEN   !nodal forces desired
       iw = iw + 3                    !update pointer
       strsg(iw-2:iw)= str          !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process moments
   IF( bst_momen /= 0 )THEN
     s   = elvar(4:6)                 ! Membrane stresses
     str(1) = l11**2*s(1)+2d0*l12*l11*s(3)+l12**2*s(2)
     str(2) = l21**2*s(1)+2d0*l21*l11*s(3)+l22**2*s(2)
     str(3) = l11*l21*s(1)+(l22*l11+l21*l12)*s(3)+l22*l12*s(2)
     IF( bst_momen > 1 )THEN   !Gauss points moments desired?
       iv = iv+3               !update pointer
       eset%elvar(iv-2:iv,1,ielem) = str  !assign
     END IF
     IF( MOD(bst_momen,2) == 1)THEN   !nodal moments desired
       iw = iw + 3                    !update pointer
       strsg(iw-2:iw)= str          !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process shear forces
   IF( bst_shear /= 0 )THEN
     s(1:2)   = elvar(7:8)                 ! Shear stresses (???)
     str(1) = l11*s(1)+l12*s(2)
     str(2) = l21*s(1)+l22*s(2)
     IF( bst_shear > 1 )THEN   !Gauss points forces desired?
       iv = iv+2               !update pointer
       eset%elvar(iv-1:iv,1,ielem) = str(1:2)  !assign
     END IF
     IF( MOD(bst_shear,2) == 1)THEN   !nodal forces desired
       iw = iw + 2                    !update pointer
       strsg(iw-1:iw)= str          !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process log strains
   needst = (bst_fldma /= 0) .OR. bst_wfldFZ .OR. bst_wfldSZ
   IF( (bst_logst /= 0) .OR. (bst_curva /= 0) .OR. needst )THEN
     IF ( etype == 25 )THEN !Total Lagrangean elements
       stran(1) = DOT_PRODUCT(tt(1:3,1),tt(1:3,1))
       stran(2) = DOT_PRODUCT(tt(1:3,2),tt(1:3,2))
       stran(3) = DOT_PRODUCT(tt(1:3,1),tt(1:3,2))
       stran(4:6) = 0d0 !do not compute curvatures
       stra0(4:6) = 0d0 !do not compute curvatures
     ELSE
       CALL stra14(a,b,x0,sides,stra0) !only original curvatures are important
       CALL stra14(a,b,x ,sides,stran) !these is what we expect in TLF
     END IF

     CALL lgst2d(stran(1:3),r1,r2,lb(1:2))   !Log strains STRAN = three components
     !lb(3) = LOG(elvar(9))               !
     lb(3) = -lb(1) - lb(2)             !isochoric approximation

     ! Curvatures in the common system
     stran(4:6) = stran(4:6) - stra0(4:6) ! change in curvatures
     stran(1) = stran(4)*l11**2+stran(5)*l12**2+stran(6)*l11*l12
     stran(2) = stran(5)*l11**2+stran(4)*l12**2-stran(6)*l11*l12
     stran(3) = (-stran(4)+stran(5))*l12*l11+stran(6)*(l11**2-l12**2)/2d0

     IF( bst_logst > 1)THEN            !Gauss points strains desired?
       iv = iv+3                         !update pointer
       eset%elvar(iv-2:iv,1,ielem) = lb   !assign
     END IF
     IF( MOD(bst_logst,2) == 1)THEN      !nodal strains desired
       iw = iw + 3                       !update pointer
       strsg(iw-2:iw)= lb              !assign to auxiliar array for smoothing
     END IF
     IF( bst_curva > 1)THEN              !Gauss points curvatures desired?
       iv = iv+3                         !update pointer
       eset%elvar(iv-2:iv,1,ielem) = stran(1:3)    !assign
     END IF
     IF( MOD(bst_curva,2) == 1)THEN      !nodal curvatures desired
       iw = iw + 3                       !update pointer
       strsg(iw-2:iw)= stran(1:3)      !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process thickness ratio
   IF( bst_thrat /= 0 )THEN
     IF( bst_thrat > 1 )THEN            !Gauss points Thickness ratio desired?
       iv = iv+1                        !update pointer
       eset%elvar(iv,1,ielem) = elvar(9)   !assign
     END IF
     IF( MOD(bst_thrat,2) == 1)THEN     !nodal Thickness ratio desired
       iw = iw + 1                      !update pointer
       strsg(iw)= elvar(9)          !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process Equivalent Plastic Strain
   IF( bst_eqpst /= 0 )THEN
     IF( bst_eqpst > 1 )THEN          !Gauss points Eq. Pl strain desired?
       iv = iv+1                        !update pointer
       eset%elvar(iv,1,ielem) = elvar(10)   !assign
       iv = iv+1                        !update pointer
       eset%elvar(iv,1,ielem) = elvar(11)   !assign
     END IF
     IF( MOD(bst_eqpst,2) == 1)THEN   !nodal Eq. Pl strain desired
       iw = iw + 1                      !update pointer
       strsg(iw)= elvar(10)          !assign to auxiliar array for smoothing
       iw = iw + 1                      !update pointer
       strsg(iw)= elvar(11)          !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process Equivalent Von Mises Stress
   IF( bst_vmise /= 0 )THEN
     IF( bst_vmise > 1 )THEN           !Gauss points Von Mises stress desired?
       iv = iv+1                       !update pointer
       eset%elvar(iv,1,ielem) = elvar(12)   !assign
       iv = iv+1                       !update pointer
       eset%elvar(iv,1,ielem) = elvar(13)   !assign
     END IF
     IF( MOD(bst_vmise,2) == 1)THEN     !nodal Eq. Mises stress desired
       iw = iw + 1                      !update pointer
       strsg(iw)= elvar(12)         !assign to auxiliar array for smoothing
       iw = iw + 1                      !update pointer
       strsg(iw)= elvar(13)         !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process FLD Map
   IF( bst_fldma /= 0 )THEN
     CALL CalDisFLD(flc%npt,flc%cv,lb(2),lb(1),fld)
     IF( bst_fldma > 1 )THEN            !Gauss points Fld Map desired
       iv = iv+1                        !update pointer
       eset%elvar(iv,1,ielem) = fld     !assign
     END IF
     IF( MOD(bst_fldma,2) == 1)THEN     !nodal Fld Maps desired
       iw = iw + 1                      !update pointer
       strsg(iw)= fld                   !assign to auxiliar array for smoothing
     END IF
   END IF
   !                            Process Forming Limit Diagram
   IF( bst_wfldFZ .OR. bst_wfldSZ )THEN
     iv = iv + 2                       !update pointer
     eset%elvar(iv-1:iv,1,ielem) = fld     !assign
     iw = iw + 2                       !update pointer
     strsg(iw-1:iw)= lb(1:2)           !assign to auxiliar array for smoothing
   END IF

   IF( naux > 0 )THEN                  !if additionao variables
     n  = iv+1                         !auxiliar to first value
     iv = iv+naux                      !update pointer (last value)
     eset%elvar(n:iv,1,ielem) = elvar(14:nstre)    !assign
   END IF

   IF( nvarn == 0 )CYCLE                  !no variables to smooth, cycle

   !     smoothing computations

   !     assemble global smoothing matrix & rhs for stresses

   DO i=1,nnode                           !for each nodal vertex
     n = nodes(lnods(i),1)                  !associated node
     accpn(n) = accpn(n) + 1./area2       !add nodal factor (area)
     DO j=1,nvarn                         !for each smoothed variable
       vargs(j,n) = vargs(j,n) + strsg(j)/area2  !sum
     END DO
   END DO

 END DO

 !           release auxiliar arrays
 DEALLOCATE( elvar )
 IF( nvarn > 0 ) DEALLOCATE( strsg )

 RETURN
 100 fin = .TRUE.      !abnormal end of file detected
 RETURN
 END SUBROUTINE gaus12
