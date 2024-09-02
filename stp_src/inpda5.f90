 SUBROUTINE inpda5( nel, sname, etype )
 !
 !  read data for a 3-D SOLID set
 !
 USE data_db
 IMPLICIT NONE
 INTEGER, INTENT(IN) ::  nel,  & !number of elements
                         etype   !element type
 CHARACTER (len=30), INTENT(IN) :: sname

 INTEGER (kind=4) ielem,n,nvar,nvarn,i,j,locax,matno
 REAL (kind=8) :: angle(3),rs(3,3),re(3,3),alpha,beta,gamma,factor,t(3),l,t1(3),t2(3)
 TYPE (sol3d), POINTER :: eset
 LOGICAL :: sang,lang,locsy
 TYPE (udmat), POINTER :: postd

 ALLOCATE (eset)               !get memory for this set
 NULLIFY (eset%next)           !nullify pointer to next set
 eset%etype = etype

 IF( sol3d_sets == 0 )THEN   !for the first 3-D SOLID set
   sol3d_head => eset
   sol3d_nvarg = 0           !initializes
   IF( sol3d_stres > 1 ) sol3d_nvarg = sol3d_nvarg + 6
   IF( sol3d_logst > 1 ) sol3d_nvarg = sol3d_nvarg + 3
   IF( sol3d_thrat > 1 ) sol3d_nvarg = sol3d_nvarg + 1
   IF( sol3d_eqpst > 1 ) sol3d_nvarg = sol3d_nvarg + 1
   IF( sol3d_vmise > 1 ) sol3d_nvarg = sol3d_nvarg + 1
   IF( sol3d_fldma > 1 ) sol3d_nvarg = sol3d_nvarg + 1
   sol3d_nvarn = 0           !initializes
   IF( MOD(sol3d_stres,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 6
   IF( MOD(sol3d_logst,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 3
   IF( MOD(sol3d_thrat,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 1
   IF( MOD(sol3d_eqpst,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 1
   IF( MOD(sol3d_vmise,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 1
   IF( MOD(sol3d_fldma,2) == 1 ) sol3d_nvarn = sol3d_nvarn + 1
   sol3d_auxiliar_nodes = 0           !initializes
 ELSE                        !for subsequent sets
   sol3d_tail%next => eset
 END IF
 sol3d_tail => eset          !last set position

 sol3d_sets = sol3d_sets + 1 !increase number of sets
 eset%set = nsets            !set position (possibly unnecessary)

 eset%sname = sname            !set name
 eset%nelem = nel              !number of elements in the set

 READ(17) eset%ngaus,   &   !number of gauss point in each space direction
          eset%nnode,   &   !number of nodes per element and
          eset%nstre,   &   !number of stress variables to read
          locax        !local system option

  ! read user defined data
 IF( eset%nstre > 8 )THEN
   nvarn = eset%nstre - 8
   READ (17) nvar,matno
   CALL new_post_data(postd,nvar)
   postd%matno = matno
   postd%nvarv = nvarn
   DO i=1,nvar
     READ (17) postd%type(i),postd%dim(i),postd%name(:,i)
   END DO
   sol3d_nvarg = sol3d_nvarg + nvarn
 END IF
 IF( etype == 16 .OR. etype == 27 .OR. etype == 29) READ(17) type27gp,gppos27
 ALLOCATE ( eset%lnods(eset%nnode,nel),eset%matno(nel) )               !Connectivities & materials
 IF(sol3d_nvarg > 0) ALLOCATE(eset%elvar(sol3d_nvarg,eset%ngaus,nel))  !Gauss point variables

 ! read connectivities (first element material only)
 DO ielem=1,nel
   READ(17) eset%matno(ielem),(eset%lnods(n,ielem),n=1,eset%nnode)
 END DO
 IF( eset%nstre > 8 )postd%matno = eset%matno(1)

 locsy = (sol3d_thrat /= 0 .OR. sol3d_logst /= 0 )
 IF( locsy ) THEN
   ALLOCATE(eset%dirt(ndime,ndime,nel))   !local systems
   DO ielem=1,nel                         !default
     eset%dirt(:,1,ielem) = (/ 1d0,0d0,0d0 /)
     eset%dirt(:,2,ielem) = (/ 0d0,1d0,0d0 /)
     eset%dirt(:,3,ielem) = (/ 0d0,0d0,1d0 /)
   END DO
 END IF
 ! read Euler angles for local system definition

 READ(17) angle
 factor = ASIN(1d0)/90d0         !pi/180
 sang = ANY( angle /= 0d0 )  !set angles
 IF( sang )THEN
   alpha = angle(1)*factor
   beta  = angle(2)*factor
   gamma = angle(3)*factor
   CALL inrotm(alpha,beta,gamma,rs(1,1))   !set rotation matrix
 END IF
 j = eset%nnode/2                  !number of nodes in each face
 IF( etype == 27 .OR. etype == 29 ) j =3
 DO ielem=1,nel
   READ(17) angle(:)               !read element angles
   IF( locsy )THEN                 !If local systems required
     lang = ANY(angle /= 0d0)      !.TRUE. if Euler angles exists for this element
     IF( lang )THEN
       alpha = angle(1)*factor
       beta  = angle(2)*factor
       gamma = angle(3)*factor
       CALL inrotm(alpha,beta,gamma,re(1,1))  !element rotation matrix
       IF( sang ) re = MATMUL(rs,re)          !compound rotation matrix
     ELSE IF( sang ) THEN
       re = rs
     END IF
     IF( locax /= 0 )THEN    !use connectivities to define shell normal
       t = 0d0                    !initializes
       DO i=1,j                   !sums aristas
         t = t + coord(:,eset%lnods(i+j,ielem)) - coord(:,eset%lnods(i,ielem))
       END DO
       CALL vecuni(3,t,l)         !unit vector
       eset%dirt(:,3,ielem) = t   !assign
       SELECT CASE(locax)
       CASE (1)
         l = (t(2)*t(2)+t(3)*t(3)) !component in th Y-Z plane
         IF( l < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  Y-Z plane
           t2 = (/ -t(3), 0d0, t(1) /) !choose t2 orthogonal to global Y direction
           CALL vecuni(3,t2,l)
           CALL vecpro(t2,t,t1)
         ELSE       !         SELECT local y=t(:,1) in the global YZ plane
           t1 = (/ 0d0, -t(3), t(2)  /)
           t2 = (/ l, -t(2)*t(1), -t(3)*t(1) /)
           CALL vecuni(3,t1,l)   !     normalizes t1 & t2
           CALL vecuni(3,t2,l)
         END IF
       CASE (2)
         l = (t(3)*t(3)+t(1)*t(1)) !component in th Z-X plane
         IF( l  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  Z-Y plane
           t2 = (/ t(2), -t(1), 0d0 /) !choose t2 orthogonal to global Z direction
           CALL vecuni(3,t2,l)
           CALL vecpro(t2,t,t1)
         ELSE       !         SELECT local z=t(:,1) in the global ZX plane
           t1 = (/ t(3), 0d0, -t(1)  /)
           t2 = (/  -t(1)*t(2), l, -t(3)*t(2) /)
           CALL vecuni(3,t1,l)   !     normalizes t1 & t2
           CALL vecuni(3,t2,l)
         END IF
       CASE (3)
         l = (t(1)*t(1)+t(2)*t(2)) !component in th X-Y plane
         IF( l  < 1.0d-5) THEN         !If t(:,3) is almost orthogonal to  X-Y plane
           t2 = (/ 0d0, t(3), -t(2) /) !choose t(:,2) orthogonal to global X direction
           CALL vecuni(3,t2,l)
           CALL vecpro(t2,t,t1)
         ELSE       !         SELECT local x=t(:,1) in the global xy plane
           t1 = (/ -t(2), t(1) , 0d0 /)
           t2 = (/ -t(1)*t(3), -t(2)*t(3), l /)
           CALL vecuni(3,t1,l)   !     normalizes t1 & t2
           CALL vecuni(3,t2,l)
         END IF
       END SELECT
       IF( lang .OR. sang )THEN
         alpha= DOT_PRODUCT(t1,re(:,1))
         beta = DOT_PRODUCT(t2,re(:,1))
         l = SQRT(alpha**2+beta**2)
         alpha= alpha/l
         beta = beta/l
         t  = t1
         t1 =  alpha*t1+ beta*t2
         t2 = -beta*t  + alpha*t2
       END IF
       eset%dirt(:,1,ielem) = t1    !assign
       eset%dirt(:,2,ielem) = t2    !assign
     ELSE
       IF( sang .OR. lang )  eset%dirt(:,:,ielem) = re
     END IF
   END IF
 END DO
 IF( etype == 27 .AND. eset%nnode == 15  )THEN
   eset%first_l = last_label
   DO ielem=1,nel
     DO j=1,3
       last_label = last_label + 1
       eset%lnods(9+j,ielem) = last_label
     END DO
   END DO
   sol3d_auxiliar_nodes = sol3d_auxiliar_nodes + 3*nel
 END IF

 RETURN
 END SUBROUTINE inpda5
