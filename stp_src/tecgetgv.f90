 SUBROUTINE tecgetgv(flag,gv,nodes,nps,ngv)
 !
 !  transfer global node information to elemental nodes
 !
 USE data_db, ONLY : npoin, ndime, ndofn, ndoft, addof, &
                     t_ngv, lc_var, sh_var, ps_var, lc_vara, sh_vara, ps_vara, &
                     wtdisp, wsdisp, widisp, waddof, wtempe, wveloc, waccel, weuler, wangve, wangac, &
                     coors, coord, coorf, euler, displ, dispi, veloc, anvel, accel, anacc, tempe, psi
 USE cont_db, ONLY : wwear, wwrink, wpress, wwork, wrinkn, presn
 IMPLICIT NONE
 ! dummy arguments
 LOGICAL, INTENT(IN) :: flag  !.FALSE. coordinates  .TRUE. all nodal variables
 INTEGER(kind=4), INTENT(IN) :: nps  !number of active nodes
 INTEGER(kind=4), POINTER :: nodes(:,:)  !number of active nodes
 INTEGER(kind=4), INTENT(OUT) :: ngv  !number of global variables
 REAL(kind=4), POINTER :: gv(:,:)     !global variables
 ! local variables
 INTEGER(kind=4) :: i,n,iv,nd
 REAL(kind=8), POINTER :: x(:,:)      !coordinates


 IF( flag )THEN  !Solution
   x => coorf
   ngv = t_ngv           !number of nodal global variables
   ALLOCATE( gv(ngv,nps) )  !get memory for global variables
   nd = ndofn - ndime
   !pass all the nodal values
   DO i=1,npoin
     n= nodes(i,1)                       !local node number
     IF( n == 0 )CYCLE                   !if node does not belong to this mesh
     iv = 1                              !initializes pointer
     !deformed geometry
     gv(iv:ndime,n) = x(:,i)                      ! coordinates
     iv = iv + ndime                                !update pointer
     !total displacements
     IF( wtdisp )THEN
       gv(iv:iv+ndime-1,n) = displ(:,i)           !displacements
       iv = iv + ndime                              !update pointer
     END IF
     !stage displacements
     IF( wsdisp )THEN
       gv(iv:iv+ndime-1,n) = displ(:,i) + coord(:,n) - coors(:,n)
       iv = iv + ndime                              !update pointer
     END IF
     !initial displacements
     IF( widisp )THEN
       gv(iv:iv+ndime-1,n) = dispi(:,i)
       iv = iv + ndime                              !update pointer
     END IF
     !additional displacements
     IF( addof .AND. waddof )THEN
       gv(iv:iv+1,n) = psi(:,i)
       iv = iv + 2                                  !update pointer
     END IF
     !Velocities
     IF( wveloc )THEN
       gv(iv:iv+ndime-1,n) = veloc(:,i)           !velocities
       iv = iv + ndime                              !update pointer
     END IF
     !Accelerations
     IF( waccel )THEN
       gv(iv:iv+ndime-1,n) = accel(:,i)           !accelerations
       iv = iv + ndime                              !update pointer
     END IF
     !Euler Angles
     IF( weuler )THEN
       IF( ndime == 2 )THEN
         gv(iv,n) = euler(1,i)
         iv = iv + 1
       ELSE
         gv(iv:iv+2,n) = euler(:,i)               !Euler angles
         iv = iv + 3                                !update pointer
       END IF
     END IF
     !Angular velocities
     IF( wangve )THEN
       gv(iv:iv+nd-1,n) = anvel(:,i)              !Angular velocities
       iv = iv + nd                                 !update pointer
     END IF
     !Angular accelerations
     IF( wangac )THEN
       gv(iv:iv+nd-1,n) = anacc(:,i)              !Angular accelerations
       iv = iv + nd                                 !update pointer
     END IF
     ! Temperatures
     IF( wtempe )THEN
       gv(iv:iv+ndoft-1,n) = tempe(:,i)           !temperatures
       iv = iv + ndoft                              !update pointer
     END IF
     ! Wearing Work
     IF( wwear  )THEN
       gv(iv,n) = wwork(i)                        !friction work
       iv = iv + 1                                  !update pointer
     END IF
     ! Gap-information
     IF( wwrink )THEN
       gv(iv:iv+1,n) = wrinkn(1:2,i)              !Gap actual and minimum
       gv(iv+2,n) = wrinkn(1,i)-wrinkn(2,i)         !Gap difference
       iv = iv + 3                                  !update pointer
     END IF
     ! binder-Pressure
     IF( wpress )THEN
       gv(iv:iv+1,n) = presn(:,i)                 !binder pressure
       iv = iv + 2                                  !update pointer
     END IF
   END DO

 ELSE   ! Grid
   IF(wsdisp) THEN          !if stage displacements
     x => coors                 !point to stage displacements
   ELSE                     !if initial displacements
     x => coord                 !point to initial displacements
   END IF
   ngv = ndime           !number of global variables
   ALLOCATE( gv(ngv,nps) )  !get memory for global variables
   DO n=1,npoin                !for each global node
     IF( nodes(n,1) > 0 )gv(:,nodes(n,1)) = x(:,n)  !pass required information
   END DO
   IF( ASSOCIATED(sh_vara) )DEALLOCATE(sh_vara,ps_vara,lc_vara)
   ALLOCATE(sh_vara(ndime),ps_vara(ndime),lc_vara(ndime))
   sh_vara = 1  !shared
   ps_vara = 0  !active
   lc_vara = 1  !nodal
 END IF

 RETURN
 END SUBROUTINE tecgetgv


 SUBROUTINE tecgenvarshls(iset) !Generate Variables Share List and Passive Variables list
 ! generate string for shared variables  ==> sh_var
 ! generate string for passive variables ==> ps_var
 ! generate string for variable location ==> lc_var
 USE data_db, ONLY : sh_var, ps_var, lc_var, t_ngv, t_ev, t_nsv, t_npv
 IMPLICIT NONE
 !Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: iset
 !local variables
 CHARACTER(len=3 ) :: ch       !function to convert a 3 digit integer into a string
 INTEGER (kind=4) :: i,sv,pv,j,k,l,lv


 ! initializes shared list including global nodal variables
 sv = 15             !initializes
 sh_var='VARSHARELIST=(['
 DO i=1,t_ngv          !for each nodal variable
   sv = sv + 4         !update pointer
   sh_var(sv-3:sv) = ch(i)//','    !
 END DO

 k = t_ngv             !last variable included

 ! initializes passive list including previous sets variables
 pv = 16             !initializes
 ps_var='PASSIVEVARLIST=['
 DO i=1,iset-1
   DO j=1,2                       !for nodal and cell-centered variables
     DO l=1,t_ev(j,i)               !for each elemental variables
       k = k + 1                       !update number of variables
       pv = pv + 4                     !update pointer
       ps_var(pv-3:pv) = ch(k)//','    !
     END DO
   END DO
 END DO

 ! include list of variables of this set
 !for nodal
 DO l=1,t_ev(1,iset)               !for each nodal variables
   k = k + 1                       !update number of variables
   sv = sv + 4                     !update pointer
   sh_var(sv-3:sv) = ch(k)//','    !
 END DO
 !for cell-centered variables
 ! initializes list including variables location
 IF(t_ev(2,iset) > 0 )THEN
   lv = 14             !initializes
   lc_var='VARLOCATION=(['
   DO l=1,t_ev(2,iset)               !for each nodal variables
     k = k + 1                       !update number of variables
     lv = lv + 4                     !update pointer
     lc_var(lv-3:lv) = ch(k)//','    !
   END DO
   lc_var(lv:lv+15) =']=CELLCENTERED )'  !close (overwrite last comma)
   lv = lv + 15
 END IF
 sh_var(sv:sv+1) = '])'  !close (overwrite last comma)

 ! close passive list including next sets variables
 DO i=iset+1,10
   DO j=1,2                       !for nodal and cell-centered variables
     DO l=1,t_ev(j,i)               !for each elemental variables
       k = k + 1                       !update number of variables
       pv = pv + 4                     !update pointer
       ps_var(pv-3:pv) = ch(k)//','    !
     END DO
   END DO
 END DO
 ps_var(pv:pv) = ']' !close (overwrite last comma if passive variables exist)
 t_nsv = (sv-15)/4
 t_npv = (pv-16)/4

 RETURN
 END SUBROUTINE tecgenvarshls

 SUBROUTINE tecgenvarshlsbin(iset) !Generate Variables Share List and Passive Variables list
 ! generate array for shared variables  ==> sh_vara
 ! generate array for passive variables ==> ps_vara
 ! generate array for variable location ==> lc_vara
 USE data_db, ONLY : sh_vara, ps_vara, lc_vara, t_ngv, t_ev, t_nsv, t_npv
 IMPLICIT NONE
 !Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: iset
 !local variables
 INTEGER (kind=4) :: i,j,k,l


 ! initializes shared list including global nodal variables
 k = t_ngv             !last global variable included
 sh_vara = 0           !initializes shared variables as NOT SHARED
 sh_vara(1:k) = 1      !global variables are shared variables
 ps_vara = 0           !initializes passive variables as ACTIVE
 lc_vara = 1           !initializes variables as NODAL

 ! initializes passive list including previous sets variables
 DO i=1,iset-1
   DO j=1,2                       !for nodal and cell-centered variables
     DO l=1,t_ev(j,i)             !for each elemental variables
       k = k + 1                  !update number of variables
       ps_vara(k) = 1             !PASSIVE
     END DO
   END DO
 END DO

 ! include list of variables of this set
 !for nodal
 DO l=1,t_ev(1,iset)               !for each nodal variables
   k = k + 1                       !update number of variables
   sh_vara(k) = 1                  !SHARED
 END DO
 !for cell-centered variables
 ! initializes list including variables location
 IF(t_ev(2,iset) > 0 )THEN
   DO l=1,t_ev(2,iset)               !for each nodal variables
     k = k + 1                       !update number of variables
     lc_vara(k) = 0                  !CELL-CENTERED
   END DO
 END IF

 ! close passive list including next sets variables
 DO i=iset+1,10
   DO j=1,2                       !for nodal and cell-centered variables
     DO l=1,t_ev(j,i)             !for each elemental variables
       k = k + 1                  !update number of variables
       ps_vara(k) = 1             !PASSIVE
     END DO
   END DO
 END DO
 t_nsv = COUNT(sh_vara == 1)
 t_npv = COUNT(ps_vara == 1)
 RETURN
 END SUBROUTINE tecgenvarshlsbin
