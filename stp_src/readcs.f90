 SUBROUTINE readcs ( )
 !
 !   read surface data to compute surface contact press
 !
 USE cont_db
 USE data_db, ONLY : ndime
 IMPLICIT NONE

 CHARACTER (len=30) :: sname
 LOGICAL :: spress,swrink,cpress
 INTEGER (kind=4) :: ncnod,nsegm
 TYPE (surf_db), POINTER :: surf

 INTEGER (kind=4) :: i,j

 nsurf = 0                      !initializes number of surfaced read
 CALL ini_srf(shead,stail)      !initializes surface data_base

 DO    !for each stored surface
   READ(44) sname,spress,swrink,cpress  !name.
   IF( .NOT.spress .AND. .NOT.swrink ) EXIT  !why not CYCLE (because a special line is read)
   READ(44) ncnod,nsegm                 !number of nodes and segments
   CALL new_surf(surf)                  !get memory for surface
   ALLOCATE( surf%lcnod(ncnod), surf%lcseg(ndime,nsegm) )  !get memory for nodes and connectivities
   READ(44) (surf%lcnod(i),i=1,ncnod)                      !read nodes
   READ(44)((surf%lcseg(j,i),j=1,ndime),i=1,nsegm)         !read connectivities
   surf%sname = sname                   !name
   surf%press = spress                  !write press
   surf%cpress= cpress                  !never know
   surf%wrink = swrink                  !write wrinkles (gap information)
   surf%ncnod = ncnod                   !number of nodes
   surf%nsegm = nsegm                   !number of segments
   nsurf = nsurf+1                      !increase number of surfaces
   CALL add_srf(surf,shead,stail)       !add to list of surfaces
   IF( spress )press = .TRUE.           !global flag (press)
   IF( swrink )wrink = .TRUE.           !global flag (wrink)
 END DO

 RETURN
 END SUBROUTINE readcs
