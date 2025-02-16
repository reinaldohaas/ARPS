!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE LLXY_MODULE
!
! This module handles transformations between model grid coordinates and 
!   latitude-longitude coordinates. The actual transformations are done through
!   the map_utils module. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wrf_llxy_module

!   use gridinfo_module
!   use list_module
   use map_utils
!   use module_debug
!   use misc_definitions_module
 
   ! Parameters
   integer, parameter :: MAX_SOURCE_LEVELS = 1
 
   ! Variables
   integer :: current_nest_number
   integer :: SOURCE_PROJ = 0
   ! The following arrays hold values for all available domains 
   ! NOTE: The entries in the arrays for "domain 0" are used for projection
   !       information of user-specified source data
   type (proj_info), dimension(-MAX_SOURCE_LEVELS:MAXDOM) :: proj_stack
 
   ! The projection and domain that we have computed constants for
   integer :: computed_proj = INVALID
   integer :: computed_domain = INVALID
 
   contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: push_source_projection
   !
   ! Purpose: 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine push_source_projection(iprojection, user_stand_lon, user_truelat1, user_truelat2, &
                                  user_dxkm, user_dykm, user_dlat, user_dlon, user_known_x, &
                                  user_known_y, user_known_lat, user_known_lon, earth_radius)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: iprojection
      real, intent(in) :: user_stand_lon, user_truelat1, user_truelat2, user_dxkm, user_dykm, &
                          user_dlat, user_dlon, &
                          user_known_x, user_known_y, user_known_lat, user_known_lon
      real, intent(in), optional :: earth_radius

      SOURCE_PROJ = SOURCE_PROJ-1
      if (SOURCE_PROJ < -MAX_SOURCE_LEVELS) then
!         call mprintf(.true.,ERROR,'In push_user_projection(), too many levels of user projections.')
         call arpsstop('ERROR: In push_user_projection(), too many levels of user projections.',1)
      end if
  
      call map_init(proj_stack(SOURCE_PROJ))

      if (iprojection == PROJ_LATLON) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_MERC) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      truelat1=user_truelat1, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
!      else if (iprojection == PROJ_CYL) then
!         call mprintf(.true.,ERROR,'Should not have PROJ_CYL as projection for ' &
!                          //'source data in push_source_projection()')
!  
!      else if (iprojection == PROJ_CASSINI) then
!         call mprintf(.true.,ERROR,'Should not have PROJ_CASSINI as projection for ' &
!                          //'source data in push_source_projection()')
!  
      else if (iprojection == PROJ_LC) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)

      else if (iprojection == PROJ_ALBERS_NAD83) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_PS) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)

      else if (iprojection == PROJ_PS_WGS84) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_GAUSS) then
         call map_set(iprojection, proj_stack(SOURCE_PROJ), &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      nlat=nint(user_dlat), &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_ROTLL) then
  ! BUG: Implement this projection.
  
      end if
     
   end subroutine push_source_projection
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: pop_source_projection
   !
   ! Purpose: 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine pop_source_projection()
 
      implicit none
  
      SOURCE_PROJ = SOURCE_PROJ+1
      
!      call mprintf((SOURCE_PROJ > 0), ERROR, &
!                   'In pop_user_projection(), projection stack has overflowed.')
      IF (SOURCE_PROJ > 0) THEN
      call arpsstop('ERROR: In pop_user_projection(), projection stack has overflowed.',1)
      END IF
 
   end subroutine pop_source_projection
 
 
!#ifdef _METGRID
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: set_domain_projection
   !
   ! Purpose: 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine set_domain_projection(iprojection, user_stand_lon, user_truelat1, user_truelat2, &
                                  user_dxkm, user_dykm, user_dlat, user_dlon, &
                                  user_xdim, user_ydim, user_known_x, &
                                  user_known_y, user_known_lat, user_known_lon, &
                                  user_pole_lat, user_pole_lon, earth_radius)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: iprojection
      integer, intent(in) :: user_xdim, user_ydim
      real, intent(in) :: user_stand_lon, user_truelat1, user_truelat2, &
                          user_dxkm, user_dykm, user_dlat, user_dlon, &
                          user_known_x, user_known_y, user_known_lat, user_known_lon, &
                          user_pole_lat, user_pole_lon
      real, intent(in), optional :: earth_radius
  
      current_nest_number = 1

      call map_init(proj_stack(current_nest_number))
  
      if (iprojection == PROJ_LATLON) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_MERC) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      truelat1=user_truelat1, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_CYL) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      stdlon=user_stand_lon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_CASSINI) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      latinc=user_dlat, &
                      loninc=user_dlon, &
                      dx=user_dxkm,        &
                      dy=user_dykm,        &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      lat0=user_pole_lat, &
                      lon0=user_pole_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_LC) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_ALBERS_NAD83) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      truelat1=user_truelat1, &
                      truelat2=user_truelat2, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_PS) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm, &
                      r_earth=earth_radius)

      else if (iprojection == PROJ_PS_WGS84) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      truelat1=user_truelat1, &
                      stdlon=user_stand_lon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      knowni=user_known_x, &
                      knownj=user_known_y, &
                      dx=user_dxkm)
  
      else if (iprojection == PROJ_GAUSS) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      nlat=nint(user_dlat), &
                      loninc=user_dlon, &
                      r_earth=earth_radius)
  
      else if (iprojection == PROJ_ROTLL) then
         call map_set(iprojection, proj_stack(current_nest_number), &
                      ixdim=user_xdim, &
                      jydim=user_ydim, &
                      phi=user_dlat, &
                      lambda=user_dlon, &
                      lat1=user_known_lat, &
                      lon1=user_known_lon, &
                      stagger=HH, &
                      latinc=user_dykm, &
                      loninc=user_dxkm, &
                      r_earth=earth_radius)
  
      end if
     
   end subroutine set_domain_projection
!#endif


!#ifdef _GEOGRID
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   ! Name: compute_nest_locations
!   !
!   ! Purpose: This routine computes the variables necessary in determining the 
!   !   location of all nests without reference to the parent or coarse domains.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   subroutine compute_nest_locations()
! 
!      implicit none
!  
!      ! Local variables
!      integer :: i
!      real :: temp_known_x, temp_known_y, temp_known_lat, temp_known_lon, &
!              temp_dxkm, temp_dykm, temp_dlat, temp_dlon
!  
!      ! Set location of coarse/mother domain
!      call map_init(proj_stack(1))
!  
!      if (iproj_type == PROJ_LATLON) then
!         call map_set(iproj_type, proj_stack(1), &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      latinc=dykm, &
!                      loninc=dxkm)
!   
!      else if (iproj_type == PROJ_MERC) then
!         call map_set(iproj_type, proj_stack(1), &
!                      truelat1=truelat1, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      dx=dxkm)
!  
!      else if (iproj_type == PROJ_CYL) then
!         call map_set(iproj_type, proj_stack(1), &
!                      latinc=dlatdeg, &
!                      loninc=dlondeg, &
!                      stdlon=stand_lon)
!  
!      else if (iproj_type == PROJ_CASSINI) then
!         call map_set(iproj_type, proj_stack(1), &
!                      latinc=dlatdeg, &
!                      loninc=dlondeg, &
!                      dx=dxkm,       &
!                      dy=dykm,       &
!                      stdlon=stand_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      lat0=pole_lat, &
!                      lon0=pole_lon, &
!                      lat1=known_lat, &
!                      lon1=known_lon)
!  
!      else if (iproj_type == PROJ_LC) then
!         call map_set(iproj_type, proj_stack(1), &
!                      truelat1=truelat1, &
!                      truelat2=truelat2, &
!                      stdlon=stand_lon, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      dx=dxkm)
!  
!      else if (iproj_type == PROJ_ALBERS_NAD83) then
!         call map_set(iproj_type, proj_stack(1), &
!                      truelat1=truelat1, &
!                      truelat2=truelat2, &
!                      stdlon=stand_lon, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      dx=dxkm)
!  
!      else if (iproj_type == PROJ_PS) then
!         call map_set(iproj_type, proj_stack(1), &
!                      truelat1=truelat1, &
!                      stdlon=stand_lon, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      dx=dxkm)
!
!      else if (iproj_type == PROJ_PS_WGS84) then
!         call map_set(iproj_type, proj_stack(1), &
!                      truelat1=truelat1, &
!                      stdlon=stand_lon, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      knowni=known_x, &
!                      knownj=known_y, &
!                      dx=dxkm)
!  
!      else if (iproj_type == PROJ_GAUSS) then
!         call map_set(iproj_type, proj_stack(current_nest_number), &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      nlat=nint(dykm), &
!                      loninc=dxkm)
!  
!      else if (iproj_type == PROJ_ROTLL) then
!         call map_set(iproj_type, proj_stack(1), &
!                      ixdim=ixdim(1), &
!                      jydim=jydim(1), &
!                      phi=phi, &
!                      lambda=lambda, &
!                      lat1=known_lat, &
!                      lon1=known_lon, &
!                      latinc=dykm, &
!                      loninc=dxkm, &
!                      stagger=HH)
!   
!      end if
!  
!      ! Now we can compute lat/lon <-> x/y for coarse domain
!      call select_domain(1)
!  
!      ! Call a recursive procedure to find the lat/lon of the centerpoint for 
!      !   each domain
!      do i=2,n_domains
!  
!         temp_known_x = real(ixdim(i))/2.
!         temp_known_y = real(jydim(i))/2.
!
!         call find_known_latlon(i, temp_known_x, temp_known_y, &
!                                temp_known_lat, temp_known_lon, &
!                                temp_dxkm, temp_dykm, temp_dlat, temp_dlon)
!   
!         if (iproj_type == PROJ_LATLON) then
!            call map_set(iproj_type, proj_stack(i), &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         latinc=temp_dlat, &
!                         loninc=temp_dlon)
!   
!         else if (iproj_type == PROJ_MERC) then
!            call map_set(iproj_type, proj_stack(i), &
!                         truelat1=truelat1, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         dx=temp_dxkm)
!    
!         else if (iproj_type == PROJ_CYL) then
!            call mprintf(.true.,ERROR,'Don''t know how to do nesting with PROJ_CYL ' &
!                                      //'in compute_nest_locations()')
!  
!         else if (iproj_type == PROJ_CASSINI) then
!            call map_set(iproj_type, proj_stack(i), &
!                         latinc=temp_dlat, &
!                         loninc=temp_dlon, &
!                         dx=temp_dxkm,  &
!                         dy=temp_dykm,  &
!                         stdlon=stand_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         lat0=pole_lat, &
!                         lon0=pole_lon, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon)
!    
!         else if (iproj_type == PROJ_LC) then
!            call map_set(iproj_type, proj_stack(i), &
!                         truelat1=truelat1, &
!                         truelat2=truelat2, &
!                         stdlon=stand_lon, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         dx=temp_dxkm)
!    
!         else if (iproj_type == PROJ_ALBERS_NAD83) then
!            call map_set(iproj_type, proj_stack(i), &
!                         truelat1=truelat1, &
!                         truelat2=truelat2, &
!                         stdlon=stand_lon, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         dx=temp_dxkm)
!   
!         else if (iproj_type == PROJ_PS) then
!            call map_set(iproj_type, proj_stack(i), &
!                         truelat1=truelat1, &
!                         stdlon=stand_lon, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         dx=temp_dxkm)
!
!         else if (iproj_type == PROJ_PS_WGS84) then
!            call map_set(iproj_type, proj_stack(i), &
!                         truelat1=truelat1, &
!                         stdlon=stand_lon, &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         knowni=temp_known_x, &
!                         knownj=temp_known_y, &
!                         dx=temp_dxkm)
!   
!         else if (iproj_type == PROJ_GAUSS) then
!            call map_set(iproj_type, proj_stack(current_nest_number), &
!                         lat1=temp_known_lat, &
!                         lon1=temp_known_lon, &
!                         nlat=nint(temp_dykm), &
!                         loninc=temp_dxkm)
!   
!         else if (iproj_type == PROJ_ROTLL) then
!   ! BUG: Implement this projection.
!   
!         end if
!  
!      end do
! 
!   end subroutine compute_nest_locations
! 
! 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   ! Name: find_known_latlon
!   !
!   ! Purpose: This recursive routine computes the latitude and longitude for a 
!   !   specified x/y location in the given nest number, and also computes the
!   !   grid spacing
!   !
!   ! NOTE: This routine assumes that xytoll will work correctly for the 
!   !       coarse domain.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   recursive subroutine find_known_latlon(n, rx, ry, rlat, rlon, dx, dy, dlat, dlon)
! 
!      implicit none
!  
!      ! Arguments
!      integer, intent(in) :: n
!      real, intent(in) :: rx, ry
!      real, intent(out) :: rlat, rlon, dx, dy, dlat, dlon
!  
!      ! Local variables
!      real :: x_in_parent, y_in_parent
!  
!      if (n == 1) then   ! Stopping case for the recursion
!  
!         dx = dxkm 
!         dy = dykm 
!         dlat = dlatdeg 
!         dlon = dlondeg 
!         call ij_to_latlon(proj_stack(current_nest_number), rx, ry, rlat, rlon)
!  
!         return
!  
!      else               ! Recursive case
!   
!         x_in_parent = (rx - ((parent_grid_ratio(n)+1.)/2.)) &
!                      / parent_grid_ratio(n) + parent_ll_x(n)
!         y_in_parent = (ry - ((parent_grid_ratio(n)+1.)/2.)) &
!                      / parent_grid_ratio(n) + parent_ll_y(n)
!   
!         call find_known_latlon(parent_id(n), x_in_parent, y_in_parent, rlat, rlon, dx, dy, dlat, dlon)
!   
!         dx = dx / parent_grid_ratio(n)
!         dy = dy / parent_grid_ratio(n)
!         dlat = dlat / parent_grid_ratio(n)
!         dlon = dlon / parent_grid_ratio(n)
!      end if 
! 
!   end subroutine find_known_latlon
!
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   ! Name: compute_nest_level_info
!   !
!   ! Purpose: This routine computes the parameters describing a nesting level for 
!   !          NMM grids.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   subroutine compute_nest_level_info()
!
!      implicit none
!
!      ! Local variables
!      integer :: i, nest_level, temp
!      type (list) :: level_list 
!
!      call list_init(level_list)
!
!      ! Set location of coarse/mother domain
!      call map_init(proj_stack(1))
!
!      call map_set(PROJ_ROTLL, proj_stack(1), &
!                   ixdim=ixdim(1), &
!                   jydim=jydim(1), &
!                   phi=phi, &
!                   lambda=lambda, &
!                   lat1=known_lat, &
!                   lon1=known_lon, &
!                   latinc=dykm, &
!                   loninc=dxkm, &
!                   stagger=HH)
!
!      parent_ur_x(1) = real(ixdim(1))
!      parent_ur_y(1) = real(jydim(1))
!
!      do i=2,n_domains
!
!         nest_level = get_nest_level(i)
!
!         if (.not. list_search(level_list, ikey=nest_level, ivalue=temp)) then
!
!            call list_insert(level_list, ikey=nest_level, ivalue=nest_level)
!
!            ixdim(nest_level) = ixdim(1)*(3**(nest_level-1))-(3**(nest_level-1)-1)
!            jydim(nest_level) = jydim(1)*(3**(nest_level-1))-(3**(nest_level-1)-1)
!
!            parent_ur_x(nest_level) = ixdim(nest_level)
!            parent_ur_y(nest_level) = jydim(nest_level)
!
!            call map_set(PROJ_ROTLL, proj_stack(nest_level), &
!                         ixdim = ixdim(nest_level), &
!                         jydim = jydim(nest_level), &
!                         phi    = phi, &
!                         lambda = lambda, &
!                         lat1=known_lat, &
!                         lon1=known_lon, &
!                         latinc=(dykm/real((3**(nest_level-1)))), &
!                         loninc=(dxkm/real((3**(nest_level-1)))), &
!                         stagger=HH)
!         end if
!
!      end do
!
!      call list_destroy(level_list)
!
!   end subroutine compute_nest_level_info
!
!   
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   ! Name: get_domain_resolution
!   !
!   ! Purpose:
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   subroutine get_domain_resolution(dom_dx, dom_dy)
!
!      implicit none
!
!      ! Arguments
!      real, intent(out) :: dom_dx, dom_dy
!
!      ! The proj_info structure only stores dx, so set both dom_dx and dom_dy to dx
!      dom_dx = proj_stack(current_nest_number)%dx
!      dom_dy = proj_stack(current_nest_number)%dx
!
!   end subroutine get_domain_resolution
!
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   ! Name: get_nest_level
!   !
!   ! Purpose: This function returns, given a grid ID number, the nesting level of
!   !   that domain; the coarse domain is taken to have nesting level 1.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   function get_nest_level(i)
!      
!      implicit none
!
!      ! Arguments
!      integer, intent(in) :: i
!
!      ! Local variables
!      integer :: j
!
!      ! Return value
!      integer :: get_nest_level
!
!      ! If argument is the coarse domain, return
!      if (i == 1) then
!         get_nest_level = 1
!         return
!      end if
!
!      if (i > MAX_DOMAINS) then
!         call mprintf(.true., ERROR, &
!                      'get_nest_level() called with invalid grid ID of %i.',i1=i)
!      end if
!
!      ! If not the coarse domain, then nesting level is at least 2
!      ! Yes, this looks silly. But we do not have a grid_id array, so
!      !    we must check on parent_id
!      get_nest_level = 2
!
!      j = i
!      do while (parent_id(j) /= 1)
!         j = parent_id(j)
!         get_nest_level = get_nest_level + 1
!         
!         ! Sanity check
!         if (get_nest_level > MAX_DOMAINS) then
!            call mprintf(.true., ERROR, &
!                         'Spooky nesting setup encountered in get_nest_level().')
!         end if
!      end do
!
!   end function get_nest_level
!#endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: select_domain
   !
   ! Purpose: This routine is used to select which nest x/y <-> lat/lon 
   !   conversions will be with respect to. For example, selecting domain 2 will
   !   cause the llxy routine to compute x/y locations with respect to domain 2
   !   given a lat/lon.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine select_domain(domain_num)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: domain_num
  
!#ifdef _GEOGRID
!      if (domain_num > n_domains) then
!         call mprintf(.true.,ERROR,'In select_domain(), selected domain is greater than n_domains.')
!      end if
!#endif
!#ifdef _METGRID
      if (domain_num > 1) then
!         call mprintf(.true.,ERROR,'In select_domain(), selected domain is greater than 1.')
         call arpsstop('ERROR: In select_domain(), selected domain is greater than 1.',1)
      end if
!#endif
  
      current_nest_number = domain_num
 
   end subroutine select_domain
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: iget_selected_domain
   !
   ! Purpose: This function returns the number of the currently selected nest. 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   function iget_selected_domain()
 
      implicit none
  
      ! Return value
      integer :: iget_selected_domain
      
      iget_selected_domain = current_nest_number
 
   end function iget_selected_domain 
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: lltoxy
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine wps_lltoxy(xlat, xlon, x, y, stagger, comp_ll)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: stagger
      real, intent(in) :: xlat, xlon
      real, intent(out) :: x, y
      logical, optional, intent(in) :: comp_ll

      ! Local variables
      logical :: save_comp_ll
  
      ! Account for grid staggering
      if (stagger == HH) then
         proj_stack(current_nest_number)%stagger = HH
      else if (stagger == VV) then
         proj_stack(current_nest_number)%stagger = VV
      end if
  
      if (present(comp_ll)) then
         save_comp_ll = proj_stack(current_nest_number)%comp_ll
         proj_stack(current_nest_number)%comp_ll = comp_ll
      end if

      call latlon_to_ij(proj_stack(current_nest_number), xlat, xlon, x, y)

      if (present(comp_ll)) then
         proj_stack(current_nest_number)%comp_ll = save_comp_ll
      end if
  
      ! Account for grid staggering
      if (stagger == U) then
         x = x + 0.5
      else if (stagger == V) then
         y = y + 0.5
      end if
 
   end subroutine wps_lltoxy

   subroutine wps_lltoxy_nmm(xlat, xlon, x, y, stagger, comp_ll)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: stagger
      real, intent(in) :: xlat, xlon
      real, intent(out) :: x, y
      logical, optional, intent(in) :: comp_ll

      ! Local variables
      logical :: save_comp_ll
  
      ! Account for grid staggering
      if (stagger == HH) then
         proj_stack(current_nest_number)%stagger = HH
      else if (stagger == VV) then
         proj_stack(current_nest_number)%stagger = VV
      end if
  
      if (present(comp_ll)) then
         save_comp_ll = proj_stack(current_nest_number)%comp_ll
         proj_stack(current_nest_number)%comp_ll = comp_ll
      end if

      !call latlon_to_ij(proj_stack(current_nest_number), xlat, xlon, x, y)
      CALL llij_rotlatlon_org(xlat,xlon,proj_stack(current_nest_number),x,y)

      if (present(comp_ll)) then
         proj_stack(current_nest_number)%comp_ll = save_comp_ll
      end if
  
      ! Account for grid staggering
      if (stagger == U) then
         x = x + 0.5
      else if (stagger == V) then
         y = y + 0.5
      end if
 
   end subroutine wps_lltoxy_nmm
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   ! Name: lltoxy
   !
   ! Purpose:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   subroutine wps_xytoll(x, y, xlat, xlon, stagger, comp_ll)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: stagger
      real, intent(in) :: x, y
      real, intent(out) :: xlat, xlon
      logical, optional, intent(in) :: comp_ll

      ! Local variables
      real :: rx, ry
      logical :: save_comp_ll
  
      ! Account for grid staggering; we cannot modify x and y, so modify local
      !   copies of them
      if (stagger == U) then
         rx = x - 0.5
         ry = y
      else if (stagger == V) then
         rx = x
         ry = y - 0.5
      else if (stagger == HH) then
         proj_stack(current_nest_number)%stagger = HH
         rx = x
         ry = y
      else if (stagger == VV) then
         proj_stack(current_nest_number)%stagger = VV
         rx = x
         ry = y
      else
         rx = x
         ry = y
      end if

      if (present(comp_ll)) then
         save_comp_ll = proj_stack(current_nest_number)%comp_ll
         proj_stack(current_nest_number)%comp_ll = comp_ll
      end if
  
      call ij_to_latlon(proj_stack(current_nest_number), rx, ry, xlat, xlon)

      if (present(comp_ll)) then
         proj_stack(current_nest_number)%comp_ll = save_comp_ll
      end if
 
   end subroutine wps_xytoll

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: met_to_map                                                             !
   !                                                                              !
   ! Purpose: Rotate Earth-relative winds to grid-relative winds                  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine met_to_map(u, u_mask, v, v_mask, &
                         us1, us2, ue1, ue2, &
                         vs1, vs2, ve1, ve2, &
                         xlon_u, xlon_v, xlat_u, xlat_v)

      implicit none

      ! Arguments
      integer, intent(in) :: us1, us2, ue1, ue2, vs1, vs2, ve1, ve2
      real, dimension(us1:ue1,us2:ue2) :: u, xlon_u, xlat_u
      real, dimension(vs1:ve1,vs2:ve2) :: v, xlon_v, xlat_v
      REAL, intent(in) :: u_mask, v_mask

      INTEGER :: orig_selected_projection

      orig_selected_projection = iget_selected_domain()
      call select_domain(1)
      call metmap_xform(u, u_mask, v, v_mask, &
                        us1, us2, ue1, ue2, &
                        vs1, vs2, ve1, ve2, &
                        xlon_u, xlon_v, xlat_u, xlat_v, -1)
      call select_domain(orig_selected_projection)

   end subroutine met_to_map

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: map_to_met                                                             !
   !                                                                              !
   ! Purpose: Rotate grid-relative winds to Earth-relative winds                  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine map_to_met(u, u_mask, v, v_mask, &
                         us1, us2, ue1, ue2, &
                         vs1, vs2, ve1, ve2, &
                         xlon_u, xlon_v, xlat_u, xlat_v)

      implicit none

      ! Arguments
      integer, intent(in) :: us1, us2, ue1, ue2, vs1, vs2, ve1, ve2
!      real, pointer, dimension(:,:) :: u, v, xlon_u, xlon_v, xlat_u, xlat_v
!      type (bitarray), intent(in) :: u_mask, v_mask
      real, dimension(us1:ue1,us2:ue2) :: u, xlon_u, xlat_u
      real, dimension(vs1:ve1,vs2:ve2) :: v, xlon_v, xlat_v
      REAL, intent(in) :: u_mask, v_mask

      INTEGER :: orig_selected_projection

      orig_selected_projection = iget_selected_domain()
      call select_domain(SOURCE_PROJ)
      call metmap_xform(u, u_mask, v, v_mask, &
                        us1, us2, ue1, ue2, &
                        vs1, vs2, ve1, ve2, &
                        xlon_u, xlon_v, xlat_u, xlat_v, 1)
      call select_domain(orig_selected_projection)

   end subroutine map_to_met


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: metmap_xform                                                           !
   !                                                                              !
   ! Purpose: Do the actual work of rotating winds for C grid.                    !
   !          If idir= 1, rotate grid-relative winds to Earth-relative winds      !
   !          If idir=-1, rotate Earth-relative winds to grid-relative winds      !
   !                                                                              !
   ! ASSUMPTIONS: 1) MEMORY ORDER IS XY.                                          !
   !              2) U ARRAY HAS ONE MORE COLUMN THAN THE V ARRAY, AND V ARRAY    !
   !                 HAS ONE MORE ROW THAN U ARRAY.                               !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine metmap_xform(u, u_mask, v, v_mask, &
                           us1, us2, ue1, ue2, vs1, vs2, ve1, ve2, &
                           xlon_u, xlon_v, xlat_u, xlat_v, idir)

      implicit none

      ! Arguments
      integer, intent(in) :: us1, us2, ue1, ue2, vs1, vs2, ve1, ve2, idir
!      real, pointer, dimension(:,:) :: u, v, xlon_u, xlon_v, xlat_u, xlat_v
!      type (bitarray), intent(in) :: u_mask, v_mask
      real, dimension(us1:ue1,us2:ue2) :: u, xlon_u, xlat_u
      real, dimension(vs1:ve1,vs2:ve2) :: v, xlon_v, xlat_v
      REAL, intent(in) :: u_mask, v_mask

      ! Local variables
      integer :: i, j
      real :: u_weight, v_weight
      real :: u_map, v_map, alpha, diff
      real, pointer, dimension(:,:) :: u_new, v_new, u_mult, v_mult
      logical :: do_last_col_u, do_last_row_u, do_last_col_v, do_last_row_v

      ! If the proj_info structure has not been initialized, we don't have
      !   information about the projection and standard longitude.
      if (proj_stack(current_nest_number)%init) then

         ! Only rotate winds for Lambert conformal, polar stereographic, or Cassini
         if ((proj_stack(current_nest_number)%code == PROJ_LC) .or. &
             (proj_stack(current_nest_number)%code == PROJ_PS) .or. &
             (proj_stack(current_nest_number)%code == PROJ_CASSINI)) then
!            call mprintf((idir ==  1),LOGFILE,'Rotating map winds to earth winds.')
!            call mprintf((idir == -1),LOGFILE,'Rotating earth winds to grid winds')
!            WRITE(6,*) 'Rotating earth winds to grid winds'

            allocate(u_mult(us1:ue1,us2:ue2))
            allocate(v_mult(vs1:ve1,vs2:ve2))

            do j=us2,ue2
               do i=us1,ue1
!                  if (bitarray_test(u_mask, i-us1+1, j-us2+1)) then
                     u_mult(i,j) = 1.
!                  else
!                     u_mult(i,j) = 0.
!                  end if
               end do
            end do

            do j=vs2,ve2
               do i=vs1,ve1
!                  if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                     v_mult(i,j) = 1.
!                  else
!                     v_mult(i,j) = 0.
!                  end if
               end do
            end do

            if (ue1-us1 == ve1-vs1) then
               do_last_col_u = .false.
               do_last_col_v = .true.
            else
               do_last_col_u = .true.
               do_last_col_v = .false.
            end if

            if (ue2-us2 == ve2-vs2) then
               do_last_row_u = .true.
               do_last_row_v = .false.
            else
               do_last_row_u = .false.
               do_last_row_v = .true.
            end if

            ! Create arrays to hold rotated winds
            allocate(u_new(us1:ue1, us2:ue2))
            allocate(v_new(vs1:ve1, vs2:ve2))

            ! Rotate U field
            do j=us2,ue2
               do i=us1,ue1

                  diff = idir * (xlon_u(i,j) - proj_stack(current_nest_number)%stdlon)
                  if (diff > 180.) then
                     diff = diff - 360.
                  else if (diff < -180.) then
                     diff = diff + 360.
                  end if

                  ! Calculate the rotation angle, alpha, in radians
                  if (proj_stack(current_nest_number)%code == PROJ_LC) then
                     alpha = diff * proj_stack(current_nest_number)%cone * rad_per_deg * proj_stack(current_nest_number)%hemi 
                  else if (proj_stack(current_nest_number)%code == PROJ_CASSINI) then
                     if (j == ue2) then
                        diff = xlon_u(i,j)-xlon_u(i,j-1)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_u(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_u(i,j)-xlat_u(i,j-1))*rad_per_deg    &
                                     )
                     else if (j == us2) then
                        diff = xlon_u(i,j+1)-xlon_u(i,j)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_u(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_u(i,j+1)-xlat_u(i,j))*rad_per_deg    &
                                     )
                     else
                        diff = xlon_u(i,j+1)-xlon_u(i,j-1)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_u(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_u(i,j+1)-xlat_u(i,j-1))*rad_per_deg    &
                                     )
                     end if
                  else
                     alpha = diff * rad_per_deg * proj_stack(current_nest_number)%hemi 
                  end if
                 
                  v_weight = 0.

                  ! On C grid, take U_ij, and get V value at the same lat/lon
                  !   by averaging the four surrounding V points
!                  if (bitarray_test(u_mask, i-us1+1, j-us2+1)) then
                     u_map = u(i,j)
                     if (i == us1) then
                        if (j == ue2 .and. do_last_row_u) then
                           v_weight = v_mult(i,j)
                           v_map = v(i,j)*v_mult(i,j)
                        else
                           v_weight = v_mult(i,j) + v_mult(i,j+1)
                           v_map = v(i,j)*v_mult(i,j) + v(i,j+1)*v_mult(i,j+1)
                        end if 
                     else if (i == ue1 .and. do_last_col_u) then
                        if (j == ue2 .and. do_last_row_u) then
                           v_weight = v_mult(i-1,j)
                           v_map = v(i-1,j)
                        else
                           v_weight = v_mult(i-1,j) + v_mult(i-1,j+1) 
                           v_map = v(i-1,j)*v_mult(i-1,j) + v(i-1,j+1)*v_mult(i-1,j+1)
                        end if 
                     else if (j == ue2 .and. do_last_row_u) then
                        v_weight = v_mult(i-1,j-1) + v_mult(i,j-1)
                        v_map = v(i-1,j-1)*v_mult(i-1,j-1) + v(i,j-1)*v_mult(i,j-1)
                     else
                        v_weight = v_mult(i-1,j) + v_mult(i-1,j+1) + v_mult(i,j) + v_mult(i,j+1)
                        v_map = v(i-1,j)*v_mult(i-1,j) + v(i-1,j+1)*v_mult(i-1,j+1) + v(i,j)*v_mult(i,j) + v(i,j+1)*v_mult(i,j+1)
                     end if
                     if (v_weight > 0.) then
                        u_new(i,j) = cos(alpha)*u_map + sin(alpha)*v_map/v_weight
                     else
                        u_new(i,j) = u(i,j)
                     end if
!                  else
!                     u_new(i,j) = u(i,j)
!                  end if

               end do
            end do

            ! Rotate V field
            do j=vs2,ve2
               do i=vs1,ve1

                  diff = idir * (xlon_v(i,j) - proj_stack(current_nest_number)%stdlon)
                  if (diff > 180.) then
                     diff = diff - 360.
                  else if (diff < -180.) then
                     diff = diff + 360.
                  end if

                  if (proj_stack(current_nest_number)%code == PROJ_LC) then
                     alpha = diff * proj_stack(current_nest_number)%cone * rad_per_deg * proj_stack(current_nest_number)%hemi 
                  else if (proj_stack(current_nest_number)%code == PROJ_CASSINI) then
                     if (j == ve2) then
                        diff = xlon_v(i,j)-xlon_v(i,j-1)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_v(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_v(i,j)-xlat_v(i,j-1))*rad_per_deg    &
                                     )
                     else if (j == vs2) then
                        diff = xlon_v(i,j+1)-xlon_v(i,j)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_v(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_v(i,j+1)-xlat_v(i,j))*rad_per_deg    &
                                     )
                     else
                        diff = xlon_v(i,j+1)-xlon_v(i,j-1)
                        if (diff > 180.) then
                           diff = diff - 360.
                        else if (diff < -180.) then
                           diff = diff + 360.
                        end if
                        alpha = atan2(   (-cos(xlat_v(i,j)*rad_per_deg) * diff*rad_per_deg),   &
                                                        (xlat_v(i,j+1)-xlat_v(i,j-1))*rad_per_deg    &
                                     )
                     end if
                  else
                     alpha = diff * rad_per_deg * proj_stack(current_nest_number)%hemi 
                  end if

                  u_weight = 0.

!                  if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                     v_map = v(i,j)
                     if (j == vs2) then
                        if (i == ve1 .and. do_last_col_v) then
                           u_weight = u_mult(i,j)
                           u_map = u(i,j)*u_mult(i,j)
                        else
                           u_weight = u_mult(i,j) + u_mult(i+1,j)
                           u_map = u(i,j)*u_mult(i,j) + u(i+1,j)*u_mult(i+1,j)
                        end if 
                     else if (j == ve2 .and. do_last_row_v) then
                        if (i == ve1 .and. do_last_col_v) then
                           u_weight = u_mult(i,j-1)
                           u_map = u(i,j-1)*u_mult(i,j-1)
                        else
                           u_weight = u_mult(i,j-1) + u_mult(i+1,j-1)
                           u_map = u(i,j-1)*u_mult(i,j-1) + u(i+1,j-1)*u_mult(i+1,j-1)
                        end if 
                     else if (i == ve1 .and. do_last_col_v) then
                        u_weight = u_mult(i,j) + u_mult(i,j-1)
                        u_map = u(i,j)*u_mult(i,j) + u(i,j-1)*u_mult(i,j-1)
                     else
                        u_weight = u_mult(i,j-1) + u_mult(i,j) + u_mult(i+1,j-1) + u_mult(i+1,j)
                        u_map = u(i,j-1)*u_mult(i,j-1) + u(i,j)*u_mult(i,j) + u(i+1,j-1)*u_mult(i+1,j-1) + u(i+1,j)*u_mult(i+1,j)
                     end if
                     if (u_weight > 0.) then
                        v_new(i,j) = -sin(alpha)*u_map/u_weight + cos(alpha)*v_map
                     else
                        v_new(i,j) = v(i,j)
                     end if
!                  else
!                     v_new(i,j) = v(i,j)
!                  end if

               end do
            end do

            ! Copy rotated winds back into argument arrays
            u = u_new 
            v = v_new 

            deallocate(u_new)
            deallocate(v_new)
            deallocate(u_mult)
            deallocate(v_mult)
         end if

      else
!         call mprintf(.true.,ERROR,'In metmap_xform(), uninitialized proj_info structure.')
         WRITE(6,*) 'ERROR,In metmap_xform(), uninitialized proj_info structure.'
      end if
 
   end subroutine metmap_xform


   !
   ! NMM Wind Rotation Code
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: map_to_met_nmm                                                         !
   !                                                                              !
   ! Purpose: Rotate grid-relative winds to Earth-relative winds                  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine map_to_met_nmm(u, u_mask, v, v_mask, &
                             vs1, vs2, ve1, ve2, &
                             xlat_v, xlon_v)

      implicit none

      ! Arguments
      integer, intent(in) :: vs1, vs2, ve1, ve2
!      real, pointer, dimension(:,:) :: u, v, xlat_v, xlon_v
!      type (bitarray), intent(in) :: u_mask, v_mask
      real, dimension(vs1:ve1,vs2:ve2) :: u, v, xlat_v, xlon_v
      REAL, intent(in) :: u_mask, v_mask

      INTEGER :: orig_selected_projection

      orig_selected_projection = iget_selected_domain()
      call select_domain(1)
      call metmap_xform_nmm(u, u_mask, v, v_mask, &
                            vs1, vs2, ve1, ve2, &
                            xlat_v, xlon_v, 1)
      call select_domain(orig_selected_projection)

   end subroutine map_to_met_nmm
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Name: met_to_map_nmm                                                         !
   !                                                                              !
   ! Purpose: Rotate Earth-relative winds to grid-relative winds                  !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine met_to_map_nmm(u, u_mask, v, v_mask, &
                             vs1, vs2, ve1, ve2, &
                             xlat_v, xlon_v)

      implicit none

      ! Arguments
      integer, intent(in) :: vs1, vs2, ve1, ve2
!      real, pointer, dimension(:,:) :: u, v, xlat_v, xlon_v
!      type (bitarray), intent(in) :: u_mask, v_mask
      real, dimension(vs1:ve1,vs2:ve2) :: u, v, xlat_v, xlon_v
      REAL, intent(in) :: u_mask, v_mask

      INTEGER :: orig_selected_projection

      orig_selected_projection = iget_selected_domain()
      call select_domain(1)
      call metmap_xform_nmm(u, u_mask, v, v_mask, &
                            vs1, vs2, ve1, ve2, &
                            xlat_v, xlon_v, -1)
      call select_domain(orig_selected_projection)

   end subroutine met_to_map_nmm


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: metmap_xform_nmm                                                       !
   !                                                                              !
   ! Purpose: Do the actual work of rotating winds for E grid.                    !
   !          If idir= 1, rotate grid-relative winds to Earth-relative winds      !
   !          If idir=-1, rotate Earth-relative winds to grid-relative winds      !
   !                                                                              !
   ! ASSUMPTIONS: 1) MEMORY ORDER IS XY.                                          !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine metmap_xform_nmm(u, u_mask, v, v_mask, &
                               vs1, vs2, ve1, ve2, &
                               xlat_v, xlon_v, idir)

      implicit none

      ! Arguments
      integer, intent(in) :: vs1, vs2, ve1, ve2, idir
!      real, pointer, dimension(:,:) :: u, v, xlat_v, xlon_v
!      type (bitarray), intent(in) :: u_mask, v_mask
      real, dimension(vs1:ve1,vs2:ve2) :: u, v, xlat_v, xlon_v
      REAL, intent(in) :: u_mask, v_mask

      ! Local variables
      integer :: i, j
      real :: u_map, v_map, diff, alpha
      real :: phi0, lmbd0, big_denominator, relm, rlat_v,rlon_v, clontemp
      real :: sin_phi0, cos_phi0,  cos_alpha, sin_alpha
      real, pointer, dimension(:,:) :: u_new, v_new


      ! If the proj_info structure has not been initialized, we don't have
      !   information about the projection and standard longitude.
      if (proj_stack(current_nest_number)%init) then

         if (proj_stack(current_nest_number)%code == PROJ_ROTLL) then

!            call mprintf((idir ==  1),LOGFILE,'Rotating map winds to earth winds.')
!            call mprintf((idir == -1),LOGFILE,'Rotating earth winds to grid winds')
!            WRITE(6,*) 'Rotating earth winds to grid winds for PROJ_ROTIL.'
   
            ! Create arrays to hold rotated winds
            allocate(u_new(vs1:ve1, vs2:ve2))
            allocate(v_new(vs1:ve1, vs2:ve2))

            phi0  = proj_stack(current_nest_number)%lat1 * rad_per_deg

            clontemp= proj_stack(current_nest_number)%lon1

            if (clontemp .lt. 0.) then
               lmbd0 = (clontemp + 360.) * rad_per_deg
            else
               lmbd0 = (clontemp) * rad_per_deg
            endif

            sin_phi0  = sin(phi0)
            cos_phi0  = cos(phi0)

            do j=vs2,ve2
               do i=vs1,ve1

                  ! Calculate the sine and cosine of rotation angle
                  rlat_v = xlat_v(i,j) * rad_per_deg
                  rlon_v = xlon_v(i,j) * rad_per_deg
                  relm = rlon_v - lmbd0
                  big_denominator = cos(asin( &
                                       cos_phi0 * sin(rlat_v) - &
                                       sin_phi0 * cos(rlat_v) * cos(relm) &
                                        )   )

                  sin_alpha = sin_phi0 * sin(relm)  /  &
                                         big_denominator

                  cos_alpha = (cos_phi0 * cos(rlat_v) + &
                               sin_phi0 * sin(rlat_v) * cos(relm))  /  &
                                            big_denominator
   
                  ! Rotate U field
!                  if (bitarray_test(u_mask, i-vs1+1, j-vs2+1)) then
                     u_map = u(i,j)
!                     if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                        v_map = v(i,j)
!                     else
!                        v_map = 0.
!                     end if
                     
                     u_new(i,j) = cos_alpha*u_map + idir*sin_alpha*v_map
!                  else
!                     u_new(i,j) = u(i,j)
!                  end if
                       
                  ! Rotate V field
!                  if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                     v_map = v(i,j)
!                     if (bitarray_test(u_mask, i-vs1+1, j-vs2+1)) then
                        u_map = u(i,j)
!                     else
!                        u_map = 0.
!                     end if
                     
                     v_new(i,j) = -idir*sin_alpha*u_map + cos_alpha*v_map
!                  else
!                     v_new(i,j) = v(i,j)
!                  end if
   
               end do
            end do

            ! Copy rotated winds back into argument arrays
            u = u_new 
            v = v_new 
   
            deallocate(u_new)
            deallocate(v_new)
   
         ! Only rotate winds for Lambert conformal, polar stereographic, or Cassini
         else if ((proj_stack(current_nest_number)%code == PROJ_LC) .or. &
                  (proj_stack(current_nest_number)%code == PROJ_PS) .or. &
                  (proj_stack(current_nest_number)%code == PROJ_CASSINI)) then

!            call mprintf((idir ==  1),LOGFILE,'Rotating map winds to earth winds.')
!            call mprintf((idir == -1),LOGFILE,'Rotating earth winds to grid winds')
!            WRITE(6,*) 'Rotating earth winds to grid winds for PROJ_LC/PS/CASSINI'

            ! Create arrays to hold rotated winds
            allocate(u_new(vs1:ve1, vs2:ve2))
            allocate(v_new(vs1:ve1, vs2:ve2))

            do j=vs2,ve2
               do i=vs1,ve1

                  diff = idir * (xlon_v(i,j) - proj_stack(current_nest_number)%stdlon)
                  if (diff > 180.) then
                     diff = diff - 360.
                  else if (diff < -180.) then
                     diff = diff + 360.
                  end if

                  ! Calculate the rotation angle, alpha, in radians
                  if (proj_stack(current_nest_number)%code == PROJ_LC) then
                     alpha = diff * proj_stack(current_nest_number)%cone * &
                             rad_per_deg * proj_stack(current_nest_number)%hemi
                  else
                     alpha = diff * rad_per_deg * proj_stack(current_nest_number)%hemi 
                  end if

                  ! Rotate U field
!                  if (bitarray_test(u_mask, i-vs1+1, j-vs2+1)) then
                     u_map = u(i,j)
!                     if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                        v_map = v(i,j)
!                     else
!                        v_map = 0.
!                     end if
                     
                     u_new(i,j) = cos(alpha)*u_map + idir*sin(alpha)*v_map
!                  else
!                     u_new(i,j) = u(i,j)
!                  end if
                       
                  ! Rotate V field
!                  if (bitarray_test(v_mask, i-vs1+1, j-vs2+1)) then
                     v_map = v(i,j)
!                     if (bitarray_test(u_mask, i-vs1+1, j-vs2+1)) then
                        u_map = u(i,j)
!                     else
!                        u_map = 0.
!                     end if
                     
                     v_new(i,j) = -idir*sin(alpha)*u_map + cos(alpha)*v_map
!                  else
!                     v_new(i,j) = v(i,j)
!                  end if

               end do
            end do

            ! Copy rotated winds back into argument arrays
            u = u_new 
            v = v_new 
   
            deallocate(u_new)
            deallocate(v_new)

         end if

      else
!         call mprintf(.true.,ERROR,'In metmap_xform_nmm(), uninitialized proj_info structure.')
         WRITE(6,*) 'ERROR,In metmap_xform_nmm(), uninitialized proj_info structure.'
      end if

   end subroutine metmap_xform_nmm

end module wrf_llxy_module
