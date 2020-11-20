!##########################################################
!##########################################################
!##########################################################
subroutine get_ind_from_ray_position(iray,ilevel,igrid,ind)
  use amr_commons
  use poisson_commons
  use ray_commons
  implicit none

  integer,intent(in)  :: iray,igrid,ilevel 
  integer,intent(out) :: ind 
  !--------------------------------------------------------------
  ! This subroutine finds on which of the eight cells of the grid 
  ! the ray is in
  !--------------------------------------------------------------

  integer  :: ix,iy,iz,i,imax
  real(dp) :: dx,dmax,x_g,y_g,z_g,x_r,y_r,z_r,x_c,y_c,z_c
  real(dp),dimension(1:3) :: ray_dir

  ! Cartesian coordinates of the grid centre wrt the observer
  x_g = xg(igrid,1)-ray_x_obs
  y_g = xg(igrid,2)-ray_y_obs
  z_g = xg(igrid,3)-ray_z_obs

  x_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
  y_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
  z_r = ray_coord(iray,1)*dcos(ray_coord(iray,2))

  ! Cell size
  dx  = 0.5D0**ilevel

  ! Deal with the case with a ray sitting at cell interfaces
  if(x_r.eq.x_g .or. y_r.eq.y_g .or. z_r.eq.z_g) then

     ! Directions of the rays
     if(ray_no_bending) then
        ! Ray moving towards the observer
        ray_dir(1) = dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
        ray_dir(2) = dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
        ray_dir(3) = dcos(ray_coord(iray,2))
     else
        ! Ray not moving towards the obsever
     end if

     if(x_r.eq.x_g) x_r=x_r-dsign(0.5D0,ray_dir(1))*dx
     if(y_r.eq.y_g) y_r=y_r-dsign(0.5D0,ray_dir(2))*dx
     if(z_r.eq.z_g) z_r=z_r-dsign(0.5D0,ray_dir(3))*dx

  end if

  ! Cartesian coordinates of the grid lower-left corner wrt the observer
  x_c = x_g-dx
  y_c = y_g-dx
  z_c = z_g-dx

  ! Find local position of the ray-containing cell inside the grid
  ix  = int((x_r-x_c)/dx)
  iy  = int((y_r-y_c)/dx)
  iz  = int((z_r-z_c)/dx)
  ! Fix cases where ix,iy,iz = 2, set it to 1
  if(ix.eq.2) ix = 1
  if(iy.eq.2) iy = 1
  if(iz.eq.2) iz = 1

  ind = 4*iz+2*iy+ix+1

  ! Sanity check
  if(ind.lt.1 .or. ind.gt.8) then
     write(*,*) 'get_ind_from_ray_position: ind exceeds bound; please check.', ix,iy,iz, ind, dx
     write(*,*) 'get_ind_from_ray_position: ind exceeds bound; please check.', x_r, y_r, z_r, x_c, y_c, z_c
     stop
  end if

end subroutine get_ind_from_ray_position

!##########################################################
!##########################################################
!##########################################################

subroutine get_cell_center(ilevel,iray,ind,delta_xc)
  use amr_commons
  use poisson_commons
  use ray_commons
  implicit none
  integer,                intent(in ) :: ilevel,iray,ind
  real(dp),dimension(1:3),intent(out) :: delta_xc 
  !-----------------------------------------------------------------------
  ! Get the Cartesian coordinates of the center of the ray_containing cell
  ! relative to grid centre
  !-----------------------------------------------------------------------

  integer :: ix, iy, iz
  real(dp) :: dx
  
  dx  = 0.5D0**ilevel
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  delta_xc(1)=(dble(ix)-0.5D0)*dx
  delta_xc(2)=(dble(iy)-0.5D0)*dx
  delta_xc(3)=(dble(iz)-0.5D0)*dx

end subroutine get_cell_center

!##########################################################
!##########################################################
!##########################################################

subroutine find_direction(iray,ilevel,ind,idim,inbor)
  use pm_commons
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  implicit none

  integer,intent(in ) :: iray,ilevel,ind
  integer,intent(out) :: idim,inbor

  !------------------------------------------------------------------
  ! Determines the face of a cubic cell via which the ray will leave.
  ! idim: direction x,y,z 
  ! inbor: positive/negative side of direction
  !   E.g., idim=1 & inbor=1 means leaving via face with smaller 
  !         x-coordinate in the x-direction (left-going ray)
  !------------------------------------------------------------------

  integer  :: imax,i
  real(dp) :: dx,dmax,tmp,dray
  real(dp) :: x_target
  logical  :: cross_edge
  integer,  dimension(1:3) :: edge_faces,xx
  real(dp), dimension(1:3) :: xo        ! Cartesian coordinates of observer
  real(dp), dimension(1:3) :: xgo       ! Cartesian coordinates of grid centre wrt observer
  real(dp), dimension(1:3) :: xco       ! Cartesian coordinates of cell centre wrt observer
  real(dp), dimension(1:3) :: delta_xc  ! Cartesian coordinates of cell centre wrt grid centre
  real(dp), dimension(1:3) :: ray_dir   ! direction of ray
  real(dp), dimension(1:3) :: ray_pos   ! Cartesian coordinates of ray wrt the observer

  ! Numerical factors
  dx = 0.5D0**ilevel                             ! cell size at level #ilevel

  ! Rates at which the ray goes in each direction (note the -ve signs)
  if(ray_no_bending) then
     ! Ray moving towards the observer
     ray_dir(1) = -dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
     ray_dir(2) = -dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
     ray_dir(3) = -dcos(ray_coord(iray,2))
  else
     ! Ray not moving towards the observer
  end if

  ! Cartesian coordinates of the rays wrt the observer
  ray_pos(1) = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
  ray_pos(2) = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
  ray_pos(3) = ray_coord(iray,1)*dcos(ray_coord(iray,2))

  ! Cartesian coordinates of the observer
  xo(1) = ray_x_obs
  xo(2) = ray_y_obs
  xo(3) = ray_z_obs

  ! Find the direction of fastest progagation
  dmax = 0.0D0                                   ! distance along the fastest direction if ray travels by 1.0
  imax = 0                                       ! ray travels fastest in the imax-th direction
  do i=1,ndim
     if(dabs(ray_dir(i)).gt.dmax) then
        dmax = dabs(ray_dir(i))
        imax = i
     end if
  end do

  ! Cartesian coordinates of the grid centre wrt the observer
  do i=1,ndim
     xgo(i) = xg(ray_grid(iray,1),i)-xo(i)
  end do

  ! Cartesian coordinates of the cell centre wrt the grid centre
  call get_cell_center(ilevel,iray,ind,delta_xc)

  ! Cartesian coordinates of the cell centre wrt the observer   
  do i=1,ndim
     xco(i) = xgo(i)+delta_xc(i)
  end do

  ! Target cell surface coordinate in the direction of fastest propagation
  ! Note: abs(ray_dir(imax)) is always larger than 0.0D0
  x_target = xco(imax)+dsign(0.5D0,ray_dir(imax))*dx  

  ! Ray distance to said target cell surface
  dray = dabs(ray_pos(imax)-x_target)/dmax

  idim = imax                                    ! assumed direction of ray crossing
  cross_edge = .false.                           ! assume ray does not cross an edge
  edge_faces = 0                                 ! assume ray does not cross an edge

  do i=1,ndim
     if(i.eq.imax) cycle                         ! this direction has been considered - cycle
     if(dabs(ray_dir(i))+1.0D0.eq.1.0D0) cycle   ! ray doesn't move in this direction - cycle

     tmp = xco(i)+dsign(0.5D0,ray_dir(i))*dx     ! target cell surface coordinte in i-th direction
     tmp = dabs(ray_pos(i)-tmp)                  ! ray distance to said target cell surface

     ! ray crosses an edge/vertex of the cell 
     if(dray*dabs(ray_dir(i))+1.0D0.eq.tmp+1.0D0) then
        cross_edge       = .true.
        edge_faces(i   ) = 1
        edge_faces(idim) = 1
     end if

     ! ray will cross a face in i-th direction earlier!
     if(dray*dabs(ray_dir(i))+1.0D0.gt.tmp+1.0D0) then 
        dray = tmp/dabs(ray_dir(i))              ! update dray 
        idim = i                                 ! update idim
     end if
  end do

  ! Find which neighbouring cell the rays is going into
  if(ray_dir(idim).lt.0.0D0) inbor = 1           ! ray travels in -x/-y/-z direction
  if(ray_dir(idim).gt.0.0D0) inbor = 2           ! ray travels in +x/+y/+z direction
 
  ! Update ray coordinates
  if(ray_no_bending) then
     ! Ray moving towards the observer: only need to update ray_coord(iray,1)
     ray_coord(iray,1) = ray_coord(iray,1)-dray
  else
     ! Ray not moving towards the observer: need to update ray_coord(iray,1:3)
  end if

  ! Deal with the cases of ray ending in current cell
  if(ray_no_bending) then
     ! Ray moving towards the observer
     if(ray_coord(iray,1).le.ray_epsilon) then
        ray_coord(iray,1) = 0.0D0       
        inbor = 0                                ! flag that ray ends in current cell (end of ray)
        return
     end if
  else
     ! Ray not moving towards the observer
  end if
  
  ! Deal with the cases of ray ending in current cell
  if(ray_no_bending) then
     ! Ray moving towards the observer
!    do i=1,ndim
!       if(ray_pos(i)+xo(i)+ray_dir(i)*dray+1.0D0.le.1.0D0) then
        if(ray_coord(iray,1)*dcos(ray_coord(iray,2))-dabs(ray_z_obs).lt.ray_epsilon) then
           inbor = 0                                ! flag that ray ends in current cell (end of box)
           ray_coord(iray,1) = dabs(ray_z_obs)/dabs(ray_dir(3))
           return
        end if
!    end do
  else
     ! Ray not moving towards the observer
  end if

  if(cross_edge) then                            ! ray has crossed an edge
     if(edge_faces(idim).eq.0) return            ! ray crosses another face before crossing the edge!

     ! flag that ray goes through an edge for warning!
     xx    = -1
     inbor = 14                                  ! the index of current cell in its 27 neighbours
     do i=1,ndim
        if(edge_faces(i).eq.0)  cycle
        if(ray_dir(i).ge.0.0D0) xx(i) = 1
        inbor = inbor+3**(i-1)*xx(i)             ! the index of the cell which the ray is going into
     end do
     inbor = -inbor                              ! flag as negative to distinguish
  end if

end subroutine find_direction

