!===============================
! ray_write_cell
!===============================
subroutine ray_write_fields()
  use amr_commons
  use poisson_commons
  use ray_commons
  use pm_commons
  implicit none

  ! For testing porpoises
  ! Arguments:
  !            iray: to identify grid
  !            ind:  to identify cell
  !            istep: number of step

  real (dp), dimension(1:RAY_NFIELDS, 1:8) :: ray_fields

  integer :: ngrid
  integer :: i_interp, ilevel, ind, igrid, icell, iskip, igrid_mg
  integer :: ix, iy, iz
  real(dp), dimension(1:3) :: xc
  real(dp), dimension(1:3) :: delta_target
  integer :: itarget
  real(dp), dimension(1:3) :: x_cell
  real(dp) :: dx
  real(dp), dimension(1:RAY_NFIELDS) :: mean_fields

  integer :: i, j

  character(LEN=5)::nchar
  character*1000 file_name

  integer :: cell_corner

  real (dp), dimension(1:6) :: tidal

  integer, save :: first = 1


  ! Check that this is the first call.  If not, do nothing.
  if(first==1) then
     first = 0
  else
     return
  endif
  
  call title(myid,nchar)
  file_name = 'ray_fields_'//trim(nchar)//'.dat'
  open(unit=1001, file=file_name, form='formatted')

  write(*,*) "ray_write_fields:  aexp = ", aexp

  do ilevel=levelmin, nlevelmax
     call ray_analytic_rho(ilevel)
     call force_fine(ilevel,1)
  enddo

  do i_interp = 1, 2
     
     do ilevel=levelmin, nlevelmax

        !call ray_analytic_rho(ilevel)

        ngrid=active(ilevel)%ngrid
        dx  = 0.5d0**ilevel

        if(myid==1) write(*,*) "output_ray: dx = ", ilevel, dx

        ! Loop over cells
        do ind=1,twotondim
           iskip = ncoarse+(ind-1)*ngridmax

           ! Loop over active grids
           do igrid_mg=1,ngrid
              igrid = active(ilevel)%igrid(igrid_mg)
              icell = iskip + igrid

              ! Write only the highest refinement
              if(son(icell).ne.0) cycle  

              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              xc(1)=(dble(ix)-0.5D0)*dx
              xc(2)=(dble(iy)-0.5D0)*dx
              xc(3)=(dble(iz)-0.5D0)*dx

              x_cell(1) = xg(igrid,1) + xc(1)
              x_cell(2) = xg(igrid,2) + xc(2)
              x_cell(3) = xg(igrid,3) + xc(3)

              if( ( x_cell(1).gt.0.0D0 .and. x_cell(1).lt.(0.5D0+dx) ) .or. &
                   ( x_cell(2).gt.0.0D0 .and. x_cell(2).lt.(0.5D0+dx) ) .or. &
                   ( x_cell(3).gt.0.0D0 .and. x_cell(3).lt.(0.5D0+dx) ) ) then

                 if(i_interp==1) then
                    call ray_fields_in_corners(ilevel, igrid, ind, ray_fields)
                 else if(i_interp==2) then
                    call ray_interpolate(ilevel, igrid, ind, ray_fields)
                 else
                    !call ray_interpolate2(ilevel, igrid, ind, ray_fields)
                 endif
                 
                 call tidal_tensor(ilevel, igrid, icell, ind, tidal)

                 do i=1,RAY_NFIELDS
                    mean_fields(i) = 0.0D0
                    do j=1,8
                       mean_fields(i) = mean_fields(i) + ray_fields(i,j)
                    enddo
                    mean_fields(i) = mean_fields(i) / 8.0D0
                 enddo
              endif


              ! This is the whole grid
              write(1001,*) &
                   1, &
                   x_cell(1), &
                   x_cell(2), &
                   x_cell(3), &
                   mean_fields(1), &
                   mean_fields(2), &
                   mean_fields(3), &
                   mean_fields(4), &
                   mean_fields(5), &
                   mean_fields(6), &
                   mean_fields(7), &
                   mean_fields(8), &
                   mean_fields(9), &
                   mean_fields(10), &
                   mean_fields(11), &
                   rho(icell), &
                   tidal(1), &
                   tidal(2), &
                   tidal(3), &
                   tidal(4), &
                   tidal(5), &
                   tidal(6), &
                   f(icell,1), &
                   f(icell,2), &
                   f(icell,3), &
                   2, &
                   ilevel, &
                   i_interp

              ! This are the values in the centre of the nodes.
              cell_corner = 0
              if( x_cell(1).gt.0.5D0 .and. x_cell(1).lt.(0.5D0+dx) ) then
                 write(1001,*) &
                      1, &
                      x_cell(1), &
                      x_cell(2), &
                      x_cell(3), &
                      mean_fields(1), &
                      mean_fields(2), &
                      mean_fields(3), &
                      mean_fields(4), &
                      mean_fields(5), &
                      mean_fields(6), &
                      mean_fields(7), &
                      mean_fields(8), &
                      mean_fields(9), &
                      mean_fields(10), &
                      mean_fields(11), &
                      rho(icell), &
                      tidal(1), &
                      tidal(2), &
                      tidal(3), &
                      tidal(4), &
                      tidal(5), &
                      tidal(6), &
                      f(icell,1), &
                      f(icell,2), &
                      f(icell,3), &
                      cell_corner, &
                      ilevel, &
                      i_interp
              endif
              if( x_cell(2).gt.0.5D0 .and. x_cell(2).lt.(0.5D0+dx) ) then
                 write(1001,*) &
                      2, &
                      x_cell(1), &
                      x_cell(2), &
                      x_cell(3), &
                      mean_fields(1), &
                      mean_fields(2), &
                      mean_fields(3), &
                      mean_fields(4), &
                      mean_fields(5), &
                      mean_fields(6), &
                      mean_fields(7), &
                      mean_fields(8), &
                      mean_fields(9), &
                      mean_fields(10), &
                      mean_fields(11), &
                      rho(icell), &
                      tidal(1), &
                      tidal(2), &
                      tidal(3), &
                      tidal(4), &
                      tidal(5), &
                      tidal(6), &
                      f(icell,1), &
                      f(icell,2), &
                      f(icell,3), &
                      cell_corner, &
                      ilevel, &
                      i_interp
              endif
              if( x_cell(3).gt.0.5D0 .and. x_cell(3).lt.(0.5D0+dx) ) then
                 write(1001,*) &
                      3, &
                      x_cell(1), &
                      x_cell(2), &
                      x_cell(3), &
                      mean_fields(1), &
                      mean_fields(2), &
                      mean_fields(3), &
                      mean_fields(4), &
                      mean_fields(5), &
                      mean_fields(6), &
                      mean_fields(7), &
                      mean_fields(8), &
                      mean_fields(9), &
                      mean_fields(10), &
                      mean_fields(11), &
                      rho(icell), &
                      tidal(1), &
                      tidal(2), &
                      tidal(3), &
                      tidal(4), &
                      tidal(5), &
                      tidal(6), &
                      f(icell,1), &
                      f(icell,2), &
                      f(icell,3), &
                      cell_corner, &
                      ilevel, &
                      i_interp

              endif

              ! These are the values in the corners.
              cell_corner = 1
              do itarget = 1, 8

                 iz=(itarget-1)/4
                 iy=(itarget-1-4*iz)/2
                 ix=(itarget-1-2*iy-4*iz)
                 delta_target(1)=(dble(ix)-0.5D0)*dx
                 delta_target(2)=(dble(iy)-0.5D0)*dx
                 delta_target(3)=(dble(iz)-0.5D0)*dx

                 x_cell(1) = xg(igrid,1) + xc(1) + delta_target(1)
                 x_cell(2) = xg(igrid,2) + xc(2) + delta_target(2)
                 x_cell(3) = xg(igrid,3) + xc(3) + delta_target(3)

                 if( x_cell(1).gt.(0.5D0-dx/2.0D0) .and. x_cell(1).lt.(0.5D0+dx/2.0D0) ) then
                    write(1001,*) &
                         1, &
                         x_cell(1), &
                         x_cell(2), &
                         x_cell(3), &
                         ray_fields(1,itarget), &
                         ray_fields(2,itarget), &
                         ray_fields(3,itarget), &
                         ray_fields(4,itarget), &
                         ray_fields(5,itarget), &
                         ray_fields(6,itarget), &
                         ray_fields(7,itarget), &
                         ray_fields(8,itarget), &
                         ray_fields(9,itarget), &
                         ray_fields(10,itarget), &
                         ray_fields(11,itarget), &
                         rho(icell), &
                         tidal(1), &
                         tidal(2), &
                         tidal(3), &
                         tidal(4), &
                         tidal(5), &
                         tidal(6), &
                         f(icell,1), &
                         f(icell,2), &
                         f(icell,3), &
                         cell_corner, &
                         ilevel, &
                         i_interp
                 endif
                 if( x_cell(2).gt.(0.5D0-dx/2.0D0) .and. x_cell(2).lt.(0.5D0+dx/2.0D0) ) then
                    write(1001,*) &
                         2, &
                         x_cell(1), &
                         x_cell(2), &
                         x_cell(3), &
                         ray_fields(1,itarget), &
                         ray_fields(2,itarget), &
                         ray_fields(3,itarget), &
                         ray_fields(4,itarget), &
                         ray_fields(5,itarget), &
                         ray_fields(6,itarget), &
                         ray_fields(7,itarget), &
                         ray_fields(8,itarget), &
                         ray_fields(9,itarget), &
                         ray_fields(10,itarget), &
                         ray_fields(11,itarget), &
                         rho(icell), &
                         tidal(1), &
                         tidal(2), &
                         tidal(3), &
                         tidal(4), &
                         tidal(5), &
                         tidal(6), &
                         f(icell,1), &
                         f(icell,2), &
                         f(icell,3), &
                         cell_corner, &
                         ilevel, &
                         i_interp
                 endif
                 if( x_cell(3).gt.(0.5D0-dx/2.0D0) .and. x_cell(3).lt.(0.5D0+dx/2.0D0) ) then
                    write(1001,*) &
                         3, &
                         x_cell(1), &
                         x_cell(2), &
                         x_cell(3), &
                         ray_fields(1,itarget), &
                         ray_fields(2,itarget), &
                         ray_fields(3,itarget), &
                         ray_fields(4,itarget), &
                         ray_fields(5,itarget), &
                         ray_fields(6,itarget), &
                         ray_fields(7,itarget), &
                         ray_fields(8,itarget), &
                         ray_fields(9,itarget), &
                         ray_fields(10,itarget), &
                         ray_fields(11,itarget), &
                         rho(icell), &
                         tidal(1), &
                         tidal(2), &
                         tidal(3), &
                         tidal(4), &
                         tidal(5), &
                         tidal(6), &
                         f(icell,1), &
                         f(icell,2), &
                         f(icell,3), &
                         cell_corner, &
                         ilevel, &
                         i_interp
                 endif

              enddo  ! itarget

           end do
        end do

     enddo  ! ilevel
  enddo  ! i_interp

  close(1001)


end subroutine ray_write_fields
