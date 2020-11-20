!===============================
! ray_write_cell
!===============================
subroutine ray_analytic_rho(ilevel)
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

  real (dp), dimension(1:100, 1:8) :: ray_fields

  integer :: ngrid
  integer :: ilevel, ind, igrid, icell, iskip, igrid_mg
  integer :: ix, iy, iz
  real(dp), dimension(1:3) :: xc
  real(dp), dimension(1:3) :: x_cell
  real(dp) :: dx
  real(dp), dimension(1:100) :: mean_fields

  integer :: i, j

  logical :: found_boundary
  integer :: found

  real(dp) :: r2
  real(dp) :: b2, c2  ! Semiaxis of the thing

  if(myid==1) write(*,*) "Entering ray_analytic_rho..."

  ngrid=active(ilevel)%ngrid
  dx  = 0.5d0**ilevel

  b2 = 1.0  !0.6**2
  c2 = 1.0  !0.7**2

  ! Loop over cells
  do ind=1,twotondim
     iskip = ncoarse+(ind-1)*ngridmax

     ! Loop over active grids
     do igrid_mg=1,ngrid
        igrid = active(ilevel)%igrid(igrid_mg)
        icell = iskip + igrid

        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(1)=(dble(ix)-0.5D0)*dx
        xc(2)=(dble(iy)-0.5D0)*dx
        xc(3)=(dble(iz)-0.5D0)*dx

        x_cell(1) = xg(igrid,1) + xc(1)
        x_cell(2) = xg(igrid,2) + xc(2)
        x_cell(3) = xg(igrid,3) + xc(3)

        r2 = (x_cell(1)-0.5D0)**2 + (x_cell(2)-0.5D0)**2/b2 + (x_cell(3)-0.5D0)**2/c2

        ! This is in the code presentation paper:
        rho(icell) = 10000.0D0 * exp(-r2/((1.5D0/ray_Lbox)**2));
        rho(icell) = 10000.0D0 * exp(-r2/((0.7D0/ray_Lbox)**2));
        !rho(icell) = 10000.0D0 * exp(-r2/((0.1D0/ray_Lbox)**2));
        !rho(icell) = 0.0001*x_cell(1)
        ! This is for ISW:
        

        !rho(icell) = 10000.0D0 * (x_cell(3)-0.5D0)**3*(x_cell(1)-0.5D0)**2
        !rho(icell) = 10000.0D0 * (x_cell(3)-0.5D0)**2*(x_cell(1)-0.5D0)**2
        
        phi(icell) = -rho(icell)

        !phi(icell) = 100.0D0*x_cell(2)**3+10000.0D0
        
     end do
  end do

  ! Communicate
  call make_virtual_fine_dp(rho(1),ilevel)
  call make_virtual_fine_dp(phi(1),ilevel)
  
end subroutine ray_analytic_rho
