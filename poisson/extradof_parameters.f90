module extradof_parameters
  use amr_parameters

  ! cubic vector Galileon parameters
  ! TODO tune beta to be similar to DGP (look at python routine)
  real(dp) :: param_b3=

  ! speed of light in standard units
  real(dp) :: sol=299792458.0d0

  integer  :: ngs_fine_extradof_pre   = 30  ! pre-smoothing factor
  integer  :: ngs_fine_extradof_pst   = 3   ! post-smoothing factor
  integer  :: ngs_coarse_extradof_pre = 30  ! pre-smoothing factor
  integer  :: ngs_coarse_extradof_pst = 3   ! post-smoothing factor

end module extradof_parameters
