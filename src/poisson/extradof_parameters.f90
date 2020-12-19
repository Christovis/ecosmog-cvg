module extradof_parameters
  use amr_parameters

  ! free cubic vector Galileon parameter
  real(dp) :: param_b3=0.000001

  ! DGP model parameters  Chr 30/03/20
  real(dp) :: param_o=0.25D0    ! Omega_rc=1/(4*H0^2*rc^2) where rc is DGP scale
  real(dp) :: param_l=1.0d0      ! Omega_Lambda (only used if nDGP is true and nDGP2 is false)
  logical  :: nDGP=.true.        ! set to true to do nDGP, otherwise do sDGP
  logical  :: nDGP2=.false.       ! set to true to do nDGP with LCDM background, false for nDGP with a cosmological constant

  ! speed of light [meter/sec]
  real(dp) :: sol=299792458.0d0

  integer  :: ngs_fine_extradof_pre   = 30  ! pre-smoothing factor
  integer  :: ngs_fine_extradof_pst   = 3   ! post-smoothing factor
  integer  :: ngs_coarse_extradof_pre = 30  ! pre-smoothing factor
  integer  :: ngs_coarse_extradof_pst = 3   ! post-smoothing factor

end module extradof_parameters
