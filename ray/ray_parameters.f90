module ray_parameters
  use amr_parameters
 
  ! Number for fuzzy comparison
  real(dp) :: ray_epsilon = 1.0d-13

  ! <><><><><><><><><><><><><><><><
  ! Set parameters for ray geometry
  ! <><><><><><><><><><><><><><><><

  ! Number of rays to follow in each direction
  integer  :: ray_Nx = 256 
  integer  :: ray_Ny = 256

  ! Opening of the field of view in x & y directions, in degrees
  real(dp) :: ray_opening_x = 5.0D0 
  real(dp) :: ray_opening_y = 5.0D0

  ! Angles of the centre of field of view, in degress
  real(dp) :: ray_cof_x = 0.0D0
  real(dp) :: ray_cof_y = 0.0D0

  ! Coordinate of the observer (in cartesian coordinates in units of box size)
  ! By defaul it is the center of the z = 0 face
  real(dp) :: ray_x_obs = 0.5D0
  real(dp) :: ray_y_obs = 0.5D0
  real(dp) :: ray_z_obs = 0.0D0 

  ! Box size in Mpc/h
  real(dp) :: ray_Lbox = 250.0D0                                   ! This is the same as boxlen_ini but with double precision

  ! source redshift
  real(dp) :: ray_z_s = 1.0D0

  ! Tests:
  logical  :: ray_step_test = .false.
 
  ! <><><><><><><><><><><><><><><><><><><
  ! Set which quantities to be calculated
  ! <><><><><><><><><><><><><><><><><><><

  logical  :: ray_do_kappa = .true.
  logical  :: ray_do_shear = .false.
  logical  :: ray_do_iswrs = .false.                               ! Baojiu-Feb-2019 
  logical  :: ray_do_deflt = .false.                               ! Baojiu-Feb-2019 
  ! +++ add more variables here to calculate other observables +++
  ! +++ also update amr/read_params.f90 and namelist/cosmo.nml +++

  ! <><><><><><><><><><><><><><><><><>
  ! General parameters for methodology
  ! <><><><><><><><><><><><><><><><><>

  ! Method to compute kappa (1: method A, 2: method B, 3: method A & B)
  integer  :: ray_kappa_method = 3

  ! Set to integrate fields at cell centers
  logical  :: ray_do_ngp = .false.

  ! Flag indicating whether bending of light rays is required
  logical  :: ray_no_bending = .true.

  ! Flag indicating whether want whole simulation to finish after all rays are finished 
  logical  :: ray_stop_after_raytracing = .true.

  ! <><><><><><><><><><><><>
  ! Ray output-related stuff                                       
  ! <><><><><><><><><><><><>

  logical           :: ray_multi_out = .false.                     
  logical           :: ray_afresh    = .true. 
  integer           :: ray_nout      = 0                                          

  integer ,parameter               :: RAY_MAXOUT = 1000                        
  real(dp),dimension(1:RAY_MAXOUT) :: ray_zout   = -9.9D0            

end module ray_parameters
