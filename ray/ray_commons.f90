module ray_commons 
  use amr_commons
  use ray_parameters

  ! <><><><><><><><><><><><><><><><><
  ! Ray arrays and related quantities
  ! <><><><><><><><><><><><><><><><><
  real(dp), allocatable, dimension(:,:) :: ray_coord               ! Ray_coordinates (chi,theta,phi)
  real(dp), allocatable, dimension(:,:) :: ray_quants              ! ray tracing quantities (Baojiu-Feb-2019)

  ! Number of ray quantities 
  integer :: RAY_NQUANTS                                           ! number of ray quantities to calculate (Baojiu-Feb-2019)

  ! Number of physical fields
  integer :: RAY_NFIELDS                                           ! number of physical fields to be used  (Baojiu-Feb-2019)

  ! Indices of ray quantities 
  integer :: i_ray_kappa                                           ! Baojiu-Feb-2019
  integer :: i_ray_shear                                           ! Baojiu-Feb-2019
  integer :: i_ray_deflt                                           ! Baojiu-Feb-2019
  integer :: i_ray_iswrs                                           ! Baojiu-Feb-2019
  ! +++ add new stuff needed to calculate other observables +++

  ! Indices of physical fields
  integer :: i_ray_tidal                                           ! Baojiu-Feb-2019
  integer :: i_ray_grads                                           ! Baojiu-Feb-2019
  integer :: i_ray_phidt                                           ! Baojiu-Feb-2019
  ! +++ add new stuff needed to calculate other observables +++

  ! Arrays related to ray identity
  integer                               :: ray_nrays               ! number of rays on my CPU
  integer,  allocatable, dimension(:  ) :: ray_id                  ! ray ID
  integer,  allocatable, dimension(:,:) :: ray_grid                ! (grid ind, info of cell_ind,icpu & ilevel)
  integer,  allocatable, dimension(:  ) :: ex_father               ! father cell before kill_grid is called

  ! Rays to store corner variables
  integer                                 :: ray_ncells            ! number of cells containing rays
  integer,  allocatable, dimension(:    ) :: ray_in_cell           ! number of rays in a given cell 
  real(dp), allocatable, dimension(:,:,:) :: ray_stored_fields     ! store corner force values (Baojiu-Feb-2019)

  ! <><><><><><><><><><><><><><><><>
  ! Communication structure for rays
  ! <><><><><><><><><><><><><><><><>
  type ray_communicator
    integer                             :: ngrid                   ! number of rays to be communicated
    integer ,dimension(:),pointer       :: igrid                   ! local index of the grid
    integer ,dimension(:),pointer       :: igrid2      
    integer ,dimension(:),pointer       :: ivar                    ! container for integer  variables to be communicated
    real(dp),dimension(:),pointer       :: fvar                    ! container for floating variables to be communicated
  end type ray_communicator

  ! Ray emission and reception communicators
  type(ray_communicator),allocatable,dimension(:) :: ray_emission
  type(ray_communicator),allocatable,dimension(:) :: ray_reception

  ! To store the value of aexp at the beggining
  ! of a particle time step (used to determine the 
  ! maximum distance a ray can travel in a time step)
  real(dp),allocatable,dimension(:) :: aexp_old_ray,aexp_new_ray   ! Sownak

  ! <><><><><><><><><><>
  ! Searching facilities
  ! <><><><><><><><><><>
  ! Arrays for identification of neighbouring nodes
  integer, dimension(1:27,    1:8) :: ray_kkk,ray_lll
  integer, dimension(1:8, 1:8,1:8) :: ray_kkk_mmm,ray_lll_mmm
  integer, dimension(1:8,     1:8) :: ray_ccc
  real(dp),dimension(1:8         ) :: ray_bbbb
  integer, dimension(1:3, 1:2,1:8) :: ray_iii,ray_jjj
  integer, dimension(1:8,     1:8) :: ray_ooo
  
  ! <><><><><><><><><><><
  ! Light cone quantities
  ! <><><><><><><><><><><
  real(dp) :: ray_chi_s                                            ! source comoving distance

  ! Ray integration status flags
  logical :: ray_all_started                                       ! flag all rays have started  integration
  logical :: ray_all_finished                                      ! flag all rays have finished integration
  logical :: ray_end_of_sim                                        ! flag that we are ready to end the simulation

  ! Box ID in the tiling
  logical :: ray_first_box

  ! Flag to be used to decide ray status
  logical :: ray_initialized = .false.

  ! Ray info output flags 
  integer :: ray_iout

end module ray_commons
