module extradof_commons 
  use amr_parameters
  
  ! cvG fields
  ! sf : longitudinal/chi-mode of proca field
  ! sf_src_mean : mean of momentarily solved 5th-force field, meaning sf or cbf
  !               (subroutine in multigrid_fine_commons_extradof.f90)
  ! sf_src_mean2 : the same as sf_src_mean but it is not overwrittin within a
  !                time-step
  ! sf_grad : gradiant of sf
  ! sf_lp : laplacian of longitudinal mode, \nabla^2\chi
  ! cbf : transverse/B-mode of proca field
  !       it interchangeably contains the different components of B_i
  !       bf_src1 : nabla_i\nabla_j\chi\nabdla_j\nabla^2\chi
  !       bf_src2 : nabla_i\chi\nabdla^2\nabla^2\chi
  !       bf_src3 : nabla^2\chi\nabdla_i\nabla^2\chi
  !       bf_src4 : nabla_j\chi\nabla_i\nabdla_j\nabla^2\chi
  real(dp),allocatable,dimension(:)   :: sf
  real(dp)                            :: sf_src_mean
  real(dp),dimension(1:4)             :: sf_src_mean2
  real(dp),allocatable,dimension(:)   :: sf_src
  real(dp),allocatable,dimension(:,:) :: sf_grad
  real(dp),allocatable,dimension(:)   :: sf_lp
  real(dp),allocatable,dimension(:,:) :: cbf

  ! cvG, These are time-dependent background functions
  ! (are defined in ./amr/adaptive_loop.f90)
  real(dp) :: alpha_cvg  ! alpha function, needed for 5th force
  real(dp) :: beta_cvg   ! beta function, needed for 5th force
  real(dp) :: rc_cvg     ! cross-over function Rc^2, needed for Vainshtein radius

  ! DGP, background related functions
  real(dp) :: beta_dgp                             ! beta(a), Chr 30/03/20
  real(dp) :: gamma_dgp                            ! Chr 30/03/20

end module extradof_commons
