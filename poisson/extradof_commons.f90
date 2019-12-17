module extradof_commons 
  use amr_parameters

  ! sf : longitudinal vector field mode, \chi
  ! sf_grad : gradiant of sf
  ! cbf : field of the cubic transverse mode, B_i
  real(dp)                            :: sf_src_mean
  real(dp),allocatable,dimension(:)   :: sf
  real(dp),allocatable,dimension(:)   :: sf_src
  real(dp),allocatable,dimension(:,:) :: sf_grad
  real(dp),allocatable,dimension(:,:) :: cbf
  ! sf_lp : laplacian of longitudinal mode, \nabla^2\chi
  ! The following terms are the source-terms for B_i
  ! bf_src1 : nabla_i\nabla_j\chi\nabdla_j\nabla^2\chi
  ! bf_src2 : nabla_i\chi\nabdla^2\nabla^2\chi
  ! bf_src3 : nabla^2\chi\nabdla_i\nabla^2\chi
  ! bf_src4 : nabla_j\chi\nabla_i\nabdla_j\nabla^2\chi
  real(dp),allocatable,dimension(:)   :: sf_lp     ! gradient of ps
  !real(dp),allocatable,dimension(:)   :: bf_src1   ! gradient of ps
  !real(dp),allocatable,dimension(:)   :: bf_src2   ! gradient of ps
  !real(dp),allocatable,dimension(:)   :: bf_src3   ! gradient of ps
  !real(dp),allocatable,dimension(:)   :: bf_src4   ! gradient of ps

  ! These are time-dependent background functions
  real(dp) :: alpha_cvg  ! alpha function, needed for 5th force
  real(dp) :: beta_cvg   ! beta function, needed for 5th force
  real(dp) :: rc_cvg     ! cross-over function, needed for Vainshtein radius

end module extradof_commons
