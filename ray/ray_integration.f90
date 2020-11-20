! <><><><><><><><><><><><><><><><><><><><><><><><><><><
! This is the main subroutine that does ray integration
! <><><><><><><><><><><><><><><><><><><><><><><><><><><
subroutine ray_integrate(iray,ilevel,ind,chi_A,ray_fields)
  !------------------------------------------------------------------
  ! Integrates [g(chi,chi_s)*quantity] along the ray path in a cell 
  ! 
  ! g(chi, chi_S) = chi(chi_S-chi)/chi_S is the lensing kernel 
  ! quantity is any quantity we want to integrate (density, nabla^2Phi,
  !      nabla_inabla_i Phi, etc.
  ! 
  ! The notation used here tries to follow that of the paper
  ! Point A --> point in cell where the integration ends
  ! Point B --> point in cell where the integration starts
  ! Point 1 --> vertice 1 of cell, ie, that with lowest x,y,z coords
  ! 
  ! a_cell  --> diff. in x-dir between Point A and Point 1
  ! b_cell  --> diff. in y-dir between Point A and Point 1
  ! c_cell  --> diff. in z-dir between Point A and Point 1
  ! 
  ! alpha_* --> coefficients of the trilinear interpolation of int_qua 
  ! d_1234  --> d_N coefficients in the integral
  !
  ! Equations are references to paper: https://arxiv.org/pdf/1601.02012.pdf
  !
  !------------------------------------------------------------------
  use pm_commons
  use poisson_commons
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  implicit none

  integer, intent(in) :: iray,ilevel,ind
  real(dp),intent(in) :: chi_A
! integer, intent(in) :: icell
  real (dp), dimension(1:RAY_NFIELDS,1:8),intent(in) :: ray_fields
 
  real(dp) :: dx
  integer  :: i
  real(dp) :: a_cell,b_cell,c_cell                                 ! Difference in x,y,z of point A and vertice 1 of cell 
  real(dp) :: R_cell                                               ! chi_B - chi_A, ending point of the integral (recall paper notation)
  real(dp),dimension(1:3) :: xo                                    ! Cartesian coordinates of the observer
  real(dp),dimension(1:3) :: chiA_pos                              ! Cartesian coordinates of point A wrt observer
  real(dp),dimension(1:3) :: xgo                                   ! Cartesian coordinates of grid centre wrt observer
  real(dp),dimension(1:3) :: xco                                   ! Cartesian coordinates of cell centre wrt observer
  real(dp),dimension(1:3) :: delta_xc                              ! Cartesian coordinates of cell centre wrt grid centre
  real(dp),dimension(1:4) :: d_1234_corr                           ! coefficients d_N (N=1..4) in the integral formula of the correction

  real(dp) :: sinthe,costhe,sinphi,cosphi,sin2phi,sin2the,cos2phi
  real(dp) :: Dx_A,Dy_A,Dz_A,Dx_B,Dy_B,Dz_B
  real(dp) :: nablachiphi_at_A,nablachiphi_at_B,g_at_A,g_at_B,correction_term
  real(dp) :: integral_cell_corr                                   ! integration result correction
  real(dp) :: aexp_mean

  ! Integrand (which will be a linear combination of ray_fields)
  real(dp),dimension(1:8) :: nablachiphi                           ! correction force term in Method A
  real(dp),dimension(1:8) :: alpha_nablachiphi                     ! trilinear interpol. coeff.

  real(dp),dimension(1:RAY_NQUANTS,1:8) :: alpha                   ! trilinear interpolation coeff (Baojiu-Feb-2019; united to one array)
  real(dp),dimension(1:RAY_NQUANTS,1:8) :: int_qua                 ! quantities to be integrated   (Baojiu-Feb-2019; united to one array)
  real(dp),dimension(1:RAY_NQUANTS,1:4) :: d_1234                  ! coefficients d_N (N=1..4) in the integral formula (Baojiu-Feb-2019; united to one array) 
  real(dp),dimension(1:RAY_NQUANTS    ) :: integral_cell           ! integration result (Baojiu-Feb-2019; united to one array)
  
  ! Useful trigonometric quantities
  sinthe  = dsin(      ray_coord(iray,2))
  costhe  = dcos(      ray_coord(iray,2))
  sinphi  = dsin(      ray_coord(iray,3))
  cosphi  = dcos(      ray_coord(iray,3))
  sin2the = dsin(2.0D0*ray_coord(iray,2))
  sin2phi = dsin(2.0D0*ray_coord(iray,3))
  cos2phi = dcos(2.0D0*ray_coord(iray,3))

  ! Initialise int_qua
  int_qua = 0.0D0                                                  ! Baojiu-Feb-2019

  ! Define quantity to be integrated:
  if(ray_do_kappa) then 
     if(ray_kappa_method.ne.2) then 
        ! Loop over corners of this cell.
        do i=1,8
           ! Quantity is \delta/a Eq.(43)
           int_qua(i_ray_kappa,i) = (ray_fields(1,i)-1.0D0)/aexp

           ! Compute the corrections to the density integral for Method A
           ! nablachiphi is \nabla_\chi\Phi; see text above Eq.(44)
           nablachiphi(i) = -(sinthe*cosphi)*ray_fields(i_ray_grads  ,i) - &                          ! x term
                             (sinthe*sinphi)*ray_fields(i_ray_grads+1,i) - &                          ! y term
                             (costhe       )*ray_fields(i_ray_grads+2,i)                              ! z term
        end do
     end if
     if(ray_kappa_method.ne.1) then
        ! Loop over corners of this cell.
        do i=1,8
           ! Quantity is \nabla^2_\xi\Phi 
           ! Eq.(39)+Eq.(40)
           int_qua(i_ray_kappa+1,i) = (sinphi**2+costhe**2*cosphi**2)*ray_fields(i_ray_tidal  ,i) + & ! xx term
                                      (cosphi**2+costhe**2*sinphi**2)*ray_fields(i_ray_tidal+1,i) + & ! yy term
                                      (sinthe**2                    )*ray_fields(i_ray_tidal+2,i) + & ! zz term
                                      (-sinthe**2*sin2phi           )*ray_fields(i_ray_tidal+3,i) + & ! xy term
                                      (-cosphi*sin2the              )*ray_fields(i_ray_tidal+4,i) + & ! xz term
                                      (-sinphi*sin2the              )*ray_fields(i_ray_tidal+5,i)     ! yz term
        enddo
     end if
  end if

  if(ray_do_shear) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nabla^1\nabla_1_\Phi-\nabla^2\nabla_2_\Phi
        ! Eq.(39)-Eq.(40)
        int_qua(i_ray_shear  ,i) = (sinphi**2-costhe**2*cosphi**2)*ray_fields(i_ray_tidal  ,i) + & !xx term
                                   (cosphi**2-costhe**2*sinphi**2)*ray_fields(i_ray_tidal+1,i) - & !yy term
                                   (sinthe**2                    )*ray_fields(i_ray_tidal+2,i) - & !zz term
                                   (sin2phi*(costhe**2+1.0D0)    )*ray_fields(i_ray_tidal+3,i) + & !xy term
                                   (cosphi*sin2the               )*ray_fields(i_ray_tidal+4,i) + & !xz term
                                   (sinphi*sin2the               )*ray_fields(i_ray_tidal+5,i)     !yz term
        ! Quantity is 2*\nabla^1\nabla_2_\Phi
        ! Eq.(41)
        int_qua(i_ray_shear+1,i) = 2.0*(                                                         &
                                   (-costhe*cosphi*sinphi        )*ray_fields(i_ray_tidal  ,i) + & !xx term
                                   ( costhe*cosphi*sinphi        )*ray_fields(i_ray_tidal+1,i) + & !yy term
!                                  ( 0.0D0                       )*ray_fields(i_ray_tidal+2,i) + & !zz term
                                   ( costhe*(cosphi**2-sinphi**2))*ray_fields(i_ray_tidal+3,i) + & !xy term
                                   ( sinphi*sinthe               )*ray_fields(i_ray_tidal+4,i) + & !xz term
                                   (-cosphi*sinthe               )*ray_fields(i_ray_tidal+5,i))    !yz term
     end do
  end if

  if(ray_do_deflt) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is \nable_1\Phi
        int_qua(i_ray_deflt  ,i) = (costhe*cosphi)*ray_fields(i_ray_grads  ,i) + &       !x term
                                   (costhe*sinphi)*ray_fields(i_ray_grads+1,i) + &       !y term
                                   (sinthe)       *ray_fields(i_ray_grads+2,i)           !z term
        ! Quantity is \nable_2\Phi
        int_qua(i_ray_deflt+1,i) = (sinthe*cosphi)*ray_fields(i_ray_grads+1,i) - &       !y term
                                   (sinthe*sinphi)*ray_fields(i_ray_grads  ,i)           !x term
     end do

  end if

  if(ray_do_iswrs) then
     ! Loop over corners of this cell.
     do i=1,8
        ! Quantity is phidot
        int_qua(i_ray_iswrs  ,i) = 2.0D0*ray_fields(i_ray_phidt+1,i)
     end do
  end if

  ! +++ add new stuff here to calculate other observables +++

  ! Ending point of the integral: chi_B - chi_A (recall paper variable)
  R_cell = ray_coord(iray,1)-chi_A

  ! cell size of ilevel
  dx = 0.5D0**ilevel

  ! Cartesian coordinates of Point A wrt the observer
  chiA_pos(1) = chi_A*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
  chiA_pos(2) = chi_A*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
  chiA_pos(3) = chi_A*dcos(ray_coord(iray,2))

  ! Cartesian coordinates of the observer
  xo(1) = ray_x_obs
  xo(2) = ray_y_obs
  xo(3) = ray_z_obs

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

  ! Differences in x,y,z direction between Point A and Point 1  
  a_cell = chiA_pos(1)-(xco(1)-0.5D0*dx)                           ! diff in x-dir
  b_cell = chiA_pos(2)-(xco(2)-0.5D0*dx)                           ! diff in y-dir
  c_cell = chiA_pos(3)-(xco(3)-0.5D0*dx)                           ! diff in z-dir

  ! alpha coefficients of the trilinear interpolation; Eq.(15)
  alpha(1:RAY_NQUANTS,1) = int_qua(1:RAY_NQUANTS,1)
  alpha(1:RAY_NQUANTS,2) = int_qua(1:RAY_NQUANTS,2)-int_qua(1:RAY_NQUANTS,1)
  alpha(1:RAY_NQUANTS,3) = int_qua(1:RAY_NQUANTS,3)-int_qua(1:RAY_NQUANTS,1)
  alpha(1:RAY_NQUANTS,4) = int_qua(1:RAY_NQUANTS,5)-int_qua(1:RAY_NQUANTS,1)
  alpha(1:RAY_NQUANTS,5) = int_qua(1:RAY_NQUANTS,4)-int_qua(1:RAY_NQUANTS,3)-int_qua(1:RAY_NQUANTS,2)+int_qua(1:RAY_NQUANTS,1) 
  alpha(1:RAY_NQUANTS,6) = int_qua(1:RAY_NQUANTS,7)-int_qua(1:RAY_NQUANTS,5)-int_qua(1:RAY_NQUANTS,3)+int_qua(1:RAY_NQUANTS,1) 
  alpha(1:RAY_NQUANTS,7) = int_qua(1:RAY_NQUANTS,6)-int_qua(1:RAY_NQUANTS,5)-int_qua(1:RAY_NQUANTS,2)+int_qua(1:RAY_NQUANTS,1) 
  alpha(1:RAY_NQUANTS,8) = int_qua(1:RAY_NQUANTS,8)-int_qua(1:RAY_NQUANTS,7)-int_qua(1:RAY_NQUANTS,6)-int_qua(1:RAY_NQUANTS,4) + & 
                           int_qua(1:RAY_NQUANTS,2)+int_qua(1:RAY_NQUANTS,5)+int_qua(1:RAY_NQUANTS,3)-int_qua(1:RAY_NQUANTS,1)

  ! --------------------------------------------------------------------------- !
  ! alexandre_block_alexandre_block_alexandre_block  ==
  ! This is to compute the corrections to the density integral for Method A 
  alpha_nablachiphi(1) = nablachiphi(1)
  alpha_nablachiphi(2) = nablachiphi(2)-nablachiphi(1)
  alpha_nablachiphi(3) = nablachiphi(3)-nablachiphi(1)
  alpha_nablachiphi(4) = nablachiphi(5)-nablachiphi(1)
  alpha_nablachiphi(5) = nablachiphi(4)-nablachiphi(3) - &
                         nablachiphi(2)+nablachiphi(1) 
  alpha_nablachiphi(6) = nablachiphi(7)-nablachiphi(5) - &
                         nablachiphi(3)+nablachiphi(1) 
  alpha_nablachiphi(7) = nablachiphi(6)-nablachiphi(5) - &
                         nablachiphi(2)+nablachiphi(1) 
  alpha_nablachiphi(8) = nablachiphi(8)-nablachiphi(7) - &
                         nablachiphi(6)-nablachiphi(4) + & 
                         nablachiphi(2)+nablachiphi(5) + &
                         nablachiphi(3)-nablachiphi(1)
  ! alexandre_block_alexandre_block_alexandre_block
  ! --------------------------------------------------------------------------- !

  ! d_N coefficients in the integral formula; Eq.(B1)-(B4)
  d_1234(1:RAY_NQUANTS,1) =                                    alpha(1:RAY_NQUANTS,1)                               &
                            + (                                alpha(1:RAY_NQUANTS,2)*a_cell                      + &
                                                               alpha(1:RAY_NQUANTS,3)       *b_cell               + &
                                                               alpha(1:RAY_NQUANTS,4)              *c_cell)/dx      &
                            + (                                alpha(1:RAY_NQUANTS,5)*a_cell*b_cell               + &
                                                               alpha(1:RAY_NQUANTS,6)       *b_cell*c_cell        + &
                                                               alpha(1:RAY_NQUANTS,7)*a_cell       *c_cell)/dx**2   &
                            + (                                alpha(1:RAY_NQUANTS,8)*a_cell*b_cell*c_cell)/dx**3
  d_1234(1:RAY_NQUANTS,2) =   (sinthe*cosphi*                  alpha(1:RAY_NQUANTS,2)                             + & 
                               sinthe*sinphi*                  alpha(1:RAY_NQUANTS,3)                             + &
                               costhe       *                  alpha(1:RAY_NQUANTS,4)                     )/dx      & 
                            + (sinthe*cosphi*(                 alpha(1:RAY_NQUANTS,7)              *c_cell        + &
                                                               alpha(1:RAY_NQUANTS,5)       *b_cell)              + &
                               sinthe*sinphi*(                 alpha(1:RAY_NQUANTS,6)              *c_cell        + &
                                                               alpha(1:RAY_NQUANTS,5)*a_cell)                     + &
                               costhe       *(                 alpha(1:RAY_NQUANTS,7)*a_cell                      + &
                                                               alpha(1:RAY_NQUANTS,6)       *b_cell)      )/dx**2   &
                            + (sinthe*cosphi                  *alpha(1:RAY_NQUANTS,8)       *b_cell*c_cell        + &
                               sinthe*sinphi                  *alpha(1:RAY_NQUANTS,8)*a_cell       *c_cell        + &
                               costhe                         *alpha(1:RAY_NQUANTS,8)*a_cell*b_cell       )/dx**3
  d_1234(1:RAY_NQUANTS,3) =   (sinthe**2*cosphi*sinphi        *alpha(1:RAY_NQUANTS,5)                             + &
                               sinthe   *cosphi*costhe        *alpha(1:RAY_NQUANTS,7)                             + &
                               sinthe   *sinphi*costhe        *alpha(1:RAY_NQUANTS,6)                     )/dx**2   & 
                            + (sinthe**2*cosphi*sinphi        *alpha(1:RAY_NQUANTS,8)               *c_cell       + &
                               sinthe   *cosphi*costhe        *alpha(1:RAY_NQUANTS,8)       *b_cell               + &
                               sinthe   *sinphi*costhe        *alpha(1:RAY_NQUANTS,8)*a_cell              )/dx**3
  d_1234(1:RAY_NQUANTS,4) =   (sinthe**2*costhe*cosphi*sinphi)*alpha(1:RAY_NQUANTS,8              )/dx**3 
  
  ! Compute the integral in the cell; Eq.(17)
  integral_cell = 0.0D0
  do i=1,4
     ! WL convergence
     if(ray_do_kappa) then
        integral_cell(i_ray_kappa:i_ray_kappa+1) = integral_cell(i_ray_kappa:i_ray_kappa+1)                       + &
                                                          d_1234(i_ray_kappa:i_ray_kappa+1,i)*(                     &
                                                      R_cell**(i+2)                              /(dble(i)+2.0D0) + & 
                                                      R_cell**(i+1)*(2.0D0*chi_A-ray_chi_s)      /(dble(i)+1.0D0) + &
                                                      R_cell**(i  )*(      chi_A-ray_chi_s)*chi_A/(dble(i)      ))/ &
                                                   ray_chi_s
     end if
     ! WL shear
     if(ray_do_shear) then
        integral_cell(i_ray_shear:i_ray_shear+1) = integral_cell(i_ray_shear:i_ray_shear+1)                       + &
                                                          d_1234(i_ray_shear:i_ray_shear+1,i)*(                     &
                                                      R_cell**(i+2)                              /(dble(i)+2.0D0) + & 
                                                      R_cell**(i+1)*(2.0D0*chi_A-ray_chi_s)      /(dble(i)+1.0D0) + &
                                                      R_cell**(i  )*(      chi_A-ray_chi_s)*chi_A/(dble(i)      ))/ &
                                                   ray_chi_s
     end if
     ! WL deflection angles
     if(ray_do_deflt) then
        integral_cell(i_ray_deflt:i_ray_deflt+1) = integral_cell(i_ray_deflt:i_ray_deflt+1)                       + &
                                                          d_1234(i_ray_deflt:i_ray_deflt+1,i)*(                     &
                                                      R_cell**(i+2)                              /(dble(i)+2.0D0) + & 
                                                      R_cell**(i+1)*(2.0D0*chi_A-ray_chi_s)      /(dble(i)+1.0D0) + &
                                                      R_cell**(i  )*(      chi_A-ray_chi_s)*chi_A/(dble(i)      ))/ &
                                                   ray_chi_s
     end if
     ! ISW-RS effect
     if(ray_do_iswrs) then
        integral_cell(i_ray_iswrs)               = integral_cell(i_ray_iswrs)                                     - &
                                                          d_1234(i_ray_iswrs,i)*R_cell**i/(dble(i)) 
     end if 
    ! +++ add new stuff here to calculate other observables +++
  end do

  if(ray_do_kappa.and.ray_kappa_method.ne.2) then
     ! This is to compute the corrections to the density integral for Method A
     ! This modification block computes the term integral([g'\nabla_\chi\Phi]dchi)
     ! d_N coefficients in the integral formula
     d_1234_corr(1) =                                    alpha_nablachiphi(1)                               &
                      + (                                alpha_nablachiphi(2)*a_cell                      + &
                                                         alpha_nablachiphi(3)       *b_cell               + &
                                                         alpha_nablachiphi(4)              *c_cell)/dx      &
                      + (                                alpha_nablachiphi(5)*a_cell*b_cell               + &
                                                         alpha_nablachiphi(6)       *b_cell*c_cell        + &
                                                         alpha_nablachiphi(7)*a_cell       *c_cell)/dx**2   &
                      +                                  alpha_nablachiphi(8)*a_cell*b_cell*c_cell /dx**3
     d_1234_corr(2) =   (sinthe*cosphi                  *alpha_nablachiphi(2)                             + &
                         sinthe*sinphi                  *alpha_nablachiphi(3)                             + &
                         costhe                         *alpha_nablachiphi(4)                     )/dx      & 
                      + (sinthe*cosphi*(                 alpha_nablachiphi(7)              *c_cell        + &
                                                         alpha_nablachiphi(5)       *b_cell       )       + &
                         sinthe*sinphi*(                 alpha_nablachiphi(6)              *c_cell        + &
                                                         alpha_nablachiphi(5)*a_cell              )       + &
                         costhe*       (                 alpha_nablachiphi(7)*a_cell                      + &
                                                         alpha_nablachiphi(6)       *b_cell      ))/dx**2   &
                      + (sinthe*cosphi                  *alpha_nablachiphi(8)       *b_cell*c_cell        + &
                         sinthe*sinphi                  *alpha_nablachiphi(8)*a_cell       *c_cell        + &
                         costhe                         *alpha_nablachiphi(8)*a_cell*b_cell       )/dx**3
     d_1234_corr(3) =   (sinthe**2*cosphi*sinphi        *alpha_nablachiphi(5)                             + &
                         sinthe   *cosphi*costhe        *alpha_nablachiphi(7)                             + &
                         sinthe   *sinphi*costhe        *alpha_nablachiphi(6)                     )/dx**2   & 
                      + (sinthe**2*cosphi*sinphi        *alpha_nablachiphi(8)              *c_cell        + &
                         sinthe   *cosphi*costhe        *alpha_nablachiphi(8)       *b_cell               + &
                         sinthe*   sinphi*costhe        *alpha_nablachiphi(8)*a_cell              )/dx**3
     d_1234_corr(4) =   (sinthe**2*costhe*cosphi*sinphi)*alpha_nablachiphi(8)/dx**3 
  
     ! Compute the integral in the cell
     integral_cell_corr = 0.0D0
     correction_term    = 0.0D0 
     do i=1,4
        integral_cell_corr = integral_cell_corr                                                           + &
                             d_1234_corr(i)*(R_cell**(i+1)* 2.0D0                 /(dble(i)+1.0D0)        + &
                                             R_cell**(i  )*(2.0D0*chi_A-ray_chi_s)/(dble(i)      ))
     end do
     integral_cell_corr = integral_cell_corr/ray_chi_s
  
     ! This is to compute the corrections to the density integral for Method A
     !  This modification block computes the term integral([g\nabla_\chi\Phi]dchi)
     ! Compute the displacements Dx, Dy, Dz for point A
     Dx_A = a_cell/dx
     Dy_A = b_cell/dx
     Dz_A = c_cell/dx
     ! Compute the displacements Dx, Dy, Dz for point B
     Dx_B = (a_cell+R_cell*sinthe*cosphi)/dx
     Dy_B = (b_cell+R_cell*sinthe*sinphi)/dx
     Dz_B = (c_cell+R_cell*costhe       )/dx
     ! Compute \nabla_\chi\Phi at point A
     nablachiphi_at_A = alpha_nablachiphi(1)                + &
                        alpha_nablachiphi(2)*Dx_A           + &
                        alpha_nablachiphi(3)*Dy_A           + &
                        alpha_nablachiphi(4)*Dz_A           + &
                        alpha_nablachiphi(5)*Dx_A*Dy_A      + &
                        alpha_nablachiphi(6)*Dy_A*Dz_A      + &
                        alpha_nablachiphi(7)*Dx_A*Dz_A      + &
                        alpha_nablachiphi(8)*Dx_A*Dy_A*Dz_A
     ! Compute \nabla_\chi\Phi at point B
     nablachiphi_at_B = alpha_nablachiphi(1)                + & 
                        alpha_nablachiphi(2)*Dx_B           + &
                        alpha_nablachiphi(3)*Dy_B           + &
                        alpha_nablachiphi(4)*Dz_B           + &
                        alpha_nablachiphi(5)*Dx_B*Dy_B      + &
                        alpha_nablachiphi(6)*Dy_B*Dz_B      + &
                        alpha_nablachiphi(7)*Dx_B*Dz_B      + &
                        alpha_nablachiphi(8)*Dx_B*Dy_B*Dz_B
     ! Compute f(chi,chi_s) at point A and point B
     g_at_A = (ray_chi_s-chi_A            )*chi_A            /ray_chi_s
     g_at_B = (ray_chi_s-ray_coord(iray,1))*ray_coord(iray,1)/ray_chi_s

     correction_term = g_at_B*nablachiphi_at_B-g_at_A*nablachiphi_at_A
  end if

  if(ray_do_kappa.and.ray_kappa_method.ne.2) then
     ! These ray_kappa are being computed in units of (km/s)^2 
     ! density integral + both corrections
     ray_quants(iray,i_ray_kappa  ) = ray_quants(iray,i_ray_kappa)                                        + &
                        (1.5D0*omega_m)*integral_cell(i_ray_kappa)      *ray_Lbox**2*100.0D0**2           + &
                                    (correction_term+integral_cell_corr)*ray_Lbox**2*100.0D0**2/aexp**2
  end if
  if(ray_do_kappa.and.ray_kappa_method.ne.1) then
     ! This should have units of (km/s)^2
     ray_quants(iray,i_ray_kappa+1) = ray_quants(iray,i_ray_kappa+1)                                      + &
                                        integral_cell(i_ray_kappa+1)    *ray_Lbox**2*100.0D0**2/aexp**2
  end if
  if(ray_do_shear) then
     !This should have units of (km/s)^2
     ray_quants(iray,i_ray_shear:i_ray_shear+1) =                                                           &
                              ray_quants(iray,i_ray_shear:i_ray_shear+1)                                  + &
                                integral_cell(i_ray_shear:i_ray_shear+1)*ray_Lbox**2*100.0D0**2/aexp**2
  end if
  if(ray_do_deflt) then
     !This should have units of (km/s)^2
     ray_quants(iray,i_ray_deflt:i_ray_deflt+1) =                                                           &
                              ray_quants(iray,i_ray_deflt:i_ray_deflt+1)                                  + &
                                integral_cell(i_ray_deflt:i_ray_deflt+1)*ray_Lbox**2*100.0D0**2/aexp**2
  end if
  if(ray_do_iswrs) then
     aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
     ray_quants(iray,i_ray_iswrs)        = ray_quants(iray,i_ray_iswrs)                                   + &
                                             integral_cell(i_ray_iswrs) *ray_Lbox**3*100.0D0**3*aexp_mean
  end if
  ! +++ add new stuff here to calculate other observables +++

end subroutine ray_integrate
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 
subroutine ray_integrate_ngp(iray,ilevel,ind,chi_A,ray_fields_c)
  use pm_commons
  use poisson_commons
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters

  implicit none

  integer, intent(in) :: iray,ilevel,ind
  real(dp),intent(in) :: chi_A
  real(dp),dimension(1:RAY_NFIELDS),intent(in) :: ray_fields_c
  !------------------------------------------------------------------
  ! Integrates [g(chi, chi_S) * quantity] along the ray path in a cell, but
  ! assuming quantity is constant inside the cell, so in practice in integrates:
  ! g(chi, chi_S) and multiplies by the quantity
  ! 
  ! g(chi, chi_S) = chi(chi_S-chi)/chi_S is the lensing kernel 
  ! quantity is any quantity we want (density, nabla^2Phi,
  !      nabla_inabla_i Phi, etc.
  ! 
  ! The notation used here tries to follow that of the paper
  ! Point A --> point in cell where the integration ends
  ! Point B --> point in cell where the integration starts
  ! 
  ! Note that for ISW, the kernel is a constant
  ! 
  !------------------------------------------------------------------
 
  real(dp)                :: dx
  integer                 :: i
  real(dp)                :: chi_B                                 ! chi_B
  real(dp),dimension(1:3) :: xo                                    ! Cartesian coordinates of the observer
  real(dp),dimension(1:3) :: chiA_pos                              ! Cartesian coordinates of point A     wrt observer
  real(dp),dimension(1:3) :: xgo                                   ! Cartesian coordinates of grid centre wrt observer
  real(dp),dimension(1:3) :: xco                                   ! Cartesian coordinates of cell centre wrt observer
  real(dp),dimension(1:3) :: delta_xc                              ! Cartesian coordinates of cell centre wrt grid centre
  real(dp) :: sinthe,costhe,sinphi,cosphi,sin2phi,sin2the,cos2phi
  real(dp) :: integral_cell_corr
  real(dp) :: aexp_mean

  real(dp),dimension(1:RAY_NQUANTS) :: int_qua                     ! quantities to be integrated (Baojiu-Feb-2019; united to one array)
  real(dp),dimension(1:RAY_NQUANTS) :: integral_cell
  
  ! Useful trigonometric quantities
  sinthe  = dsin(      ray_coord(iray,2))
  costhe  = dcos(      ray_coord(iray,2))
  sinphi  = dsin(      ray_coord(iray,3))
  cosphi  = dcos(      ray_coord(iray,3))
  sin2the = dsin(2.0D0*ray_coord(iray,2))
  sin2phi = dsin(2.0D0*ray_coord(iray,3))
  cos2phi = dcos(2.0D0*ray_coord(iray,3))

  ! Initialise arrays
  int_qua       = 0.0D0
  integral_cell = 0.0D0
  
  ! Ending point of the integral: chi_B (paper variable)
  chi_B = ray_coord(iray,1)

  ! Define quantity to be integrated:
  if(ray_do_kappa) then 
     if(ray_kappa_method.ne.2) then 
        int_qua(i_ray_kappa  ) =                                  (ray_fields_c(1)-1.0D0)/aexp
     end if
     if(ray_kappa_method.ne.1) then
        ! Quantity is \nabla^2_\xi\Phi
        int_qua(i_ray_kappa+1) = (sinphi**2 + costhe**2*cosphi**2)*ray_fields_c(i_ray_tidal  ) + & !xx term
                                 (cosphi**2 + costhe**2*sinphi**2)*ray_fields_c(i_ray_tidal+1) + & !yy term
                                 (sinthe**2                      )*ray_fields_c(i_ray_tidal+2) + & !zz term
                                 (-sinthe**2*sin2phi             )*ray_fields_c(i_ray_tidal+3) + & !xy term
                                 (-cosphi*sin2the                )*ray_fields_c(i_ray_tidal+4) + & !xz term
                                 (-sinphi*sin2the                )*ray_fields_c(i_ray_tidal+5)     !yz term
     end if
     integral_cell(i_ray_kappa:i_ray_kappa+1)   = -((chi_B**2-chi_A**2)/2.0D0-(chi_B**3-chi_A**3)/(3.0D0*ray_chi_s)) * &
                                                  int_qua(i_ray_kappa:i_ray_kappa+1)
     ray_quants(iray,i_ray_kappa:i_ray_kappa+1) = ray_quants(iray,i_ray_kappa:i_ray_kappa+1)                               + &
                                                    integral_cell(i_ray_kappa:i_ray_kappa+1)*ray_Lbox**2*100.D0**2/aexp**2
  endif
  if(ray_do_shear) then
     ! Quantity is \nabla^1\nabla_1_\Phi - \nabla^2\nabla_2_\Phi
     int_qua(i_ray_shear  ) =    (sinphi**2 - costhe**2*cosphi**2)*ray_fields_c(i_ray_tidal  ) + & !xx term
                                 (cosphi**2 - costhe**2*sinphi**2)*ray_fields_c(i_ray_tidal+1) + & !yy term
                                 (-sinthe**2                     )*ray_fields_c(i_ray_tidal+2) + & !zz term
                                 (-(costhe**2+1.0D0)*sin2phi     )*ray_fields_c(i_ray_tidal+3) + & !xy term
                                 (cosphi*sin2the                 )*ray_fields_c(i_ray_tidal+4) + & !xz term
                                 (sinphi*sin2the                 )*ray_fields_c(i_ray_tidal+5)     !yz term
     ! Quantity is \nabla^1\nabla_2_\Phi
     int_qua(i_ray_shear+1) =    (-costhe* cosphi*sinphi  /sinthe)*ray_fields_c(i_ray_tidal  ) + & !xx term
                                 ( costhe* cosphi*sinphi  /sinthe)*ray_fields_c(i_ray_tidal+1) + & !yy term
!                                ( 0.0D0                         )*ray_fields_c(i_ray_tidal+2) + & !zz term
                                 ( costhe*cos2phi         /sinthe)*ray_fields_c(i_ray_tidal+3) + & !xy term
                                 ( sinphi                        )*ray_fields_c(i_ray_tidal+4) + & !xz term
                                 (-cosphi                        )*ray_fields_c(i_ray_tidal+5)     !yz term
     integral_cell(i_ray_shear:i_ray_shear+1)   = -((chi_B**2-chi_A**2)/2.0D0-(chi_B**3-chi_A**3)/(3.0D0*ray_chi_s)) * &
                                                  int_qua(i_ray_shear:i_ray_shear+1)
     ray_quants(iray,i_ray_shear:i_ray_shear+1) = ray_quants(iray,i_ray_shear:i_ray_shear+1)                               + &
                                                    integral_cell(i_ray_shear:i_ray_shear+1)*ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_deflt) then
     ! Quantity is \nable_1\Phi
     int_qua(i_ray_deflt  ) =    (costhe*cosphi)                  *ray_fields_c(i_ray_grads  ) + &  !x term
                                 (costhe*sinphi)                  *ray_fields_c(i_ray_grads+1) + &  !y term
                                 (sinthe)                         *ray_fields_c(i_ray_grads+2)      !z term
     ! Quantity is \nable_2\Phi
     int_qua(i_ray_deflt+1) =    (sinthe*cosphi)                  *ray_fields_c(i_ray_grads+1) - &  !y term
                                 (sinthe*sinphi)                  *ray_fields_c(i_ray_grads  )      !x term

     integral_cell(i_ray_deflt:i_ray_deflt+1)   = -((chi_B**2-chi_A**2)/2.0D0-(chi_B**3-chi_A**3)/(3.0D0*ray_chi_s)) * &
                                                  int_qua(i_ray_deflt:i_ray_deflt+1)
     ray_quants(iray,i_ray_deflt:i_ray_deflt+1) = ray_quants(iray,i_ray_deflt:i_ray_deflt+1)                               + &
                                                    integral_cell(i_ray_deflt:i_ray_deflt+1)*ray_Lbox**2*100.D0**2/aexp**2
  end if
  if(ray_do_iswrs) then
     int_qua(i_ray_iswrs  ) = 2.0D0*ray_fields_c(i_ray_phidt+1) 
     integral_cell(i_ray_iswrs) = int_qua(i_ray_iswrs)*(chi_A-chi_B)
     aexp_mean              = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
     ray_quants(iray,i_ray_iswrs) = ray_quants(iray,i_ray_iswrs)+integral_cell(i_ray_iswrs)*aexp_mean*ray_Lbox**3*100.0D0**3
  end if
  ! +++ add new stuff here to calculate other observables +++
  
end subroutine ray_integrate_ngp
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< 
