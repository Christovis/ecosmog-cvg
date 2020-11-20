! ------------------------------------------------------------------------------
! Multigrid extradof solver for refined AMR levels
! ------------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Note: We solve fine_fine-level explicitly and need to distinguish between
! the vector components using isf.
!   if isf > 1 : transverse/B-mode components
!   if isf = 1 : longitudinal/chi-mode
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -------------------------------------------------------------------
!     potential             sf            active_mg(myid,ilevel)%u(:,1)
!     truncation error      N/A           active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                      N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!     restricted sf         N/A           active_mg(myid,ilevel)%u(:,5)
!     restricted dens       N/A           active_mg(myid,ilevel)%u(:,6)
!     -------------------------------------------------------------------

! ------------------------------------------------------------------------------
! Computate residual of the fine (AMR level) and store into f(:,1)
! ------------------------------------------------------------------------------
subroutine cmp_residual_mg_fine_extradof(ilevel,isf)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel,isf

   real(dp) :: dx,dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,idim,inbor
   integer  :: iskip_amr
   integer  :: igshift,igrid_nbor_amr,icell_nbor_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: eta,op,dop                       ! extradof Poisson equ.
   real(dp) :: sfc,sf_lp_lp                     ! chi-mode
   real(dp) :: bf_src1,bf_src2,bf_src3,bf_src4  ! B-mode source terms
   real(dp), dimension(1:3) :: nb_sum_sfi,nb_diff_sfi,nb_sum_sfij
   real(dp), dimension(1:3) :: nb_sum_sf_lpi,nb_diff_sf_lpi,nb_sum_sf_lpij
   real(dp), dimension(1:3) :: nb_sum_bi
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,jnbor
   logical  :: bdy

   if(.not.(extradof4).and.(isf > 1)) then
      write(*,*) "In cmp_residual_mg_fine_extradof: extradof4 is not set to true"
      call clean_stop
   endif

   dx  = 0.5d0**ilevel
   dx2 = dx*dx

   ngrid=active(ilevel)%ngrid

   if(isf.eq.1) then
      ! longitudinal/chi-mode of proca field
      do ind=1,twotondim
         iskip_amr = ncoarse+(ind-1)*ngridmax
         do igrid_mg=1,ngrid,nvector
            nbatch=MIN(nvector,ngrid-igrid_mg+1)
            do i=1,nbatch
               igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
               icell_amr(i) = iskip_amr+igrid_amr(i)
            end do
            call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
            do i=1,nbatch
               sfc = sf(icell_amr(i))
               if(flag2(icell_amr(i))/ngridmax==0) then
                  bdy = .false.
                  do jnbor=1,27
                     if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                        bdy = .true.
                     end if
                  end do
                  if(bdy==.true.) cycle
                  nb_sum_sfi(1)  = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))
                  nb_sum_sfi(2)  = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))
                  nb_sum_sfi(3)  = sf(nbors_cells(i, 5))+sf(nbors_cells(i,23))

                  nb_diff_sf_lpi(1)  = sf_lp(nbors_cells(i,13))-sf_lp(nbors_cells(i,15))
                  nb_diff_sf_lpi(2)  = sf_lp(nbors_cells(i,11))-sf_lp(nbors_cells(i,17))
                  nb_diff_sf_lpi(3)  = sf_lp(nbors_cells(i, 5))-sf_lp(nbors_cells(i,23))

                  nb_sum_sfij(1) = sf(nbors_cells(i,18))+sf(nbors_cells(i,10)) - &
                                 & sf(nbors_cells(i,12))-sf(nbors_cells(i,16)) ! sf_{,xy}
                  nb_sum_sfij(2) = sf(nbors_cells(i,24))+sf(nbors_cells(i, 4)) - &
                                 & sf(nbors_cells(i, 6))-sf(nbors_cells(i,22)) ! sf_{,xz}
                  nb_sum_sfij(3) = sf(nbors_cells(i,26))+sf(nbors_cells(i, 2)) - &
                                 & sf(nbors_cells(i, 8))-sf(nbors_cells(i,20)) ! sf_{,yz}

                  nb_sum_sf_lpij(1) = sf_lp(nbors_cells(i,18))+sf_lp(nbors_cells(i,10)) - &
                                    & sf_lp(nbors_cells(i,12))-sf_lp(nbors_cells(i,16)) ! sl_lp_{,xy}
                  nb_sum_sf_lpij(2) = sf_lp(nbors_cells(i,24))+sf_lp(nbors_cells(i, 4)) - &
                                    & sf_lp(nbors_cells(i, 6))-sf_lp(nbors_cells(i,22)) ! sl_lp_{,xz}
                  nb_sum_sf_lpij(3) = sf_lp(nbors_cells(i,26))+sf_lp(nbors_cells(i, 2)) - &
                                    & sf_lp(nbors_cells(i, 8))-sf_lp(nbors_cells(i,20)) ! sl_lp_{,yz}

                  eta = 1.0d0/3.0d0/dx2**2*nb_sum_sfi(1) *(2.0d0*nb_sum_sfi(1)-nb_sum_sfi(2)-nb_sum_sfi(3))    + &
                        & 1.0d0/3.0d0/dx2**2*nb_sum_sfi(2) *(2.0d0*nb_sum_sfi(2)-nb_sum_sfi(1)-nb_sum_sfi(3))    + &
                        & 1.0d0/3.0d0/dx2**2*nb_sum_sfi(3) *(2.0d0*nb_sum_sfi(3)-nb_sum_sfi(1)-nb_sum_sfi(2))    + &
                        &     0.125d0/dx2**2*nb_sum_sfij(1)*       nb_sum_sfij(1)                                + &
                        &     0.125d0/dx2**2*nb_sum_sfij(2)*       nb_sum_sfij(2)                                + &
                        &     0.125d0/dx2**2*nb_sum_sfij(3)*       nb_sum_sfij(3)

                  eta = eta + 1.5d0*(alpha_cvg/beta_cvg)*omega_m*aexp*(rho(icell_amr(i))-rho_tot)
                  eta = alpha_cvg**2 + 8d0/3d0 * eta

                  if(eta<0.0d0) then
                     write(*,*) 'In cmp_residual_mg_fine_extradof: imaginary square root!',eta
                     stop
                  end if
                  if (eta.eq.0.0d0) then
                     eta = 1.0d-9
                  end if
                  if(alpha_cvg>=0.d0) then
                     op  = op-0.75d0*(-alpha_cvg + sqrt(eta))*dx2
                  else
                     op  = op-0.75d0*(-alpha_cvg - sqrt(eta))*dx2
                  end if

                  f(icell_amr(i),1) = -op/dx2
               else
                  f(icell_amr(i),1) = 0.0d0
               end if
            end do
         end do
      end do
   else
      ! transverse/B-mode of proca field
      do ind=1,twotondim
         iskip_amr = ncoarse+(ind-1)*ngridmax
         do igrid_mg=1,ngrid,nvector
            nbatch=MIN(nvector,ngrid-igrid_mg+1)
            do i=1,nbatch
               igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
               icell_amr(i) = iskip_amr+igrid_amr(i)
            end do
            call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
            do i=1,nbatch
               sfc = cbf(icell_amr(i),isf-1)
               if(flag2(icell_amr(i))/ngridmax==0) then
                  bdy = .false.
                  do jnbor=1,27
                     if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                        bdy = .true.
                     end if
                  end do
                  if(bdy==.true.) cycle
                  nb_sum_sfi(1)  = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))-2.0d0*sf(icell_amr(i)) ! sf_{,xx}
                  nb_sum_sfi(2)  = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))-2.0d0*sf(icell_amr(i)) ! sf_{,yy}
                  nb_sum_sfi(3)  = sf(nbors_cells(i,5))+sf(nbors_cells(i,23))-2.0d0*sf(icell_amr(i)) ! sf_{,zz}

                  nb_sum_sf_lpi(1)  = sf_lp(nbors_cells(i,13))+sf_lp(nbors_cells(i,15)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,xx}
                  nb_sum_sf_lpi(2)  = sf_lp(nbors_cells(i,11))+sf_lp(nbors_cells(i,17)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,yy}
                  nb_sum_sf_lpi(3)  = sf_lp(nbors_cells(i,5))+sf_lp(nbors_cells(i,23)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,zz}

                  sf_lp_lp = nb_sum_sf_lpi(1) + nb_sum_sf_lpi(2) + nb_sum_sf_lpi(3)

                  nb_diff_sfi(1)  = sf(nbors_cells(i,15))-sf(nbors_cells(i,13))
                  nb_diff_sfi(2)  = sf(nbors_cells(i,17))-sf(nbors_cells(i,11))
                  nb_diff_sfi(3)  = sf(nbors_cells(i, 23))-sf(nbors_cells(i,5))

                  nb_diff_sf_lpi(1)  = sf_lp(nbors_cells(i,15))-sf_lp(nbors_cells(i,13))
                  nb_diff_sf_lpi(2)  = sf_lp(nbors_cells(i,17))-sf_lp(nbors_cells(i,11))
                  nb_diff_sf_lpi(3)  = sf_lp(nbors_cells(i, 23))-sf_lp(nbors_cells(i,5))

                  nb_sum_sfij(1) = sf(nbors_cells(i,18))+sf(nbors_cells(i,10)) - &
                                 & sf(nbors_cells(i,12))-sf(nbors_cells(i,16)) ! sf_{,xy}
                  nb_sum_sfij(2) = sf(nbors_cells(i,24))+sf(nbors_cells(i, 4)) - &
                                 & sf(nbors_cells(i, 6))-sf(nbors_cells(i,22)) ! sf_{,xz}
                  nb_sum_sfij(3) = sf(nbors_cells(i,26))+sf(nbors_cells(i, 2)) - &
                                 & sf(nbors_cells(i, 8))-sf(nbors_cells(i,20)) ! sf_{,yz}

                  nb_sum_sf_lpij(1) = sf_lp(nbors_cells(i,18))+sf_lp(nbors_cells(i,10)) - &
                                    & sf_lp(nbors_cells(i,12))-sf_lp(nbors_cells(i,16)) ! sf_lp_{,xy}
                  nb_sum_sf_lpij(2) = sf_lp(nbors_cells(i,24))+sf_lp(nbors_cells(i, 4)) - &
                                    & sf_lp(nbors_cells(i, 6))-sf_lp(nbors_cells(i,22)) ! sf_lp_{,xz}
                  nb_sum_sf_lpij(3) = sf_lp(nbors_cells(i,26))+sf_lp(nbors_cells(i, 2)) - &
                                    & sf_lp(nbors_cells(i, 8))-sf_lp(nbors_cells(i,20)) ! sf_lp_{,yz}

                  if(isf.eq.1) then
                     bf_src1 = nb_sum_sfij(1)*nb_diff_sf_lpi(2) &
                             + nb_sum_sfij(2)*nb_diff_sf_lpi(3) &
                             + nb_sum_sfi(1)*nb_diff_sf_lpi(1)
                     bf_src2 = nb_diff_sfi(1) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(1)
                     bf_src4 = -nb_diff_sfi(2)*nb_sum_sf_lpij(1) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpij(2) &
                             - nb_diff_sfi(1)*nb_sum_sf_lpi(1)
                  endif

                  if(isf.eq.2) then
                     bf_src1 = nb_sum_sfij(1)*nb_diff_sf_lpi(1) &
                             + nb_sum_sfij(3)*nb_diff_sf_lpi(3) &
                             + nb_sum_sfi(2)*nb_diff_sf_lpi(2)
                     bf_src2 = nb_diff_sfi(2) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(2)
                     bf_src4 = -nb_diff_sfi(1)*nb_sum_sf_lpij(1) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpij(3) &
                             - nb_diff_sfi(2)*nb_sum_sf_lpi(2)
                  endif

                  if(isf.eq.3) then
                     bf_src1 = nb_sum_sfij(2)*nb_diff_sf_lpi(1) &
                             + nb_sum_sfij(3)*nb_diff_sf_lpi(2) &
                             + nb_sum_sfi(3)*nb_diff_sf_lpi(3)
                     bf_src2 = nb_diff_sfi(3) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(3)
                     bf_src4 = -nb_diff_sfi(1)*nb_sum_sf_lpij(2) &
                             - nb_diff_sfi(2)*nb_sum_sf_lpij(3) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpi(3)
                  endif

                  nb_sum_bi(1) = cbf(nbors_cells(i,13),isf-1)+cbf(nbors_cells(i,15),isf-1)
                  nb_sum_bi(2) = cbf(nbors_cells(i,11),isf-1)+cbf(nbors_cells(i,17),isf-1)
                  nb_sum_bi(3) = cbf(nbors_cells(i,5),isf-1)+cbf(nbors_cells(i,23),isf-1)

                  ! op = lhs of scalar field Poission (Eq. 59)
                  op = nb_sum_bi(1)+nb_sum_bi(2)+nb_sum_bi(3) - 6.0d0*sfc
                  ! op = lhs - rhs of scalar field Poission (Eq. 59)
                  op = op - param_b3*(bf_src1 + bf_src2 + bf_src3 + bf_src4)/dx**3
                  dop = -6.0d0

                  f(icell_amr(i),1) = -op/dx2
               else
                  f(icell_amr(i),1) = 0.0d0
               end if
            end do
         end do
      end do

   endif
end subroutine cmp_residual_mg_fine_extradof

! ##################################################################
! ##################################################################

subroutine cmp_residual_norm2_fine_extradof(ilevel,norm2,n_cell_f)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel
   real(dp), intent(out) :: norm2,n_cell_f

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg
   integer  :: igrid_amr,icell_amr,iskip_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**2
   ngrid=active(ilevel)%ngrid

   norm2    = 0.0d0
   n_cell_f = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         norm2    = norm2+abs(f(icell_amr,1))**2
         n_cell_f = n_cell_f+1.0d0
      end do
   end do
!  norm2 = dx2*norm2

end subroutine cmp_residual_norm2_fine_extradof

! ------------------------------------------------------------------------------
! Computate sf_lp of the AMR level
! ------------------------------------------------------------------------------
subroutine cmp_sf_lp_fine_extradof(ilevel)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel

   real(dp) :: dx,dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,idim,inbor
   integer  :: iskip_amr
   integer  :: igshift,igrid_nbor_amr,icell_nbor_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: sfc,op0
   real(dp), dimension(1:3) :: nb_sum_sfp
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,jnbor
   logical  :: bdy
   integer  :: ix,iy,iz
   real(dp) :: xx,yy,zz

   dx  = 0.5d0**ilevel
   dx2 = dx**2

   ngrid=active(ilevel)%ngrid

   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      do igrid_mg=1,ngrid,nvector
         nbatch=MIN(nvector,ngrid-igrid_mg+1)
         do i=1,nbatch
            igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
            icell_amr(i) = iskip_amr+igrid_amr(i)
         end do
         call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
         do i=1,nbatch
            sfc = sf(icell_amr(i))
            if(flag2(icell_amr(i))/ngridmax==0) then
               bdy = .false.
               do jnbor=1,27
                  if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                     bdy = .true.
                  end if
               end do
               if(bdy==.true.) cycle
               nb_sum_sfp(1) = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))
               nb_sum_sfp(2) = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))
               nb_sum_sfp(3) = sf(nbors_cells(i, 5))+sf(nbors_cells(i,23))

               op0 = nb_sum_sfp(1)+nb_sum_sfp(2)+nb_sum_sfp(3) - 6.0d0*sfc

               sf_lp(icell_amr(i)) = op0

            end if
         end do
      end do
   end do

end subroutine cmp_sf_lp_fine_extradof

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------
subroutine gauss_seidel_mg_fine_extradof(ilevel,isf,redstep)
   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters
   implicit none
   integer, intent(in) :: ilevel,isf
   logical, intent(in) :: redstep

   integer, dimension(1:3,1:2,1:8) :: iii,jjj
   integer, dimension(1:3,1:8)     :: ired,iblack

   real(dp) :: dx,dx2,dx3
   integer  :: ngrid
   integer  :: ind,ind0,igrid_mg,idim,inbor
   integer  :: iskip_amr
   integer  :: igshift,igrid_nbor_amr,icell_nbor_amr

   real(dp) :: dtwondim = (twondim)

   real(dp) :: eta,sfc,op,dop,sf_lp_lp
   real(dp) :: bf_src1,bf_src2,bf_src3,bf_src4  ! B-mode source terms
   real(dp), dimension(1:3) :: nb_sum_sfi,nb_sum_sfij,nb_diff_sfi
   real(dp), dimension(1:3) :: nb_sum_sf_lpi,nb_sum_sf_lpij,nb_diff_sf_lpi
   real(dp), dimension(1:3) :: nb_sum_bi
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,j,jnbor
   logical  :: bdy
   real(dp) :: sf0,sor

   ! Set constants
   dx  = 0.5d0**ilevel
   dx2 = dx*dx
   dx3 = dx2*dx

   ired  (1,1:8)=(/1,0,0,0,0,0,0,0/)
   iblack(1,1:8)=(/2,0,0,0,0,0,0,0/)
   ired  (2,1:8)=(/1,4,0,0,0,0,0,0/)
   iblack(2,1:8)=(/2,3,0,0,0,0,0,0/)
   ired  (3,1:8)=(/1,4,6,7,2,3,5,8/)
   iblack(3,1:8)=(/1,2,3,4,5,6,7,8/)

   ngrid=active(ilevel)%ngrid ! How many active grid at level #ilevel on myid

   if(isf.eq.1) then
      ! Loop over cells, with red/black ordering
      do ind0=1,twotondim      ! Only half of the cells for a red or black sweep
         if(redstep) then
            ind = ired  (ndim,ind0)
         else
            ind = iblack(ndim,ind0)
         end if

         iskip_amr = ncoarse+(ind-1)*ngridmax

         ! chi-term; Loop over active grids
         do igrid_mg=1,ngrid,nvector
            nbatch=MIN(nvector,ngrid-igrid_mg+1)
            do i=1,nbatch
               igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
               icell_amr(i) = iskip_amr+igrid_amr(i)
            end do
            call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
            do i=1,nbatch
               sfc = sf(icell_amr(i))

               if(flag2(icell_amr(i))/ngridmax==0) then
                  bdy = .false.
                  do jnbor=1,27
                     if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                        bdy = .true.
                     end if
                  end do
                  if(bdy==.true.) cycle
                  nb_sum_sfi(1)  = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))
                  nb_sum_sfi(2)  = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))
                  nb_sum_sfi(3)  = sf(nbors_cells(i, 5))+sf(nbors_cells(i,23))

                  nb_diff_sf_lpi(1)  = sf_lp(nbors_cells(i,13))-sf_lp(nbors_cells(i,15))
                  nb_diff_sf_lpi(2)  = sf_lp(nbors_cells(i,11))-sf_lp(nbors_cells(i,17))
                  nb_diff_sf_lpi(3)  = sf_lp(nbors_cells(i, 5))-sf_lp(nbors_cells(i,23))

                  nb_sum_sfij(1) = sf(nbors_cells(i,18))+sf(nbors_cells(i,10)) - &
                                 & sf(nbors_cells(i,12))-sf(nbors_cells(i,16)) ! sf_{,xy}
                  nb_sum_sfij(2) = sf(nbors_cells(i,24))+sf(nbors_cells(i, 4)) - &
                                 & sf(nbors_cells(i, 6))-sf(nbors_cells(i,22)) ! sf_{,xz}
                  nb_sum_sfij(3) = sf(nbors_cells(i,26))+sf(nbors_cells(i, 2)) - &
                                 & sf(nbors_cells(i, 8))-sf(nbors_cells(i,20)) ! sf_{,yz}

                  nb_sum_sf_lpij(1) = sf_lp(nbors_cells(i,18))+sf_lp(nbors_cells(i,10)) - &
                                    & sf_lp(nbors_cells(i,12))-sf_lp(nbors_cells(i,16)) ! sf_lp_{,xy}
                  nb_sum_sf_lpij(2) = sf_lp(nbors_cells(i,24))+sf_lp(nbors_cells(i, 4)) - &
                                    & sf_lp(nbors_cells(i, 6))-sf_lp(nbors_cells(i,22)) ! sf_lp_{,xz}
                  nb_sum_sf_lpij(3) = sf_lp(nbors_cells(i,26))+sf_lp(nbors_cells(i, 2)) - &
                                    & sf_lp(nbors_cells(i, 8))-sf_lp(nbors_cells(i,20)) ! sf_lp_{,yz}

                  eta = 1.0d0/3.0d0/dx2**2*nb_sum_sfi(1) *(2.0d0*nb_sum_sfi(1)-nb_sum_sfi(2)-nb_sum_sfi(3))    + &
                        & 1.0d0/3.0d0/dx2**2*nb_sum_sfi(2) *(2.0d0*nb_sum_sfi(2)-nb_sum_sfi(1)-nb_sum_sfi(3))    + &
                        & 1.0d0/3.0d0/dx2**2*nb_sum_sfi(3) *(2.0d0*nb_sum_sfi(3)-nb_sum_sfi(1)-nb_sum_sfi(2))    + &
                        &     0.125d0/dx2**2*nb_sum_sfij(1)*       nb_sum_sfij(1)                                + &
                        &     0.125d0/dx2**2*nb_sum_sfij(2)*       nb_sum_sfij(2)                                + &
                        &     0.125d0/dx2**2*nb_sum_sfij(3)*       nb_sum_sfij(3)

                  ! derived form Eq. 69
                  eta = eta + 1.5d0*(alpha_cvg/beta_cvg)*omega_m*aexp*(rho(icell_amr(i))-rho_tot)
                  eta = alpha_cvg**2 + 8.0d0/3.0d0 * eta

                  if(eta<0.0d0) then
                     write(*,*) 'In gauss_seidel_mg_fine_extradof: imaginary square root!',eta
                     stop
                  end if

                  if (eta.eq.0.0d0) then
                     eta = 1.0d-9
                  end if

                  op = nb_sum_sfi(1)+nb_sum_sfi(2)+nb_sum_sfi(3) - 6.0d0*sfc

                  if(alpha_cvg>=0.d0) then
                     ! derived form Eq. 69
                     op  = op - (0.75d0*(-alpha_cvg+sqrt(eta)) - sf_src_mean)*dx2
                     dop = -6.0d0
                     sf(icell_amr(i)) = sf(icell_amr(i))-op/dop
                     sf_src(icell_amr(i)) = 0.75d0*(-alpha_cvg + sqrt(eta))
                  else
                     ! derived form Eq. 69
                     op  = op - (0.75d0*(-alpha_cvg-sqrt(eta)) - sf_src_mean)*dx2
                     dop = -6.0d0
                     sf(icell_amr(i)) = sf(icell_amr(i))-op/dop
                     sf_src(icell_amr(i)) = 0.75d0*(-alpha_cvg - sqrt(eta))
                  end if

               end if
            end do
         end do
      end do

   else
      ! B-terms; Loop over cells, with red/black ordering
      do ind0=1,twotondim      ! Only half of the cells for a red or black sweep
         if(redstep) then
            ind = ired  (ndim,ind0)
         else
            ind = iblack(ndim,ind0)
         end if

         iskip_amr = ncoarse+(ind-1)*ngridmax

         ! Loop over active grids
         do igrid_mg=1,ngrid,nvector
            nbatch=MIN(nvector,ngrid-igrid_mg+1)
            do i=1,nbatch
               igrid_amr(i) = active(ilevel)%igrid(i+igrid_mg-1)
               icell_amr(i) = iskip_amr+igrid_amr(i)
            end do
            call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
            do i=1,nbatch
               sfc = cbf(icell_amr(i),isf-1)

               if(flag2(icell_amr(i))/ngridmax==0) then
                  bdy = .false.
                  do jnbor=1,27
                     if(flag2(nbors_cells(i,jnbor))/ngridmax.ne.0) then
                        bdy = .true.
                     end if
                  end do
                  if(bdy==.true.) cycle
                  nb_sum_sfi(1)  = sf(nbors_cells(i,13))+sf(nbors_cells(i,15))-2.0d0*sf(icell_amr(i)) ! sf_{,xx}
                  nb_sum_sfi(2)  = sf(nbors_cells(i,11))+sf(nbors_cells(i,17))-2.0d0*sf(icell_amr(i)) ! sf_{,yy}
                  nb_sum_sfi(3)  = sf(nbors_cells(i,5))+sf(nbors_cells(i,23))-2.0d0*sf(icell_amr(i)) ! sf_{,zz}
                  nb_sum_bi(1) = cbf(nbors_cells(i,13),isf-1)+cbf(nbors_cells(i,15),isf-1)
                  nb_sum_bi(2) = cbf(nbors_cells(i,11),isf-1)+cbf(nbors_cells(i,17),isf-1)
                  nb_sum_bi(3) = cbf(nbors_cells(i,5),isf-1)+cbf(nbors_cells(i,23),isf-1)


                  nb_sum_sf_lpi(1)  = sf_lp(nbors_cells(i,13))+sf_lp(nbors_cells(i,15)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,xx}
                  nb_sum_sf_lpi(2)  = sf_lp(nbors_cells(i,11))+sf_lp(nbors_cells(i,17)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,yy}
                  nb_sum_sf_lpi(3)  = sf_lp(nbors_cells(i,5))+sf_lp(nbors_cells(i,23)) &
                                    - 2.0d0*sf_lp(icell_amr(i)) ! sf_lp_{,zz}

                  sf_lp_lp = nb_sum_sf_lpi(1) + nb_sum_sf_lpi(2) + nb_sum_sf_lpi(3)

                  nb_diff_sfi(1)  = sf(nbors_cells(i,15))-sf(nbors_cells(i,13))
                  nb_diff_sfi(2)  = sf(nbors_cells(i,17))-sf(nbors_cells(i,11))
                  nb_diff_sfi(3)  = sf(nbors_cells(i, 23))-sf(nbors_cells(i,5))

                  nb_diff_sf_lpi(1)  = sf_lp(nbors_cells(i,15))-sf_lp(nbors_cells(i,13))
                  nb_diff_sf_lpi(2)  = sf_lp(nbors_cells(i,17))-sf_lp(nbors_cells(i,11))
                  nb_diff_sf_lpi(3)  = sf_lp(nbors_cells(i, 23))-sf_lp(nbors_cells(i,5))

                  nb_sum_sfij(1) = sf(nbors_cells(i,18))+sf(nbors_cells(i,10)) - &
                                 & sf(nbors_cells(i,12))-sf(nbors_cells(i,16)) ! sf_{,xy}
                  nb_sum_sfij(2) = sf(nbors_cells(i,24))+sf(nbors_cells(i, 4)) - &
                                 & sf(nbors_cells(i, 6))-sf(nbors_cells(i,22)) ! sf_{,xz}
                  nb_sum_sfij(3) = sf(nbors_cells(i,26))+sf(nbors_cells(i, 2)) - &
                                 & sf(nbors_cells(i, 8))-sf(nbors_cells(i,20)) ! sf_{,yz}

                  nb_sum_sf_lpij(1) = sf_lp(nbors_cells(i,18))+sf_lp(nbors_cells(i,10)) - &
                                    & sf_lp(nbors_cells(i,12))-sf_lp(nbors_cells(i,16)) ! sf_lp_{,xy}
                  nb_sum_sf_lpij(2) = sf_lp(nbors_cells(i,24))+sf_lp(nbors_cells(i, 4)) - &
                                    & sf_lp(nbors_cells(i, 6))-sf_lp(nbors_cells(i,22)) ! sf_lp_{,xz}
                  nb_sum_sf_lpij(3) = sf_lp(nbors_cells(i,26))+sf_lp(nbors_cells(i, 2)) - &
                                    & sf_lp(nbors_cells(i, 8))-sf_lp(nbors_cells(i,20)) ! sf_lp_{,yz}

                  if(isf.eq.1) then
                     bf_src1 = nb_sum_sfij(1)*nb_diff_sf_lpi(2) &
                             + nb_sum_sfij(2)*nb_diff_sf_lpi(3) &
                             + nb_sum_sfi(1)*nb_diff_sf_lpi(1)
                     bf_src2 = nb_diff_sfi(1) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(1)
                     bf_src4 = -nb_diff_sfi(2)*nb_sum_sf_lpij(1) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpij(2) &
                             - nb_diff_sfi(1)*nb_sum_sf_lpi(1)
                  endif

                  if(isf.eq.2) then
                     bf_src1 = nb_sum_sfij(1)*nb_diff_sf_lpi(1) &
                             + nb_sum_sfij(3)*nb_diff_sf_lpi(3) &
                             + nb_sum_sfi(2)*nb_diff_sf_lpi(2)
                     bf_src2 = nb_diff_sfi(2) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(2)
                     bf_src4 = -nb_diff_sfi(1)*nb_sum_sf_lpij(1) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpij(3) &
                             - nb_diff_sfi(2)*nb_sum_sf_lpi(2)
                  endif

                  if(isf.eq.3) then
                     bf_src1 = nb_sum_sfij(2)*nb_diff_sf_lpi(1) &
                             + nb_sum_sfij(3)*nb_diff_sf_lpi(2) &
                             + nb_sum_sfi(3)*nb_diff_sf_lpi(3)
                     bf_src2 = nb_diff_sfi(3) * sf_lp_lp
                     bf_src3 = -sf_lp(icell_amr(i)) * nb_diff_sf_lpi(3)
                     bf_src4 = -nb_diff_sfi(1)*nb_sum_sf_lpij(2) &
                             - nb_diff_sfi(2)*nb_sum_sf_lpij(3) &
                             - nb_diff_sfi(3)*nb_sum_sf_lpi(3)
                  endif

                  ! op = lhs - rhs Eq. 59
                  op = nb_sum_bi(1)+nb_sum_bi(2)+nb_sum_bi(3)-6.0d0*sfc
                  op = op - param_b3*(bf_src1 + bf_src2 + bf_src3 + bf_src4)/dx3
                  dop = -6.0d0
                  cbf(icell_amr(i),isf-1) = cbf(icell_amr(i),isf-1)-op/dop

               end if
            end do
         end do
      end do
   endif
end subroutine gauss_seidel_mg_fine_extradof

! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine_reverse_extradof(ifinelevel)

   use amr_commons
   use poisson_commons
   use extradof_commons

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg) ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr              ! amr fine-cell index
         ! Is fine cell masked?
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)                         ! amr coarse-cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax ! amr coarse-grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                    ! cpu for coarse cell

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                       ! mg coarse-grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                        ! mg coarse-cell index

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=f(icell_f_amr,1)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
   end do
end subroutine restrict_residual_fine_reverse_extradof


! ------------------------------------------------------------------------
! extradof restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_extradof_fine_reverse_extradof(ifinelevel,isf)
   use amr_commons
   use poisson_commons
   use extradof_commons
   implicit none
   integer, intent(in) :: ifinelevel,isf

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: sf1
   real(dp) :: dtwotondim = (twotondim)
   integer,dimension(:),allocatable::n_masked

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   allocate(n_masked(1:active(ifinelevel)%ngrid))
   n_masked(1:active(ifinelevel)%ngrid)=0
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index
         if(f(icell_f_amr,3)<=0d0) n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)                                ! amr coarse-cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax        ! amr coarse-grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                           ! cpu for coarse cell

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                              ! mg coarse-grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                               ! mg coarse-cell index

         ! If coarse cell masked, it is boundary and R\tilde{u} is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Restriction to compute the sf value in the coarse cell
         if(isf.eq.1) then
            sf1=sf(icell_f_amr)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         else
            sf1=cbf(icell_f_amr,isf-1)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         endif
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)+sf1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_extradof_fine_reverse_extradof

! ------------------------------------------------------------------------
! density restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_density_fine_reverse_extradof(ifinelevel)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: rho1
   real(dp) :: dtwotondim = (twotondim)

   integer  :: icoarselevel
   integer  :: ix,iy,iz
   real(dp) :: ctilde,xc,dx
   integer,dimension(:),allocatable::n_masked

   icoarselevel=ifinelevel-1

   dx = (0.5d0)**ifinelevel
   ctilde = sol/boxlen_ini/100000.0d0

   allocate(n_masked(1:active(ifinelevel)%ngrid))
   n_masked(1:active(ifinelevel)%ngrid)=0
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr
         if(f(icell_f_amr,3)<=0d0) n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)               ! amr fine-grid index
         icell_f_amr=igrid_f_amr+iskip_f_amr                            ! amr fine-cell index
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)                                ! amr coarse-cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax        ! amr coarse-grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                           ! cpu for coarse cell

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                              ! mg coarse-grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                               ! mg coarse-cell index

         ! If coarse cell masked, it is boundary and R rho is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Restriction to compute the rho value in the coarse cell
         rho1=(sf_src(icell_f_amr)-sf_src_mean)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)+rho1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_density_fine_reverse_extradof


! ------------------------------------------------------------------------
! Interpolation (prolongation) and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine_extradof(ifinelevel)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i,ind_father,ind_average,ind_f,iskip_f_amr
   integer  :: ngrid_f,istart,nbatch
   integer  :: icell_c_amr,igrid_c_amr,igrid_c_mg,icell_c_mg
   integer  :: icoarselevel,ind_c,cpu_amr

   real(dp) :: a,b,c,d,coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save               :: igrid_f_amr, icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_father_grids
   real(dp), dimension(1:nvector), save               :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active(ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active(ifinelevel)%igrid(istart+i-1) ! amr fine grid index
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))                 ! amr coarse cell index
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids,nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax

         do i=1,nbatch
            ! note that the same icell_amr stores different things here and above
            icell_amr(i) = iskip_f_amr+igrid_f_amr(i)        ! amr fine cell index
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               if(f(icell_amr(i),3)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
!              ind_c       = (icell_c_amr-ncoarse-1)/ngridmax + 1
               ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               cpu_amr     = cpu_map(father(igrid_c_amr))
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               if(igrid_c_mg<=0) cycle

               icell_c_mg=(ind_c-1)*active_mg(cpu_amr,icoarselevel)%ngrid+igrid_c_mg
               ! only unmasked coarse cells contribute to the fine-cell correction
               if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0.0) cycle
               corr(i)=corr(i)+coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,1) &
                      &       -coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)
            end do
         end do
         ! Correct the extradof on the fine cells
         do i=1,nbatch
            sf(icell_amr(i))=sf(icell_amr(i))+corr(i)
         end do
      end do
   end do
   ! End loop over grids
end subroutine interpolate_and_correct_fine_extradof
