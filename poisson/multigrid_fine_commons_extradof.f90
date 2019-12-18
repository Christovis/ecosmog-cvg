! ------------------------------------------------------------------------------
! Multigrid extradof solver for refined AMR levels
! ------------------------------------------------------------------------------
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     extradof              sf            active_mg(myid,ilevel)%u(:,1)
!     truncation error      N/A           active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)         N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!     restricted sf         N/A           active_mg(myid,ilevel)%u(:,5)
!     restricted dens       N/A           active_mg(myid,ilevel)%u(:,6)
!     -----------------------------------------------------------------


! ------------------------------------------------------------------------------
! Main multigrid routine for the extradof, called by amr_step
! ------------------------------------------------------------------------------
subroutine multigrid_fine_extradof(ilevel,icount,isf)

   use amr_commons
   use poisson_commons
   use poisson_parameters
   use extradof_commons
   use extradof_parameters
   use pm_commons

   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel,icount,isf

   integer, parameter  :: MAXITER  = 50

   integer  :: ifine,i,iter,idim,info,icpu,ii,jj
   real(dp) :: res_norm2,res_norm2_tot,residual,residual_old
   real(dp) :: n_cell_f,n_cell_f_tot,n_cell_c,n_cell_c_tot
   integer  :: ind,igrid_mg,icell_mg,igrid_amr,icell_amr,ind_cell
   real(dp) :: trunc_norm2,trunc_norm2_tot,trunc_err

   logical  :: allmasked,allmasked_tot
   integer  :: ix,iy,iz
   real(dp) :: xx,yy,zz,dx
   integer  :: oscillation

   if(numbtot(1,ilevel)==0) return

   if(verbose) print '(A,I2)','V-cycle: entering fine multigrid at level ',ilevel

   if(.not.(extradof4).and.(isf > 1)) then
      write(*,*) "In multigrid_fine_extradof: isf can't be larger than 1 if extradof4 is False"
      call clean_stop
   endif


   residual_old = 1.0d0
   ! ---------------------------------------------------------------------
   ! Prepare first guess, mask and BCs at finest level
   ! ---------------------------------------------------------------------
   if(ilevel>levelmin)then
      ! if on the domain level, initialise sf and sf_lp to previous step values;
      !   (1) this is automatically done without calling any subroutines.
      ! if on refinements, initialise by interpolation from coarser level:
      !   (1) for sf_lp, real interpolation using old scheme;
      !   (2) for sf, real interpolation using new or old scheme (-DOLDINTERP);
      call make_initial_extradof(ilevel,icount)
   else
      ! If on domain grid, do nothing because we can use the result of last
      ! step as the initial guess of the current step
   endif

   ! communicate the cv-Galileon related arrays on virtual boundaries
   call make_virtual_fine_dp(sf(1),ilevel)
   if(extradof4) then
      call make_virtual_fine_dp(sf_lp(1),ilevel)
      do idim=1,ndim
         call make_virtual_fine_dp(cbf(1,idim),ilevel)
      enddo
   endif

   ! call make_fine_bc_rhs_extradof(ilevel)

   ! ---------------------------------------------------------------------
   ! Restrict relevant quantities up
   ! ---------------------------------------------------------------------
   if(ilevel>1) then
      ! The restricted density field
      call restrict_density_fine_reverse_extradof(ilevel)
      call make_reverse_mg_dp(6,ilevel-1)
      call make_virtual_mg_dp(6,ilevel-1)

      ! Restrict relevant quantities up for coarser levels
      do ifine=(ilevel-1),levelmin_mg,-1
         ! The restricted density field
         call restrict_density_coarse_reverse_extradof(ifine)
         call make_reverse_mg_dp(6,ifine-1)
         call make_virtual_mg_dp(6,ifine-1)
      end do
   end if

   if(verbose) then
      dx  = 0.5d0**ilevel
      do ind=1,twotondim
          do igrid_mg=1,active(ilevel)%ngrid
             igrid_amr = active(ilevel)%igrid(igrid_mg)     ! Grid amr index
             icell_amr = ncoarse+(ind-1)*ngridmax+igrid_amr ! Cell amr index
             iz = (ind-1)/4
             iy = (ind-1-4*iz)/2
             ix = (ind-1-4*iz-2*iy)
             xx = (dble(ix)-0.5d0)*dx+xg(igrid_amr,1)
             yy = (dble(iy)-0.5d0)*dx+xg(igrid_amr,2)
             zz = (dble(iz)-0.5d0)*dx+xg(igrid_amr,3)
             if(yy-0.5d0<0.002 .and. yy>0.5d0 .and. zz-0.5d0<0.002 .and. &
                zz>0.5d0 .and. xx-0.5d0<0.002 .and. xx>0.5d0 .and. &
                ilevel==levelmin) then
                write(*,*) 'last step:','test ps = ',ps(icell_amr),'test sf = ',sf(icell_amr)
             end if
          end do
       end do
    endif

   ! ---------------------------------------------------------------------
   ! Initiate solve at fine level
   ! ---------------------------------------------------------------------
   iter = 0
   oscillation = 0
   sf_src_mean = 0.0d0
   main_iteration_loop: do
      iter=iter+1
      ! Pre-smoothing
      do i=1,ngs_fine_extradof_pre
         call gauss_seidel_mg_fine_extradof(ilevel,isf,.true.)
         if(isf.eq.1) then
            call make_virtual_fine_dp(sf(1),ilevel)
            call make_virtual_fine_dp(sf_src(1),ilevel)
            call gauss_seidel_sf_src_mean_extradof(ilevel)
         else
            call make_virtual_fine_dp(cbf(1,isf-1),ilevel)
         endif
         call gauss_seidel_mg_fine_extradof(ilevel,isf,.false.)
         if(isf.eq.1) then
            call make_virtual_fine_dp(sf(1),ilevel)
            call make_virtual_fine_dp(sf_src(1),ilevel)
            call gauss_seidel_sf_src_mean_extradof(ilevel)
         else
            call make_virtual_fine_dp(cbf(1,isf-1),ilevel)
         endif
      end do

      ! Compute residual and restrict into upper level RHS
      call cmp_residual_mg_fine_extradof(ilevel,isf)
      call make_virtual_fine_dp(f(1,1),ilevel)

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
         active_mg(icpu,ilevel-1)%u(:,2)=0.0d0
         active_mg(icpu,ilevel-1)%u(:,5)=0.0d0
      end do

      ! Restrict and do communications
      call restrict_residual_fine_reverse_extradof(ilevel)
      call make_reverse_mg_dp(2,ilevel-1)
      call make_virtual_mg_dp(2,ilevel-1)
      call restrict_extradof_fine_reverse_extradof(ilevel,isf)
      call make_reverse_mg_dp(5,ilevel-1)
      call make_virtual_mg_dp(5,ilevel-1)
      call make_physical_rhs_coarse_extradof(ilevel-1)

      if(ilevel>1 .and. levelmin_mg<ilevel) then
         ! Make initial guess at upper level before solve
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,1)=active_mg(icpu,ilevel-1)%u(:,5)
         end do

         ! Multigrid-solve the upper level
         call recursive_multigrid_coarse_extradof(ilevel-1,safe_mode(ilevel))

         ! Interpolate coarse solution and correct fine solution
         call interpolate_and_correct_fine_extradof(ilevel)
         if(isf.eq.1) then
            call make_virtual_fine_dp(sf(1),ilevel)
         else
            call make_virtual_fine_dp(cbf(1,isf-1),ilevel)
         endif
      end if

      ! Post-smoothing
      do i=1,ngs_fine_extradof_pst
         call gauss_seidel_mg_fine_extradof(ilevel,isf,.true.)
         if(isf.eq.1) then
            call make_virtual_fine_dp(sf(1),ilevel)
            call make_virtual_fine_dp(sf_src(1),ilevel)
            call gauss_seidel_sf_src_mean_extradof(ilevel)
         else
            call make_virtual_fine_dp(cbf(1,isf-1),ilevel)
         endif
         call gauss_seidel_mg_fine_extradof(ilevel,isf,.false.)
         if(isf.eq.1) then
            call make_virtual_fine_dp(sf(1),ilevel)
            call make_virtual_fine_dp(sf_src(1),ilevel)
            call gauss_seidel_sf_src_mean_extradof(ilevel)
         else
            call make_virtual_fine_dp(cbf(1,isf-1),ilevel)
         endif
      end do

      ! Update fine residual
      call cmp_residual_mg_fine_extradof(ilevel,isf)
      call make_virtual_fine_dp(f(1,1),ilevel)
      call cmp_residual_norm2_fine_extradof(ilevel,res_norm2,n_cell_f)
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(n_cell_f ,n_cell_f_tot ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      res_norm2 = res_norm2_tot
      n_cell_f  = n_cell_f_tot
#endif
      residual = sqrt(res_norm2)/n_cell_f

      call cmp_uvar_norm2_coarse_extradof(2,ilevel-1,trunc_norm2,n_cell_c)
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(trunc_norm2,trunc_norm2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(n_cell_c   ,n_cell_c_tot   ,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      trunc_norm2 = trunc_norm2_tot
      n_cell_c    = n_cell_c_tot
#endif
      trunc_err = sqrt(trunc_norm2)/n_cell_c

      ! Verbosity
      if(verbose) print '(A,I4,A,1e13.6,A,1e13.6)','   ==> Step=',iter,', Residual =',residual,', Truncation error= ',trunc_err

      ! Converged?
      ! start of TODO: think about different convergence criteria for sf and cbf
      if(ilevel==levelmin) then
         if(residual<1.0d-8 .or. iter>=MAXITER .or. &
            abs(residual-residual_old)<1.0d-9  .or. &
            (residual<=0.33*trunc_err .and. iter>=20) .or. &
            (residual>residual_old) .or. aexp>0.38) exit
      else
         if(residual<1.0d-8 .or. iter>=MAXITER .or. &
            abs(residual-residual_old)<1.0d-9  .or. &
            (residual<=0.33*trunc_err .and. iter>=20)) exit
      end if

      if(residual<=residual_old) then
         if(oscillation>=1) then
            if(residual<0.33*trunc_err) exit
         end if
      end if

      if(residual>residual_old) then
         oscillation = oscillation+1
         if(residual<0.01*trunc_err .or. (oscillation>1 .and. residual<0.1*trunc_err)) exit
      end if
      ! end of TODO

      residual_old = residual
   end do main_iteration_loop

   if(extradof4) then
      call cmp_sf_lp_fine_extradof (      ilevel)
      call make_virtual_fine_dp(sf_lp (1),ilevel)
   endif

   if(myid==1) print '(A,I5,A,I5,A,1E15.6, A, 1e15.6)','   ==> Level=',ilevel,' Step=',&
               iter,' Residual=',residual,' Truncation error=',trunc_err

   if(myid==1 .and. iter==MAXITER) print *,'Warning: Fine multigrid extradof eqn fails to converge.'
   if(residual>1.0d-5 .and. residual>0.25*trunc_err) print *,'Warning2: Fine multigrid extradof eqn fails to converge.'

   if(verbose) then
      do ind=1,twotondim
         do igrid_mg=1,active(ilevel)%ngrid
            igrid_amr = active(ilevel)%igrid(igrid_mg)     ! Grid amr index
            icell_amr = ncoarse+(ind-1)*ngridmax+igrid_amr ! Cell amr index
            iz = (ind-1)/4
            iy = (ind-1-4*iz)/2
            ix = (ind-1-4*iz-2*iy)
            xx = (dble(ix)-0.5d0)*dx+xg(igrid_amr,1)
            yy = (dble(iy)-0.5d0)*dx+xg(igrid_amr,2)
            zz = (dble(iz)-0.5d0)*dx+xg(igrid_amr,3)
            if(yy-0.5d0<0.002 .and. yy>0.5d0 .and. zz-0.5d0<0.002 .and. &
               zz>0.5d0 .and. xx-0.5d0<0.002 .and. xx>0.5d0 .and. &
               ilevel==levelmin) then
               write(*,*) 'this step:','test ps = ',ps(icell_amr),'test sf = ',sf(icell_amr)
            end if
         end do
      end do
   endif

   ! ---------------------------------------------------------------------
   ! Cleanup MG levels after solve complete
   ! ---------------------------------------------------------------------
   do ifine=1,ilevel-1
      call cleanup_mg_level(ifine)
   end do

end subroutine multigrid_fine_extradof


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------------
recursive subroutine recursive_multigrid_coarse_extradof(ifinelevel,safe)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ifinelevel
   logical, intent(in) :: safe

   integer :: i,icpu,info,icycle,ncycle
   integer :: ind,igrid,icell_mg

   if(verbose) then
      if(ifinelevel>levelmin_mg) then
         print '(A,I2)','V-cycle: entering coarse multigrid at level ',ifinelevel
      else
         print '(A,I2)','V-cycle: entering coarsest multigrid level ',ifinelevel
      end if
   end if

   if(ifinelevel<=levelmin_mg) then
      if(ifinelevel<levelmin_mg) return
      do i=1,ngs_coarse_extradof_pst
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.true.)
         call make_virtual_mg_dp(1,ifinelevel)
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.false.)
         call make_virtual_mg_dp(1,ifinelevel)
      end do
      if(verbose) print '(A,I2)','V-cycle: leaving coarsest multigrid level ',ifinelevel
      return
   end if

   if(safe) then
      ncycle=ncycles_coarse_safe
   else
      ncycle=1
   endif

   do icycle=1,ncycle
      ! Pre-smoothing
      do i=1,ngs_coarse_extradof_pre
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.true.)
         call make_virtual_mg_dp(1,ifinelevel)
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.false.)
         call make_virtual_mg_dp(1,ifinelevel)
      end do

      ! Compute residual and restrict into upper level
      call cmp_residual_mg_coarse_extradof(ifinelevel)
      call make_virtual_mg_dp(3,ifinelevel)

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,2)=0.0d0
         active_mg(icpu,ifinelevel-1)%u(:,5)=0.0d0
      end do
      ! Restrict and do communications
      call restrict_residual_coarse_reverse_extradof(ifinelevel)
      call make_reverse_mg_dp(2,ifinelevel-1)
      call make_virtual_mg_dp(2,ifinelevel-1)
      call restrict_extradof_coarse_reverse_extradof(ifinelevel)
      call make_reverse_mg_dp(5,ifinelevel-1)
      call make_virtual_mg_dp(5,ifinelevel-1)
      call make_physical_rhs_coarse_extradof(ifinelevel-1)

      ! Initial guess from upper level before solve
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,1)=active_mg(icpu,ifinelevel-1)%u(:,5)
      end do

      ! Multigrid-solve the upper level
      call recursive_multigrid_coarse_extradof(ifinelevel-1,safe)

      ! Interpolate coarse solution and correct back into fine solution
      call interpolate_and_correct_coarse_extradof(ifinelevel)
      call make_virtual_mg_dp(1,ifinelevel)

      ! Post-smoothing
      do i=1,ngs_coarse_extradof_pst
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.true.)
         call make_virtual_mg_dp(1,ifinelevel)
         call gauss_seidel_mg_coarse_extradof(ifinelevel,safe,.false.)
         call make_virtual_mg_dp(1,ifinelevel)
      end do

   end do

   if(verbose) print '(A,I2)','V-cycle: leaving coarse multigrid at level ',ifinelevel

end subroutine recursive_multigrid_coarse_extradof

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------------
! Preprocess the fine (AMR) level vars to account for boundary conditions
! for the extradof solver
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
! Sets BC-modified RHS for the extradof solver into f(:,2), this is simply
! zero due to the form of the extradof equation (see paper)
!
! This subroutine is essentially useless
! ------------------------------------------------------------------------------
subroutine make_fine_bc_rhs_extradof(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters
   implicit none
   integer, intent(in) :: ilevel

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx,oneoverdx2,sf_b,nb_mask,nb_sf,w

   ! Arrays for vectorized interpol_extradof
   real(dp), dimension(1:nvector,1:twotondim) :: sf_int
   integer,  dimension(1:nvector) :: ind_cell

   integer  :: ngrid
   integer  :: ind,igrid_mg,idim,inbor
   integer  :: igrid_amr,icell_amr,iskip_amr
   integer  :: igshift,igrid_nbor_amr,icell_nbor_amr
   integer  :: ifathercell_nbor_amr

   integer  :: nx_loc
   real(dp) :: scale, fourpi
   integer  :: iz,iy,ix
   real(dp) :: xc,drho

   ! Set constants
   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) fourpi = 1.5D0*omega_m*aexp*scale

   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg) ! Grid amr index
         icell_amr = iskip_amr+igrid_amr ! Cell amr index

         f(icell_amr,2) = 0.0d0
         if(f(icell_amr,3)<=0.0) cycle ! Do not process masked cells
      end do
   end do

end subroutine make_fine_bc_rhs_extradof

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------------
! Multigrid level cleanup
! ------------------------------------------------------------------------------
subroutine cleanup_mg_level(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons

   implicit none

   integer, intent(in) :: ilevel
   integer :: igrid,icpu,cur_grid,cur_cpu

   ! ---------------------------------------------------------------------
   ! Cleanup lookup table
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      do igrid=1,active_mg(icpu,ilevel)%ngrid
         cur_grid=active_mg(icpu,ilevel)%igrid(igrid)
         cur_cpu=cpu_map(father(cur_grid))
         if(cur_cpu==myid) then
            lookup_mg(cur_grid)=0
         else
            lookup_mg(cur_grid)=-mod(flag2(cur_grid),ngridmax)
         end if
      end do
   end do

   ! ---------------------------------------------------------------------
   ! Deallocate communicators
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0)then
         deallocate(active_mg(icpu,ilevel)%igrid)
         deallocate(active_mg(icpu,ilevel)%u)
         deallocate(active_mg(icpu,ilevel)%f)
      endif
      active_mg(icpu,ilevel)%ngrid=0
      if(emission_mg(icpu,ilevel)%ngrid>0)then
         deallocate(emission_mg(icpu,ilevel)%igrid)
         deallocate(emission_mg(icpu,ilevel)%u)
         deallocate(emission_mg(icpu,ilevel)%f)
      endif
      emission_mg(icpu,ilevel)%ngrid=0
   end do

end subroutine cleanup_mg_level

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------------
! Multigrid level restoration for the Poisson equation solver next
! ------------------------------------------------------------------------------
subroutine restore_mg_level_extradof(ilevel)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: igrid, icpu, cur_grid, cur_cpu

   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0) then
         active_mg(icpu,ilevel)%u(:,1)=0.0d0
         active_mg(icpu,ilevel)%u(:,2)=0.0d0
         active_mg(icpu,ilevel)%u(:,3)=0.0d0
         active_mg(icpu,ilevel)%u(:,5)=0.0d0
         active_mg(icpu,ilevel)%u(:,6)=0.0d0
      endif

      if(emission_mg(icpu,ilevel)%ngrid>0) then
         emission_mg(icpu,ilevel)%u(:,1)=0.0d0
         emission_mg(icpu,ilevel)%u(:,2)=0.0d0
         emission_mg(icpu,ilevel)%u(:,3)=0.0d0
         emission_mg(icpu,ilevel)%u(:,5)=0.0d0
         emission_mg(icpu,ilevel)%u(:,6)=0.0d0
      endif

   end do

end subroutine restore_mg_level_extradof

subroutine restore_amr_level_extradof(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons

   implicit none

   integer,intent(in) :: ilevel

   integer :: ngrid,igrid,ind,ind_cell,ind_grid,iskip,i

   ngrid = active(ilevel)%ngrid

   do ind=1,twotondim
      iskip = ncoarse+(ind-1)*ngridmax
      do igrid=1,ngrid
           ind_grid = active(ilevel)%igrid(igrid)
           ind_cell = ind_grid+iskip
           do i=1,ndim-1
              f(ind_cell,i) = 0.0d0
           end do
      end do
   end do

end subroutine restore_amr_level_extradof

! ------------------------------------------------------------------------------
! Computate mean value of sf_src
! ------------------------------------------------------------------------------
subroutine gauss_seidel_sf_src_mean_extradof(ilevel)
   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters
   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel
   integer  :: ind,iskip,i,info
   real(dp) :: test_src,test_srcbar

   if(ilevel.ne.levelmin) return
   test_src=0.0D0

   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,active(ilevel)%ngrid
         test_src=test_src+sf_src(active(ilevel)%igrid(i)+iskip)
      end do
   end do

#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(test_src,test_srcbar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif

#ifdef WITHOUTMPI
   test_srcbar=test_src
#endif

   test_srcbar=test_srcbar/dble((2**levelmin)**3)

   if(verbose) write(*,*) 'The average srouce is',test_srcbar

   sf_src_mean = test_srcbar

end subroutine gauss_seidel_sf_src_mean_extradof
