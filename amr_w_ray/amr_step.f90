!-------------------------------------------------------------------------------
! This file contains 1 subroutine:
!    amr_step:
!-------------------------------------------------------------------------------

recursive subroutine amr_step(ilevel,icount)
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, do not change it,  !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  use amr_commons
  use amr_parameters
  use pm_commons
  use hydro_commons
  use poisson_commons
  use poisson_parameters
  use extradof_commons
  use extradof_parameters
  use ray_commons
  use ray_parameters
#ifdef RT
  use rt_hydro_commons
  use SED_module
#endif
  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer :: ilevel,icount
  integer  :: i,idim,ivar,info
  integer  :: icell
  logical  :: ok_defrag,ok_output
  real(dp) :: test_rho,test_rhobar
  integer  :: ind,iskip
  logical,save :: first_step=.true.

  !----------------------exact output scale factor-----------------------!
  integer  :: ii 
  !----------------------exact output scale factor-----------------------!
  real(dp) :: ray_a_s                                                                   ! RAY_RAMSES
  real(dp) :: ray_aout                                                                  ! RAY_RAMSES
  real(dp) :: dcom_ray                                                                  ! RAY_RAMSES
  real(dp) :: zexp,chi_wave,coverH0_ray,chi_start                                       ! RAY_RAMSES
  real(dp) :: temp_dist_x,temp_dist_y,temp_dist_z                                       ! RAY_RAMSES
  real(dp) :: pi                                                                        ! RAY_RAMSES
  integer  :: ilevel2                                                                   ! RAY_RAMSES

  if(numbtot(1,ilevel)==0) return

  !----------------------exact output scale factor-----------------------!
  if(cosmo.and.dabs(t_next).lt.1.0D-8)then
     ! Find neighbouring scale factors
     ii=1
     do while(aexp_frw(ii)>aout(iout).and.ii<n_frw)
        ii=ii+1
     end do
     ! Interpolate expansion factor for the next step
     t_next = tau_frw(ii  )*(aout(iout)-aexp_frw(ii-1))/(aexp_frw(ii  )-aexp_frw(ii-1)) + &
            & tau_frw(ii-1)*(aout(iout)-aexp_frw(ii  ))/(aexp_frw(ii-1)-aexp_frw(ii  ))
  end if
  !----------------------exact output scale factor-----------------------!

  if(ray) then                                                                          ! RAY_RAMSES
     ray_a_s = 1.0d0/(1.0d0+ray_z_s)                                                    ! RAY_RAMSES
     if(.not.allocated(aexp_old_ray)) then                                              ! RAY_RAMSES
        allocate(aexp_old_ray(1:nlevelmax))                                             ! RAY_RAMSES
        aexp_old_ray = -1.0D0 ! abs_for_sync                                            ! RAY_RAMSES
        if(ray_multi_out.and.ray_nout.gt.0) then                                        ! RAY_RAMSES       
           ray_iout = 1                                                                 ! RAY_RAMSES       
        end if                                                                          ! RAY_RAMSES       
        allocate(aexp_new_ray(1:nlevelmax))                                             ! RAY_RAMSES
        aexp_new_ray = aexp                                                             ! RAY_RAMSES
     end if                                                                             ! RAY_RAMSES
     pi          = 4.0d0*datan(1.0d0)                                                   ! RAY_RAMSES
     coverH0_ray = 299792.458d0/100.d0  !c/100/h [Mpc/h]                                ! RAY_RAMSES
  end if                                                                                ! RAY_RAMSES

#ifdef OUTPUT_PHIDOT
  if(.not.allocated(aexp_old_ray)) then
     allocate(aexp_old_ray(1:nlevelmax))
     aexp_old_ray = -1.0D0 ! abs_for_sync
     if(ray_multi_out.and.ray_nout.gt.0) then
        ray_iout = 1
     end if
     allocate(aexp_new_ray(1:nlevelmax))
     aexp_new_ray = aexp
  end if
#endif

  if(verbose) write(*,999) icount,ilevel

  !--------------------------------
  ! Initialise the 5th-force arrays
  !--------------------------------
  ! Baojiu-----------------------------18-09-14---------------------------------
  !---- Changed the original "levelmin.lt.nlevelmax" to "levelmin.le.nlevelmax"
  !---- because levelmin=levelmax=nlevelmax is used to do no-refinement simus,
  !---- and in this case we need to make sure initial guesses are properly set.
  !------------------------------------18-09-14---------------------------------
  if(levelmin.le.nlevelmax .and..not. static)then
     if(extradof) then
        ! on domain level: use zero as initial guess
        if(ilevel.eq.levelmin .and. aexp.le.2.0d0) then
           do icell=1,ncoarse+twotondim*ngridmax
              sf(icell)=0.0d0
           end do
           if(extradof4) then
              do icell=1,ncoarse+twotondim*ngridmax
                 cbf(icell,1:3)=0.0d0
              end do
           end if
        end if
     end if
  end if
  !------------------------------------18-09-14---------------------------------

  !----------------------------------------------------------
  ! Make new refinements and update boundaries for everything
  !----------------------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin .or. icount>1) then
        do i=ilevel,nlevelmax
           if(i>levelmin) then

              ! Build communicators
              call build_comm(i)

              ! Update boundaries
              call make_virtual_fine_int(cpu_map(1),i)

              if(hydro) then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif

                 if(simple_boundary) call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if

              if(extradof)then
                 call make_virtual_fine_dp(sf(1),i)
                 call make_virtual_fine_dp(sf_src(1),i)
                 !----------------------------------------------
                 ! Update chi- & B-values on the virtual boundaries:
                 ! Can only be done after communicator is built.
                 ! Need to be done after refinement is created.
                 !----------------------------------------------
                 do idim=1,ndim
                    call make_virtual_fine_dp(sf_grad(1,idim),i)
                 end do

                 if(extradof4) then
                    call make_virtual_fine_dp(sf_lp(1),i)
                    do idim=1,ndim
                       call make_virtual_fine_dp(cbf(1,idim),i)
                    end do
                 end if

              end if

           end if

           ! Refine grids
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin .lt. nlevelmax) then
     if(ilevel==levelmin) then
        if(nremap>0) then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0 .and. first_step) then
              first_step = .false.
           else
              if(MOD(nstep_coarse,nremap)==0) then
                 call load_balance
                 call defrag
                 ok_defrag = .true.
              end if
           end if
        end if
     endif
  end if

  !--------------------------------------
  ! Update sink cloud particle properties
  !--------------------------------------
  if(sink)call update_cloud(ilevel,.false.)

  !-----------------
  ! Particle leakage
  !-----------------
  if(pic)call make_tree_fine(ilevel)

  !------------------------
  ! Output results to files
  !------------------------
  ok_output=.false.
  if(ilevel==levelmin)then
     !----------------------exact output scale factor-----------------------!
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout)-1.0D-8.or.t>=tout(iout))then
     !----------------------exact output scale factor-----------------------!
!    if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout)       .or.t>=tout(iout))then
        if(.not.ok_defrag)then
           call defrag
        endif
        ok_output=.true.
        call dump_all

        !----------------------exact output scale factor-----------------------!
        if(cosmo)then
           ! Find neighbouring scale factors
           ii=1
           do while(aexp_frw(ii)>aout(iout).and.ii<n_frw)
              ii=ii+1
           end do
           ! Interpolate expansion factor for the next step
           t_next = tau_frw(ii  )*(aout(iout)-aexp_frw(ii-1))/(aexp_frw(ii  )-aexp_frw(ii-1)) + &
                  & tau_frw(ii-1)*(aout(iout)-aexp_frw(ii  ))/(aexp_frw(ii-1)-aexp_frw(ii  ))
        end if
        !----------------------exact output scale factor-----------------------!

        if(gas_analytics) call gas_ana

        ! Run the clumpfinder
        ! (produce output, don't keep arrays alive on output)
        ! Note: create_output is used to destinguish between the case where
        ! the clumpfinder is called from create_sink or directly from amr_step.
        if(clumpfind .and. ndim==3) call clump_finder(.true.)

        ! Dump lightcone
        if(lightcone) call output_cone()

     endif

     ! Important can't be done in sink routines because it must be done after dump all
     if(sink)acc_rate=0.

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(aexp>=amovout(imov).or.t>=tmovout(imov))then
        call output_frame()
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
     if(hydro.and.star.and.eta_sn>0.and.f_w>0)call kinetic_feedback

  endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson) then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)

#ifdef OUTPUT_PHIDOT
     if(aexp_old_ray(ilevel).lt.0.0d0 .and. ilevel.ne.levelmin) then                  
        aexp_old_ray(ilevel) = (aexp_old_ray(ilevel-1)+aexp_new_ray(ilevel-1))/2.0d0
     else                                                   
        aexp_old_ray(ilevel) = aexp_new_ray(ilevel)    
     end if                     
#else 
     if(ray) then                                                                       ! RAY_RAMSES
        if(aexp_old_ray(ilevel).lt.0.0d0 .and. ilevel.ne.levelmin) then                 ! RAY_RAMSES
           aexp_old_ray(ilevel) = (aexp_old_ray(ilevel-1)+aexp_new_ray(ilevel-1))/2.0d0 ! RAY_RAMSES
        else                                                                            ! RAY_RAMSES
           aexp_old_ray(ilevel) = aexp_new_ray(ilevel)                                  ! RAY_RAMSES
        end if                                                                          ! RAY_RAMSES
     end if                                                                             ! RAY_RAMSES
#endif

     call rho_fine(ilevel,icount)

     if(extradof) call save_extradof_old(ilevel)
  endif

  if(ilevel==levelmin) then
     test_rho=0.0D0
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           test_rho=test_rho+rho(active(ilevel)%igrid(i)+iskip)
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(test_rho,test_rhobar,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     test_rhobar=test_rho
#endif
     test_rhobar=test_rhobar/dble((2**levelmin)**3)
     if(verbose) write(*,*) 'The average density at a=',aexp,' is',test_rhobar
  end if

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic) then
     ! Remove particles to finer levels
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson) then

     ! Remove gravity source term with half time step and old force
     if(hydro) then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

     ! Compute gravitational potential
     if(ilevel>levelmin) then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if

#ifdef OUTPUT_PHIDOT
     aexp_new_ray(ilevel) = aexp                         
#else
     if(ray) then                                                                       ! RAY_RAMSES
        aexp_new_ray(ilevel) = aexp                                                     ! RAY_RAMSES
     end if                                                                             ! RAY_RAMSES
#endif

     !when there is no old potential...
     if (nstep==0) call save_phi_old(ilevel)

     !--------------------------------------------------------------------------
     ! Compute the chi & B and its gradient of chi
     !--------------------------------------------------------------------------
     if(extradof) then  ! if extradof is True
        call multigrid_fine_extradof(ilevel,icount,1)
        ! get gradient of sf for:
        !   (1) fifth force calculation
        !   (2) interpolation to find sf for finer level
        call sf_grad_fine_extradof(ilevel,icount)

        if(nstep==0) call save_extradof_old(ilevel)

        if(extradof4) then
           ! compute transverse/B-mode of cv-Galileon based on longitudinal/chi-mode
           call multigrid_fine_extradof(ilevel,icount,2)
           call multigrid_fine_extradof(ilevel,icount,3)
           call multigrid_fine_extradof(ilevel,icount,4)
        end if
     end if
     !--------------------------------------------------------------------------
     ! End of Test chi-distribution
     !--------------------------------------------------------------------------

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Thermal feedback from stars
     if(hydro.and.star.and.eta_sn>0)call thermal_feedback(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
        call synchro_fine(ilevel)
     end if

     if(hydro)then
        ! Compute Bondi-Hoyle accretion parameters
        if(sink.and.bondi) call bondi_hoyle(ilevel)

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end if

  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)

  !----------------------exact output scale factor-----------------------!
  if(ilevel.eq.levelmin .and. cosmo) then
     if(t+dtnew(ilevel).gt.t_next .and. t.le.t_next) then
        dtnew(ilevel) = t_next-t
     end if
     if(myid.eq.1) write(*,'(100e15.6)') 1.0,t,dtnew(ilevel),t+dtnew(ilevel),t_next,aout(iout)
  end if
  !----------------------exact output scale factor-----------------------!

  if(ilevel>levelmin) then
     dtnew(ilevel) = MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
  if(hydro) call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
  if(rt) call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax) then
     if(numbtot(1,ilevel+1)>0) then
        if(nsubcycle(ilevel)==2) then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1) = dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1) = dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)

        if(ray.and..not.ray_initialized) then                                           ! RAY_RAMSES
           zexp        = 1.0d0/aexp-1.0d0                                               ! RAY_RAMSES
           chi_wave    = dcom_ray(0.0d0,zexp,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox ! RAY_RAMSES
           temp_dist_z = dabs(ray_z_obs)/ray_Lbox+1.0d0                                 ! RAY_RAMSES
           temp_dist_x = dtan(ray_opening_x/360.0d0*pi)*temp_dist_z                     ! RAY_RAMSES
           temp_dist_y = dtan(ray_opening_y/360.0d0*pi)*temp_dist_z                     ! RAY_RAMSES
           chi_start   = dsqrt(temp_dist_x**2+temp_dist_y**2+temp_dist_z**2)            ! RAY_RAMSES
           ray_chi_s   = dcom_ray(0.0d0,ray_z_s,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox ! RAY_RAMSES
           if(chi_wave.lt.dmin1(chi_start,ray_chi_s)*1.01d0) then                       ! RAY_RAMSES
              ilevel2= levelmin                                                         ! RAY_RAMSES
              call init_ray(ilevel2) ! Initialise the rays                              ! RAY_RAMSES
           end if                                                                       ! RAY_RAMSES
        end if                                                                          ! RAY_RAMSES
        if(ray.and.ray_initialized) then                                                ! RAY_RAMSES
           if(ray_multi_out.and.ray_nout.gt.0) then                                     ! RAY_RAMSES
              ray_aout = 1.0D0/(1.0D0+ray_zout(ray_iout))                               ! RAY_RAMSES
              if(ray_aout.lt.aexp.and.ray_aout.gt.aexp_old_ray(ilevel)) then            ! RAY_RAMSES
                 call ray_step(ilevel,ray_aout)                                         ! RAY_RAMSES
                 ray_iout = ray_iout+1                                                  ! RAY_RAMSES
              end if                                                                    ! RAY_RAMSES
           end if                                                                       ! RAY_RAMSES
           call ray_step(ilevel,-1.0D0)                                                 ! RAY_RAMSES
        end if                                                                          ! RAY_RAMSES

        if(sink) call update_sink(ilevel)
     end if
  else
     call update_time(ilevel)

     if(ray.and..not.ray_initialized) then                                              ! RAY_RAMSES
        zexp        = 1.0d0/aexp-1.0d0                                                  ! RAY_RAMSES
        chi_wave    = dcom_ray(0.0d0,zexp,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox  ! RAY_RAMSES
        temp_dist_z = dabs(ray_z_obs)/ray_Lbox+1.0d0  ! units in box size               ! RAY_RAMSES
        temp_dist_x = dtan(ray_opening_x/360.0d0*pi)*temp_dist_z                        ! RAY_RAMSES
        temp_dist_y = dtan(ray_opening_y/360.0d0*pi)*temp_dist_z                        ! RAY_RAMSES
        chi_start   = dsqrt(temp_dist_x**2+temp_dist_y**2+temp_dist_z**2)               ! RAY_RAMSES
        ray_chi_s   = dcom_ray(0.0d0,ray_z_s,omega_m,omega_l,coverH0_ray,1.0d-6)/ray_Lbox ! RAY_RAMSES
        if(chi_wave.lt.dmin1(chi_start,ray_chi_s)*1.01d0) then                          ! RAY_RAMSES
           ilevel2= levelmin                                                            ! RAY_RAMSES
           call init_ray(ilevel2) ! Initialise the rays                                 ! RAY_RAMSES
        end if                                                                          ! RAY_RAMSES
     end if                                                                             ! RAY_RAMSES
     if(ray.and.ray_initialized) then                                                   ! RAY_RAMSES
        if(ray_multi_out.and.ray_nout.gt.0) then                                        ! RAY_RAMSES
           ray_aout = 1.0D0/(1.0D0+ray_zout(ray_iout))                                  ! RAY_RAMSES
!          write(*,*) 'xxxxx',aexp,aexp_old_ray(ilevel),ray_aout                        ! RAY_RAMSES
           if(ray_aout.lt.aexp.and.ray_aout.gt.aexp_old_ray(ilevel)) then               ! RAY_RAMSES
              call ray_step(ilevel,ray_aout)                                            ! RAY_RAMSES
              ray_iout = ray_iout+1                                                     ! RAY_RAMSES
           end if                                                                       ! RAY_RAMSES
        end if                                                                          ! RAY_RAMSES
        call ray_step(ilevel,-1.0D0)                                                    ! RAY_RAMSES
     end if                                                                             ! RAY_RAMSES

     if(sink) call update_sink(ilevel)
  end if

#ifdef RT
  ! Add stellar radiation sources
  if(rt.and.rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
#endif

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call move_fine(ilevel) ! Only remaining particles
  end if

  !-----------
  ! Hydro step
  !-----------
  if(hydro) then

     ! Hyperbolic solver
     call godunov_fine(ilevel)

     ! Reverse update boundaries
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call set_uold(ilevel)

     ! Density threshold or Bondi accretion onto sink particle
     if(sink) call grow_sink(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step
     if(poisson) call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
     call upload_fine(ilevel)

  endif

#ifdef RT
  !---------------
  ! Radiation step
  !---------------
  if(rt)then
     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     ! Restriction operator
     call rt_upload_fine(ilevel)
  endif
#endif

  !-------------------------------
  ! Source term in leaf cells only
  !-------------------------------
  if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
  if(hydro.and.star)call star_formation(ilevel)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif
#ifdef RT
  if(rt)then
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)
  end if
#endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.static) call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic) call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin) then
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin) then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
     if(ilevel==levelmin) call create_sink
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin) then
     if(nsubcycle(ilevel-1)==1) dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2) dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step
