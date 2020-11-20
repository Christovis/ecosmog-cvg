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
  logical,save :: first_step=.true.

  ! Declare test parameters
  !    general
  !real(dp) :: test_rho,test_rhobar
  integer  :: ind,iskip,ngrid,igrid_amr,icell_amr
  integer  :: ix,iy,iz,countc1,countc2             ! BH_test
  real(dp),dimension(1:twotondim,1:ndim) :: xc     ! BH_test
  real(dp) :: dx,Amp                               ! BH_test
  !    one-dimensional
  real(dp) :: pii                                  ! BH_test_1D
  !    spherical overdensity
  real(dp) :: rr,xs,ys,zs,R0                       ! BH_spherical
  real(dp) :: ctilde, twoac2, aomega               ! BH_spherical

  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,999) icount,ilevel

  ! Parameters to set size of spherical overdensity
  ctilde = sol/boxlen_ini/100000.0d0    !GR_test
  twoac2 = 2.0d0*(aexp*ctilde)**2       !GR_test
  aomega = aexp*omega_m

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
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout))then
        if(.not.ok_defrag)then
           call defrag
        endif
        ok_output=.true.
        call dump_all

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
#ifndef OLDINTERP
     call save_phi_old(ilevel)
#endif
     call rho_fine(ilevel,icount)
#ifndef OLDINTERP
     if(extradof) call save_extradof_old(ilevel)
#endif
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

     !when there is no old potential...
#ifndef OLDINTERP
     if (nstep==0) call save_phi_old(ilevel)
#endif

     !--------------------------------------------------------------------------
     ! Start of Test chi-distribution
     !--------------------------------------------------------------------------
     ! Compute the chi & B and its gradient of chi
     if(extradof) then  ! if extradof is True
        !   call multigrid_fine_extradof(ilevel,icount,1)
        !   ! get gradient of sf for:
        !   !   (1) fifth force calculation
        !   !   (2) interpolation to find sf for finer level
        !   call sf_grad_fine_extradof(ilevel,icount)

        ! initialise test distribution of longitudinal/chi-mode cv-Galileon
        ngrid=active(ilevel)%ngrid
        do ind=1,twotondim  ! Loop over cells
           iskip = ncoarse + (ind-1)*ngridmax

           ! Set position of cell centers relative to grid center
           iz = (ind-1)/4         ! BH_test
           iy = (ind-1-4*iz)/2    ! BH_test
           ix =  ind-1-4*iz-2*iy  ! BH_test
           dx=0.5d0**ilevel       ! BH_test
           if(ndim>0) xc(ind,1) = (dble(ix)-0.5d0)*dx  ! BH_test
           if(ndim>1) xc(ind,2) = (dble(iy)-0.5d0)*dx  ! BH_test
           if(ndim>2) xc(ind,3) = (dble(iz)-0.5d0)*dx  ! BH_test
           pii   =  4.0d0*datan(1.0d0)  ! BH_1D

           do i=1,ngrid  ! Loop over active grids
              igrid_amr = active(ilevel)%igrid(i)
              icell_amr = igrid_amr + iskip

              do idim=1,ndim    ! BH_test
                 if (idim.eq.1)  xs = xg(igrid_amr,idim) + xc(ind,idim)
                 if (idim.eq.2)  ys = xg(igrid_amr,idim) + xc(ind,idim)
                 if (idim.eq.3)  zs = xg(igrid_amr,idim) + xc(ind,idim)
              end do            ! BH_test

              ! Homogeneous field test of longitudinal/chi-mode cv-Galileon
              !sf(icell_amr) = 1.0d0

              ! one-dimensional cosine test of longitudinal/chi-mode cv-Galileon
              sf(icell_amr) = Amp*dcos(2.0d0*pii*xs)

              ! spherically symmetric top-hat test of longitudinal/chi-mode cv-Galileon
              rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
              R0=0.01d0
              Amp=1.0d0
              if(rr<0.004D0) then
                 sf(icell_amr) = Amp
              end if

              ! spherically symmetric gaussian test of longitudinal/chi-mode cv-Galileon
              rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
              Amp=1.0d0  ! sigma
              sf(icell_amr) = 1/(Amp*dsqrt(2*pii))*dexp(-0.5*rr/Amp)**2

           end do ! Loop over grids
        end do ! Loop over cells

        ! share initial distribution of longitudinal/chi-mode with all cpu's
        call make_virtual_fine_dp(sf(1)  ,ilevel)

#ifnd   ef OLDINTERP
           if(nstep==0) call save_extradof_old(ilevel)
#endi   f

           if(extradof4) then
              ! compute transverse/B-mode of cv-Galileon based on longitudinal/chi-mode
              call multigrid_fine_extradof(ilevel,icount,2)
              call multigrid_fine_extradof(ilevel,icount,3)
              call multigrid_fine_extradof(ilevel,icount,4)

              ! TODO for the b-mode, use write(*,*) to check b-mode
              !if(verbose) then
              !   write(*,*) "The transverse/B-mode components are: "
              !   write(*,'(3(f16.10,2x))') cbf(1), cbf(2), cbf(3)
              !end if
           endif


        ! Write transverse/B-mode to files to compare with analytic solution
        ngrid = active(ilevel)%ngrid

        do ind=1,twotondim  ! Loop over cells
           iskip = ncoarse + (ind-1)*ngridmax

           ! Set position of cell centers relative to grid center
           iz = (ind-1)/4        ! BH_test
           iy = (ind-1-4*iz)/2   ! BH_test
           ix =  ind-1-4*iz-2*iy ! BH_test
           dx=0.5d0**ilevel      ! BH_test
           if(ndim>0) xc(ind,1) = (dble(ix)-0.5d0)*dx ! BH_test
           if(ndim>1) xc(ind,2) = (dble(iy)-0.5d0)*dx ! BH_test
           if(ndim>2) xc(ind,3) = (dble(iz)-0.5d0)*dx ! BH_test

           do i=1,ngrid  ! Loop over active grids
              igrid_amr = active(ilevel)%igrid(i)
              icell_amr = igrid_amr + iskip

              do idim=1,ndim    ! BH_test
                 if (idim.eq.1)  xs = xg(igrid_amr,idim) + xc(ind,idim)
                 if (idim.eq.2)  ys = xg(igrid_amr,idim) + xc(ind,idim)
                 if (idim.eq.3)  zs = xg(igrid_amr,idim) + xc(ind,idim)
              end do            ! BH_test

              ! Circular plots
              !rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
              !if(rr<0.004d0) then
              !   write(*,'(6(f12.6,2x))') myid, xs, ys, zs, rr, gr_pot(icell_amr,4)
              !end if

              ! Spherically symmetric test
              if(xs>0.5d0.and.zs>0.5D0.and.zs<0.5D0+1.0D0/256.0D0.and.ys>0.5D0.and.ys<0.5D0+1.0D0/256.0D0) then
                 write(*,'(8(f16.10,2x))') xs, ys, zs, sf(icell_amr), src_mean, cbf(icell_amr,2), cbf(icell_amr,3), cbf(icell_amr,4)
              end if

              !! 1D test
              !if(zs>0.5D0.and.zs<0.5D0+1.0D0/256.0D0.and.ys>0.5D0.and.ys<0.5D0+1.0D0/256.0D0) then
              !   write(*,'(6(f16.10,2x))') xs, ys, zs, gr_pot(icell_amr,6), src_mean, gr_pot(icell_amr,2)
              !end if

           end do ! Loop over grids
        end do ! Loop over cells

     endif  ! if extradof is True
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
        if(sink) call update_sink(ilevel)
     end if
  else
     call update_time(ilevel)
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
