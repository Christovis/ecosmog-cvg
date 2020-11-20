recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use gr_commons
  use gr_parameters
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::mpi_err
#endif
  integer::ilevel,icount
  integer::igr,igrp,igrm
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!

  integer::i,idim,ivar
  logical::ok_defrag,output_now_all
  logical,save::first_step=.true.
  integer,dimension(1:5) :: gr_ord1
  integer,dimension(1:6) :: gr_ord2

  ! Declare test parameters
!  real(dp) :: test_rho,test_rhobar
  integer  :: ind,iskip,ngrid,igrid_amr,icell_amr

  integer  :: ix,iy,iz,countc1,countc2             ! BH_test
  real(dp),dimension(1:twotondim,1:ndim) :: xc    ! BH_test
  real(dp) :: dx                                  ! BH_test
!  real(dp) :: yy,zz                               ! BH_test
  real(dp) :: rr,xs,ys,zs,R0,Amp,Amp2                   ! BH_spherical
  real(dp) :: pii                             ! BH_test_1D
  real(dp) :: ctilde, twoac2, aomega
  real(dp) :: rng

  ctilde = sol/boxlen_ini/100000.0d0    !GR_test
  twoac2 = 2.0d0*(aexp*ctilde)**2       !GR_test
  aomega = aexp*omega_m

  ! Specific ordering of gr_pots for synchro and move steps
  gr_ord1(1:5)=(/9,8,7,5,6   /) ! Smallest to largest (sync)
  gr_ord2(1:6)=(/6,5,7,8,9,10/) ! Largest to smallest (move)

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
                               call timer('refine','start')
  if(levelmin.lt.nlevelmax .and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
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
                 if(momentum_feedback)call make_virtual_fine_dp(pstarold(1),i)
                 if(simple_boundary)call make_boundary_hydro(i)
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

              if(gr)then
                 do igrp=1,10
                    call make_virtual_fine_dp(gr_pot(1,igrp),i)
                 end do
                 do igrm=1,4
                    call make_virtual_fine_dp(gr_mat(1,igrm),i)
                 end do
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
                               call timer('load balance','start')
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              if(nrestart.eq.nrestart_quad) restart_remap=.true.
              if(restart_remap) then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
#if NDIM==3
                               call timer('sinks','start')
  if(sink)call update_cloud(ilevel)
#endif
  !-----------------
  ! Particle leakage
  !-----------------
                               call timer('particles','start')
  if(pic)call make_tree_fine(ilevel)

  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then

#ifdef WITHOUTMPI
     output_now_all = output_now
#else
     ! check if any of the processes received a signal for output
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(output_now,output_now_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_err)
#endif
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now_all.EQV..true.)then
                               call timer('io','start')
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        ! CAREFUL: create_output is used to destinguish between the case where
        ! the clumpfinder is called from create_sink or directly from amr_step.
#if NDIM==3
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)
#endif

        ! Dump lightcone
        if(lightcone .and. ndim==3) call output_cone()

        if (output_now_all.EQV..true.) then
          output_now=.false.
        endif

     endif

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
                               call timer('movie','start')
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
                               call timer('feedback','start')
     if(hydro.and.star.and.eta_sn>0.and.f_w>0)call kinetic_feedback

  endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('poisson','start')
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
                               call timer('rho','start')
     if(gr) gr2=.true.
     if(gr.and.gr_newtonian) gr2=.false.
     call rho_fine(ilevel,icount)

     ! Newtonian synchro
     if(gr.and.gr_newtonian)then
        ! Compute gravitational potential
        if(ilevel>levelmin)then
           if(ilevel .ge. cg_levelmin) then
              call phi_fine_cg(ilevel,icount)
           else
              call multigrid_fine(ilevel,icount)
           end if
        else
           call multigrid_fine(levelmin,icount)
        end if

        ! Compute gravitational acceleration
        call force_fine(ilevel,icount)

        ! Synchronize remaining particles for gravity
        if(pic)then
                                  call timer('particles','start')
           if(static_dm.or.static_stars)then
              call synchro_fine_static(ilevel)
           else
              call synchro_fine(ilevel)
           end if
        end if
        gr2=.true.
        call rho_fine(ilevel,icount)
        ! Synchronize remaining particles for gravity
        if(pic)then
                                  call timer('particles','start')
           if(static_dm.or.static_stars)then
              call synchro_fine_static(ilevel)
           else
              call synchro_fine(ilevel)
           end if
        end if
     end if
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
                               call timer('particles','start')
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
                               call timer('poisson','start')

     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

     if(.not.gr)then
        ! Compute gravitational potential
        if(ilevel>levelmin)then
           if(ilevel .ge. cg_levelmin) then
              call phi_fine_cg(ilevel,icount)
           else
              call multigrid_fine(ilevel,icount)
           end if
        else
           call multigrid_fine(levelmin,icount)
        end if
        !when there is no old potential...
        if (nstep==0)call save_phi_old(ilevel)

        ! Compute gravitational acceleration
        call force_fine(ilevel,icount)

        ! Synchronize remaining particles for gravity
        if(pic)then
                                  call timer('particles','start')
           if(static_dm.or.static_stars)then
              call synchro_fine_static(ilevel)
           else call synchro_fine(ilevel)
           end if
        end if

     else  ! if gr
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!! GR Test block below !!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do igr=1,6 ! don't need
           if(igr>1.and.igr<6) cycle
           ngrid=active(ilevel)%ngrid

!           countc1=0
!           countc2=0
           ! Loop over cells
           do ind=1,twotondim
              iskip = ncoarse + (ind-1)*ngridmax

              ! Set position of cell centers relative to grid center
              iz = (ind-1)/4        ! BH_test
              iy = (ind-1-4*iz)/2   ! BH_test
              ix =  ind-1-4*iz-2*iy ! BH_test
              dx=0.5d0**ilevel      ! BH_test
              if(ndim>0) xc(ind,1) = (dble(ix)-0.5d0)*dx ! BH_test
              if(ndim>1) xc(ind,2) = (dble(iy)-0.5d0)*dx ! BH_test
              if(ndim>2) xc(ind,3) = (dble(iz)-0.5d0)*dx ! BH_test

              pii   =  4.0d0*datan(1.0d0)         ! BH_1D

              ! Loop over active grids
              do i=1,ngrid
                 igrid_amr = active(ilevel)%igrid(i)
                 icell_amr = igrid_amr + iskip

                 do idim=1,ndim    ! BH_test
                    if (idim.eq.1)  xs = xg(igrid_amr,idim) + xc(ind,idim)
                    if (idim.eq.2)  ys = xg(igrid_amr,idim) + xc(ind,idim)
                    if (idim.eq.3)  zs = xg(igrid_amr,idim) + xc(ind,idim)
                 end do            ! BH_test

                 ! Fix Psi (appearing as source)
!                 gr_pot(icell_amr,5)=0.0d0

!                 ! Homogeneous density field test
!                 rho(icell_amr) = 1.0d0

                 ! Sine density field test  (in my case this is chi=Adcos(2pix))
                 rho(icell_amr) = dcos(2.0d0*pii*xs)            ! Aij test

                 ! INVERSE APPROACH SOURCE
!                 Amp=1.0d0
!                 rho(icell_amr)=(3.0d0*twoac2*aomega-Amp*(8.0d0*pii**2*twoac2+15.0d0*aomega)*dsin(2.0d0*pii*xs))/3.0d0/aomega/(twoac2+Amp*dsin(2.0d0*pii*xs))

                 ! SPHERICALLY SYMMETRIC TESTS
                 rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
                 R0=0.01d0
                 Amp=1.0d0

!                 ! Set imperfect spherical tophat source
!                 rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
!                 if(rr<0.004D0) then
!                    gr_mat(icell_amr,1)= 256.0d0**3 ! Mind the grid size
!                 end if

                  ! Spherically symmetric sources -- INVERSE APPROACH!
!                 rho(icell_amr)=(8.0d0*twoac2*Amp*rr**2-12.0d0*twoac2*Amp*R0+3.0d0*(twoac2*dexp(rr**2/R0)-5.0d0*Amp)*aomega*R0**2)/(3.0d0*aomega*R0**2*(twoac2*dexp(rr**2/R0)+Amp))
!

                 ! Store initial value for latter comparison
!                 gr_pot(icell_amr,2)=gr_pot(icell_amr,6)

              end do ! Loop over grids
           end do ! Loop over cells

           !call make_virtual_fine_dp(gr_pot(1,6),ilevel)
           !call make_virtual_fine_dp(gr_pot(1,5),ilevel)
           !call make_virtual_fine_dp(gr_pot(1,2),ilevel)
           !call make_virtual_fine_dp(gr_mat(1,4),ilevel)
           !call make_virtual_fine_dp(gr_mat(1,1),ilevel)
           !call make_virtual_fine_dp(rho   (1)  ,ilevel)
           call make_virtual_fine_dp(sf(1)  ,ilevel)

           ! Solve one of the GR fields using the previous info
           call multigrid_fine_gr(ilevel,icount,igr)
           call multigrid_fine_extradof(ilevel,icount,2)  ! last place for isf

           ! Write & Exit condition
           if(igr==6) then
              ngrid=active(ilevel)%ngrid
              ! Loop over cells
              do ind=1,twotondim
                 iskip = ncoarse + (ind-1)*ngridmax

                 iz = (ind-1)/4        ! BH_test
                 iy = (ind-1-4*iz)/2   ! BH_test
                 ix =  ind-1-4*iz-2*iy ! BH_test
                 dx=0.5d0**ilevel      ! BH_test
                 if(ndim>0) xc(ind,1) = (dble(ix)-0.5d0)*dx ! BH_test
                 if(ndim>1) xc(ind,2) = (dble(iy)-0.5d0)*dx ! BH_test
                 if(ndim>2) xc(ind,3) = (dble(iz)-0.5d0)*dx ! BH_test

                 ! Loop over active grids
                 do i=1,ngrid
                    igrid_amr = active(ilevel)%igrid(i)
                    icell_amr = igrid_amr + iskip

                    do idim=1,ndim    ! BH_test
                       if (idim.eq.1)  xs = xg(igrid_amr,idim) + xc(ind,idim)
                       if (idim.eq.2)  ys = xg(igrid_amr,idim) + xc(ind,idim)
                       if (idim.eq.3)  zs = xg(igrid_amr,idim) + xc(ind,idim)
                    end do            ! BH_test

                    ! Circular plots
!                    rr = dsqrt((xs-0.5d0)**2+(ys-0.5d0)**2+(zs-0.5d0)**2)
!                    if(rr<0.004d0) then
!                       write(*,'(6(f12.6,2x))') myid, xs, ys, zs, rr, gr_pot(icell_amr,4)
!                    end if

                    ! Spherically symmetric test
                    if(xs>0.5d0.and.zs>0.5D0.and.zs<0.5D0+1.0D0/256.0D0.and.ys>0.5D0.and.ys<0.5D0+1.0D0/256.0D0) then
                       write(*,'(6(f16.10,2x))') xs, ys, zs, gr_pot(icell_amr,6), src_mean, gr_pot(icell_amr,2)
                    end if

!                    ! 1D test
!                    if(zs>0.5D0.and.zs<0.5D0+1.0D0/256.0D0.and.ys>0.5D0.and.ys<0.5D0+1.0D0/256.0D0) then
!                       write(*,'(6(f16.10,2x))') xs, ys, zs, gr_pot(icell_amr,6), src_mean, gr_pot(icell_amr,2)
!                    end if

                 end do ! Loop over grids
              end do ! Loop over cells

!              write(*,*) 'myid =',myid, 'countc1 =',countc1, 'src_mean =',src_mean
              if (myid==1) write(*,*) 'GR test finished for non-linear GR Solver. Now exiting.'
              call clean_stop
           end if ! Write condition

        end do ! igr loop

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!! End GR test block !!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Use GR fields from smallest to largest for synchro
        do i=1,5
           igrp = gr_ord1(i)
           ! Compute force contribution from gr_pot
           call force_fine_gr(ilevel,icount,igrp)

           ! Synchronize remaining particles for gravity
           if(pic)then
                                     call timer('particles','start')
              if(static_dm.or.static_stars)then
                 call synchro_fine_static(ilevel)
              else
                 call synchro_fine_gr(ilevel,igrp)
              end if
           end if
        end do
     end if
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!! End GR test block !!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(hydro)then
                                  call timer('poisson','start')

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

        ! Compute Bondi-Hoyle accretion parameters
#if NDIM==3
                               call timer('sinks','start')
        if(sink)call collect_acczone_avg(ilevel)
#endif
     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
                               call timer('radiative transfer','start')
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
                               call timer('radiative transfer','start')
  if(rt)call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
#if NDIM==3
        if(sink)call update_sink(ilevel)
#endif
     end if
  else
     call update_time(ilevel)
#if NDIM==3
     if(sink)call update_sink(ilevel)
#endif
  end if

  ! Thermal feedback from stars
#if NDIM==3
                               call timer('feedback','start')
  if(hydro.and.star.and.eta_sn>0)call thermal_feedback(ilevel)
#endif

  ! Density threshold or Bondi accretion onto sink particle
#if NDIM==3
  if(sink)then
                               call timer('sinks','start')
     call grow_sink(ilevel,.false.)
  end if
#endif
  !-----------
  ! Hydro step
  !-----------
  if((hydro).and.(.not.static_gas))then

     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
     call godunov_fine(ilevel)

     ! Reverse update boundaries
                               call timer('hydro - rev ghostzones','start')
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
     if(momentum_feedback)then
        call make_virtual_reverse_dp(pstarnew(1),ilevel)
     endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step
                               call timer('poisson','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
                               call timer('hydro upload fine','start')
     call upload_fine(ilevel)

  endif

  !---------------------
  ! Do RT/Chemistry step
  !---------------------
#ifdef RT
  if(rt .and. rt_advect) then
                               call timer('radiative transfer','start')
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
                               call timer('cooling','start')
     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
                               call timer('radiative transfer','start')
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(aexp)
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
                               call timer('cooling','start')
     if(cosmo) call update_coolrates_tables(dble(aexp))
                               call timer('radiative transfer','start')
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
                               call timer('cooling','start')
  if((hydro).and.(.not.static_gas)) then
    if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
  endif
#endif

  !---------------
  ! Move particles
  !---------------
  if(.not.gr)then
     if(pic)then
                                  call timer('particles','start')
        if(static_dm.or.static_stars)then
           call move_fine_static(ilevel) ! Only remaining particles
        else
           call move_fine(ilevel)        ! Only remaining particles
        end if
     end if
  else
     ! Use GR fields from largest to smallest for move
     do i=1,6
        igrp = gr_ord2(i)
        ! Compute force contribution from gr_pot
        call force_fine_gr(ilevel,icount,igrp)
        if(pic)then
                                  call timer('particles','start')
           if(static_dm.or.static_stars)then
              call move_fine_static(ilevel)      ! Only remaining particles
           else
              call move_fine_gr(ilevel,igrp)     ! Only remaining particles
           end if
        end if
     end do
  end if

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
#if NDIM==3
                               call timer('feedback','start')
  if(hydro.and.star.and.(.not.static_gas))call star_formation(ilevel)
#endif
  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if((hydro).and.(.not.static_gas))then
                               call timer('hydro - ghostzones','start')
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
     if(momentum_feedback)call make_virtual_fine_dp(pstarold(1),ilevel)
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
  if((hydro).and.(.not.static_gas))then
     if(eta_mag>0d0.and.ilevel==levelmin)then
                               call timer('hydro - diffusion','start')
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)) call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
                               call timer('particles','start')
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
                               call timer('aton','start')
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
                               call timer('sinks','start')
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
#if NDIM==3
     if(ilevel==levelmin)call create_sink
#endif
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

#ifdef RT
subroutine rt_step(ilevel)
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, t, dtnew, myid
  use rt_cooling_module, only: update_UVrates
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel

!--------------------------------------------------------------------------
!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries.
!--------------------------------------------------------------------------

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=t ; t=t-t_left

  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2.0**(ilevel-levelmin))
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     ! If (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel
     if (i_substep > 1) call rt_set_unew(ilevel)

     if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

                               call timer('cooling','start')
     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
                               call timer('radiative transfer','start')

     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

  !900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

end subroutine rt_step
#endif
