subroutine adaptive_loop
  use amr_commons
  use amr_parameters
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use extradof_commons
  use extradof_parameters
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer      :: ilevel,idim,ivar,info
  real(kind=8) :: tt1,tt2
  real(kind=4) :: real_mem,real_mem_tot

  real(dp)     :: EE                       ! used if extradof2=T

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif


  call init_amr                            ! Initialise AMR variables
  call init_time                           ! Initialise time variables
  if(hydro) call init_hydro                ! Initialise hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       call rt_init_hydro                  ! Initialise radiation variables
#endif
  if(poisson) call init_poisson_extradof   ! Initialise poisson and extradof variables
#ifdef ATON
  if(aton)call init_radiation              ! Initialise radiation variables
#endif
  if(nrestart==0) call init_refine         ! Build initial AMR grid
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))          ! Initialise cooling look up table
  if(pic) call init_part                   ! Initialise particle variables
  if(pic) call init_tree                   ! Initialise particle tree
  if(nrestart==0) call init_refine_2       ! Build initial AMR grid again

#ifndef WITHOUTMPI
  tt2 = MPI_WTIME(info)
  if(myid==1) write(*,*) 'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*) 'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0) write(*,999) ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1) write(*,*) 'Starting time integration' 

  if(extradof2) then
     open(unit=1234, file='background.txt')
  end if

  do

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif
     
     if(extradof2) then
        ! cv-Galileon background (analytical solution on the tracker)
        EE       = dsqrt(0.5d0*(omega_m/aexp**3+dsqrt(omega_m**2/aexp**6+4.0d0*(1.0d0-omega_m))))
        rc_cvg   = EE**2*(param_b3/(2.0d0*(1.0d0-omega_m)))**(2.0d0/3.0d0)

        alpha_cvg  = 0.5d0*(0.5d0*param_b3/(1.0d0-omega_m))**(1.0d0/3.0d0) &
                   * (dsqrt(omega_m**2/aexp**6+4.0d0*(1.0d0-omega_m))-omega_m/aexp**3)

        beta_cvg   = 0.5d0*(0.5d0*param_b3/(1.0d0-omega_m))**(1.0d0/3.0d0) &
                   * (5.0d0*omega_m/aexp**3+3.0d0*omega_m**2/aexp**6/dsqrt(omega_m**2/aexp**6+4.0d0*(1.0d0-omega_m))) &
                   + param_b3
        
        write(1234,'(100e15.6)') aexp,rc_cvg,alpha_cvg,beta_cvg
     end if

     if(verbose) write(*,*) 'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy

#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax .and..not.static)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro) then
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
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if

           if(extradof)then
              call make_virtual_fine_dp(sf(1),ilevel)
              call make_virtual_fine_dp(sf_src(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(sf_grad(1,idim),ilevel)
              end do
              if(extradof4)then
                 call make_virtual_fine_dp(sf_lp(1),ilevel)
                 do idim=1,ndim
                    call make_virtual_fine_dp(cbf(1,idim),ilevel)
                 end do
              end if
           end if
           if(ilevel<levelmin) call refine_fine(ilevel)
        end do
     end if

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax .and..not. static) then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro) then
              call upload_fine(ilevel)
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
              if(simple_boundary) call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt) then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary) call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson) then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
           ! vector-field book-keeping
           if(extradof) then
              call make_virtual_fine_dp(sf(1),ilevel)
              call make_virtual_fine_dp(sf_src(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(sf_grad(1,idim),ilevel)
              end do
              if(extradof4)then
                 call make_virtual_fine_dp(sf_lp(1),ilevel)
                 do idim=1,ndim
                    call make_virtual_fine_dp(cbf(1,idim),ilevel)
                 end do
              end if
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     end if

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0) then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1) then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#endif

  end do
  
  if(extradof2) then
     close(1234)
  endif

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
