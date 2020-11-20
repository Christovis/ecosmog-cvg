
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_ray_grid_refinement(i_arg)
  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  use pm_commons ! for debuging
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in) :: i_arg

  integer  :: old_ind,old_grid,old_cell,old_level
  integer  :: new_ind,new_grid,new_cell,new_level
  integer  ::iray
  integer  :: ix,iy,iz
  integer  :: new_cpu,old_cpu
  real(dp) :: dx,x_g,y_g,z_g,x_r,y_r,z_r
  logical  :: out_of_bound                                         ! for test use only

#ifndef WITHOUTMPI
  integer,dimension(ncpu)                 :: nrays_comm
#endif

  !============================================
  ! determine the size of communication commons
  !============================================

  nrays_comm(1:ncpu) = 0
  
  do iray=1,ray_nrays

     old_grid  =     ray_grid(iray,1)                       ! index of old grid (local)
     old_cpu   =     ray_grid(iray,2)/1000
     old_ind   =    (ray_grid(iray,2)-old_cpu*1000)/100     ! index of old cell relative to old grid
     old_level = mod(ray_grid(iray,2)-old_cpu*1000 ,100)    ! old level #
     old_cell  = ncoarse+(old_ind-1)*ngridmax+old_grid      ! index of old cell (local)

     if(old_grid.eq.0) then
        write(*,*) 'check_ray: old_grid index incorrect; please check.'
        stop
     end if

     ! Grid has been destroyed
     if(father(old_grid).eq.0) then

!       if(ray_step_test) write(*,*) "check_ray: grid has been destroyed"
        
        new_level = old_level-1                             ! coarse level #
        new_cell  = ex_father(old_grid)                     ! index of carse cell (local);
                                                            ! ex_father stores ex father cell index of the killed grid
        new_ind   = (new_cell-ncoarse-1)/ngridmax+1         ! index of coarse cell relative to coarse grid (local)
        new_grid  = new_cell-(ncoarse+(new_ind-1)*ngridmax) ! index of coarse grid (local)

        new_cpu   = cpu_map(new_cell)                       ! coarse grid CPU #

!       write(*,*) "check_ray1",ray_id(iray),myid,new_cpu,old_level,new_level,old_grid,new_grid
        
        if(new_grid.eq.0) then
           write(*,*) 'check_ray: new_grid index out of bound; please check (after destroying grid).'
           stop
        end if
        if(new_level.le.0 .or. new_cpu.le.0) then
           write(*,*) 'check_ray: new_level out of bound (1); please check.',new_level,new_cpu
           stop
        end if

     ! Grid has been further refined
     else if(son(old_cell).gt.0) then
        
!       if(ray_step_test) write(*,*) "check_ray: grid has been refined"
        
        new_grid  = son(old_cell)       ! index of find grid (local)
        new_level = old_level+1         ! fine level #
        new_cpu   = cpu_map(old_cell)   ! fine grid CPU #

        call get_ind_from_ray_position(iray,new_level,new_grid,new_ind)

        if(new_grid.eq.0) then
           write(*,*) 'check_ray: new_grid index out of bound (2); please check.'
           stop
        end if

     ! Grid remains unchanged
     else
        new_grid  = old_grid
        new_level = old_level
        new_cpu   = old_cpu
        if(old_cpu.ne.cpu_map(father(old_grid))) new_cpu = cpu_map(father(old_grid))
        new_ind   = old_ind

!       write(*,*) "check_ray3",ray_id(iray),myid,new_cpu,old_level,new_level,old_grid,new_grid
     end if

     if(new_level.le.0 .or. new_cpu.le.0) then
        write(*,*) 'check_ray: new_level out of bound(3); please check.',new_level,new_cpu
        stop
     end if

     ! Treat myid's grids and others' grids differently
     if(new_cpu.eq.myid) then
        ray_grid(iray,1) =  new_grid 
     else
        ray_grid(iray,1) = -new_grid                  ! mark that the new grid belongs to a different CPU
        nrays_comm(new_cpu) = nrays_comm(new_cpu)+1   ! number of rays to be sent to that CPU
     end if
     ray_grid(iray,2) = 1000*new_cpu+100*new_ind+new_level

     ! Test if new grid actually contains the ray
     if(ray_step_test .and. .false.) then
        x_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dcos(ray_coord(iray,3))
        y_r = ray_coord(iray,1)*dsin(ray_coord(iray,2))*dsin(ray_coord(iray,3))
        z_r = ray_coord(iray,1)*dcos(ray_coord(iray,2))
        dx  = 0.5D0**new_level
        out_of_bound = .false.
        if(x_r.lt.xg(new_grid,1)-ray_x_obs-dx .or. x_r.gt.xg(new_grid,1)-ray_x_obs+dx) then
           write(*,*) xg(new_grid,1)-ray_x_obs-dx, x_r, xg(new_grid,1)-ray_x_obs+dx
           out_of_bound = .true.
        endif
        if(y_r.lt.xg(new_grid,2)-ray_y_obs-dx .or. y_r.gt.xg(new_grid,2)-ray_y_obs+dx) then
           write(*,*) xg(new_grid,2)-ray_y_obs-dx, y_r, xg(new_grid,2)-ray_y_obs+dx
           out_of_bound = .true.
        endif
        if(z_r.lt.xg(new_grid,3)-ray_z_obs-dx .or. z_r.gt.xg(new_grid,3)-ray_z_obs+dx) then
           write(*,*) xg(new_grid,3)-ray_z_obs-dx, z_r, xg(new_grid,3)-ray_z_obs+dx
           out_of_bound = .true.
        endif
        if(out_of_bound) then
           write(*,*) 'check_ray: ray outside boundary of new grid; please check.', new_level, new_grid, "(", xg(new_grid,1), xg(new_grid,2), xg(new_grid,3), ")  (", x_r, y_r, z_r, ") (", &
                xg(new_grid,1)-ray_x_obs-dx, xg(new_grid,2)-ray_y_obs-dx, xg(new_grid,3)-ray_z_obs-dx, ") (", &
                xg(new_grid,1)-ray_x_obs+dx, xg(new_grid,2)-ray_y_obs+dx, xg(new_grid,3)-ray_z_obs+dx, ")"
           write(*,*) ray_coord(iray,1), dsin(ray_coord(iray,2)), dcos(ray_coord(iray,3))
           !stop
        end if
     endif

  end do ! loop over rays

  !-----------------
  ! Communicate rays
  !-----------------
  call ray_communicate(i_arg,0,nrays_comm)
  

end subroutine check_ray_grid_refinement

!===========================================
! ray_communicate
!===========================================
subroutine ray_communicate(i_arg,j_arg,nrays_comm)

  use amr_commons
  use amr_parameters
  use ray_commons
  use ray_parameters
  use pm_commons ! for debuging
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer,intent(in) :: i_arg,j_arg
  integer,intent(in),dimension(ncpu)                 :: nrays_comm
  
  integer  :: new_grid,new_level
  integer  :: iray,jray,kray,nrays1,nrays2,nrays,nrays_tot
  integer  :: icpu,icpu2,ig
  integer  :: test_count
  integer  :: iloop  

  integer::info,tag=101
  integer,dimension(ncpu) :: reqsend,reqrecv
#ifndef WITHOUTMPI
  integer,dimension(ncpu)                 :: ig_loc
  integer,dimension(ncpu)                 :: sendbuf,recvbuf
  integer,dimension(MPI_STATUS_SIZE,ncpu) :: statuses
  integer                                 :: countsend,countrecv,ncache
  integer, allocatable,dimension(:,:)     :: tmparr1
  real(dp),allocatable,dimension(:,:)     :: tmparr2
  integer                                 :: cur_size,new_size 
  integer                                 :: chunk_size 
#endif

  !-----------------------------------------------------------------------------
  ! - This routine communicate rays among processors.
  ! - The rays to be communicated must be flagged with negative ray_grid(iray,1)
  ! - The number of rays to be communicated to each cpu should be given in
  !    nrays_comm(new_cpu)
  ! - This routine updates nrays_tot.
  !-----------------------------------------------------------------------------
  
#ifndef WITHOUTMPI

  !============================
  ! build communication commons
  !============================

  sendbuf = 0
  recvbuf = 0  
  ! Allocate emission commons to contain grid indices myid is going to send to other CPUs
  do icpu=1,ncpu
     ncache = nrays_comm(icpu)
     if(ncache>0) then
        allocate(ray_emission(icpu)%igrid (1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%igrid2(1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%ivar  (1:(nrays_comm(icpu))))
        allocate(ray_emission(icpu)%fvar  (1:(nrays_comm(icpu))))
        ray_emission(icpu)%ngrid = nrays_comm(icpu)                   ! number of rays to be sent to icpu
        sendbuf(icpu) = ray_emission(icpu)%ngrid
        ig_loc(icpu)  = 1                                            
     end if
  end do

  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate grid index
  do icpu=1,ncpu
     ray_reception(icpu)%ngrid=recvbuf(icpu)
     ncache = ray_reception(icpu)%ngrid
     if(ncache>0) then
        allocate(ray_reception(icpu)%igrid (1:ncache))
        allocate(ray_reception(icpu)%igrid2(1:ncache))
        allocate(ray_reception(icpu)%ivar  (1:ncache))
        allocate(ray_reception(icpu)%fvar  (1:ncache))
     end if
  end do

  ! Receive grid list    
  countrecv=0
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache>0) then
        countrecv = countrecv+1
        call MPI_IRECV(ray_reception(icpu)%igrid,ncache, &
                       & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission arrays
  nrays2 = ray_nrays                                               ! temporarily store number of rays on myid
  do iray=1,ray_nrays
     ! Consider only rays no longer belonging to myid grids
     if(ray_grid(iray,1).gt.0) cycle

     ! Find the new CPU to which the ray belongs
     icpu = ray_grid(iray,2)/1000
     ! Find the new level to which the ray belongs
     new_level = mod(ray_grid(iray,2)-icpu*1000,100)

     if(reception(icpu,new_level)%ngrid.eq.0) cycle

     do ig=1,reception(icpu,new_level)%ngrid
        if(reception(icpu,new_level)%igrid(ig).eq.(abs(ray_grid(iray,1)))) then
           ray_emission(icpu)%igrid(ig_loc(icpu)) = ig
!          if(myid.eq.1) write(*,*) 'emission:',myid,icpu,ig_loc(icpu),new_level,ig,ray_grid(iray,2),ray_id(iray),ray_emission(2)%ngrid,ray_grid(iray,1)
           ig_loc(icpu) = ig_loc(icpu)+1
           exit
        end if
     end do

     nrays2 = nrays2-1
  end do

  ! Send grid list   
  countsend=0
  do icpu=1,ncpu
     ncache = ray_emission(icpu)%ngrid
     if(ncache>0) then
        countsend = countsend+1
        call MPI_ISEND(ray_emission(icpu)%igrid,ncache,  &
                       & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Update ray number on myid. At the end of this: 
  ! ray_nrays remains as the old number of rays;
  ! nrays1 is the old number of rays plus rays to be received from other CPUs;
  ! nrays2 is the old number of rays minus rays to be sent to other CPUs.
  nrays1 = ray_nrays                                               ! temporarily store number of rays on myid
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid                            ! number of rays received from icpu
     if(ncache.le.0) cycle

     if(icpu.ne.myid) then 
        nrays1 = nrays1+ncache                                     ! update number of rays on myid
     end if
  end do

  !==================
  ! do communications
  !==================

  cur_size = size(ray_id)                                          ! stores the current size of ray-related arrays

  kray = ray_nrays
  call make_virtual_ray_grid_int(ray_grid(:,2))

  chunk_size = 0                                                   
  do icpu=1,ncpu                                                 
     chunk_size = chunk_size+ray_reception(icpu)%ngrid             
  end do                                                          
  chunk_size = chunk_size-(cur_size-ray_nrays)+5                   
  if(chunk_size.lt.100) chunk_size = 100         

  do icpu=1,ncpu

     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray  = kray+1
        new_level = mod(ray_reception(icpu)%ivar(ig),100)
!       if(myid.eq.2) write(*,*) 'reception:',myid,icpu,ig,new_level,ray_reception(icpu)%igrid(ig),ray_reception(icpu)%ivar(ig),ray_reception(1)%ngrid
        new_grid  = emission(icpu,new_level)%igrid(ray_reception(icpu)%igrid(ig))
        if(new_grid.le.0) then
           write(*,*) 'ray_communicate: ray grid incorrect.',myid,icpu,new_level
           stop
        end if
        ! Extend array size if necessary
        if(kray.gt.cur_size) then
           new_size = cur_size+chunk_size

           allocate   (tmparr1(1:cur_size,1:2          ))
           ! <><><><><><><><><><><><><><><><><><><><
           ! United array for ray-tracing quantities               ! Baojiu-Feb-2019
           ! <><><><><><><><><><><><><><><><><><><><
           if(RAY_NQUANTS.lt.3) then
              allocate(tmparr2(1:cur_size,1:3          ))
           else
              allocate(tmparr2(1:cur_size,1:RAY_NQUANTS))
           end if

           ! Array for grid index info
           tmparr1(1:cur_size,1:2) = ray_grid(1:cur_size,1:2)
           deallocate(ray_grid                )
           allocate  (ray_grid(1:new_size,1:2))
           ray_grid = 0
           ray_grid(1:cur_size,1:2) = tmparr1(1:cur_size,1:2)

           ! Array for ray global id
           tmparr1(1:cur_size,1) = ray_id(1:cur_size)
           deallocate(ray_id            )
           allocate  (ray_id(1:new_size))
           ray_id = 0
           ray_id(1:cur_size) = tmparr1(1:cur_size,1)

           ! Array for ray coordinates
           tmparr2(1:cur_size,1:3) = ray_coord(1:cur_size,1:3)
           deallocate(ray_coord                )
           allocate  (ray_coord(1:new_size,1:3))
           ray_coord = 0.0D0
           ray_coord(1:cur_size,1:3) = tmparr2(1:cur_size,1:3)

           ! <><><><><><><><><><><><><><><><><><><><
           ! United array for ray tracing quantities               ! Baojiu-Feb-2019 
           ! <><><><><><><><><><><><><><><><><><><><
           tmparr2(1:cur_size,1:RAY_NQUANTS) = ray_quants(1:cur_size,1:RAY_NQUANTS)     
           deallocate(ray_quants                          )      
           allocate  (ray_quants(1:new_size,1:RAY_NQUANTS))   
           ray_quants = 0.0D0     
           ray_quants(1:cur_size,1:RAY_NQUANTS) = tmparr2(1:cur_size,1:RAY_NQUANTS) 
 
           deallocate(tmparr1)
           deallocate(tmparr2)

           cur_size = new_size
        end if
        ray_grid(kray,1) = new_grid
        ray_grid(kray,2) = ray_reception(icpu)%ivar(ig)
     end do 
  end do

  kray = ray_nrays
  call make_virtual_ray_grid_int(ray_id(:))
  do icpu=1,ncpu
     ncache = ray_reception(icpu)%ngrid
     if(ncache.le.0)  cycle
     if(icpu.eq.myid) cycle
   
     do ig=1,ncache
        kray = kray+1
        ray_id(kray) = ray_reception(icpu)%ivar(ig)
     end do 
  end do

  ! <><><><><><><><><><><><><><><><><><><><
  ! ray_coord components united to one loop                        ! Baojiu-Feb-2019
  ! <><><><><><><><><><><><><><><><><><><><
  do iloop=1,3
     kray = ray_nrays
     call make_virtual_ray_grid_dp(ray_coord(:,iloop))
     do icpu=1,ncpu
        ncache = ray_reception(icpu)%ngrid
        if(ncache.le.0)  cycle
        if(icpu.eq.myid) cycle
   
        do ig=1,ncache
           kray = kray+1
           ray_coord(kray,iloop) = ray_reception(icpu)%fvar(ig)
        end do 
     end do
  end do

  ! <><><><><><><><><><><><><><><><><><><
  ! All ray quantities united to one loop                          ! Baojiu-Feb-2019
  ! <><><><><><><><><><><><><><><><><><><
  if(RAY_NQUANTS>0) then                                         
     do iloop=1,RAY_NQUANTS
        kray = ray_nrays
        call make_virtual_ray_grid_dp(ray_quants(:,iloop))              
        do icpu=1,ncpu                       
           ncache = ray_reception(icpu)%ngrid       
           if(ncache.le.0)  cycle                 
           if(icpu.eq.myid) cycle           
   
           do ig=1,ncache            
              kray = kray+1
              ray_quants(kray,iloop) = ray_reception(icpu)%fvar(ig) 
           end do 
        end do
     end do
  else
     write(*,*) 'check_ray: RAY_NQUANTS should be larger than 0!'
     call clean_stop
  end if                           

  ! Defrag arrays so that old myid rays are listed first
  jray = 0
  do iray=1,nrays1
     if(ray_grid (iray,1).gt.0) then
        jray                = jray+1
        ray_grid (jray,1:2) = ray_grid (iray,1:2)
        ray_id   (jray    ) = ray_id   (iray    )
        ray_coord(jray,1:3) = ray_coord(iray,1:3)
        ! <><><><><><><><><><><><><><><><><><><><
        ! United array for ray-tracing quantities                  ! Baojiu-Feb-2019
        ! <><><><><><><><><><><><><><><><><><><><
        ray_quants(jray,1:RAY_NQUANTS) = ray_quants(iray,1:RAY_NQUANTS)
     endif
  end do

  ray_nrays = jray

#endif

  nrays = ray_nrays

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nrays,nrays_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  nrays_tot=nrays
#endif
  if(nrays_tot.ne.ray_Nx*ray_Ny) then
     write(*,*) 'ray_communicate: incorrect ray number.'
     write(*,*) 'after communication myid ',myid,' has ', nrays1,' rays.'
     if(myid.eq.1) write(*,*) 'there are ',nrays_tot,' rays in total, while it should be ',ray_Nx*ray_Ny,'!'
     stop
  end if
  if(j_arg.eq.0) then
!    write(*,888) i_arg,nrays_tot,myid,nrays
  else
!    write(*,999) nrays_tot,myid,nrays
  endif

  ! Deallocate the ray communication commons
  do icpu=1,ncpu
     if(ray_emission (icpu)%ngrid>0) then
        ray_emission (icpu)%ngrid = 0
        deallocate(ray_emission (icpu)%igrid )
        deallocate(ray_emission (icpu)%igrid2)
        deallocate(ray_emission (icpu)%ivar  )
        deallocate(ray_emission (icpu)%fvar  )
     end if
     if(ray_reception(icpu)%ngrid>0) then
        ray_reception(icpu)%ngrid = 0
        deallocate(ray_reception(icpu)%igrid )
        deallocate(ray_reception(icpu)%igrid2)
        deallocate(ray_reception(icpu)%ivar  )
        deallocate(ray_reception(icpu)%fvar  )
     end if
  end do

888 format('Level',I3,' refined, nrays_tot =',I8,', nrays (myid =',I4,') =',I6)
999 format('Rays have changed CPUs, nrays_tot =',I8,', nrays (myid =',I4,') =',I6)

end subroutine ray_communicate
