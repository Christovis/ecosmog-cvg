!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_ray_grid_int(xx)
  use amr_commons
  use ray_commons
  use amr_parameters
  use ray_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu) :: statuses
#endif
  ! -----------------------------------------
  ! This routine communicates rays among CPUs
  ! -----------------------------------------

  integer,dimension(1:ray_nrays) :: xx
  integer                        :: icpu,i,ncache,ind,iray,ilevel
  integer                        :: countsend,countrecv
  integer                        :: info,buf_count,tag=101
  integer,dimension(ncpu)        :: reqsend,reqrecv,ig_loc

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=ray_reception(icpu)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(ray_reception(icpu)%ivar,ncache, &
                     & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  ig_loc(1:ncpu) = 0
  do iray=1,ray_nrays
     ! Considier only rays no longer blonging to myid grids
     if(ray_grid(iray,1).gt.0) cycle
     ! Find the new CPU to which the ray belongs
     icpu = ray_grid(iray,2)/1000 ! Sownak
     ! Fill in ray emission commons
     ig_loc(icpu) = ig_loc(icpu)+1 
     ray_emission(icpu)%ivar(ig_loc(icpu)) = xx(iray)        
!    if(myid.eq.1.and.icpu.eq.2) write(*,*) 'emission:',myid,icpu,ig_loc(icpu),ray_emission(icpu)%ivar(ig_loc(icpu))
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=ray_emission(icpu)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(ray_emission(icpu)%ivar,ncache, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends   
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_ray_grid_int')

end subroutine make_virtual_ray_grid_int


!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_virtual_ray_grid_dp(xx)
  use amr_commons
  use ray_commons
  use amr_parameters
  use ray_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu) :: statuses
#endif
  ! -----------------------------------------
  ! This routine communicates rays among CPUs
  ! -----------------------------------------
  real(dp),dimension(1:ray_nrays) :: xx
  integer                         :: icpu,i,ncache,ind,iray,ilevel
  integer                         :: countsend,countrecv
  integer                         :: info,buf_count,tag=101
  integer,dimension(ncpu)         :: reqsend,reqrecv,ig_loc

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=ray_reception(icpu)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(ray_reception(icpu)%fvar,ncache, &
                     & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  ig_loc(1:ncpu) = 0
  do iray=1,ray_nrays
     ! Considier only rays no longer blonging to myid grids
     if(ray_grid(iray,1).gt.0) cycle
     ! Find the new CPU to which the ray belongs
     icpu = ray_grid(iray,2)/1000 ! Sownak
     ! Fill in ray emission commons
     ig_loc(icpu) = ig_loc(icpu)+1
     ray_emission(icpu)%fvar(ig_loc(icpu)) = xx(iray)
!    if(myid.eq.1) write(*,*) 'emission:',myid,icpu,ray_emission(icpu)%fvar(ig_loc(icpu))
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=ray_emission(icpu)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(ray_emission(icpu)%fvar,ncache, &
                     & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends   
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_ray_grid_dp')

end subroutine make_virtual_ray_grid_dp
