subroutine backup_poisson(filename)
  use amr_commons
  use poisson_commons
  use extradof_commons

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,info
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc

  integer::ix,iy,iz                                  !for the coordinates
  real(dp)::dx                                       !for the coordinates
  real(dp),dimension(1:twotondim,1:3)::xc            !cell coordinates relative to grids

  if(verbose)write(*,*)'Entering backup_poisson'

  ilun=ncpu+myid+10
     
  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  !do ilevel=1,nlevelmax
  do ilevel=levelmin,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              
              ! Write coordinates
              dx=0.5d0**ilevel
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              if(ndim>0)xc(ind,1)=(dble(ix)-0.5d0)*dx
              if(ndim>1)xc(ind,2)=(dble(iy)-0.5d0)*dx
              if(ndim>2)xc(ind,3)=(dble(iz)-0.5d0)*dx
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=xg(ind_grid(i),ivar)+xc(ind,ivar)
                 end do
                 write(ilun)xdp
              end do

              ! Write vector Galileon fields
              do i=1,ncache
                 xdp(i)=sf(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
              !do ivar=1,ndim
              !   do i=1,ncache
              !      xdp(i)=cbf(ind_grid(i)+iskip,ivar)
              !   end do
              !   write(ilun)xdp
              !end do

              ! Write potential
              do i=1,ncache
                 xdp(i)=phi(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
              
              ! Write force
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do

           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)
    
  ! Send the token
!#ifndef WITHOUTMPI
!  if(IOGROUPSIZE>0) then
!     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
!        dummy_io=1
!        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
!             & MPI_COMM_WORLD,info2)
!     end if
!  endif
!#endif

end subroutine backup_poisson


