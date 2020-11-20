!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine sf_grad_fine_extradof(ilevel,icount)

  use amr_commons
  use pm_commons
  use poisson_commons
  use extradof_commons
  use extradof_parameters

  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer :: ilevel,icount

  !----------------------------------------------------------
  ! Covering routine to compute the gradient of the extradof
  !----------------------------------------------------------

  integer  :: igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer  :: idim
  real(dp) :: dx
  real(dp),dimension(1:twotondim,1:3) :: xc

  integer ,dimension(1:nvector),save :: ind_grid,ind_cell,ind_cell_father
 
  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,111) ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     call gradient_extradof(ind_grid,ngrid,ilevel,icount)
  end do

  do idim=1,ndim
     call make_virtual_fine_dp(sf_grad(1,idim),ilevel)
  end do
 
111 format('   Entering sf_grad_fine for level ',I2)

end subroutine sf_grad_fine_extradof

!#########################################################
!#########################################################
!#########################################################
!#########################################################

subroutine gradient_extradof(ind_grid,ngrid,ilevel,icount)

  use amr_commons
  use pm_commons
  use poisson_commons
  use extradof_commons
  use extradof_parameters

  implicit none

  integer::ngrid,ilevel,icount           ! number of grids, level #
  integer,dimension(1:nvector)::ind_grid ! grid index

  !----------------------------------------------------------------
  ! This routine compute the 3-grad of the extradof for all cells
  ! in grids ind_grid(:) at level ilevel, using five-points FDA
  !----------------------------------------------------------------

  integer::i,idim,ind,iskip
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::sf1,sf2,sf3,sf4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::sf_left,sf_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)  = nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)  = nbor(ind_grid(i),2*idim)
        igridn(i,2*idim-1) = son(ind_left (i,idim))
        igridn(i,2*idim)   = son(ind_right(i,idim))
     end do
  end do

  ! Interpolate extradof from upper level
  do idim=1,ndim
     call interpol_extradof(ind_left (1,idim),sf_left (1,1,idim),ngrid,ilevel,icount)
     call interpol_extradof(ind_right(1,idim),sf_right(1,1,idim),ngrid,ilevel,icount)
  end do

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     do idim=1,ndim
        ! Loop over nodes
        ! ig1 is the relative position of the neighbouring grid which contains the
        !        first left cell whose extradof will be used in the 5-point FDA
        ! ig2 is the relative position of the neighbouring grid which contains the
        !        first right cell whose extradof will be used in the 5-point FDA
        ! ig3 is the relative position of the neighbouring grid which contains the
        !        second left cell whose extradof will be used in the 5-point FDA
        ! ig4 is the relative position of the neighbouring grid which contains the
        !        second right cell whose extradof will be used in the 5-point FDA
        ! id1 is the position (1-8) of the first left cell in its father grid
        ! id2 is the position (1-8) of the first right cell in its father grid
        ! id3 is the position (1-8) of the second left cell in its father grid
        ! id4 is the position (1-8) of the second right cell in its father grid
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax
        
        ! Gather extradof
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              sf1(i)=sf(igridn(i,ig1)+ih1)
           else
              sf1(i)=sf_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              sf2(i)=sf(igridn(i,ig2)+ih2)
           else
              sf2(i)=sf_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              sf3(i)=sf(igridn(i,ig3)+ih3)
           else
              sf3(i)=sf_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              sf4(i)=sf(igridn(i,ig4)+ih4)
           else
              sf4(i)=sf_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           sf_grad(ind_cell(i),idim)=a*(sf1(i)-sf2(i))-b*(sf3(i)-sf4(i)) ! 5-point FDA
        end do
     end do
  end do

end subroutine gradient_extradof
