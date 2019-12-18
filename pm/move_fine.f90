!-------------------------------------------------------------------------------
! This file contains 3 subroutine:
!    1) move_fine: use move1 or move2 to compute particles position
!    2) move1: use CIC to compute force on particle
!    3) move2: use TIC to compute force on particle
!
! Note: the contribution of the fifth-force is in
! synchro_fine.f90 equivalent to move_fine.f90
!-------------------------------------------------------------------------------

subroutine move_fine(ilevel)
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart)
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ip=ip+1
           ind_part(ip)=ipart
           ind_grid_part(ip)=ig
           if(ip==nvector)then
#ifdef TSC
              call move2(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
#else
              call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
#endif
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
#ifdef TSC
  if(ip>0)call move2(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
#else
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
#endif

111 format('   Entering move_fine for level ',I2)

end subroutine move_fine


subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  use amr_commons
  use amr_parameters
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  use extradof_commons
  use extradof_parameters

  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,length,dx_loc,scale,vol_loc,r2
  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: f0,f1,f2

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,&
       & ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=dd(j,idim)
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=dd(j,idim)
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+uold(indp(j,ind),idim+1)*vol(j,ind)
           end do
        end do
     end do
  endif
  if(poisson)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              ! cv-Galileon (add fifth force to the total force if appropriate)
              if(extradof .and. extradof2 .and. .not.extradof3) then
                 ! full case
                 f1         = f(indp(j,ind),idim) + alpha_cvg*sf_grad(indp(j,ind),idim)
                 ff(j,idim) = ff(j,idim) + f1*vol(j,ind)
              else if(.not.extradof .and. extradof2 .and. extradof3) then
                 ! linearized case
                 ! \nabla^2\Phi = \Omega_m*a*\rho*(1 + \alpha/\beta)
                 ! \alpha/\beta is the 5th-force to Newtonian force ratio
                 ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*(1.0d0 + alpha_cvg/beta_cvg)*vol(j,ind)
              else
                 ! LambdaCDM case
                 ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
              end if
           end do
        end do
     end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static.or.tracer)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
        end do
     endif
  end do

  ! For sink cloud particle only
  if(sink)then
     ! Overwrite cloud particle velocity with sink velocity
     do idim=1,ndim
        do j=1,np
           isink=-idp(ind_part(j))
           if(isink>0)then
              new_vp(j,idim)=vsnew(isink,idim,ilevel)
           end if
        end do
     end do
  end if

  ! Store velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Update position
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+new_vp(j,idim)*dtnew(ilevel)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine move1


subroutine move2(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse TSC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then TSC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  use amr_commons
  use amr_parameters
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  use extradof_commons
  use extradof_parameters

  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,length,dx_loc,scale,vol_loc,r2
  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::ff,new_xp,new_vp
  real(dp),dimension(1:nvector,1:ndim),save::x,cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: f0,f1,f2

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0) error=.true.
     end do
  end do
  if(error) then
     write(*,*) 'problem in move2'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0) then
              write(*,*) xp(ind_part(j),1:ndim),x0(ind_grid_part(j),1:ndim),x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! TSC at level ilevel; a particle contributes
  !     to three cells in each dimension
  ! cl: position of leftmost cell centre
  ! cc: position of central cell centre
  ! cr: position of rightmost cell centre
  ! wl: weighting function for leftmost cell
  ! wc: weighting function for central cell
  ! wr: weighting function for rightmost cell
  do idim=1,ndim
     do j=1,np
        cl(j,idim)=dble(int(x(j,idim)))-0.5D0
        cc(j,idim)=dble(int(x(j,idim)))+0.5D0
        cr(j,idim)=dble(int(x(j,idim)))+1.5D0
        wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
        wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
        wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
     end do
  end do

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igl(j,idim)=(int(cl(j,idim)))/2
        igc(j,idim)=(int(cc(j,idim)))/2
        igr(j,idim)=(int(cr(j,idim)))/2
     end do
  end do
  do j=1,np
     kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
     kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
     kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
     kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
     kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
     kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
     kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
     kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
     kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
     kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
     kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
     kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
     kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
     kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
     kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
     kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
     kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
     kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
     kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
     kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
     kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
     kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
     kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
     kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
     kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
     kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
     kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
  end do
  do ind=1,threetondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,threetondim
     do j=1,np
        ok(j)=ok(j).and.(igrid(j,ind)>0)
     end do
  end do

  ! If particle not entirely sits in the refinement,
  ! rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j)) then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! In this case, has to redo TSC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j)) then
           cl(j,idim)=dble(int(x(j,idim)))-0.5D0
           cc(j,idim)=dble(int(x(j,idim)))+0.5D0
           cr(j,idim)=dble(int(x(j,idim)))+1.5D0
           wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
           wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
           wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j)) then
           icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
           icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
           icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
        else
           icl(j,idim)=0
           icc(j,idim)=1
           icr(j,idim)=2
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
        icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
        icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
     else
        icell(j,1 )=1+icl(j,1)+3*icl(j,2)+9*icl(j,3)
        icell(j,2 )=1+icc(j,1)+3*icl(j,2)+9*icl(j,3)
        icell(j,3 )=1+icr(j,1)+3*icl(j,2)+9*icl(j,3)
        icell(j,4 )=1+icl(j,1)+3*icc(j,2)+9*icl(j,3)
        icell(j,5 )=1+icc(j,1)+3*icc(j,2)+9*icl(j,3)
        icell(j,6 )=1+icr(j,1)+3*icc(j,2)+9*icl(j,3)
        icell(j,7 )=1+icl(j,1)+3*icr(j,2)+9*icl(j,3)
        icell(j,8 )=1+icc(j,1)+3*icr(j,2)+9*icl(j,3)
        icell(j,9 )=1+icr(j,1)+3*icr(j,2)+9*icl(j,3)
        icell(j,10)=1+icl(j,1)+3*icl(j,2)+9*icc(j,3)
        icell(j,11)=1+icc(j,1)+3*icl(j,2)+9*icc(j,3)
        icell(j,12)=1+icr(j,1)+3*icl(j,2)+9*icc(j,3)
        icell(j,13)=1+icl(j,1)+3*icc(j,2)+9*icc(j,3)
        icell(j,14)=1+icc(j,1)+3*icc(j,2)+9*icc(j,3)
        icell(j,15)=1+icr(j,1)+3*icc(j,2)+9*icc(j,3)
        icell(j,16)=1+icl(j,1)+3*icr(j,2)+9*icc(j,3)
        icell(j,17)=1+icc(j,1)+3*icr(j,2)+9*icc(j,3)
        icell(j,18)=1+icr(j,1)+3*icr(j,2)+9*icc(j,3)
        icell(j,19)=1+icl(j,1)+3*icl(j,2)+9*icr(j,3)
        icell(j,20)=1+icc(j,1)+3*icl(j,2)+9*icr(j,3)
        icell(j,21)=1+icr(j,1)+3*icl(j,2)+9*icr(j,3)
        icell(j,22)=1+icl(j,1)+3*icc(j,2)+9*icr(j,3)
        icell(j,23)=1+icc(j,1)+3*icc(j,2)+9*icr(j,3)
        icell(j,24)=1+icr(j,1)+3*icc(j,2)+9*icr(j,3)
        icell(j,25)=1+icl(j,1)+3*icr(j,2)+9*icr(j,3)
        icell(j,26)=1+icc(j,1)+3*icr(j,2)+9*icr(j,3)
        icell(j,27)=1+icr(j,1)+3*icr(j,2)+9*icr(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do ind=1,threetondim
     do j=1,np
        if(ok(j)) then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
  do j=1,np
     vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
     vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
     vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
     vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
     vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
     vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
     vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
     vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
     vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
     vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
     vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
     vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
     vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
     vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
     vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
     vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
     vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
     vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
     vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
     vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
     vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
     vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
     vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
     vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
     vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
     vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
     vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
  end do

  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,threetondim
        do idim=1,ndim
           do j=1,np
              ff(j,idim)=ff(j,idim)+uold(indp(j,ind),idim+1)*vol(j,ind)
           end do
        end do
     end do
  endif
  if(poisson)then
     do ind=1,threetondim
        do idim=1,ndim
           do j=1,np
              ! cv-Galileon (add fifth force to the total force if appropriate)
              if(extradof .and. extradof2 .and. .not.extradof3) then
                 ! full
                 f1         = f(indp(j,ind),idim) + alpha_cvg*sf_grad(indp(j,ind),idim)
                 ff(j,idim) = ff(j,idim) + f1*vol(j,ind)
              else if(.not.extradof .and. extradof2 .and. extradof3) then
                 ! linearized
                 ! \nabla^2\Phi = \Omega_m*a*\rho*(1 + \alpha/\beta)
                 ! \alpha/\beta is the 5th-force to Newtonian force ratio
                 ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*(1.0d0 + alpha_cvg/beta_cvg)*vol(j,ind)
              else
                 ! LambdaCDM
                 ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
              end if
           end do
        end do
     end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static.or.tracer)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
        end do
     endif
  end do

  ! For sink cloud particle only
  if(sink)then
     ! Overwrite cloud particle velocity with sink velocity
     do idim=1,ndim
        do j=1,np
           isink=-idp(ind_part(j))
           if(isink>0)then
              new_vp(j,idim)=vsnew(isink,idim,ilevel)
           end if
        end do
     end do
  end if

  ! Store velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Update position
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+new_vp(j,idim)*dtnew(ilevel)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine move2
