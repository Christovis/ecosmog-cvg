! ------------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------------
! This file contains all MG-coarse-level related routines.
!
! Note: We use interpolation from fine_fine-level onto coarse-level.
! Therefore we don't have to differentiate between the vector-components
! of the Proca field. (Meaning we can neglect isf)
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     scalar field          sf            active_mg(myid,ilevel)%u(:,1)
!     truncation error      N/A           active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                      N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!     restricted sf         N/A           active_mg(myid,ilevel)%u(:,5)
!     restricted dens       N/A           active_mg(myid,ilevel)%u(:,6)
!     -----------------------------------------------------------------

! ------------------------------------------------------------------------------
! Compute residual, and stores it into active_mg(myid,ilevel)%u(:,3)
! ------------------------------------------------------------------------------
subroutine cmp_residual_mg_coarse_extradof(ilevel)
   ! On the coardse level we simplify by using interpolating amr-levels.
   ! Therefore this subroutine does not need to differentiate between the
   ! vector-components. (Meaning we can neglect isf)

   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,iskip_mg,iskip_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: eta,sfc,op
   real(dp), dimension(1:3) :: nb_sum_sfi,nb_sum_sfj
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr,icell_mg
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,j,ind2
   integer,  dimension(1:threetondim) :: jgrid_amr,jgrid_mg,jcell_mg,cpu_nbor_amr

   dx2 = (0.5d0**ilevel)**2.0

   ngrid=active_mg(myid,ilevel)%ngrid

   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = ncoarse+(ind-1)*ngridmax

      do igrid_mg=1,ngrid,nvector
         nbatch=MIN(nvector,ngrid-igrid_mg+1)
         do i=1,nbatch
            igrid_amr(i) = active_mg(myid,ilevel)%igrid(i+igrid_mg-1)
            icell_amr(i) = iskip_amr+igrid_amr(i)
            icell_mg(i)  = iskip_mg+igrid_mg+i-1
         end do
         call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
         do i=1,nbatch
            sfc = active_mg(myid,ilevel)%u(icell_mg(i),1)
            if(.not. btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0)) then
               do j=2,26,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               do j=11,17,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               do j=5,23,18
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do

               nb_sum_sfi(1) = active_mg(cpu_nbor_amr(13),ilevel)%u(jcell_mg(13),1) + &
                             & active_mg(cpu_nbor_amr(15),ilevel)%u(jcell_mg(15),1)
               nb_sum_sfi(2) = active_mg(cpu_nbor_amr(11),ilevel)%u(jcell_mg(11),1) + &
                             & active_mg(cpu_nbor_amr(17),ilevel)%u(jcell_mg(17),1)
               nb_sum_sfi(3) = active_mg(cpu_nbor_amr(5 ),ilevel)%u(jcell_mg(5 ),1) + &
                             & active_mg(cpu_nbor_amr(23),ilevel)%u(jcell_mg(23),1)

               ! op = lhs - rhs of scalar field Poission (Eq. 59)
               op  = nb_sum_sfi(1)+nb_sum_sfi(2)+nb_sum_sfi(3)-6.0d0*sfc &
                 - active_mg(myid,ilevel)%u(icell_mg(i),6)*dx2 &
                 - active_mg(myid,ilevel)%u(icell_mg(i),2)*dx2

            else
               op = 0.0d0
            end if
            active_mg(myid,ilevel)%u(icell_mg(i),3) = -op/dx2
         end do
      end do
   end do

end subroutine cmp_residual_mg_coarse_extradof

! ##################################################################
! ##################################################################

subroutine cmp_uvar_norm2_coarse_extradof(ivar,ilevel,norm2,n_cell_c)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer,  intent(in)  :: ilevel,ivar
   real(dp), intent(out) :: norm2,n_cell_c

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,icell_mg,iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2    = 0.0d0
   n_cell_c = 0.0d0
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg+igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0 .and. ivar/=4) cycle
         norm2    = norm2+(active_mg(myid,ilevel)%u(icell_mg,ivar))**2.0
         n_cell_c = n_cell_c+1.0d0
      end do
   end do
!  norm2 = dx2*norm2
end subroutine cmp_uvar_norm2_coarse_extradof

! ##################################################################
! ##################################################################

subroutine cmp_fvar_norm2_coarse_extradof(ivar,ilevel,norm2)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer,  intent(in)  :: ilevel,ivar
   real(dp), intent(out) :: norm2

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg,icell_mg,iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      ! Loop over active grids
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg + igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0) cycle
         norm2 = norm2+active_mg(myid,ilevel)%f(icell_mg,ivar)**2
      end do
   end do
   norm2 = dx2*norm2
end subroutine cmp_fvar_norm2_coarse_extradof

! ------------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------------
subroutine gauss_seidel_mg_coarse_extradof(ilevel,safe,redstep)
   ! On the coardse level we simplify by using interpolating amr-levels.
   ! Therefore this subroutine does not need to differentiate between the
   ! vector-components. (Meaning we can neglect isf)

   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel
   logical, intent(in) :: safe
   logical, intent(in) :: redstep
   integer, dimension(1:3,1:4)     :: ired,iblack

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind, ind0,igrid_mg
   integer  :: iskip_mg,iskip_amr
   real(dp) :: dtwondim = (twondim)
   integer  :: nac

   real(dp) :: eta,sfc,op,dop
   real(dp), dimension(1:3) :: nb_sum_sfi,nb_sum_sfj
   real(dp), dimension(1:3) :: nb_sum_sf_lpi,nb_sum_sf_lpj
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr,icell_mg
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,j,ind2
   integer,  dimension(1:threetondim) :: jgrid_amr,jgrid_mg,jcell_mg,cpu_nbor_amr

   ! Set constants
   dx2    = (0.5d0**ilevel)**2

   ired  (1,1:4)=(/1,0,0,0/)
   iblack(1,1:4)=(/2,0,0,0/)
   ired  (2,1:4)=(/1,4,0,0/)
   iblack(2,1:4)=(/2,3,0,0/)
   ired  (3,1:4)=(/1,4,6,7/)
   iblack(3,1:4)=(/2,3,5,8/)

   ngrid=active_mg(myid,ilevel)%ngrid
   ! Loop over cells, with red/black ordering
   do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_mg  = (ind-1)*ngrid
      iskip_amr = ncoarse+(ind-1)*ngridmax

      do igrid_mg=1,ngrid,nvector
         nbatch=MIN(nvector,ngrid-igrid_mg+1)
         do i=1,nbatch
            igrid_amr(i) = active_mg(myid,ilevel)%igrid(i+igrid_mg-1)
            icell_amr(i) = iskip_amr+igrid_amr(i)
            icell_mg(i)  = iskip_mg+igrid_mg+i-1
         end do
         call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
         do i=1,nbatch
            sfc = active_mg(myid,ilevel)%u(icell_mg(i),1)
            if(.not. btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0)) then
               do j=2,26,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
                  if(j.eq.14 .and. (icell_mg(i).ne.jcell_mg(j))) then
                     write(*,*) icell_mg(i),jcell_mg(j),jgrid_mg(j),igrid_mg
                     stop
                  end if
               end do
               do j=11,17,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               do j=5,23,18
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               nb_sum_sfi(1) = active_mg(cpu_nbor_amr(13),ilevel)%u(jcell_mg(13),1) + &
                             & active_mg(cpu_nbor_amr(15),ilevel)%u(jcell_mg(15),1)
               nb_sum_sfi(2) = active_mg(cpu_nbor_amr(11),ilevel)%u(jcell_mg(11),1) + &
                             & active_mg(cpu_nbor_amr(17),ilevel)%u(jcell_mg(17),1)
               nb_sum_sfi(3) = active_mg(cpu_nbor_amr(5 ),ilevel)%u(jcell_mg(5 ),1) + &
                             & active_mg(cpu_nbor_amr(23),ilevel)%u(jcell_mg(23),1)

               op  = nb_sum_sfi(1)+nb_sum_sfi(2)+nb_sum_sfi(3)-6.0d0*sfc &
                 - active_mg(myid,ilevel)%u(icell_mg(i),6)*dx2 &
                 - active_mg(myid,ilevel)%u(icell_mg(i),2)*dx2
               dop = -6.0d0

               active_mg(myid,ilevel)%u(icell_mg(i),1) = active_mg(myid,ilevel)%u(icell_mg(i),1)-op/dop
            end if
         end do
      end do
   end do

end subroutine gauss_seidel_mg_coarse_extradof

! ------------------------------------------------------------------------------
! Compute the physical rhs of the coarse level, including the corrections
! from the finer level, such as the residual etc.
! The physical rhs is stored in active_mg(myid,ilevel)%u(:,2)
! ------------------------------------------------------------------------------
subroutine make_physical_rhs_coarse_extradof(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ilevel

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind,igrid_mg
   integer  :: iskip_mg,iskip_amr
   real(dp) :: dtwondim = (twondim)

   real(dp) :: sfc,op
   real(dp), dimension(1:3) :: nb_sum_sfi,nb_sum_sfj
   integer,  dimension(1:nvector), save               :: igrid_amr,icell_amr,icell_mg
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_grids
   integer  :: nbatch,i,j,ind2
   integer,  dimension(1:threetondim) :: jgrid_amr,jgrid_mg,jcell_mg,cpu_nbor_amr

   ! Set constants
   dx2    = (0.5d0**ilevel)**2

   ngrid=active_mg(myid,ilevel)%ngrid

   ! Loop over cells
   do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = (ind-1)*ngridmax+ncoarse

      ! Loop over active grids
      do igrid_mg=1,ngrid,nvector
         nbatch=MIN(nvector,ngrid-igrid_mg+1)
         do i=1,nbatch
            igrid_amr(i) = active_mg(myid,ilevel)%igrid(i+igrid_mg-1)
            icell_amr(i) = iskip_amr+igrid_amr(i)
            icell_mg(i)  = iskip_mg+igrid_mg+i-1
         end do
         call get3cubefather(icell_amr,nbors_cells,nbors_grids,nbatch,ilevel)
         do i=1,nbatch
            sfc = active_mg(myid,ilevel)%u(icell_mg(i),5)
            if(.not. btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0)) then
               do j=2,26,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               do j=11,17,2
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do
               do j=5,23,18
                  ind2            = (nbors_cells(i,j)-ncoarse)/ngridmax+1
                  jgrid_amr(j)    = nbors_cells(i,j)-ncoarse-(ind2-1)*ngridmax
                  jgrid_mg(j)     = lookup_mg(jgrid_amr(j))
                  cpu_nbor_amr(j) = cpu_map(father(jgrid_amr(j)))
                  jcell_mg(j)     = jgrid_mg(j)+(ind2-1)*active_mg(cpu_nbor_amr(j),ilevel)%ngrid
               end do

               nb_sum_sfi(1) = active_mg(cpu_nbor_amr(13),ilevel)%u(jcell_mg(13),5) + &
                             & active_mg(cpu_nbor_amr(15),ilevel)%u(jcell_mg(15),5)
               nb_sum_sfi(2) = active_mg(cpu_nbor_amr(11),ilevel)%u(jcell_mg(11),5) + &
                             & active_mg(cpu_nbor_amr(17),ilevel)%u(jcell_mg(17),5)
               nb_sum_sfi(3) = active_mg(cpu_nbor_amr(5 ),ilevel)%u(jcell_mg(5 ),5) + &
                             & active_mg(cpu_nbor_amr(23),ilevel)%u(jcell_mg(23),5)


               op  = nb_sum_sfi(1)+nb_sum_sfi(2)+nb_sum_sfi(3)-6.0d0*sfc &
                 - active_mg(myid,ilevel)%u(icell_mg(i),6)*dx2

            else
               op = 0.0d0
            end if

            active_mg(myid,ilevel)%u(icell_mg(i),2) = active_mg(myid,ilevel)%u(icell_mg(i),2)+op/dx2
         end do
      end do
   end do
end subroutine make_physical_rhs_coarse_extradof


! ------------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------------
subroutine restrict_residual_coarse_reverse_extradof(ifinelevel)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none

   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell,ind_f_cell,cpu_amr
   integer :: iskip_c_mg
   integer :: igrid_c_amr,igrid_c_mg
   integer :: icell_c_amr,icell_c_mg
   integer :: iskip_f_mg
   integer :: igrid_f_amr,igrid_f_mg
   integer :: icell_f_mg
   integer :: icoarselevel

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Is fine cell masked?
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)
         icell_c_amr=father(igrid_f_amr)
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=active_mg(myid,ifinelevel)%u(icell_f_mg,3)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
   end do

end subroutine restrict_residual_coarse_reverse_extradof

! ------------------------------------------------------------------------------
! Restriction of the extradof to be used in the FAS (bottom-up)
! ------------------------------------------------------------------------------
subroutine restrict_extradof_coarse_reverse_extradof(ifinelevel)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: sf1
   real(dp) :: dtwotondim = (twotondim)
   integer,dimension(:),allocatable::n_masked

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allocate(n_masked(1:active_mg(myid,ifinelevel)%ngrid))
   n_masked(1:active_mg(myid,ifinelevel)%ngrid)=0

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Count the number of masked fine cells in each fine grid
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) then
            n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
         end if
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg                                  ! mg fine cell index

         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)          ! amr fine grid index
         icell_c_amr=father(igrid_f_amr)                                   ! amr coarse cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1                     ! amr coarse cell position
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1                       ! amr coarse cell position
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax           ! amr coarse grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                              ! coarse cell cpu index

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)                                 ! mg coarse grid index
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                                  ! mg coarse cell index

         ! If coarse cell is masked, it is boundary and R\tilde{u} is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0.or.n_masked(igrid_f_mg)==8) cycle

         ! Restriction to compute the sf value in the coarse cell
         sf1=active_mg(myid,ifinelevel)%u(icell_f_mg,1)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)+sf1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_extradof_coarse_reverse_extradof

subroutine restrict_density_coarse_reverse_extradof(ifinelevel)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: rho1
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   integer,dimension(:),allocatable::n_masked

   icoarselevel=ifinelevel-1

   allocate(n_masked(1:active_mg(myid,ifinelevel)%ngrid))
   n_masked(1:active_mg(myid,ifinelevel)%ngrid)=0

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Count the number of masked fine cells in each fine grid
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) then
            n_masked(igrid_f_mg)=n_masked(igrid_f_mg)+1
         end if
      end do
   end do

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg                                  ! mg fine cell index
         ! Is fine cell masked?
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)          ! amr fine grid index
         icell_c_amr=father(igrid_f_amr)                                   ! amr coarse cell index
!        ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1                     ! amr coarse cell position
         ind_c_cell=(icell_c_amr-ncoarse)/ngridmax+1                     ! amr coarse cell position
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax           ! amr coarse grid index
         cpu_amr=cpu_map(father(igrid_c_amr))                              ! coarse cell cpu index

         ! Convert to MG index, get MG coarse cell id                      ! mg coarse grid index
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg                                  ! mg coarse cell index

         ! If coarse cell is masked, it's outside boundary and R\rho is not needed
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0.or.n_masked(igrid_f_mg)==8) cycle

         ! Restriction to compute the sf value in the coarse cell
         rho1=active_mg(myid,ifinelevel)%u(icell_f_mg,6)/(dtwotondim-dble(n_masked(igrid_f_mg)))
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,6)+rho1
      end do
   end do

   deallocate(n_masked)

end subroutine restrict_density_coarse_reverse_extradof

! ------------------------------------------------------------------------------
! Interpolation (prolongation) and correction
! ------------------------------------------------------------------------------
subroutine interpolate_and_correct_coarse_extradof(ifinelevel)
   use amr_commons
   use poisson_commons
   use extradof_commons
   use extradof_parameters

   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i,ind_father,ind_average,ind_f,iskip_f_amr,ngrid_f,istart,nbatch
   integer  :: icell_c_amr,igrid_c_amr,igrid_c_mg,icell_c_mg,iskip_f_mg,icell_f_mg
   integer  :: icoarselevel,ind_c,cpu_c_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save                :: igrid_f_amr, icell_amr, cpu_amr
   integer,  dimension(1:nvector,1:threetondim), save  :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save    :: nbors_father_grids
   real(dp), dimension(1:nvector), save                :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   if(verbose) print '(A,I2)','correcting at level ',ifinelevel

   ! Loop over fine grids by vector sweeps
   ngrid_f=active_mg(myid,ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector                                        ! mg fine grid index

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active_mg(myid,ifinelevel)%igrid(istart+i-1)  ! amr fine grid index
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))                          ! amr coarse cell index
         cpu_amr(i)  =cpu_map(icell_amr(i))                           ! index of cpu containing the coarse cell
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids,nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
         iskip_f_mg  = (ind_f-1)*ngrid_f

         do i=1,nbatch
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)               ! amr fine cell index
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               icell_f_mg  = iskip_f_mg + istart+i-1
               if(active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,4)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               cpu_c_amr   = cpu_map(father(igrid_c_amr))
               if(igrid_c_mg<=0) cycle

               icell_c_mg  = (ind_c-1)*active_mg(cpu_c_amr,icoarselevel)%ngrid + igrid_c_mg
               ! only unmasked coarse cells contribute to the fine-cell correction
               if(active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,4)<=0.0) cycle
               corr(i)=corr(i)+coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,1) &
                      &       -coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,5)

            end do
         end do

         ! Correct the extradof on the fine cells
         do i=1,nbatch
            icell_f_mg  = iskip_f_mg+istart+i-1                    ! mg fine cell index
            active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) = active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1)+corr(i)
         end do

      end do
      ! End loop over cells

   end do

   if(verbose) print '(A,I2)','corrected at level ',ifinelevel

end subroutine interpolate_and_correct_coarse_extradof
