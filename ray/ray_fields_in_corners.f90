! ---------------------
! ray_fields_in_corners
! ---------------------
subroutine ray_fields_in_corners(ilevel,igrid,ind,ray_fields)
  !-----------------------------------------------------------------------------!
  ! This interpolates the values of fields from cell centers to cell vertices
  ! (this is done for ease of implementation of the ray integration in A)
  !
  ! cells) If a cell has all neighbours of the same ilevel (not at a coarse/fine
  ! interface), then the value of the vertice is given by the eight cells that
  ! share that vertice.
  !
  ! B) If a cell is at a coarse/fine interface, then the vertices on the
  ! interface side are obtained using a weighted average of the 27
  ! neighbouring father cells. This ensures the field varies continuously
  ! from the coarse to the fine side.
  ! The fine cell vertices that are further away from the coarse/fine
  ! interface are computed as in point A) above.
  !
  ! It takes as input {ilevel, igrid, ind}, and returns {int_rho} - the field at
  ! the cell vertices. The indices of the return variable int_rho(1:8), work in 
  ! the same way as the variable ind
  !-----------------------------------------------------------------------------!
  use amr_commons
  use amr_parameters
  use poisson_commons
  use ray_commons
  use ray_parameters ! for the test
  implicit none

  integer,                              intent(in   ) :: ilevel
  real(dp),dimension(1:RAY_NFIELDS,1:8),intent(inout) :: ray_fields
  integer,                              intent(in   ) :: igrid,ind

  integer  :: igrid_nbor,icell_nbor

  ! For finding 27 neighbours
  integer,dimension(1:nvector              ) :: father_cell
  integer,dimension(1:nvector,1:threetondim) :: nbors_father_cells
  integer,dimension(1:nvector,1:twotondim  ) :: nbors_father_grids

  logical,dimension(1:8) :: normal_vertice
  logical  :: neighbour_refined
  logical  :: increase_array

  integer  :: i
  integer  :: icell,old_size,new_size                               
  integer  :: itarget,ind_nbor
  integer  :: ix,iy,iz,jtarget,iind
  integer  :: xx,yy,zz,xy,yz,xz
  integer  :: igrid_nbor_fine,icell_nbor_fine,ind_nbor_fine
  integer  :: ind_coarse                                                                ! for the tidal and the case without neighbour

  real(dp) :: xtarget1,xtarget2,xtarget3
  real(dp) :: xc1,xc2,xc3
  real(dp) :: xcorner_father1,xcorner_father2,xcorner_father3
  real(dp) :: q1,q2,q3,q4,q5,q6,q7,q8
  real(dp) :: dx,dy,dz
  real(dp) :: weight,weight_tot
  real(dp) :: oneeighth

  real(dp),dimension(1:6) :: tidal
  real(dp),allocatable,dimension(:,:,:) :: temp_arr               
  real(dp),dimension(1:RAY_NFIELDS,1:8) :: ray_fields_tmp

  real(dp) :: phidot
  real(dp) :: aexp_mean,zexp_mean,E_friedmann_ray
  real(dp) :: analytic_isw  ! for testing isw
  logical  :: down_level

  !if(myid.eq.1) write(*,*) 'Entering ray_fields_in_corner'

  ! If corner fields have already calculated, use stored values  
  icell = igrid+ncoarse+(ind-1)*ngridmax                           ! ray cell index

  increase_array = .false.

  ! United ray_stored_fields array                                 ! Baojiu-Feb-2019
  ! Allocated in ray_step.f90
  if(ray_in_cell(icell).gt.0) then ! corner field values exist
     ! use stored corner field values
     ray_fields(1:RAY_NFIELDS,1:8) = ray_stored_fields(ray_in_cell(icell),1:RAY_NFIELDS,1:8) 
     return                                                        ! no need of further calculation
  end if

  father_cell(1) = father(igrid)                                   ! father cell of current grid
  oneeighth      = 1.0D0/8.0D0                                     ! numerical parameter

  ! Find 27 neighbours of the father of the current cell
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,1,ilevel)

  ! Loop over each vertice of the cell
  do itarget=1,8
     ! Initialise fields
     do i=1,RAY_NFIELDS
        ray_fields(i,itarget) = 0.0D0
     enddo
     weight_tot              = 0.0D0                               ! to be accumulated

     normal_vertice(itarget) = .true.                              ! all vertices are normal by default

     ! Loop over the 8 neighbours of the current cell
     do ind_nbor=1,8

        ! Get neighbouring grid index
        igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,ind_nbor,ind)))
        neighbour_refined = .false.

        ! Proceed separately for THREE cases 
        if(igrid_nbor.eq.0) then
           ! CASE II: neighbournig cells do not exist
           normal_vertice(itarget) = .false.                       ! mark vertice at boundary
           !write(*,*) itarget, "is not normal"
           exit 
        else
           icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,ind_nbor,ind)-1)*ngridmax !ID of neighbouring cell
           if(son(icell_nbor).eq.0) then
              ! CASE I: neighbouring cells of the same level exist
              ! add on to the sum of ray_fields for this itarget with weight unity.
              ray_fields(1,itarget) = ray_fields(1,itarget)+rho(icell_nbor)
              if(i_ray_tidal.gt.0) then                            ! Baojiu-Feb-2019
                 call tidal_tensor(ilevel,igrid_nbor,icell_nbor,ray_lll_mmm(itarget,ind_nbor,ind),tidal)
                 do i=1,6
                    ray_fields(i_ray_tidal+i-1,itarget) = ray_fields(i_ray_tidal+i-1,itarget)+tidal(i)
                 enddo
              end if                              
              if(i_ray_grads.gt.0) then                            ! Baojiu-Feb-2019 
                 do i=1,3
                    ray_fields(i_ray_grads+i-1,itarget) = ray_fields(i_ray_grads+i-1,itarget)+f(icell_nbor,i)
                 enddo
              end if 
              if(i_ray_phidt.gt.0) then                            ! Baojiu-Feb-2019
!                write(*,"(I5, E20.12, E20.12,E20.12,E20.12)") ilevel,aexp_new_ray(ilevel),aexp,aexp_old_ray(ilevel),aexp_mean
                 if(ilevel.gt.levelmin .and. dabs(phi_old(icell_nbor)/phi_old(father(igrid_nbor))-1.0D0).lt.1.0D-10) then
                    aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(father(igrid_nbor))/aexp_new_ray(ilevel-1)**2-phi_old(father(igrid_nbor))/aexp_old_ray(ilevel-1)**2)  &
                           / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1))                                                            &
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 else
                    aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2)                      &
                           / (aexp_new_ray(ilevel)-aexp_old_ray(ilevel))                                                                &
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 end if
                 ! For testing analytic solution
!                phidot = analytic_isw(ilevel,igrid,ind,itarget,aexp_mean,zexp_mean)
!                phidot = phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2
                 ray_fields(i_ray_phidt  ,itarget) = ray_fields(i_ray_phidt  ,itarget)+phi(icell_nbor)
                 ray_fields(i_ray_phidt+1,itarget) = ray_fields(i_ray_phidt+1,itarget)+phidot
              end if
              ! +++ add new stuff here to calculate other observables +++

              weight_tot = weight_tot+1.0D0
           else
              !write(*,*) "Setting neighbour_refined = true"
              neighbour_refined = .true.                           ! neighbour has been refined: dealt with separately below
              exit
           end if
        endif

     end do ! end of ind_nbor loop

     ! CASE III: neighbouring cells further refined
     if(neighbour_refined) then

        ! Initialise fields
        do i=1,RAY_NFIELDS
           ray_fields(i,itarget) = 0.0D0
        enddo
        weight_tot = 0.0D0                                         ! reset this to zero 

        ! Position of current cell in its parent grid
        iz = (ind-1)/4
        iy = (ind-1-iz*4)/2
        ix = (ind-1-iz*4-iy*2)
        !write(*,*) "ix = ", ix, iy, iz, igrid

        ! Coordinate of current cell centre wrt box
        !write(*,*) "xg = ", igrid, ix, ilevel  !, xg(igrid, 1)
        xc1 = xg(igrid,1)+(dble(ix)-0.5D0)*0.5D0**ilevel
        xc2 = xg(igrid,2)+(dble(iy)-0.5D0)*0.5D0**ilevel
        xc3 = xg(igrid,3)+(dble(iz)-0.5D0)*0.5D0**ilevel

        ! Position of target vertice wrt to cell centre
        iz = (itarget-1)/4
        iy = (itarget-1-iz*4)/2
        ix = (itarget-1-iz*4-iy*2)

        ! Coordinate of target vertice wrt to box
        xcorner_father1 = xc1+(dble(ix)-0.5D0)*0.5D0**ilevel
        xcorner_father2 = xc2+(dble(iy)-0.5D0)*0.5D0**ilevel
        xcorner_father3 = xc3+(dble(iz)-0.5D0)*0.5D0**ilevel

        ! Loop over the eight cells sharing the target vertice        
        do ind_nbor=1,8

           igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,ind_nbor,ind)))
           if(igrid_nbor.eq.0) then
              write(*,*) 'ray_fields_in_corners: grid structure problematic.'
              stop
           end if

           ! Index of neighbouring level-ilevel cell shareing the target vertice
           icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,ind_nbor,ind)-1)*ngridmax

           if(son(icell_nbor).eq.0) then                           ! neighbouring level-ilevel cell unrefined
              weight           = 1.0D0
              ! Add value to all the fields
              ray_fields(1,itarget) = ray_fields(1,itarget)+weight*rho(icell_nbor)
              if(i_ray_tidal.gt.0) then                            ! Baojiu-Feb-2019
                 call tidal_tensor(ilevel,igrid_nbor,icell_nbor,ray_lll_mmm(itarget,ind_nbor,ind),tidal)
                 do i=1,6
                    ray_fields(i_ray_tidal+i-1,itarget) = ray_fields(i_ray_tidal+i-1,itarget)+weight*tidal(i)
                 enddo
              end if                                              
              if(i_ray_grads.gt.0) then                            ! Baojiu-Feb-2019
                 do i=1,3
                    ray_fields(i_ray_grads+i-1,itarget) = ray_fields(i_ray_grads+i-1,itarget)+weight*f(icell_nbor,i)
                 enddo
              end if                               
              if(i_ray_phidt.gt.0) then
                 if(ilevel.gt.levelmin .and. dabs(phi_old(icell_nbor)/phi_old(father(igrid_nbor))-1.0D0).lt.1.0D-10) then
                    aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(father(igrid_nbor))/aexp_new_ray(ilevel-1)**2-phi_old(father(igrid_nbor))/aexp_old_ray(ilevel-1)**2) &
                           / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1)) &
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 else
                    aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2) & 
                           / (aexp_new_ray(ilevel)-aexp_old_ray(ilevel)) & 
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 end if
!                phidot = -2.0*phi(icell_nbor)/aexp_mean**2*E_friedmann_ray(zexp_mean,omega_m,omega_l)
!                phidot = analytic_isw(ilevel, igrid, ind, itarget, aexp_mean, zexp_mean)
!                phidot                 = (phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2)
                 ray_fields(i_ray_phidt  ,itarget) = ray_fields(i_ray_phidt  ,itarget)+weight*phi(icell_nbor)
                 ray_fields(i_ray_phidt+1,itarget) = ray_fields(i_ray_phidt+1,itarget)+weight*phidot
              end if
              ! +++ add new stuff here to calculate other observables +++

           else                                                    ! neighbouring level-ilevel cell refined
              weight           = 0.0D0                             ! reset to zero; to be accummulated

              igrid_nbor_fine  = son(icell_nbor)                   ! level-(ilevel+1) son grid index

              ! Find the position of the level-(ilevel+1) cell sharing the target vertice wrt its parent grid
              dx = xcorner_father1-(xg(igrid_nbor_fine,1)-0.5D0**(ilevel+1))
              dy = xcorner_father2-(xg(igrid_nbor_fine,2)-0.5D0**(ilevel+1))
              dz = xcorner_father3-(xg(igrid_nbor_fine,3)-0.5D0**(ilevel+1))
              if(dabs(dx).gt.0.5D0) dx = 1.0D0-dabs(dx)            ! periodic boundary condition
              if(dabs(dy).gt.0.5D0) dy = 1.0D0-dabs(dy)            ! periodic boundary condition
              if(dabs(dz).gt.0.5D0) dz = 1.0D0-dabs(dz)            ! periodic boundary condition
              ix = int(dx/0.5D0**(ilevel+1)*0.9999D0)
              iy = int(dy/0.5D0**(ilevel+1)*0.9999D0)
              iz = int(dz/0.5D0**(ilevel+1)*0.9999D0)
              ind_nbor_fine    = 1+ix+2*iy+4*iz

              if(ind_nbor_fine.lt.1 .or. ind_nbor_fine.gt.8) then
                 write(*,*) 'ray_fields_in_corners: ind_nbor_fine out of bound.'
                 write(*,*) ilevel, dx, dy, dz
                 call clean_stop()
              end if

              ! Find the index of the level-(ilevel+1) cell sharing the target vertice
              icell_nbor_fine  = igrid_nbor_fine+ncoarse+(ind_nbor_fine-1)*ngridmax

              weight = weight+oneeighth  ! accummulate weight 

              ! Find the position of level-ilevel cells wrt the target vertice
              iz = (ind_nbor-1)/4
              iy = (ind_nbor-1-iz*4)/2
              ix = (ind_nbor-1-iz*4-iy*2)

              ! Is the level-ilevel cell with the same y AND z coordinates refined?
              iind       = ind_nbor+(1-2*ix)
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              xx         = son(icell_nbor)

              ! Is the level-ilevel cell with the same x AND z coordinates refined?
              iind       = ind_nbor+(1-2*iy)*2
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              yy         = son(icell_nbor)

              ! Is the level-ilevel cell with the same x AND y coordinates refined?
              iind       = ind_nbor+(1-2*iz)*4
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              zz         = son(icell_nbor)

              ! Is the level-ilevel cell with the same z AND different x & y coordinates refined?
              iind       = ind_nbor+(1-2*ix)+(1-2*iy)*2
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              xy         = son(icell_nbor)

              ! Is the level-ilevel cell with the same x AND different y & z coordinates refined?
              iind       = ind_nbor+(1-2*iy)*2+(1-2*iz)*4
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              yz         = son(icell_nbor)

              ! Is the level-ilevel cell with the same y AND different x & z coordinates refined?
              iind       = ind_nbor+(1-2*ix)+(1-2*iz)*4
              igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind)))
              icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(itarget,iind,ind)-1)*ngridmax
              xz         = son(icell_nbor)

              ! Accummulate weight accordingly
              if(xx.eq.0)                       weight = weight+oneeighth  
              if(yy.eq.0)                       weight = weight+oneeighth  
              if(zz.eq.0)                       weight = weight+oneeighth  
              if(xx.eq.0.or.yy.eq.0.or.xy.eq.0) weight = weight+oneeighth   
              if(yy.eq.0.or.zz.eq.0.or.yz.eq.0) weight = weight+oneeighth   
              if(zz.eq.0.or.xx.eq.0.or.xz.eq.0) weight = weight+oneeighth   

              ! Accummulate density accordingly
              ray_fields(1,itarget) = ray_fields(1,itarget)+weight*rho(icell_nbor_fine)
              if(i_ray_tidal.gt.0) then                            ! Baojiu-Feb-2019           
                 call tidal_tensor(ilevel+1,igrid_nbor_fine,icell_nbor_fine,ind_nbor_fine,tidal)
                 do i=1,6
                    ray_fields(i_ray_tidal+i-1,itarget) = ray_fields(i_ray_tidal+i-1,itarget)+weight*tidal(i)
                 enddo
              end if                                       
              if(i_ray_grads.gt.0) then                            ! Baojiu-Feb-2019          
                 do i=1,3
                    ray_fields(i_ray_grads+i-1,itarget) = ray_fields(i_ray_grads+i-1,itarget)+weight*f(icell_nbor_fine,i)
                 enddo
              end if                  
              if(i_ray_phidt.gt.0) then                            ! Baojiu-Feb-2019 
                 if(ilevel+1.gt.levelmin .and. dabs(phi_old(icell_nbor_fine)/phi_old(father(igrid_nbor_fine))-1.0D0).lt.1.0D-10) then
                    aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(father(igrid_nbor_fine))/aexp_new_ray(ilevel)**2-phi_old(father(igrid_nbor_fine))/aexp_old_ray(ilevel)**2) &
                           / (aexp_new_ray(ilevel)-aexp_old_ray(ilevel)) &
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 else
                    aexp_mean = (aexp_new_ray(ilevel+1)+aexp_old_ray(ilevel+1))/2.0D0
                    zexp_mean = 1.0D0/aexp_mean-1.0D0
                    phidot = (phi(icell_nbor_fine)/aexp_new_ray(ilevel+1)**2-phi_old(icell_nbor_fine)/aexp_old_ray(ilevel+1)**2) & 
                           / (aexp_new_ray(ilevel+1)-aexp_old_ray(ilevel+1)) &
                           * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
                 end if
!                phidot    = (phi(icell_nbor_fine)/aexp_new_ray(ilevel+0)**2-phi_old(icell_nbor_fine)/aexp_old_ray(ilevel+0)**2)
                 ray_fields(i_ray_phidt  ,itarget) = ray_fields(i_ray_phidt  ,itarget)+weight*phi(icell_nbor_fine)
                 ray_fields(i_ray_phidt+1,itarget) = ray_fields(i_ray_phidt+1,itarget)+weight*phidot
              end if            
              ! +++ add new stuff here to calculate other observables +++
           end if

           ! Find the total weight from all cells sharing the target vertice
           weight_tot = weight_tot+weight
        end do

     end if

     ! normalise the mean value (subtract 1 to get density contrast)
     if(normal_vertice(itarget)) then
        do i=1,RAY_NFIELDS
           ray_fields(i,itarget) = ray_fields(i,itarget)/weight_tot !-1.0D0
        enddo
!       ray_fields(1,itarget) = ray_fields(1,itarget)/weight_tot   !-1.0D0
!       ray_fields(1,itarget) = ray_fields(1,itarget)/aexp_old_ray ! this is delta/a
     end if

  end do !End of itarget loop

  ! Baojiu-Feb-2019: united all arrays to a single array ray_stored_fields
  if(all(normal_vertice)) then                                     ! if not at refinement boundary, calculation finishes here

     ray_ncells = ray_ncells+1
     ray_in_cell(icell) = ray_ncells

     if(ray_ncells.gt.size(ray_stored_fields,1)) then

        increase_array = .true.

        old_size = size(ray_stored_fields,1)
!       new_size = old_size+old_size/2
        new_size = 5*old_size

        write(*,*) 'ray_fields_in_corners: enlarging corner arrays',old_size,new_size 

        allocate(temp_arr(1:old_size,1:RAY_NFIELDS,1:8))

        temp_arr = ray_stored_fields
        deallocate(ray_stored_fields                              )
        allocate  (ray_stored_fields(1:new_size,1:RAY_NFIELDS,1:8))
        ray_stored_fields                             = 0.0D0
        ray_stored_fields(1:old_size,1:RAY_NFIELDS,:) = temp_arr(1:old_size,1:RAY_NFIELDS,:)

        deallocate(temp_arr)

     end if

     ray_stored_fields(ray_ncells,1:RAY_NFIELDS,1:8) = ray_fields(1:RAY_NFIELDS,1:8)

     return

  end if

  ! CASE II: neighbournig level-ilevel cells do not exist 
  ! Loop over each vertice of the father, i.e., level-(ilevel-1), cell
  do jtarget=1,8

     weight_tot = 0.0D0                                            ! to be accumulated 
     ! Initialise fields
     do i=1,RAY_NFIELDS
        ray_fields_tmp(i,jtarget) = 0.0D0
     enddo


     ! Loop over the 8 neighbours of the father level-(ilevel-1) cell
     do ind_nbor=1,8

        ! Get the index of neighbouring level-ilevel grid
        igrid_nbor = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor,jtarget)))

        if(igrid_nbor.eq.0) then                                   ! neighbouring level-ilevel grid not exist
           icell_nbor = nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor,jtarget))
           weight     = 1.0D0
           ! grid that corresponds to this cell_nbor
           ind_coarse = (icell_nbor-ncoarse-1)/ngridmax+1
           igrid_nbor = icell_nbor-ncoarse-(ind_coarse-1)*ngridmax
           if(i_ray_tidal.gt.0) then                               ! Baojiu-Feb-2019 
              call tidal_tensor(ilevel-1,igrid_nbor,icell_nbor,ind_coarse,tidal)
           end if                 
           down_level = .true. 
        else                                                       ! neighbouring level-ilevel grid exists
           weight     = 0.0D0                                      ! to be accummulated

           ! Get the index of the neighbouring level-ilevel cell   
           icell_nbor = igrid_nbor+ncoarse+(ray_lll_mmm(jtarget,ind_nbor,jtarget)-1)*ngridmax

           weight     = weight+oneeighth                           ! accummulate weight

           ! Get the position of father cell #ind_nbor wrt the target vertice
           iz = (ind_nbor-1)/4
           iy = (ind_nbor-1-iz*4)/2
           ix = (ind_nbor-1-iz*4-iy*2)

           ! Are neighbouring cells of father cell #ind_nbor refined?
           xx = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*ix)             ,jtarget)))
           yy = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*iy)*2           ,jtarget)))
           zz = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*iz)*4           ,jtarget)))
           xy = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*ix)  +(1-2*iy)*2,jtarget)))   
           yz = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*iy)*2+(1-2*iz)*4,jtarget)))   
           xz = son(nbors_father_cells(1,ray_kkk_mmm(jtarget,ind_nbor+(1-2*ix)  +(1-2*iz)*4,jtarget)))   

           ! Accummulate weight accordingly
           if(xx.eq.0)                       weight = weight+oneeighth   
           if(yy.eq.0)                       weight = weight+oneeighth  
           if(zz.eq.0)                       weight = weight+oneeighth  
           if(xx.eq.0.or.yy.eq.0.or.xy.eq.0) weight = weight+oneeighth  
           if(yy.eq.0.or.zz.eq.0.or.yz.eq.0) weight = weight+oneeighth  
           if(zz.eq.0.or.xx.eq.0.or.xz.eq.0) weight = weight+oneeighth  

           if(i_ray_tidal.gt.0) then                               ! Baojiu-Feb-2019                   
              call tidal_tensor(ilevel,igrid_nbor,icell_nbor,ray_lll_mmm(jtarget,ind_nbor,jtarget),tidal)
           end if 
           down_level = .false.   

        end if

        weight_tot       = weight_tot+weight
        ! Add value to all the fields
        ray_fields_tmp(1,jtarget) = ray_fields_tmp(1,jtarget)+weight*rho(icell_nbor)
        if(i_ray_tidal.gt.0) then                                  ! Baojiu-Feb-2019 
           do i=1,6
              ray_fields_tmp(i_ray_tidal+i-1,jtarget) = ray_fields_tmp(i_ray_tidal+i-1,jtarget)+weight*tidal(i)
           enddo
        end if
        if(i_ray_grads.gt.0) then                                  ! Baojiu-Feb-2019
           do i=1,3
              ray_fields_tmp(i_ray_grads+i-1,jtarget) = ray_fields_tmp(i_ray_grads+i-1,jtarget)+weight*f(icell_nbor,i)
           enddo
        end if 
        if(i_ray_phidt.gt.0) then                                  ! Baojiu-Feb-2019
           if(down_level) then
              if(ilevel-1.gt.levelmin .and. dabs(phi_old(icell_nbor)/phi_old(father(igrid_nbor))-1.0D0).lt.1.0D-10) then
                 aexp_mean = (aexp_new_ray(ilevel-2)+aexp_old_ray(ilevel-2))/2.0D0
                 zexp_mean = 1.0D0/aexp_mean-1.0D0
                 phidot = (phi(father(igrid_nbor))/aexp_new_ray(ilevel-2)**2-phi_old(father(igrid_nbor))/aexp_old_ray(ilevel-2)**2) &
                        / (aexp_new_ray(ilevel-2)-aexp_old_ray(ilevel-2)) &
                        * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
              else
                 aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
                 zexp_mean = 1.0D0/aexp_mean-1.0D0
                 phidot = (phi(icell_nbor)/aexp_new_ray(ilevel-1)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel-1)**2) &
                        / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1)) &
                        * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
              end if
!             phidot = (phi(icell_nbor)/aexp_new_ray(ilevel-0)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel-0)**2)
           else
              if(ilevel.gt.levelmin .and. dabs(phi_old(icell_nbor)/phi_old(father(igrid_nbor))-1.0D0).lt.1.0D-10) then
                 aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
                 zexp_mean = 1.0D0/aexp_mean-1.0D0
                 phidot = (phi(father(igrid_nbor))/aexp_new_ray(ilevel-1)**2-phi_old(father(igrid_nbor))/aexp_old_ray(ilevel-1)**2) &
                        / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1)) &
                        * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean 
              else
                 aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
                 zexp_mean = 1.0D0/aexp_mean-1.0D0
                 phidot = (phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2) &
                        / (aexp_new_ray(ilevel)-aexp_old_ray(ilevel)) &
                        * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
              end if
!             phidot = (phi(icell_nbor)/aexp_new_ray(ilevel)**2-phi_old(icell_nbor)/aexp_old_ray(ilevel)**2)
           end if
           ray_fields_tmp(i_ray_phidt  ,jtarget) = ray_fields_tmp(i_ray_phidt  ,jtarget)+weight*phi(icell_nbor)
           ray_fields_tmp(i_ray_phidt+1,jtarget) = ray_fields_tmp(i_ray_phidt+1,jtarget)+weight*phidot
        end if   
        ! +++ add new stuff here to calculate other observables +++
     end do

     do i=1,RAY_NFIELDS
        ray_fields_tmp(i,jtarget) = ray_fields_tmp(i,jtarget)/weight_tot
     enddo

  end do

  ! Get position of current cell wrt its parent grid
  iz  = (ind-1)/4
  iy  = (ind-1-iz*4)/2
  ix  = (ind-1-iz*4-iy*2)

  ! Get coordinates of current cell centre wrt box
  xc1 = xg(igrid,1)+(dble(ix)-0.5D0)*0.5D0**ilevel
  xc2 = xg(igrid,2)+(dble(iy)-0.5D0)*0.5D0**ilevel
  xc3 = xg(igrid,3)+(dble(iz)-0.5D0)*0.5D0**ilevel

  ! Get coordinates of the vertice wrt box
  xcorner_father1 = xg(igrid,1)-0.5D0**ilevel
  xcorner_father2 = xg(igrid,2)-0.5D0**ilevel
  xcorner_father3 = xg(igrid,3)-0.5D0**ilevel

  ! For all the fields
  do i=1,RAY_NFIELDS
     ! Quantities for trilinear interpolation
     q1 = ray_fields_tmp(i,1)
     q2 = ray_fields_tmp(i,2)-ray_fields_tmp(i,1)
     q3 = ray_fields_tmp(i,3)-ray_fields_tmp(i,1)
     q4 = ray_fields_tmp(i,5)-ray_fields_tmp(i,1)
     q5 = ray_fields_tmp(i,4)-ray_fields_tmp(i,3)-ray_fields_tmp(i,2)+ray_fields_tmp(i,1)
     q6 = ray_fields_tmp(i,7)-ray_fields_tmp(i,5)-ray_fields_tmp(i,3)+ray_fields_tmp(i,1)
     q7 = ray_fields_tmp(i,6)-ray_fields_tmp(i,5)-ray_fields_tmp(i,2)+ray_fields_tmp(i,1)
     q8 = ray_fields_tmp(i,8)-ray_fields_tmp(i,7)-ray_fields_tmp(i,6)-ray_fields_tmp(i,4)+ &
          ray_fields_tmp(i,2)+ray_fields_tmp(i,5)+ray_fields_tmp(i,3)-ray_fields_tmp(i,1)

     ! Do trilinear interpolation to get interpolated value at 8 vertices of current cell
     do itarget=1,8 

        if(normal_vertice(itarget)) cycle                             ! Skip if the vertice is normal

        ! Get position of vertice wrt current cell centre
        iz = (itarget-1)/4
        iy = (itarget-1-iz*4)/2
        ix = (itarget-1-iz*4-iy*2)

        ! Get coordinates of current-cell vertice wrt box 
        xtarget1 = xc1+(dble(ix)-0.5D0)*0.5D0**ilevel
        xtarget2 = xc2+(dble(iy)-0.5D0)*0.5D0**ilevel
        xtarget3 = xc3+(dble(iz)-0.5D0)*0.5D0**ilevel

        ! Transform to unit of cell size
        dx = dabs(xtarget1-xcorner_father1)/0.5D0**(ilevel-1)
        dy = dabs(xtarget2-xcorner_father2)/0.5D0**(ilevel-1)
        dz = dabs(xtarget3-xcorner_father3)/0.5D0**(ilevel-1)

        if(dx<0.0D0 .or. dx>1.0D0) write(*,*) 'dx wrong:',dx 
        if(dy<0.0D0 .or. dy>1.0D0) write(*,*) 'dy wrong:',dy 
        if(dz<0.0D0 .or. dz>1.0D0) write(*,*) 'dz wrong:',dz 

        ! Do trilinear interpolation
        ray_fields(i,itarget) = q1+q2*dx+q3*dy+q4*dz+q5*dx*dy+q6*dy*dz+q7*dx*dz+q8*dx*dy*dz

     end do

  enddo  ! over i = fields

  ray_ncells = ray_ncells+1
  ray_in_cell(icell) = ray_ncells

  ! Baojiu-Feb-2019: united all arrays to a single array ray_stored_fields
  if(ray_ncells.gt.size(ray_stored_fields,1)) then

     increase_array = .true.

     old_size = size(ray_stored_fields,1)
     new_size = old_size+old_size/2

     write(*,*) 'ray_fields_in_corners: enlarging corner arrays',old_size,new_size 

     allocate(temp_arr(1:old_size,1:RAY_NFIELDS,1:8))

     temp_arr = ray_stored_fields
     deallocate(ray_stored_fields                              )
     allocate  (ray_stored_fields(1:new_size,1:RAY_NFIELDS,1:8))
     ray_stored_fields                             = 0.0D0
     ray_stored_fields(1:old_size,1:RAY_NFIELDS,:) = temp_arr(1:old_size,1:RAY_NFIELDS,:)

     deallocate(temp_arr)

  end if

  ray_stored_fields(ray_ncells,1:RAY_NFIELDS,1:8) = ray_fields(1:RAY_NFIELDS,1:8)

end subroutine ray_fields_in_corners

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
! Baojiu-Feb-2019: this new subroutine was added to closely mimic
!                  ray_fields_in_corners, but instead compute ray
!                  fields at the cell centres rather than corners
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
subroutine ray_fields_at_centers(ilevel,igrid,ind,ray_fields_c)
  use amr_commons
  use amr_parameters
  use poisson_commons
  use ray_commons
  use ray_parameters

  implicit none

  integer,                          intent(in   ) :: ilevel
  real(dp),dimension(1:RAY_NFIELDS),intent(inout) :: ray_fields_c
  integer,                          intent(in   ) :: igrid,ind

  integer                 :: i,icell
  real(dp)                :: phidot,aexp_mean,zexp_mean
  real(dp),dimension(1:6) :: tidal 
  real(dp) :: E_friedmann_ray

  ray_fields_c = 0.0D0
  icell = igrid+ncoarse+(ind-1)*ngridmax 

  if(i_ray_tidal.gt.0) then
     call tidal_tensor(ilevel,igrid,icell,ind,tidal)
     do i=1,6
        ray_fields_c(i_ray_tidal+i-1) = tidal(i)
     end do
  end if
  if(i_ray_grads.gt.0) then    
     do i=1,3
        ray_fields_c(i_ray_grads+i-1) = f(icell,i)
     enddo
  end if 
  if(i_ray_phidt.gt.0) then
     phidot = 0.0d0
     if(ilevel.gt.levelmin.and.dabs(phi_old(icell)/phi_old(father(igrid))-1.0D0).lt.1.0D-10) then
        ! If the cell had just been refined, then compute time derivative using coarse level for stability
        aexp_mean = (aexp_new_ray(ilevel-1)+aexp_old_ray(ilevel-1))/2.0D0
        zexp_mean = 1.0D0/aexp_mean-1.0D0
        phidot    = (phi(father(igrid))/aexp_new_ray(ilevel-1)**2-phi_old(father(igrid))/aexp_old_ray(ilevel-1)**2)  &
                  / (aexp_new_ray(ilevel-1)-aexp_old_ray(ilevel-1)) * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
     else
        ! If the cell already existed in previous steps, use current cells data
        aexp_mean = (aexp_new_ray(ilevel)+aexp_old_ray(ilevel))/2.0D0
        zexp_mean = 1.0D0/aexp_mean-1.0D0
        phidot    = (phi(icell)/aexp_new_ray(ilevel)**2-phi_old(icell)/aexp_old_ray(ilevel)**2)                      &
                  / (aexp_new_ray(ilevel)-aexp_old_ray(ilevel))                                                      &
                  * E_friedmann_ray(zexp_mean,omega_m,omega_l)*aexp_mean
     end if
     ray_fields_c(i_ray_phidt  ) = phi(icell)
     ray_fields_c(i_ray_phidt+1) = phidot
  end if

end subroutine ray_fields_at_centers


!=====================================
! analytic_isw
!=====================================
function analytic_isw(ilevel,igrid,ind,itarget,aexp_mean,zexp_mean)
  use amr_commons
  use amr_parameters
  use poisson_commons
  use ray_commons
  use ray_parameters ! for the test
  implicit none

  integer :: ilevel,igrid,ind,itarget
  real(dp) :: analytic_isw
  real(dp) :: aexp_mean, zexp_mean

  integer :: ix,iy,iz
  real(dp) :: xc1,xc2,xc3
  real(dp) :: xcorner_father1,xcorner_father2,xcorner_father3
  real(dp) :: r2
  real(dp) :: E_friedmann_ray

  ! Position of current cell in its parent grid
  iz = (ind-1)/4
  iy = (ind-1-iz*4)/2
  ix = (ind-1-iz*4-iy*2)
  !write(*,*) "ix = ", ix, iy, iz, igrid

  ! Coordinate of current cell centre wrt box
  !write(*,*) "xg = ", igrid, ix, ilevel  !, xg(igrid, 1)
  xc1 = xg(igrid,1)+(dble(ix)-0.5D0)*0.5D0**ilevel
  xc2 = xg(igrid,2)+(dble(iy)-0.5D0)*0.5D0**ilevel
  xc3 = xg(igrid,3)+(dble(iz)-0.5D0)*0.5D0**ilevel

  ! Position of target vertice wrt to cell centre
  iz = (itarget-1)/4
  iy = (itarget-1-iz*4)/2
  ix = (itarget-1-iz*4-iy*2)

  ! Coordinate of target vertice wrt to box
  xcorner_father1 = xc1+(dble(ix)-0.5D0)*0.5D0**ilevel
  xcorner_father2 = xc2+(dble(iy)-0.5D0)*0.5D0**ilevel
  xcorner_father3 = xc3+(dble(iz)-0.5D0)*0.5D0**ilevel

  r2 = (xcorner_father1-0.5D0)**2 + (xcorner_father2-0.5D0)**2 + (xcorner_father3-0.5D0)**2

  analytic_isw = -2.0*(-10000.0D0) * exp(-r2/((0.7D0/ray_Lbox)**2))/aexp_mean**2*E_friedmann_ray(zexp_mean,omega_m,omega_l);

  return 

end function analytic_isw


! -----------------------------------------------------------
! ray_interpolate
! Baojiu-Feb-2019: Note that this subroutine is not called in
!                  ray calculations but only in outputs 
! -----------------------------------------------------------
subroutine ray_interpolate(ilevel,igrid,ind,ray_fields)

  use amr_commons
  use poisson_commons
  use ray_commons
  use ray_parameters ! for the test
  implicit none

  integer :: ilevel
  real(dp),dimension(1:RAY_NFIELDS,1:8),intent(inout) :: ray_fields
  integer :: igrid, ind

  integer :: igrid_nbor,icell_nbor,ind_nbor
  real(dp) :: rho_nbor

  ! For finding 27 neighbours
  integer,dimension(1:nvector) :: father_cell
  integer,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer,dimension(1:nvector,1:twotondim)::nbors_father_grids

  integer :: itarget,iind

  integer :: ix,iy,iz,ind_coarse,igrid_coarse
  real(dp) :: xtarget,ytarget,ztarget,dx,dy,dz,xc1,xc2,xc3

  integer :: weight,weight_tot

  real(dp),dimension(1:6) :: tidal

  real (dp) :: coeff

  ! For lower level interpolation
  integer :: ind_father,indice,ind_average

  ! For loops over fields
  integer :: i

  !---------------------------------------------------------------------------
  ! This interpolates the values of fields from cell centers to cell vertices
  ! (this is done for ease of implementation of the ray integration in cells)
  !
  ! A) If a cell has all neighbours of the same ilevel (not at a coarse/fine
  ! interface), then the value of the vertice is given by the eight cells that
  ! share that vertice.
  !
  ! B) If a cell is at a coarse/fine interface, then the vertices on the interface
  ! side are obtained using a weighted average of the 27 neighbouring father
  ! cells. This ensures the field varies continuously from the coarse to the
  ! fine side. The fine cell vertices that are further away from the coarse/fine
  ! interface are computed as in point A) above.
  !
  ! It takes as input {ilevel, igrid, ind}, and returns {int_rho} - the field at
  ! the cell vertices. The indices of the return variable int_rho1:8), work in 
  ! the same way as the variable ind
  !--------------------------------------------------------------------------

  ! Find 27 father neighbours of the cell
  father_cell(1) = father(igrid)
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,1,ilevel)

  ! Loop over each vertice of the cell
  do itarget=1,8

     ! Initialise fields
     do i=1,RAY_NFIELDS
        ray_fields(i,itarget) = 0.0D0
     enddo

     ! Loop over the 8 neighbouring cells
     do iind = 1,8
        igrid_nbor = son(  nbors_father_cells(1,ray_kkk_mmm(itarget,iind,ind) )  )

        ! Proceed separately depending on whether cell has all neighbours or not
        if(igrid_nbor .eq. 0) then
           ! Cell does not have all neighbours, i.e., cell is at a coarse/fine interface

           ! Re-initialize int_rho(itarget) because for this itarget, some iind
           ! iterations on the "else" below may have been performed and changed
           ! int_rho(itarget). If we entered this if, then we need to compute
           ! int_rho in a different way and so we re-initialize it.
           do i=1,RAY_NFIELDS
              ray_fields(i,itarget) = 0.0D0
           enddo

           ! ================================================================ !
           ! ---------------------------------------------------------------- !
           ! This block determines the coordinate of the cell wrt to box

           ! Position of cell wrt grid (0, 1 values)
           iz = (ind-1          )/4
           iy = (ind-1-4*iz     )/2
           ix = (ind-1-4*iz-2*iy)
           ! Position of cell center wrt box
           xc1 = xg(igrid, 1)+(dble(ix)-0.5D0)/2.0D0**ilevel
           xc2 = xg(igrid, 2)+(dble(iy)-0.5D0)/2.0D0**ilevel
           xc3 = xg(igrid, 3)+(dble(iz)-0.5D0)/2.0D0**ilevel
           ! ---------------------------------------------------------------- !
           ! ================================================================ !
           ! ================================================================ !
           ! ---------------------------------------------------------------- !
           ! This block determines the coordinate of the vertice itarget wrt box

           ! Position of target vertex wrt cell
           ! (0,1 values -- overwrite values above)
           iz = (itarget-1          )/4
           iy = (itarget-1-4*iz     )/2
           ix = (itarget-1-4*iz-2*iy)
           ! Position of vertex itarget wrt box
           xtarget = xc1+(dble(ix)-0.5D0)/2.0D0**ilevel 
           ytarget = xc2+(dble(iy)-0.5D0)/2.0D0**ilevel 
           ztarget = xc3+(dble(iz)-0.5D0)/2.0D0**ilevel 
           ! ---------------------------------------------------------------- !
           ! ================================================================ !

           weight_tot = 0 ! This is going to be the normalization factor of the weighted average

           ! Loop over the 27 neighbouring father cells
           do ind_father=1,27
              indice = nbors_father_cells(1,ind_father) ! ID of the neighbouring father cell
              ! ================================================================ !
              ! ---------------------------------------------------------------- !
              ! This block determines the coordinate of the father neighbour
              ! ind_father wrt box

              ind_coarse   = (indice-ncoarse-1)/ngridmax+1
              igrid_coarse = indice-ncoarse-(ind_coarse-1)*ngridmax
              ! Position of the target father cell wrt to its grid
              ! (this grid is cell at a lower ilevel)
              iz = (ind_coarse-1          )/4
              iy = (ind_coarse-1-4*iz     )/2
              ix = (ind_coarse-1-4*iz-2*iy)
              ! Position of neighbouring father cells wrt box
              ! (overwrite values above)
              xc1 = xg(igrid_coarse,1)+(dble(ix)-0.5D0)/2.0D0**(ilevel-1)
              xc2 = xg(igrid_coarse,2)+(dble(iy)-0.5D0)/2.0D0**(ilevel-1)
              xc3 = xg(igrid_coarse,3)+(dble(iz)-0.5D0)/2.0D0**(ilevel-1)
              ! ---------------------------------------------------------------- !
              ! ================================================================ !
              ! ================================================================ !
              ! ---------------------------------------------------------------- !
              ! This block determines the weight of the father cell indice in
              ! the weighted average. The weighted density of this cell is added
              ! to int_rho

              ! Difference in x,y,z-dir between vertex itarget and center of neighbouring father cell indice
              dx = dabs(xc1-xtarget)*2.0D0**ilevel
              dy = dabs(xc2-ytarget)*2.0D0**ilevel
              dz = dabs(xc3-ztarget)*2.0D0**ilevel
              ! Weight of the father cell
              weight     = (2-(int(dx+1.0D-6)+1)/2)*(2-(int(dy+1.0D-6)+1)/2)*(2-(int(dz+1.0D-6)+1)/2)
              ! Add on to the total weight, which is going to be the normalization factor below 
              weight_tot = weight_tot+weight 
              ! Add on to all the fields sum for this itarget
              ray_fields(1,itarget) = ray_fields(1,itarget)+dble(weight)*rho(indice)
              if(i_ray_tidal.gt.0) then                            ! Baojiu-Feb-2019 
                 call tidal_tensor(ilevel-1,igrid_coarse,indice,ind_coarse,tidal)
                 do i=1,6
                    ray_fields(i_ray_tidal+i-1,itarget) = ray_fields(i_ray_tidal+i-1,itarget)+dble(weight)*tidal(i)
                 enddo
              end if   
              if(i_ray_grads.gt.0) then                            ! Baojiu-Feb-2019  
                 do i=1,3
                    ray_fields(i_ray_grads+i-1,itarget) = ray_fields(i_ray_grads+i-1,itarget)+dble(weight)*f(indice,i)
                 enddo
              end if   
              if(i_ray_phidt.gt.0) then                            ! Baojiu-Feb-2019
                 ray_fields(i_ray_phidt,itarget) = ray_fields(i_ray_phidt,itarget)+dble(weight)*phi(indice)
              end if
              ! +++ add new stuff here to calculate other observables +++

              ! ---------------------------------------------------------------- !
              ! ================================================================ !

              ! ================================================================ !
              ! ---------------------------------------------------------------- !
              ! This block is here for tests. It prints weights of itarget of
              ! ind of a grid. Want igrid to be one at a coarse/fine interface

              !              if (igrid == 5387) then
              !                 if (ind==1 .and. itarget == 5) then
              !                    write(*,*) 'aaasssddd', ind_father, dx+dy+dz
              !                    write(*,*)  dx, dy, dz
              !                    write(*,*) weight, weight_tot
              !                 end if
              !              end if
              ! ================================================================ !
              ! ---------------------------------------------------------------- !
           end do !End of ind_father loop

           ! Finally, normalize by weight_tot (multiplication by 8 is because
           ! this is going to be divided by 8 below -- easier for cases where neighbours exist)
           do i=1,RAY_NFIELDS
              ray_fields(i,itarget) = ray_fields(i,itarget)/dble(weight_tot)*8.0D0
           enddo

           exit !Exit because after the ind_father loop, this coarse/fine vertice is done and don't want further iind iterations to affect its value

        else
           ! Cell has all neighbours, i.e., cell is NOT at a coarse/fine interface
           ind_nbor = ray_lll_mmm(itarget,iind,ind)
           icell_nbor = igrid_nbor+(ncoarse+(ind_nbor-1)*ngridmax) !ID of neighbouring cell
           ! Add on to the sum of int_rho for this itarget with weight unity.
           ray_fields(1,itarget) = ray_fields(1,itarget)+rho(icell_nbor)
           if(i_ray_tidal.gt.0) then                               ! Baojiu-Feb-2019 
              call tidal_tensor(ilevel,igrid_nbor,icell_nbor,ind_nbor,tidal)
              do i=1,6
                 ray_fields(i_ray_tidal+i-1,itarget) = ray_fields(i_ray_tidal+i-1,itarget)+tidal(i)
              enddo
           end if
           if(i_ray_grads.gt.0) then                               ! Baojiu-Feb-2019
              do i=1,3
                 ray_fields(i_ray_grads+i-1,itarget) = ray_fields(i_ray_grads+i-1,itarget)+f(icell_nbor,i)
              enddo
           end if
           if(i_ray_phidt.gt.0) then                               ! Baojiu-Feb-2019
              ray_fields(i_ray_phidt,itarget) = ray_fields(i_ray_phidt,itarget)+phi(icell_nbor)
           end if
           ! +++ add new stuff here to calculate other observables +++
        endif

     end do !End of iind loop

     ! normalise the mean value
     do i=1,RAY_NFIELDS
        ray_fields(i,itarget) = ray_fields(i,itarget)/8.0D0
     enddo

  end do !End of itarget loop

end subroutine ray_interpolate

