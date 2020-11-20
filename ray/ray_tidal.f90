!===============================
! tidal_tensor
!===============================
subroutine tidal_tensor(ilevel,igrid,icell,ind,tidal)
  use amr_commons
  use poisson_commons
  use ray_commons
  use ray_parameters ! for the test
  use extradof_commons
  implicit none

  real(dp), dimension(1:6) :: tidal
  integer :: ilevel
  integer :: igrid, icell
  integer :: ind

  ! Calculates tidal tensor on a given cell
  ! Arguments:
  !    - nbors_father_cells_nbor: neighbouring father cells (only one of them).
  !    - tidal: to give back the tidal.
  !

  integer :: icell_minus_minus, icell_minus_plus
  integer :: icell_plus_minus, icell_plus_plus
  integer :: igrid_minus_minus, igrid_minus_plus
  integer :: igrid_plus_minus, igrid_plus_plus

  integer :: idim
  integer :: inbor
  integer :: icell_nbor
  integer :: igrid_nbor

  integer :: igshift

  real(dp) :: dx, oneoverdx2

  integer :: i

  real(dp) :: phi_nbor
  real(dp) :: sf2phi                                                                    ! RAY_RAMSES

  ! For finding 27 neighbours
  integer ,dimension(1:nvector) :: father_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids

  ! Flag to identify components that must be calculated in coarse level
  integer, dimension(1:6) :: to_do

  ! Set constants
  dx  = 0.5D0**ilevel
  oneoverdx2 = 1.0D0/(dx*dx)

  if(extradof .and. extradof2 .and..not.extradof3) then                                 ! RAY_RAMSES
     sf2phi = 3.0D0*beta_dgp/(2.0D0*beta_cvg)*alpha_cvg                                 ! RAY_RAMSES
  else if(.not.extradof .and. extradof2 .and. extradof3) then                           ! RAY_RAMSES
     sf2phi = 1.0d0+alpha_cvg/beta_cvg                                                  ! RAY_RAMSES
  endif                                                                                 ! RAY_RAMSES

  ! Initialise stuff
  do i=1,6
     tidal(i) = 0.0D0
     to_do(i) = 0  ! assume we will be able to calculate all componets in fine grid
  enddo

  ! Diagonal components of tidal tensor
  do idim=1,ndim  ! This define which component of tidal tensor we calculate
     do inbor=1,2  ! This defines two directions

        igshift = ray_iii(idim,inbor,ind)
        if(igshift==0) then
           igrid_nbor = igrid
        else
           igrid_nbor = son(nbor(igrid,igshift))
        end if

        ! If we do not have the neighbour cell, flag the component of tidal()
        if(igrid_nbor .eq. 0) then
           tidal(idim) = 0.0D0
           to_do(idim) = 1
           exit
        else
           icell_nbor = igrid_nbor+(ncoarse+(ray_jjj(idim,inbor,ind)-1)*ngridmax)
           if(extradof .and. extradof2 .and..not.extradof3) then                        ! RAY_RAMSES
              phi_nbor = phi(icell_nbor)+sf2phi*sf(icell_nbor)                          ! RAY_RAMSES
           else if(.not.extradof .and. extradof2 .and. extradof3) then                  ! RAY_RAMSES
              phi_nbor = phi(icell_nbor)*sf2phi                                         ! RAY_RAMSES
           else                                                                         ! RAY_RAMSES
              phi_nbor = phi(icell_nbor)                                                ! RAY_RAMSES
           endif                                                                        ! RAY_RAMSES
        endif

        tidal(idim) = tidal(idim)+phi_nbor
     enddo

     ! Subtract twice the central value and normalize
     if(to_do(idim).eq.0) then
        if(extradof .and. extradof2 .and..not.extradof3) then                           ! RAY_RAMSES
           tidal(idim) = (tidal(idim)-2.0D0*(phi(icell)+sf2phi*sf(icell)))*oneoverdx2   ! RAY_RAMSES
        else if(.not.extradof .and. extradof2 .and. extradof3) then                     ! RAY_RAMSES
           tidal(idim) = (tidal(idim)-2.0D0*(phi(icell)*sf2phi          ))*oneoverdx2   ! RAY_RAMSES
        else                                                                            ! RAY_RAMSES
           tidal(idim) = (tidal(idim)-2.0D0*(phi(icell)                 ))*oneoverdx2   ! RAY_RAMSES
        endif                                                                           ! RAY_RAMSES
     endif
  enddo

  ! Non diagonal components of tidal tensor
  ! Find neighbours of the cell
  father_cell(1) = father(igrid)
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,1,ilevel)
  ! ALL THIS STUFF BELOW COULD BE WRITTEN IN A LOOP.
  ! THE LOOP WILL BE A BIT OBSCURE FOR THE NORMAL USER.
  ! DO WE WANT THAT OR WE KEEP IT IN THIS MORE READABLE FORM?
  ! BOTH WAYS SHOULD HAVE SAME SPEED.        
  !dxdy:
  igrid_minus_minus = son( nbors_father_cells(1, ray_kkk(10,ind) )  )
  igrid_plus_minus  = son( nbors_father_cells(1, ray_kkk(12,ind) )  )
  igrid_plus_plus   = son( nbors_father_cells(1, ray_kkk(18,ind) )  )
  igrid_minus_plus  = son( nbors_father_cells(1, ray_kkk(16,ind) )  )
  if(igrid_minus_minus .eq.0 .or. &
       igrid_plus_minus.eq.0 .or. &
       igrid_plus_plus .eq.0 .or. &
       igrid_minus_plus.eq.0) then
     to_do(4) = 1
  else
     icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll(10,ind)-1)*ngridmax)
     icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll(12,ind)-1)*ngridmax)
     icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(18,ind)-1)*ngridmax)
     icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(16,ind)-1)*ngridmax)     
     if(extradof .and. extradof2 .and..not.extradof3) then                              ! RAY_RAMSES
        tidal(4) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0 + &                     ! RAY_RAMSES
                   (sf (icell_minus_minus) - &                                          ! RAY_RAMSES
                    sf (icell_plus_minus ) + &                                          ! RAY_RAMSES
                    sf (icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    sf (icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else if(.not.extradof .and. extradof2 .and. extradof3) then                        ! RAY_RAMSES
        tidal(4) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else                                                                               ! RAY_RAMSES
        tidal(4) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0                         ! RAY_RAMSES
     endif                                                                              ! RAY_RAMSES
  endif

  !dxdz:
  igrid_minus_minus = son(  nbors_father_cells(1, ray_kkk( 4,ind) )  )
  igrid_plus_minus  = son(  nbors_father_cells(1, ray_kkk( 6,ind) )  )
  igrid_plus_plus   = son(  nbors_father_cells(1, ray_kkk(24,ind) )  )
  igrid_minus_plus  = son(  nbors_father_cells(1, ray_kkk(22,ind) )  )
  if(igrid_minus_minus .eq.0 .or. &
       igrid_plus_minus.eq.0 .or. &
       igrid_plus_plus .eq.0 .or. &
       igrid_minus_plus.eq.0) then
     to_do(5) = 1
  else
     icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll( 4,ind)-1)*ngridmax)
     icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll( 6,ind)-1)*ngridmax)     
     icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(24,ind)-1)*ngridmax)
     icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(22,ind)-1)*ngridmax)     
     if(extradof .and. extradof2 .and..not.extradof3) then                              ! RAY_RAMSES
        tidal(5) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0 + &                     ! RAY_RAMSES
                   (sf (icell_minus_minus) - &                                          ! RAY_RAMSES
                    sf (icell_plus_minus ) + &                                          ! RAY_RAMSES
                    sf (icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    sf (icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else if(.not.extradof .and. extradof2 .and. extradof3) then                        ! RAY_RAMSES
        tidal(5) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else                                                                               ! RAY_RAMSES
        tidal(5) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0                         ! RAY_RAMSES
     endif
  endif

  !dydz:
  igrid_minus_minus = son(  nbors_father_cells(1, ray_kkk( 2,ind) )  )
  igrid_plus_minus  = son(  nbors_father_cells(1, ray_kkk( 8,ind) )  )
  igrid_plus_plus   = son(  nbors_father_cells(1, ray_kkk(26,ind) )  )
  igrid_minus_plus  = son(  nbors_father_cells(1, ray_kkk(20,ind) )  )
  if(igrid_minus_minus .eq.0 .or. &
       igrid_plus_minus.eq.0 .or. &
       igrid_plus_plus .eq.0 .or. &
       igrid_minus_plus.eq.0) then
     to_do(6) = 1
  else
     icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll( 2,ind)-1)*ngridmax)
     icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll( 8,ind)-1)*ngridmax)
     icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(26,ind)-1)*ngridmax)
     icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(20,ind)-1)*ngridmax)
     if(extradof .and. extradof2 .and..not.extradof3) then                              ! RAY_RAMSES
        tidal(6) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0 + &                     ! RAY_RAMSES
                   (sf (icell_minus_minus) - &                                          ! RAY_RAMSES
                    sf (icell_plus_minus ) + &                                          ! RAY_RAMSES
                    sf (icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    sf (icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else if(.not.extradof .and. extradof2 .and. extradof3) then                        ! RAY_RAMSES
        tidal(6) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0*sf2phi                  ! RAY_RAMSES
     else                                                                               ! RAY_RAMSES
        tidal(6) = (phi(icell_minus_minus) - &                                          ! RAY_RAMSES
                    phi(icell_plus_minus ) + &                                          ! RAY_RAMSES
                    phi(icell_plus_plus  ) - &                                          ! RAY_RAMSES
                    phi(icell_minus_plus )   )*oneoverdx2/4.0D0                         ! RAY_RAMSES
     endif                                                                              ! RAY_RAMSES
  endif

  ! Do the components that we could not calculate on the fine grid
  ! IS THERE A MORE ELEGANT WAY OF WRITING THIS LOOP?
  do i=1,6
     if(to_do(i).eq.1) then
        call tidal_tensor_from_coarse(ilevel, igrid, ind, nbors_father_cells, to_do, tidal)
        exit
     endif
  enddo

end subroutine tidal_tensor

!====================================
! tidal_tensor_from_coarse
!====================================
subroutine tidal_tensor_from_coarse(ilevel,igrid,ind,nbors_father_cells,to_do,tidal)
  use amr_commons
  use poisson_commons
  use ray_commons
  use ray_parameters ! for the test
  use extradof_commons
  implicit none

  real(dp), dimension(1:6) :: tidal
  integer :: ilevel
  integer :: igrid
  integer :: ind
  integer, dimension(1:6) :: to_do
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells

  ! Calculates tidal tensor on a given quarter of a cell by interpolating
  ! from the tidal tensor that is calculated one level below.
  !
  ! Arguments:
  !    - tidal: to return the tidal.
  !    - to_do:  flags to identity components to calculate
  !    - nbors_father_cells: we calculated this above, so no need to calculate
  !      it here again.
  !

  integer :: idim_coarse
  integer :: inbor_coarse
  integer :: icell_coarse_nbor
  integer :: igrid_coarse_nbor

  integer :: igshift

  real(dp) :: dx, oneover2dx2

  integer :: icell_minus_minus, icell_minus_plus
  integer :: icell_plus_minus, icell_plus_plus
  integer :: igrid_minus_minus, igrid_minus_plus
  integer :: igrid_plus_minus, igrid_plus_plus

  real(dp) :: phi_nbor
  real(dp) :: sf2phi                                                                    ! RAY_RAMSES

  ! For finding 27 neighbours
  !integer,dimension(1:nvector) :: ifather
  !integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids

  ! For finding 27 neighbours of neighbours
  integer,dimension(1:nvector) :: ifather_coarse
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells_coarse
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids_coarse

  real (dp) :: coeff

  integer :: icell_coarse, igrid_coarse, ind_average, ind_coarse  !, ind_father
  real (dp) :: tidal_one_cell

  ! Set constants
  dx = 0.5D0**ilevel
  oneover2dx2 = 1.0D0/((2.0D0*dx)*(2.0D0*dx))  ! Note that there is a 2!

  if(extradof .and. extradof2 .and..not.extradof3) then                                 ! RAY_RAMSES
     sf2phi = 3.0D0*beta_dgp/(2.0D0*beta_cvg)*alpha_cvg                                 ! RAY_RAMSES
  else if(.not.extradof .and. extradof2 .and. extradof3) then                           ! RAY_RAMSES
     sf2phi = 1.0d0+alpha_cvg/beta_cvg                                                  ! RAY_RAMSES
  endif                                                                                 ! RAY_RAMSES

  ! Identify <<all>> the coarse neighbours that will be needed
  do ind_average=1,twotondim
     ! Identify the cell
     icell_coarse = nbors_father_cells(1, ray_ooo(ind_average,ind))
     ind_coarse = (icell_coarse-ncoarse-1)/ngridmax+1
     igrid_coarse = icell_coarse-ncoarse-(ind_coarse-1)*ngridmax
     ! Store the father
     ifather_coarse(ind_average) = father(igrid_coarse)
  enddo
  call get3cubefather(ifather_coarse,nbors_father_cells_coarse,nbors_father_grids_coarse,8,ilevel-1)
    
  ! Third order interpolation taken from interpol_phi
  do ind_average=1,twotondim

     ! Identify the father cell (1 of 8)
     ! (This is the coarse cell where we want the Laplacian)
     icell_coarse = nbors_father_cells(1, ray_ooo(ind_average,ind))
     ind_coarse = (icell_coarse-ncoarse-1)/ngridmax+1
     igrid_coarse = icell_coarse-ncoarse-(ind_coarse-1)*ngridmax

     ! Calculate weight
     coeff=ray_bbbb(ind_average)

     ! Diagonal components of tidal tensor
     do idim_coarse=1,ndim  ! This define which component of tidal tensor we calculate
        if(to_do(idim_coarse).eq. 1) then ! Check if we need to calculate this component
           tidal_one_cell = 0.0D0
           do inbor_coarse=1,2  ! This defines two directions

              igshift = ray_iii(idim_coarse,inbor_coarse,ind_coarse)
              if(igshift==0) then
                 igrid_coarse_nbor = igrid_coarse
              else
                 igrid_coarse_nbor = son(nbor(igrid_coarse,igshift))
              end if

              ! If we do not have the neighbour cell
              if(igrid_coarse_nbor .eq. 0) then
                 write(*,*) "tidal_tensor_from_coarse:  problem is ill defined.  Aborting."
                 call clean_stop()
              else
                 icell_coarse_nbor = igrid_coarse_nbor+(ncoarse+(ray_jjj(idim_coarse,inbor_coarse,ind_coarse)-1)*ngridmax)
                 if(extradof .and. extradof2 .and..not.extradof3) then                  ! RAY_RAMSES
                    phi_nbor = phi(icell_coarse_nbor)+sf(icell_coarse_nbor)*sf2phi      ! RAY_RAMSES
                 else if(.not.extradof .and. extradof2 .and. extradof3) then            ! RAY_RAMSES
                    phi_nbor = phi(icell_coarse_nbor)                      *sf2phi      ! RAY_RAMSES
                 else                                                                   ! RAY_RAMSES
                    phi_nbor = phi(icell_coarse_nbor)                                   ! RAY_RAMSES
                 endif                                                                  ! RAY_RAMSES
              endif

              tidal_one_cell = tidal_one_cell+phi_nbor
           enddo
           ! Subtract twice the central value and normalize
           if(extradof .and. extradof2 .and..not.extradof3) then                                                  ! RAY_RAMSES
              tidal_one_cell = (tidal_one_cell-2.0D0*(phi(icell_coarse)+sf(icell_coarse)*sf2phi))*oneover2dx2     ! RAY_RAMSES
           else if(.not.extradof .and. extradof2 .and. extradof3) then                                            ! RAY_RAMSES
              tidal_one_cell = (tidal_one_cell-2.0D0*(phi(icell_coarse)                 *sf2phi))*oneover2dx2     ! RAY_RAMSES
           else                                                                                                   ! RAY_RAMSES
              tidal_one_cell = (tidal_one_cell-2.0D0*(phi(icell_coarse)                        ))*oneover2dx2     ! RAY_RAMSES
           endif                                                                                                  ! RAY_RAMSES 

           tidal(idim_coarse) = tidal(idim_coarse)+coeff*tidal_one_cell

        endif  ! to_do==1

     enddo  ! idim

     ! Non diagonal components of tidal tensor
     ! ALL THIS STUFF BELOW COULD BE WRITTEN IN A LOOP.
     ! THE LOOP WILL BE A BIT OBSCURE FOR THE NORMAL USER.
     ! DO WE WANT THAT OR WE KEEP IT IN THIS MORE READABLE FORM?
     ! BOTH WAYS SHOULD HAVE SAME SPEED.
     !dxdy:
     if(to_do(4).eq.1) then
        igrid_minus_minus = son( nbors_father_cells_coarse(ind_average, ray_kkk(10,ind_coarse) )  )
        igrid_plus_minus  = son( nbors_father_cells_coarse(ind_average, ray_kkk(12,ind_coarse) )  )
        igrid_plus_plus   = son( nbors_father_cells_coarse(ind_average, ray_kkk(18,ind_coarse) )  )
        igrid_minus_plus  = son( nbors_father_cells_coarse(ind_average, ray_kkk(16,ind_coarse) )  )

        icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll(10,ind_coarse)-1)*ngridmax)
        icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll(12,ind_coarse)-1)*ngridmax)
        icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(18,ind_coarse)-1)*ngridmax)        
        icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(16,ind_coarse)-1)*ngridmax)

        if(extradof .and. extradof2 .and..not.extradof3) then                           ! RAY_RAMSES
           tidal(4) = tidal(4) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0 + &           ! RAY_RAMSES
                      coeff*(sf (icell_minus_minus) - &                                 ! RAY_RAMSES
                             sf (icell_plus_minus ) + &                                 ! RAY_RAMSES
                             sf (icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             sf (icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else if(.not.extradof .and. extradof2 .and. extradof3) then                     ! RAY_RAMSES
           tidal(4) = tidal(4) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else                                                                            ! RAY_RAMSES
           tidal(4) = tidal(4) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0               ! RAY_RAMSES
        endif                                                                           ! RAY_RAMSES
     endif

     !dxdz:
     if(to_do(5).eq.1) then
        igrid_minus_minus = son( nbors_father_cells_coarse(ind_average, ray_kkk( 4,ind_coarse) )  )
        igrid_plus_minus  = son( nbors_father_cells_coarse(ind_average, ray_kkk( 6,ind_coarse) )  )
        igrid_plus_plus   = son( nbors_father_cells_coarse(ind_average, ray_kkk(24,ind_coarse) )  )
        igrid_minus_plus  = son( nbors_father_cells_coarse(ind_average, ray_kkk(22,ind_coarse) )  )
        
        icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll( 4,ind_coarse)-1)*ngridmax)
        icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll( 6,ind_coarse)-1)*ngridmax)
        icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(24,ind_coarse)-1)*ngridmax)        
        icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(22,ind_coarse)-1)*ngridmax)

        if(extradof .and. extradof2 .and..not.extradof3) then                           ! RAY_RAMSES
           tidal(5) = tidal(5) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0 + &           ! RAY_RAMSES
                      coeff*(sf (icell_minus_minus) - &                                 ! RAY_RAMSES
                             sf (icell_plus_minus ) + &                                 ! RAY_RAMSES
                             sf (icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             sf (icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else if(.not.extradof .and. extradof2 .and. extradof3) then                     ! RAY_RAMSES
           tidal(5) = tidal(5) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else                                                                            ! RAY_RAMSES
           tidal(5) = tidal(5) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0               ! RAY_RAMSES
        endif                                                                           ! RAY_RAMSES
     endif

     !dydz:
     if(to_do(6).eq.1) then
        igrid_minus_minus = son( nbors_father_cells_coarse(ind_average, ray_kkk( 2,ind_coarse) )  )
        igrid_plus_minus  = son( nbors_father_cells_coarse(ind_average, ray_kkk( 8,ind_coarse) )  )
        igrid_plus_plus   = son( nbors_father_cells_coarse(ind_average, ray_kkk(26,ind_coarse) )  )
        igrid_minus_plus  = son( nbors_father_cells_coarse(ind_average, ray_kkk(20,ind_coarse) )  )

        icell_minus_minus = igrid_minus_minus + (ncoarse+(ray_lll( 2,ind_coarse)-1)*ngridmax)
        icell_plus_minus  = igrid_plus_minus  + (ncoarse+(ray_lll( 8,ind_coarse)-1)*ngridmax)
        icell_plus_plus   = igrid_plus_plus   + (ncoarse+(ray_lll(26,ind_coarse)-1)*ngridmax)        
        icell_minus_plus  = igrid_minus_plus  + (ncoarse+(ray_lll(20,ind_coarse)-1)*ngridmax)

        if(extradof .and. extradof2 .and..not.extradof3) then                           ! RAY_RAMSES
           tidal(6) = tidal(6) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0 + &           ! RAY_RAMSES
                      coeff*(sf (icell_minus_minus) - &                                 ! RAY_RAMSES
                             sf (icell_plus_minus ) + &                                 ! RAY_RAMSES
                             sf (icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             sf (icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else if(.not.extradof .and. extradof2 .and. extradof3) then                     ! RAY_RAMSES
           tidal(6) = tidal(6) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0*sf2phi        ! RAY_RAMSES
        else                                                                            ! RAY_RAMSES
           tidal(6) = tidal(6) + &                                                      ! RAY_RAMSES
                      coeff*(phi(icell_minus_minus) - &                                 ! RAY_RAMSES
                             phi(icell_plus_minus ) + &                                 ! RAY_RAMSES
                             phi(icell_plus_plus  ) - &                                 ! RAY_RAMSES
                             phi(icell_minus_plus )   )*oneover2dx2/4.0D0               ! RAY_RAMSES
        endif                                                                           ! RAY_RAMSES
     endif

  enddo ! over eight cells around

end subroutine tidal_tensor_from_coarse


