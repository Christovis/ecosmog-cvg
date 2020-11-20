subroutine init_ray(ilevel)
  use pm_commons
  use amr_commons
  use amr_parameters ! for arrays of neighbours (ndim)
  use ray_commons
  use ray_parameters

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer :: ilevel

  integer :: ngrid,ncache,igrid,i,ind,ind_cell,info,n_current,icpu
  integer, dimension(1:nvector) :: ind_grid
  integer, parameter :: alloc_chunk_size = 100
  logical  :: skip
  integer  :: ix,iy,iz

  real(dp) :: pi,rt2                                               ! numerical consts.

  real(dp) :: dx                                                   ! grid size
  real(dp) :: xg1,xg2,xg3,xgo,ygo,zgo                              ! grid centre coord

  real(dp) :: dax,day                                              ! ray grid size
  integer  :: nrays_ave                                            ! ray number on each CPU on average
  integer  :: nrays,nrays_tot
  real(dp) :: opening_x,opening_y,cof_x,cof_y                      ! ray geometry
  integer  :: ii,ij,i_max,i_min,j_max,j_min                        ! ray indexing
  real(dp) :: ray_x,ray_y,ray_z                                    ! ray Cartesan coord
  real(dp) :: r_c,t_c,p_c,t_s,t_lor,t_upr,p_lor,p_upr              ! ray angular coord
  real(dp) :: xlim1,xlim2,ylim1,ylim2,x_lor,x_upr,y_lor,y_upr      ! ray boundry

  real(dp) :: a_elp,c_elp,e_elp                                    ! ellipse properties
  real(dp) :: z_foc                                                ! z-coord of ellipse

  real(dp) :: tmp,tmp1,tmp2,chi_s                                  ! intermediate quantities

  integer, allocatable,dimension(:,:) :: tmparr1                   ! temporary array
  real(dp),allocatable,dimension(:,:) :: tmparr2                   ! temporary array

  ! For testing
  character*1000 file_name
  character(LEN=5)::nchar

  ! For defining chi_s
  real(dp) :: dcom_ray  ! this is a function
  real(dp) :: coverH0_ray
  
  ! For defining arrays of neighbours
  integer  :: j,k,mmm
  integer  :: jx,jy,jz
  real(dp) :: aa,bb,cc,dd

  integer  :: ncells

  !---------------------------------------
  ! Express distances in code units
  !---------------------------------------
  ray_z_obs = ray_z_obs/ray_Lbox          
  ray_x_obs = ray_x_obs/ray_Lbox
  ray_y_obs = ray_y_obs/ray_Lbox
  
  !---------------------------------------
  ! Open file for testing
  !---------------------------------------
  if(ray_step_test) then
     call title(myid,nchar)                                        ! myid = processor ID
     file_name = 'ray_'//trim(nchar)//'.dat'
     if(myid==1) write(*,*) "init_ray:  opening file for full evolution ", file_name
     open(unit=1000, file=file_name, form='formatted')
  endif

  !-----------------------------------------------------------------------------!
  ! This subroutine finds the rays managed by CPU #myid, and does the
  ! initialisation. 
  !-----------------------------------------------------------------------------!

  if(ilevel.ne.levelmin) return
  if(verbose) write(*,*)'Entering init_ray'

  coverH0_ray = 299792.458D0/100.D0
  chi_s = dcom_ray(0.0D0, ray_z_s, omega_m, omega_l, coverH0_ray, 1.0D-6)/ray_Lbox

  if(chi_s.le.dabs(ray_z_obs)+1.0D0 .and. chi_s.gt.dabs(ray_z_obs)) then
     ray_first_box = .true.
  else
     ray_first_box = .false. 
  end if

  ! Allocate ray communication commons
  allocate(ray_emission (1:ncpu))
  allocate(ray_reception(1:ncpu))
  do icpu=1,ncpu
     ray_emission (icpu)%ngrid = 0
     ray_reception(icpu)%ngrid = 0
  end do

  allocate(ex_father(1:ngridmax))  ! father cell index for killed grids

  ncells = ngridmax*twotondim + ncoarse
  allocate(ray_in_cell(1:ncells))  ! local ID of cells containing rays
  ray_in_cell = 0

  ! Deal with straight and non-straight rays differently
  if(ray_no_bending) then

     ! Some numerical constants
     pi  = 4.0D0*datan(1.0D0)
     rt2 = dsqrt(2.0D0)

     ! Light cone parameters
     opening_x = ray_opening_x*pi/180.0D0                          ! deg -> rad
     opening_y = ray_opening_y*pi/180.0D0                          ! deg -> rad
     cof_x     = ray_cof_x    *pi/180.0D0                          ! deg -> rad
     cof_y     = ray_cof_y    *pi/180.0D0                          ! deg -> rad
     dax       = opening_x/dble(ray_Nx-1)                          ! stepsize of x grid
     day       = opening_y/dble(ray_Ny-1)                          ! stepsize of y grid

     ! Number of rays of each CPU on average
     nrays_ave = ray_Nx*ray_Ny/ncpu+1 

     ! Initial allocation - may need to extend the arrays if necessary
     allocate(ray_id   (1:nrays_ave    ))
     allocate(ray_coord(1:nrays_ave,1:3))
     allocate(ray_grid (1:nrays_ave,1:2))
     allocate(tmparr1  (1:nrays_ave,1:2))

     ! <><><><><><><><><><><><><><><><><>< 
     ! Initialise ray fields array indices                         ! Baojiu-Feb-2019
     ! <><><><><><><><><><><><><><><><><><
     RAY_NFIELDS = 1                                               ! total number of fields to use
     i_ray_tidal = 0                                               ! index of first tidal field
     i_ray_grads = 0                                               ! index of first gradient field
     i_ray_phidt = 0                                               ! index of first potential field
     ! +++ add new stuff here to calculate other observables +++
     if(ray_do_shear.or.(ray_do_kappa.and.ray_kappa_method.ne.1)) then
        i_ray_tidal = RAY_NFIELDS+1
        RAY_NFIELDS = RAY_NFIELDS+6
     end if
     if(ray_do_deflt.or.(ray_do_kappa.and.ray_kappa_method.ne.2)) then
        i_ray_grads = RAY_NFIELDS+1
        RAY_NFIELDS = RAY_NFIELDS+3
     end if
     if(ray_do_iswrs) then
        i_ray_phidt = RAY_NFIELDS+1
        RAY_NFIELDS = RAY_NFIELDS+2
     end if
     ! +++ add new stuff here to calculate other observables +++

     ! <><><><><><><><><><><><><><><><><><><
     ! Initialise ray quantity array indices                       ! Baojiu-Feb-2019
     ! <><><><><><><><><><><><><><><><><><><
     RAY_NQUANTS = 0                                               ! total number of quantities to compute
     i_ray_kappa = 0                                               ! index of kappa quantities
     i_ray_shear = 0                                               ! index of shear quantities
     i_ray_iswrs = 0                                               ! index of iswrs quantities
     i_ray_deflt = 0                                               ! index of deflt quantities
     ! +++ add new stuff here to calculate other observables +++
     if(ray_do_kappa) then
        i_ray_kappa = RAY_NQUANTS+1
        RAY_NQUANTS = RAY_NQUANTS+2
     end if
     if(ray_do_shear) then
        i_ray_shear = RAY_NQUANTS+1
        RAY_NQUANTS = RAY_NQUANTS+2
     end if
     if(ray_do_deflt) then          
        i_ray_deflt = RAY_NQUANTS+1
        RAY_NQUANTS = RAY_NQUANTS+2
     end if
     if(ray_do_iswrs) then
        i_ray_iswrs = RAY_NQUANTS+1
        RAY_NQUANTS = RAY_NQUANTS+1
     end if
     ! +++ add new stuff here to calculate other observables +++

     ! <><><><><><><><><><><><><><><><><><><><
     ! Sanity check that ray tracing is needed                     ! Baojiu-Feb-2019
     ! <><><><><><><><><><><><><><><><><><><><
     if(RAY_NQUANTS.eq.0) then
        write(*,*) 'init_ray: ray tracing is only needed if RAY_NQUANTS>0.'
        call clean_stop
     end if
     
     ! <><><><><><><><><><><><><><><><><><><><
     ! Allocations of arrays with proper sizes                     ! Baojiu-Feb-2019
     ! <><><><><><><><><><><><><><><><><><><><
     if(RAY_NQUANTS.lt.3) then
        allocate(tmparr2(1:nrays_ave,1:3          ))
     else
        allocate(tmparr2(1:nrays_ave,1:RAY_NQUANTS))
     end if
     allocate(ray_quants(1:nrays_ave,1:RAY_NQUANTS))

     ray_id     = 0
     ray_grid   = 0
     ray_coord  = 0.0D0
     ! <><><><><><><><><><><><><><><><><><><><><><><><
     ! United ray-tracing quantities to a single array             ! Baojiu-Feb-2019
     ! <><><><><><><><><><><><><><><><><><><><><><><><
     ray_quants = 0.0D0

     ! Current size of array. 
     n_current = nrays_ave

     ! Initialization of nrays to 0; nrays is the number of rays for CPU myid.
     nrays=0

     ! Loop over ilevel from levelmin to nlevelmax.
     do ilevel=levelmin,nlevelmax
 
        if(numbtot(1,ilevel).eq.0) exit                            ! this level does not have any active grids; skip

        dx     = 0.5D0**ilevel                                     ! coarse level cell size in code unit
        ncache = active(ilevel)%ngrid                              ! number of active grids on this level

        ! Loop over grids by vector sweeps find number of rays.
        do igrid=1,ncache,nvector
 
           ngrid=MIN(nvector,ncache-igrid+1)

           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)

              ! Check if the grid is at the highest refinement level;
              ! Skip this grid iff all its son cells are refined.
              skip = .true.
              do ind=1,twotondim
                 ind_cell = ind_grid(i)+ncoarse+(ind-1)*ngridmax
                 if(son(ind_cell).eq.0) skip = .false.
              end do
              if(skip) cycle

              ! Cartesan coordinates of the grid centre in local box.
              xg1 = xg(ind_grid(i),1)
              xg2 = xg(ind_grid(i),2)
              xg3 = xg(ind_grid(i),3) 

              ! Cartesan coordinates of the grid centre wrt observer.
              xgo = xg1-ray_x_obs
              ygo = xg2-ray_y_obs
              zgo = xg3-ray_z_obs

              ! Radial coordinate of the grid centre wrt observer.
              r_c = dsqrt(xgo**2+ygo**2+zgo**2)

              ! Take care of the tiling geometry.
              if(ray_first_box) then
                 ! Check if the ray lies in a sphere containing the grid.
                 if(ray_chi_s.lt.r_c-rt2*dx .or. ray_chi_s.gt.r_c+rt2*dx) cycle
              else
                 ! Check if the ray lies in a grid with extreme z-coordinate.
                 if(zgo.gt.0.0D0 .and. xg3+1.5D0*dx.lt.1.0D0) cycle
                 if(zgo.lt.0.0D0 .and. xg3-1.5D0*dx.gt.1.0D0) cycle
              end if

              ! The angle spanned by the sphere's radius to the observer.
              t_s = dasin(rt2*dx/r_c)         

              ! Angular coordintes of the grid centre wrt observer.
              if(dsqrt(xgo**2+ygo**2).le.dabs(zgo)) then           ! theta_c
                 t_c = datan(dsqrt(xgo**2+ygo**2)/zgo)  
              else
                 t_c = 0.5D0*pi-datan(zgo/dsqrt(xgo**2+ygo**2))  
              end if
              p_c = dacos(xgo/dsqrt(xgo**2+ygo**2))                ! psi_c

              if(t_c+1.0D0.lt.1.0D0) t_c = pi      +t_c
              if(ygo+1.0D0.lt.1.0D0) p_c = pi*2.0D0-p_c

              ! Angular coordinates of sphere edges wrt observer.
              if(t_c.le.t_s) then 

                 t_lor = 0.0D0     ! theta for lower sphere edge wrt observer
                 t_upr = t_c+t_s   ! theta for upper sphere edge wrt observer
                 p_lor = 0.0D0     ! psi   for lower sphere edge wrt observer
                 p_upr = 2.0D0*pi  ! psi   for upper sphere edge wrt observer

                 i_max = ray_Nx/2+int(t_upr/dax+0.5D0)
                 i_min = ray_Ny/2-int(t_upr/day+0.5D0)+1
                 i_min = min(i_max,i_min)
                 j_max = i_max
                 j_min = i_min

              else

                 t_lor = t_c-t_s ! theta for lower sphere edge wrt observer
                 t_upr = t_c+t_s ! theta for upper sphere edge wrt observer

                 ! z coordinate of the focus of the ellipse further away from centre:
                 ! The ellipse is the conic curve cut perpendicular to z direction on 
                 ! the cone containing the above-mentioned sphere; the sphere is used 
                 ! as a Dandelin spheres to find the focus of the ellipse. Note that 
                 ! the x & y coordinates of the focus are just xgo & ygo. The ellipse
                 ! can be projected to the x-y plane to assist analysis.
                 if(t_upr.lt.0.5D0*pi) z_foc = zgo-rt2*dx
                 if(t_lor.gt.0.5D0*pi) z_foc = zgo+rt2*dx
                 
                 tmp1  = dabs(dsqrt(xgo**2+ygo**2)-(z_foc*dtan(t_lor)))
                 tmp2  = dabs(dsqrt(xgo**2+ygo**2)-(z_foc*dtan(t_upr)))

                 ! Find important quantities of the ellipse
                 a_elp = 0.5D0*dabs(tmp1+tmp2) ! ellipse:length of semi-major axis
                 c_elp = 0.5D0*dabs(tmp1-tmp2) ! ellipse:linear eccentricity
                 e_elp = c_elp/a_elp           ! ellipse:eccentricity

                 ! Find the boundary containing the project ellipse in the x-y plane, 
                 ! whose four sides are parallel to the x and y axes.
                 ! angle at which the ellipse is farthest/closest from/to y-axis
                 tmp1  = dasin(-e_elp*dsin(p_c))
                 ! angle at which the ellipse is farthest/closest from/to x-axis
                 tmp2  = dacos(-e_elp*dcos(p_c))

                 xlim1 = xgo + &
                         a_elp*(1.0D0-e_elp**2) / &
                         (1.0D0+e_elp*dcos(pi-tmp1-p_c))*dcos(pi-tmp1)
                 xlim2 = xgo + &
                         a_elp*(1.0D0-e_elp**2)/ &
                         (1.0D0+e_elp*dcos(tmp1-p_c))*dcos(tmp1)
                 ylim1 = ygo + &
                         a_elp*(1.0D0-e_elp**2) / &
                         (1.0D0+e_elp*dcos(2.0D0*pi-tmp2-p_c))*dsin(2.0D0*pi-tmp2)
                 ylim2 = ygo + &
                         a_elp*(1.0D0-e_elp**2) / &
                         (1.0D0+e_elp*dcos(tmp2-p_c))*dsin(tmp2)

                 x_lor = dmin1(xlim1,xlim2)  ! lowest  x-coord on the ellipse
                 x_upr = dmax1(xlim1,xlim2)  ! highest x-coord on the ellipse
                 y_lor = dmin1(ylim1,ylim2)  ! lowest  y-coord on the ellipse
                 y_upr = dmax1(ylim1,ylim2)  ! highest y-coord on the ellipse

                 ! Find the i & j indices of the rays that fall in this boundary
                 i_max = int((datan(x_upr/dabs(z_foc))+0.5D0*opening_x)/dax)+1
                 i_min = int((datan(x_lor/dabs(z_foc))+0.5D0*opening_x)/dax)+2
                 j_max = int((datan(y_upr/dabs(z_foc))+0.5D0*opening_y)/day)+1
                 j_min = int((datan(y_lor/dabs(z_foc))+0.5D0*opening_y)/day)+2

                 if(i_min.gt.ray_Nx+1 .or. i_max.lt.1) cycle  ! outside boundary
                 if(j_min.gt.ray_Ny+1 .or. j_max.lt.1) cycle  ! outside boundary

                 if(i_max.eq.1) i_min = 1
                 if(j_max.eq.1) j_min = 1
                 if(i_min.eq.2) i_min = 1
                 if(j_min.eq.2) j_min = 1

                 if(i_min.gt.i_max .or. j_min.gt.j_max) cycle      ! grid falls between neighboring rays

              end if

              i_max = min(i_max,ray_Nx)
              j_max = min(j_max,ray_Ny)    
              i_min = max(i_min,1)
              j_min = max(j_min,1)

              ! Go through all potential rays to check if they are in this grid
              do ii=i_min,i_max
                 do ij=j_min,j_max

                    tmp = dsqrt(1.0D0+(dtan((ii-1)*dax-0.5D0*opening_x))**2+(dtan((ij-1)*day-0.5D0*opening_y))**2)
                    ! Find the Cartesan coordinates of the ray
                    if(ray_first_box) then
                       ray_z = ray_chi_s/tmp
                    else
                       ray_z = dabs(ray_z_obs) + 1.0D0
                    end if
                    ray_x = ray_z*dtan((ii-1)*dax-0.5D0*opening_x)
                    ray_y = ray_z*dtan((ij-1)*day-0.5D0*opening_y)

                    if(zgo.lt.0.0D0) ray_z = -ray_z

                    ! Check if the ray lies in the grid boundary 
                    if(ray_z.gt.zgo+dx .or. ray_z.le.zgo-dx) cycle

                    if(xgo.le.-dx .and. ray_x.ge.xgo+dx) cycle
                    if(xgo.gt.-dx .and. ray_x.gt.xgo+dx) cycle
                    if(xgo.lt. dx .and. ray_x.le.xgo-dx) cycle
                    if(xgo.ge. dx .and. ray_x.lt.xgo-dx) cycle
                    if(ygo.le.-dx .and. ray_y.ge.ygo+dx) cycle
                    if(ygo.gt.-dx .and. ray_y.gt.ygo+dx) cycle
                    if(ygo.lt. dx .and. ray_y.le.ygo-dx) cycle
                    if(ygo.ge. dx .and. ray_y.lt.ygo-dx) cycle

                    ! Check if the grid is at the highest refinement level
                    ix = int((ray_x-(xgo-dx))/dx)
                    iy = int((ray_y-(ygo-dx))/dx)
                    iz = int((ray_z-(zgo-dx))/dx)
                    if(.not.ray_first_box .and. iz.eq.2) iz = 1       ! the case of tiling: ray_z = 1.0
                    ind = 4*iz+2*iy+ix+1                              ! local index of son cell containing the ray
                    if(ind.lt.1 .or. ind.gt.twotondim) then
                       write(*,*) 'init_ray: incorrect cell index containing the ray; please check.',ix,iy,iz 
                       stop
                    end if
                    ind_cell = ind_grid(i)+ncoarse+(ind-1)*ngridmax
                    if(son(ind_cell).gt.0) cycle

                    nrays = nrays+1                                   ! Incretment ray number of myid

                    ! If nrays exceeds current array size, increase the latter by alloc_chunk_size
                    if(nrays.gt.n_current) then
                       call extend_arrays_if_too_small(alloc_chunk_size)
                       n_current = n_current+alloc_chunk_size
                    end if

                    ! Initializes ray-related arrays
                    ray_id   (nrays  )    = (ij-1)*ray_Nx+ii
                    ray_grid (nrays,1)    = ind_grid(i)
                    ray_grid (nrays,2)    = myid*1000+ind*100+ilevel
                    if(ray_first_box) then
                       ray_coord(nrays,1) = ray_chi_s
                    else
                       ray_coord(nrays,1) = dabs(ray_z)*tmp
                    end if
                    ray_coord(nrays,2)    = datan(dsqrt(ray_x**2+ray_y**2)/dabs(ray_z))
                    ray_coord(nrays,3)    = dacos(ray_x/dsqrt(ray_x**2+ray_y**2))
                    if(ray_y.lt.0.0D0) then
                       ray_coord(nrays,3) = 2.0D0*pi-ray_coord(nrays,3)
                    end if

!                   write(*,999) ray_grid(nrays,1),ind,ray_id(nrays),myid,ilevel,ray_coord(nrays,1),ray_coord(nrays,2),ray_coord(nrays,3)                        
                 end do ! loop over candidate rays
              end do ! loop over candidate rays
 
           end do ! loopo over active grids on ilevel (inner loop)
        end do ! loop over active grids on ilevel (outer loop)

     end do ! loop over level

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nrays,nrays_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     nrays_tot=nrays
#endif
     if(nrays_tot.ne.ray_Nx*ray_Ny) then
        write(*,*) 'init_ray: incorrect ray number; check if the ray bundle is all inside the box.'
        stop
     end if

     ray_nrays = nrays

     write(*,888) nrays_tot,myid,nrays

     deallocate(tmparr1)
     deallocate(tmparr2)  

  else
     ! Nonstraight ray case: needs to read ray information from restart file
     if(ray_first_box) then
       
     end if
  end if

  ray_initialized = .true.

  if(ray_first_box) then
     ray_all_started = .true.
  else
     ray_all_started = .false.
  end if

  ray_all_finished = .false. 
  ray_end_of_sim   = .false. 

  !-----------------------------------------------------------
  ! Initialise arrays for identification of neighbouring nodes
  !-----------------------------------------------------------
  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  ! THIS IS ORIGINAL RAMSES.  DO NOT USE IT HERE.
  !ray_bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)
  ! NOTE: THE ORDERING IN THIS ARRAY IS NOT THE SAME AS IN THE STANDARD CODE!!!
  ray_bbbb(:)  =(/dd, cc, cc, bb, cc, bb, bb, aa/)
  
  do i=1,8
     do j=1,8
        
        iz = (i-1)/4
        iy = (i-1-iz*4)/2
        ix = (i-1-iz*4-iy*2)

        jz = (j-1)/4
        jy = (j-1-jz*4)/2
        jx = (j-1-jz*4-jy*2)

        ! Taken from force_fine (sampling positions in the 3x3x3 father cell cube)
        ray_ccc(i,j) = 14-(ix-1)*(2*jx-1)-(iy-1)*(2*jy-1)*3-(iz-1)*(2*jz-1)*9
        
        ! Indexes are (target corner,cell).
        ! This gives eight cells surrounding a corner of a given cell.
        mmm = 1+(ix+jx)+(iy+jy)*3+(iz+jz)*9
        
        do k=1,8
        
           ! Find kkk(mmm)
           iz = (mmm-1)/9
           iy = (mmm-1-iz*9)/3
           ix = (mmm-1-iz*9-iy*3)
           
           jz = (k-1)/4
           jy = (k-1-jz*4)/2
           jx = (k-1-jz*4-jy*2)
           
           ! Indexes are (target corner, cell, ind of original cell)
           ray_kkk_mmm(i,j,k) = 14 + abs(ix-1)*((ix/2)+jx-1) + abs(iy-1)*((iy/2)+jy-1)*3 + abs(iz-1)*((iz/2)+jz-1)*9
           ray_lll_mmm(i,j,k) = k  - abs(ix-1)*(     2*jx-1) - abs(iy-1)*(     2*jy-1)*2 - abs(iz-1)*(     2*jz-1)*4
           
        enddo
     enddo
  enddo

  do i=1,27
     do j=1,8
        iz = (i-1)/9
        iy = (i-1-iz*9)/3
        ix = (i-1-iz*9-iy*3)
        
        jz = (j-1)/4
        jy = (j-1-jz*4)/2
        jx = (j-1-jz*4-jy*2)
        
        ! Neighbour father cells for each ind when we want to move in direction 1.
        ray_kkk(i,j) = 14 + abs(ix-1)*((ix/2)+jx-1) + abs(iy-1)*((iy/2)+jy-1)*3 + abs(iz-1)*((iz/2)+jz-1)*9
        ! Ind of neighbouring cells.
        ray_lll(i,j) = j  - abs(ix-1)*(     2*jx-1) - abs(iy-1)*(     2*jy-1)*2 - abs(iz-1)*(     2*jz-1)*4
        
     enddo
  enddo

  ! This comes from standard Ramses:
  ray_iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); ray_jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ray_iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); ray_jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ray_iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); ray_jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ray_iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); ray_jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ray_iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); ray_jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ray_iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); ray_jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Indexes are target cell (1 to 8), ind of cell where we want to interpolate
  ! Result moves first in x, then in y and finally in z
  ! In this way, the weights given by bbbb are Ok.
  ! To be used by tidal_tensor_from_coarse
  ray_ooo(1,1:8)=(/14,14,14,14,14,14,14,14/)
  ray_ooo(2,1:8)=(/13,15,13,15,13,15,13,15/)
  ray_ooo(3,1:8)=(/11,11,17,17,11,11,17,17/)
  ray_ooo(4,1:8)=(/10,12,16,18,10,12,16,18/)
  ray_ooo(5,1:8)=(/ 5, 5, 5, 5,23,23,23,23/)
  ray_ooo(6,1:8)=(/ 4, 6, 4, 6,22,24,22,24/)
  ray_ooo(7,1:8)=(/ 2, 2, 8, 8,20,20,26,26/)
  ray_ooo(8,1:8)=(/ 1, 3, 7, 9,19,21,25,27/)

 
888 format('init_ray: nrays_tot =',I10,', nrays (myid =',I3,') =',I10)
999 format('Grid ID',I7,', cell ind',I3,' contains ray ID',I10,' on CPU ID',I5,', at level',I3,' and coord (',3(F10.6,' '),')')

contains

  !------------------------------------------------------------------------------
  ! Extend arrays if their sizes are not big enough (deallocation and allocation)
  !------------------------------------------------------------------------------
  subroutine extend_arrays_if_too_small(chunk_size)

    ! Allocation chunk size
    integer :: chunk_size,size_kappa

    tmparr1(1:n_current,1) = ray_id(1:n_current)
    deallocate(ray_id)
    allocate  (ray_id(1:n_current+chunk_size))
    ray_id = 0
    ray_id(1:n_current) = tmparr1(1:n_current,1)

    tmparr1(1:n_current,1:2) = ray_grid(1:n_current,1:2)
    deallocate(ray_grid)
    allocate  (ray_grid(1:n_current+chunk_size,1:2))
    ray_grid(:,1:2) = 0
    ray_grid(1:n_current,1:2) = tmparr1(1:n_current,1:2)

    tmparr2(1:n_current,1:3) = ray_coord(1:n_current,1:3)
    deallocate(ray_coord)
    allocate  (ray_coord(1:n_current+chunk_size,1:3))
    ray_coord = 0.0D0
    ray_coord(1:n_current,1:3) = tmparr2(1:n_current,1:3)

    ! <><><><><><><><><><><><><><><><><><><><><><>
    ! Increase the size of united array ray_quants                 ! Baojiu-Feb-2019
    ! <><><><><><><><><><><><><><><><><><><><><><>
    tmparr2(1:n_current,1:RAY_NQUANTS) = ray_quants(1:n_current,1:RAY_NQUANTS)
    deallocate(ray_quants                                      )
    allocate  (ray_quants(1:n_current+chunk_size,1:RAY_NQUANTS)) 
    ray_quants = 0.0D0  
    ray_quants(1:n_current,1:RAY_NQUANTS) = tmparr2(1:n_current,1:RAY_NQUANTS)    
    
    deallocate(tmparr1)
    deallocate(tmparr2)  

    allocate   (tmparr1(1:n_current+chunk_size,1:2          ))
    ! <><><><><><><><><><><><><><><><><><><>
    ! Reallocate array with appropriate size                       ! Baojiu-Feb-2019
    ! <><><><><><><><><><><><><><><><><><><>
    if(RAY_NQUANTS.lt.3) then
       allocate(tmparr2(1:n_current+chunk_size,1:3          ))
    else
       allocate(tmparr2(1:n_current+chunk_size,1:RAY_NQUANTS))
    end if

  end subroutine extend_arrays_if_too_small

end subroutine init_ray

