program main
use types,only:rprec
use param
use sim_param
use io,only:openfiles,output_loop,output_final,inflow_write,avg_stats
use output_slices,only:write_slices
use fft
use immersedbc
use test_filtermodule
use topbc,only:setsponge,sponge
use bottombc,only:wt_s,wq_s
use scalars_module,only:beta_scal,obukhov,theta_all_in_one,RHS_T,RHS_Tf
use scalars_module_q,only:q_all_in_one,RHS_q,RHS_qf
use atm_thermodynamics,only:thermo
use scalars_module2,only:timestep_conditions
$if ($LVLSET)
  use level_set, only : level_set_init, level_set_cylinder_CD,  &
                        level_set_smooth_vel
$endif
use debug_mod  !--just for debugging
use messages
implicit none

! MM
!real(rprec),parameter::ug=ug_dim/u_star,vg=vg_dim/u_star
! MMi

integer, parameter :: wbase = 100  ! controls the frequency of screen diagnostics
integer, parameter :: nenergy = 10 ! frequency of writes to check_ke.dat

integer :: jx,jy,jz,i,j,k,k_global,Nzz
integer ::jt_diurnal
$if ($MPI)
  integer :: ip, np, coords(1)
$endif

real(kind=rprec),allocatable,dimension(:,:,:)::ug,vg
real(kind=rprec):: rmsdivvel,ke,kestor,testu   
real(kind=rprec)::const,tt,omega
real (rprec):: force

real(kind=rprec),dimension(4)::timestep_vars
 
integer::counter=0,plot_count=0

!---------------------------------------------------------------------


$if ($MPI)
  !--check for consistent preprocessor & param.f90 definitions of 
  !  MPI and $MPI
  if (.not. USE_MPI) then
    write (*, *) 'inconsistent use of USE_MPI and $MPI'
    stop
  end if

  call mpi_init (ierr)
  call mpi_comm_size (MPI_COMM_WORLD, np, ierr)
  call mpi_comm_rank (MPI_COMM_WORLD, global_rank, ierr)
! SKS
! mpi_init --> Initialize the MPI execution environment
! mpi_comm_size --> Determines the size of the group associated with a
! communictor
! mpi_comm_rank --> Determines the rank of the calling process in the 
! communicator
! SKS

  !--check if run-time number of processes agrees with nproc parameter
  if (np /= nproc) then
    write (*, *) 'runtime number of procs = ', np,  &
                 ' not equal to nproc = ', nproc
    stop
  end if

  !--set up a 1d cartesian topology 
  call mpi_cart_create (MPI_COMM_WORLD, 1, (/ nproc /), (/ .false. /),  &
                        .true., comm, ierr)
  !--slight problem here for ghost layers:
  !  u-node info needs to be shifted up to proc w/ rank "up",
  !  w-node info needs to be shifted down to proc w/ rank "down"
  call mpi_cart_shift (comm, 0, 1, down, up, ierr)
  call mpi_comm_rank (comm, rank, ierr)
  call mpi_cart_coords (comm, rank, 1, coords, ierr)
  coord = coords(1)  !--use coord (NOT rank) to determine global position

  write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use

  !--rank->coord and coord->rank conversions
  do ip = 0, np-1
    call mpi_cart_rank (comm, (/ ip /), rank_of_coord(ip), ierr)
    call mpi_cart_coords (comm, ip, 1, coord_of_rank(ip), ierr)
  end do

  write (*, *) 'Hello! from process with coord = ', coord

  !--set the MPI_RPREC variable
  if (rprec == kind (1.e0)) then
    MPI_RPREC = MPI_REAL
    MPI_CPREC = MPI_COMPLEX
  else if (rprec == kind (1.d0)) then
    MPI_RPREC = MPI_DOUBLE_PRECISION
    MPI_CPREC = MPI_DOUBLE_COMPLEX
  else
    write (*, *) 'error defining MPI_RPREC/MPI_CPREC'
    stop
  end if
$else
  if (nproc /= 1) then
    write (*, *) 'nproc /=1 for non-MPI run is an error'
    stop
  end if
  if (USE_MPI) then
    write (*, *) 'inconsistent use of USE_MPI and $MPI'
    stop
  end if

  !--leave this blank or put in coord
  !write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use
  chcoord = ''

$endif

!-----------------------
call read_namelist
!---------------------------

! -------- initialize time -----------
tt=0
!----------------------------

if( .not. allocated(ug)) allocate(ug(ld,ny,1:((nz-1)*nproc)))
if( .not. allocated(vg)) allocate(vg(ld,ny,1:((nz-1)*nproc)))

! Assign geostrophic winds
Nzz=(nz-1)*nproc
  do k=1,Nzz
            ug(:,:,k)=ug0/u_star
            vg(:,:,k)=vg0/u_star
  end do

! surface BCs in bottombc
call patches ()
!-----------------------

! initialize velocity and/or scalar fields (or read from file) in initial.f90
call initial()
!-----------------------------------------

!--could move this into something like initial ()
$if ($LVLSET)
  call level_set_init ()
$endif

! formulate the fft plans--may want to use FFTW_USE_WISDOM
! initialize the kx,ky arrays
call init_fft()

! Open output files      
  call openfiles()
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  print *,'Starting from time = ',jt_total
  end if
!--initialize test filter
!--this is used for lower BC, even if no dynamic model
!TSC standard dynamic or Lagrangian
  call test_filter_init (2._rprec * filter_size, G_test)

if (model == 3 .or. model == 5 .or. model == 6 .or. model == 7) then  !--scale dependent dynamic
  call test_filter_init (4._rprec * filter_size, G_test_test)
end if

if (ubc == 1) then
    call setsponge()
    print *,'sponge value calculated for damping layer'
else
    sponge=0._rprec
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  print *, 'Number of timesteps', nsteps
  print *, 'dt = , dt_dim = ', dt, dt_dim
  print *, 'Nx, Ny, Nz, Nz_total ', nx, ny, nz, nz_tot
  print *, 'Lx, Ly, Lz, Lz_total = ', L_x, L_y, L_z, L_z*nproc
  
  if (USE_MPI) print *, 'Number of processes = ', nproc
  print *, 'sampling stats every ', c_count, ' timesteps'
  print *, 'writing stats every ', p_count, ' timesteps'
  if (molec) print*, 'molecular viscosity (dimensional) ', nu_molec
end if

! MPI: u,v,w should be set at jz = 0:nz before getting here, except
! bottom process which is BOGUS (starts at 1)
! time Loop
do jt=1,nsteps   

  jt_total = jt_total + 1  !--moved from io.f90

  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
          print *, 'Coord 0 is at timestep', jt_total
  end if

  tt=tt+dt      ! advance total time

     RHSx_f = RHSx
     RHSy_f = RHSy
     RHSz_f = RHSz
  !end if

! Call obukhov to calculate the MO functions !!
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) call obukhov()
! if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) call obukhov(jt)

  !--no MPI here yet
  if (use_bldg) then
    call building_interp (u, v, w, .04_rprec, 3)
    call building_interp (dudx, dvdx, dwdx, .04_rprec, 3)
    call building_interp (dudy, dvdy, dwdy, .04_rprec, 3)
  end if

  ! kill oddballs and calculate horizontal derivatives

  ! except on bottom process (0 level set to BOGUS, starts at 1)
  call filt_da (u, dudx, dudy)
  call filt_da (v, dvdx, dvdy)
  call filt_da (w, dwdx, dwdy)

   ! finite differences
   !--MPI: on exit of ddz_uv, have dudz, dvdz at 1:nz, except
   !  bottom process has 2:nz
   call ddz_uv(dudz,u)
   call ddz_uv(dvdz,v)
   !--MPI: on exit of ddz_w, have dwdz at 0:nz-1, except top process
   !  has 0:nz, and bottom process has 1:nz-1
   call ddz_w(dwdz,w)

!TS calculate wall stress and calculate derivatives at wall
   if (dns_bc) then
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress_dns ()
     end if
   else
!TS "impose" wall stress and calculate derivatives at wall
     if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
       call wallstress ()  !--provides txz, tyz, dudz, dvdz at jz=1
                           !--MPI: bottom process only
     end if
     if(use_bldg) call walldudx_building
   end if

! compute turbulent viscosity (const.)
  if (dns_bc .and. molec) then
    call dns_stress(txx,txy,txz,tyy,tyz,tzz)
  else
    !--MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz
      ! -- MMC : All of the SGS Models are copmuted here.
    call sgs_stag()
  end if

  if(use_bldg)then
     call wallstress_building(txy,txz,tyz)
     call building_mask(u,v,w)
  endif

  if(S_FLAG.and.(jt.GE.SCAL_INIT))  then
     call theta_all_in_one
     call q_all_in_one
     call thermo()
  else
     beta_scal=0._rprec
  end if

!xx----VK -ADDED FOR SCALARS !! --xxxxxx

  $if ($MPI)
     !--exchange ghost-node information for tij
     !--send stuff up to ghost nodes
     !--move this into sgs_stag?
     call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 6,   &
                        tzz(:, :, 0), ld*ny, MPI_RPREC, down, 6,  &
                        comm, status, ierr)
  $endif

! compute divergence of SGS shear stresses     
! note: the divt's and the diagonal elements of t are equivalenced!
!--actually, they are not equivalenced in this version

  !--provides divtz 1:nz-1
  call divstress_uv(divtx, txx, txy, txz)
  call divstress_uv(divty, txy, tyy, tyz)
  !--provides divtz 1:nz-1, except 1:nz at top process
  call divstress_w(divtz, txz, tyz, tzz)

  if (VERBOSE) write (*, *) 'main about to call convec'

  !--provides RHS{x,y,z} 1:nz-1
  call convec(RHSx,RHSy,RHSz)

  if (use_bldg) call building_mask (u,v,w)

! Compute preliminary RHS matrices for pressure calculation
  RHSx(:, :, 1:nz-1) = -RHSx(:, :, 1:nz-1) - divtx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = -RHSy(:, :, 1:nz-1) - divty(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = -RHSz(:, :, 1:nz-1) - divtz(:, :, 1:nz-1)

 
  if (S_FLAG .and. (.not.passive_scalar)) then
    !--add buoyancy term...only valid for theta
    RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) + beta_scal(:, :, 1:nz-1)
  end if

  if (coriolis_forcing) then
    ! This is to put in the coriolis forcing using coriol,ug and vg as
    ! precribed in param.f90. (ug,vg) specfies the geostrophic wind vector
    ! Note that ug and vg are non-dimensional (using u_star in param.f90)
    RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) +                 &
                         coriol * v(:, :, 1:nz-1) - coriol *vg(:,:,(1+coord*(nz-1)):(nz-1)*(coord+1))

                         ! coriol * v(:, :, 1:nz-1) - ug_time_factor*coriol * vg
    RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) -                 &
                         coriol * u(:, :, 1:nz-1) + coriol *ug(:,:,(1+coord*(nz-1)):(nz-1)*(coord+1))

                         ! coriol * u(:, :, 1:nz-1) + ug_time_factor*coriol * ug
  end if

!XXXXXX%%%%%  Add damping terms to the momentum RHS %%%XXXXXXXXXXXX
  if (ubc==1 .and. damping_method==2) then !add damping terms to the momentum RHS
      do jz=1,nz-1 
      RHSx(1:nx,1:ny,jz)=RHSx(1:nx,1:ny,jz)-0.5*(sponge(jz)+sponge(jz+1))*&
                   (u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/(nx*ny))
      RHSy(1:nx,1:ny,jz)=RHSy(1:nx,1:ny,jz)-0.5*(sponge(jz)+sponge(jz+1))*&
                   (v(1:nx,1:ny,jz)-sum(v(1:nx,1:ny,jz))/(nx*ny))
      RHSz(1:nx,1:ny,jz)=RHSz(1:nx,1:ny,jz)-0.5*sponge(jz)*&
                   (w(1:nx,1:ny,jz)-sum(w(1:nx,1:ny,jz))/(nx*ny))
      end do
  elseif (ubc==1 .and. damping_method==1) then
      do jz=1,nz-1 
      RHSz(1:nx,1:ny,jz)=RHSz(1:nx,1:ny,jz)-sponge(jz)*&
                   w(1:nx,1:ny,jz)
      end do
  end if 
!XXXXXX%%%%%  Sponge/dampling block ends %%%%%%%%%%%%%%%XXXXXXXXXXXX

   !--calculate u^(*) (intermediate vel field)
   !  at this stage, p, dpdx_i are from previous time step
   !  (assumes old dpdx has NOT been added to RHSx_f, etc)
   !  we add force (mean press forcing) here so that u^(*) is as close
   !  to the final velocity as possible

   if (use_mean_p_force) then
     force = mean_p_force
   else
     force = 0._rprec
   end if

  if ((jt == 1) .and. (.not. initu)) then
    ! if initu, then this is read from the initialization file
    ! else for the first step put RHS_f=RHS
    !--i.e. at first step, take an Euler step
    RHSx_f=RHSx
    RHSy_f=RHSy
    RHSz_f=RHSz
  end if

   !--only 1:nz-1 are valid / It is actually u*
   u(:,:,1:nz-1) = u(:,:,1:nz-1) + dt*(tadv1*RHSx(:,:,1:nz-1) + &
                 tadv2 * RHSx_f(:, :, 1:nz-1) + force )
   v(:,:,1:nz-1) = v(:,:,1:nz-1) + dt*(tadv1*RHSy(:,:,1:nz-1) + &
                 tadv2 * RHSy_f(:, :, 1:nz-1) )
   w(:,:,1:nz-1) = w(:,:,1:nz-1) + dt*(tadv1*RHSz(:,:,1:nz-1) + &
                 tadv2 * RHSz_f(:, :, 1:nz-1) )

  $if ($MPI)
    !--after this point, u,v,w at jz = 0 are not useful, until updated
    u(:, :, 0) = BOGUS
    v(:, :, 0) = BOGUS
    w(:, :, 0) = BOGUS
  $endif

  !--this is an experiment
  u(:, :, nz) = BOGUS
  v(:, :, nz) = BOGUS
  w(:, :, nz) = BOGUS

  !--u, v, w at jz = nz are not useful either, except possibly w(nz), but that
  !  is supposed to zero anyway?
  !--this has to do with what bc are imposed on intermediate velocity

  !--solve Poisson equation for pressure
  !--provides p, dpdx, dpdy at 0:nz-1

  call press_stag_array (p, dpdx, dpdy)

  !--calculate dpdz here
  !--careful, p is not dimensioned the same as the others
  dpdz(1:nx,1:ny,1:nz-1) = (p(1:nx,1:ny,1:nz-1) - p(1:nx,1:ny,0:nz-2)) / dz
  dpdz(:, :, nz) = BOGUS

  !--if really wanted to, could avoid storing pressure gradients
  !  just add them directly to RHS in press_stag
  RHSx(:, :, 1:nz-1) = RHSx(:, :, 1:nz-1) - dpdx(:, :, 1:nz-1)
  RHSy(:, :, 1:nz-1) = RHSy(:, :, 1:nz-1) - dpdy(:, :, 1:nz-1)
  RHSz(:, :, 1:nz-1) = RHSz(:, :, 1:nz-1) - dpdz(:, :, 1:nz-1)

  call forcing ()
  
  !--provides u, v, w at 1:nz 
  call project ()
  
  $if ($MPI)
    !--exchange ghost-node information
    !--send stuff up to ghost layer (0) (for z-derivs)
    !--nz should already be in sync with 1 level: done in project()
    call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                       u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                       comm, status, ierr)

    call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                       v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                       comm, status, ierr)

    call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                       w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                       comm, status, ierr)
    call mpi_sendrecv (dudz(1, 1, nz-1), ld*ny, MPI_RPREC, up, 4,  &
                       dudz(1, 1, 0), ld*ny, MPI_RPREC, down, 4,   &
                       comm, status, ierr)
    call mpi_sendrecv (dvdz(1, 1, nz-1), ld*ny, MPI_RPREC, up, 5,  &
                       dvdz(1, 1, 0), ld*ny, MPI_RPREC, down, 5,   &
                       comm, status, ierr)


  $endif

  !--MPI: at this point, have u, v, w at 0:nz
!  if (modulo (jt, nenergy) == 0) call energy(ke)

!  call avg_stats ()  !--only does something once every n_avg_stats steps

  $if ($LVLSET)
    call level_set_cylinder_CD ()
  $endif

  if (modulo (jt, 100) == 0) then
    call rmsdiv (rmsdivvel)
    call timestep_conditions(timestep_vars(1),timestep_vars(2),timestep_vars(3),timestep_vars(4))
!   timestep_vars(1) is CFL and timestep_vars(2) is viscous stability limit
     if ((timestep_vars(1) .GT. 0.1_rprec) .OR. (timestep_vars(1) .LT. 0.0_rprec) ) then
    write (*, *) 'CFLu has exceeded 0.1 or became negative and is equal to = ', timestep_vars(1)
      
    stop
  end if 
 !   if (modulo (jt, 1) == 0) then

 !   call rmsdiv (rmsdivvel)
    if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
       if ((S_FLAG) .or. (coriolis_forcing)) then
         write (6, 7778) jt_total,sum(wt_s(1:nx,1:ny))/(nx*ny), sum(wq_s(1:nx,1:ny))/(nx*ny), &
                        coriolis_forcing,ug(1,1,Nzz)*u_star
         write (6, 7779) jt, dt, jt_total*(dt*z_i/u_star),rmsdivvel, &
                      timestep_vars(1),timestep_vars(2),timestep_vars(3),timestep_vars(4)
         ! SKS
        ! if ((timestep_vars(1) .OR. timestep_vars(2) .OR. timestep_vars(3)) .ge. 0.1) then
        !    print*,'** Caution: One of the CFLs > 0.1 **'
        ! end if
         ! SKS
       else
         write (6, 7777) jt_total, dt, rmsdivvel, timestep_vars(1),timestep_vars(2)
       end if
     end if
  end if
   
7777 format ('jt,dt,divu:',1x,i6.6,3(1x,e9.4))
7778 format ('jt_total,wt_s(K-m/s),wq_s(gm/kgs),coriolis, Ug(m/s):',(1x,i6.6,1x,f9.4,1x,f9.4,1x,L2,1x,f9.4))
7779 format ('jt,dt,time(s),divu,CFLu,CFLv,CFLw,visc_stab:',1x,i6.6,2(1x,E9.4),1x,E15.6,4(1x,F9.4))

!-------------------------------------------------
! SKS
 !  if (jt .ge. spectraCALC) then
 !     call spectra()
 !  end if

 !  if (mod(jt,p_count) == 0) then
 !     call calc_tke()
 !  end if
! SKS
!-------------------------------------------------


    if ((jt_total .GE. 80000) .AND. (jt_total .LE. 100000) .AND. modulo (jt, 20) == 0) then

            call write_slices()

    end if



  call output_loop ()

  if (write_inflow_file) call inflow_write () !--for creating inflow_BC file

end do  !--end time loop

if (output) call output_final (jt)

!call write_spectra()

$if ($MPI)
  call mpi_finalize (ierr)
$endif

end program main
