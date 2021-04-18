module param
use types,only:rprec

implicit none
$if ($MPI)
  include "mpif.h"
$endif

save

private rprec  !--this is dumb.
public

!--mpi stuff
$if ($MPI)
  $define $MPI_LOGICAL .true.
  $define $NPROC 384
$else
  $define $MPI_LOGICAL .false.
  $define $NPROC 1
$endif

logical, parameter :: USE_MPI = $MPI_LOGICAL

$undefine $MPI_LOGICAL

$if ($MPI)
  integer :: status(MPI_STATUS_SIZE)
$endif

character(*),parameter::path='./'



!--this stuff must be defined, even if not using MPIi
character (8) :: chcoord  !--holds character representation of coord
integer, parameter :: nproc = $NPROC  !--this must be 1 if no MPI
integer :: ierr
integer :: comm
integer :: up, down
integer :: global_rank
integer :: MPI_RPREC, MPI_CPREC
integer :: rank = -1   !--init to bogus (so its defined, even if no MPI)
integer :: coord = -1  !--same here
integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
!--end mpi stuff

logical,parameter::VERBOSE = .false. 
integer,parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
real(kind=rprec),parameter :: BOGUS = -1234567890._rprec
real(kind=rprec),parameter::pi=3.1415926535897932384626433_rprec
integer,parameter::spectraCALC=0 ! 125000
real(rprec),parameter::Pr=0.55_rprec

! -- these should be namelist parameters (default values)
integer,parameter::nx=360,ny=288,nz=384/nproc + 1

integer,parameter::nz_tot=(nz-1)*nproc + 1
integer,parameter::nx2=3*nx/2,ny2=3*ny/2
integer,parameter::lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

real(kind=rprec),parameter::L_x=5.0_rprec*pi,L_y=4.0_rprec*pi
real(kind=rprec),parameter::z_i=1000._rprec, L_z=3000._rprec/nproc
real(kind=rprec),parameter::dz=L_z/z_i/(nz-1)
real(kind=rprec),parameter::dx=L_x/nx,dy=L_y/ny

integer,protected::nsteps
logical::coriolis_forcing=.false.
real(kind=rprec),protected::u_star=0.0_rprec
real(kind=rprec),protected::dt_dim=0.1_rprec 
real(kind=rprec),protected::ug0=10.0,vg0=1.0
real(kind=rprec),protected:: dt 
real(kind=rprec),protected::coriol 
integer,protected::c_count,p_count
integer,protected::ubc=0, damping_method=1


!--- default values. if initu and initsc are false, initialize, otherwise read from previous
logical,protected:: initu=.true.,initsc=.true.,inilag=.false.,interp=.FALSE.

real(kind=rprec),protected:: inv_strength=0.0050_rprec,dTdz_top=0.0050_rprec
real(kind=rprec),protected :: T_scale=400._rprec,T_init=283.00_rprec ! unstable
real(kind=rprec),protected::inv_strength_q=-0.004_rprec,dqdz_top=-0.004_rprec
real(kind=rprec),protected :: q_init=1.00_rprec 
real(kind=rprec) :: q_scale=100._rprec

! ---- end of namelist parameters ------


real(rprec),parameter::cp=1005.0_rprec       ! specific heat capacity at constant pressure (J/kg/K)
real(rprec),parameter::Le=2.5008E06          ! Latent heat of vaporization (J/kg)
real(rprec),parameter::Rv=461.5_rprec        ! Gas constant for water vapor (J/kg/K)
real(rprec),parameter::Rd=287.04_rprec        ! Gas constant for dry air (J/kg/K)
real(rprec),parameter::rho_d=1.2923_rprec      ! density of dry air (kg/m^3)
real(rprec),parameter::rho_w=1000.0_rprec     ! water density (kg/m^3)
real(rprec),parameter::pr_surf=100000.0_rprec    ! Surface pressure (Pa)

real(rprec),parameter::vonk=0.41_rprec

integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates
logical,parameter::output=.true.
logical, parameter :: use_avgslice = .true.
integer, parameter :: average_dim_num = 1
real(rprec),parameter::nu_molec=1.14e-5_rprec
logical,parameter::use_bldg=.false.
logical,parameter::molec=.false.,sgs=.true.,dns_bc=.false.

! Model type:  1->Smagorinsky; 2->Dynamic; 3->Scale dependent
!              4->Lagrangian scale-sim     5-> Lagragian scale-dep
!              6->new Scale dependent dynamic
!              7->Lagrangian scale dependent for momentum and heat
! Models type: 1->static prandtl,          2->Dynamic
! Cs is the Smagorinsky Constant
! Co and nnn are used in the mason model for smagorisky coeff
integer,parameter::model=5,models=1,nnn=2,BETA_FLAG=1

! calc_lag_Ssim enables calculation of a Cs based on lag Ssim
! using BETA = 1 while using lag sdep for the real Cs. The Cs calc using
! calc_lag_ssim is only used for plotting purposes
! CAUTION : use only with model=5 and only when required (e.g. diurnal runs)
logical,parameter::calc_lag_Ssim=.false.

real(kind=rprec),parameter::Co=0.16_rprec
! Test filter type: 1->cut off 2->Gaussian 3->Top-hat
integer,parameter::ifilter=1
character (*), parameter :: lbc_mom = 'wall' !--'wall', 'stress free'

logical,parameter::inflow=.false.

real (rprec), parameter :: buff_end = 1._rprec    ! position of right end of buffer region, as a fraction of L_x
real (rprec), parameter :: buff_len = 0.125_rprec ! length of buffer region as a fraction of L_x
real (rprec), parameter :: face_avg = 1.0_rprec

! SKS
logical, parameter :: read_inflow_file = .false.


logical, parameter :: write_inflow_file = .false. !--records at position jx_s
integer, parameter :: jt_start_write = 15000

! forcing along top and bottom bdrys if inflow is true and force_top_bot is true,
! then the top & bottom velocities are forced to the inflow velocity
logical, parameter :: force_top_bot = .false.
logical, parameter :: use_mean_p_force = .false. ! .not.inflow
real(rprec),parameter::mean_p_force = 1._rprec * z_i/L_z/nproc

integer :: jt        ! global time-step counter
integer :: jt_total  ! used for cumulative time (see io module)

! time advance parameters (AB2)
real (rprec), parameter :: tadv1 = 1.5_rprec, tadv2 = 1._rprec - tadv1

!------xxxxxxxxx--SCALARS_PARAMETERS--xxxxxxxxx---------------
! SKS
! logical,parameter :: S_FLAG=.false.,spatial_flux_flag=.FALSE.,OB_homog_flag=.TRUE.,WALL_HOMOG_FLAG=.FALSE.
! logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.TRUE.,WALL_HOMOG_FLAG=.TRUE.
! MM Change to Netural S_FLAG= False !
logical,parameter :: S_FLAG=.TRUE.,spatial_flux_flag=.FALSE.,OB_homog_flag=.FALSE.,WALL_HOMOG_FLAG=.FALSE.

integer,parameter::DYN_init=1, SCAL_init=1

integer,parameter :: lbc=1, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0

logical,parameter :: remote_to_patch_flag=.FALSE.! create a 2 patch T_s field using the remote-sensed data
! The corresponding remote_to_patch subroutine is in bottombc.f90
integer,parameter :: diurnal_forcing_flag=0, no_days=1
logical,parameter :: jan_diurnal_run=.false.,ug_diurnal_variation=.false.
logical,parameter :: GABLS_diurnal_test=.false.


! Added a new parameter - passive_scalar for passive scalars with bldngs
logical,parameter :: passive_scalar=.false.,GABLS_test=.false.
logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
logical,parameter :: check_dt=.TRUE.
integer,parameter :: stencil_pts=4
logical,parameter :: coarse_grain_flag=.FALSE.
real(kind=rprec),parameter::g=9.81_rprec



! --- assign namelists ------

namelist /param_nml/ q_scale, T_scale, &
                     u_star, ug0, vg0, &
                     inv_strength, inv_strength_q, &
                     dTdz_top, dqdz_top, &
                     T_init, q_init, &
                     coriolis_forcing, nsteps, dt_dim, &
                     c_count, p_count, ubc, damping_method, &
                     initu, initsc, inilag 
                
!--------
contains
! ---------
subroutine read_namelist

 implicit none

!------
  integer::io1
  logical :: exst1
  character (*), parameter :: ftotal_time = path // 'total_time.dat'
  
  inquire (file=ftotal_time, exist=exst1)
  if (exst1) then
    open (1, file=ftotal_time)
    read (1, *) jt_total 
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', ftotal_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0 
  end if


  open (9873, file='input.nml')
  read (9873, nml=param_nml, iostat=io1)
  close(9873)
  dt=dt_dim*u_star/z_i
  coriol=1.3962634E-04*z_i/u_star
  
  if (jt_total .GT. 0) then
  
  if(coord==0) then
  write (*, *) 'jt_total is larger than 0, reading from previous vel_sc.out'
  end if
  
  initu = .true.
  initsc = .true.
  inilag=.false.
  
  end if
  
if(coord==0) then
write(*,  param_nml)
end if

end subroutine read_namelist

end module param
