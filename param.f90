module param
use types,only:rprec
implicit none
$if ($MPI)
  ! use mpi
  include "mpif.h"
$endif
! implicit none

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

logical,parameter::VERBOSE = .false.  !--prints small stuff to screen
                   !--use DEBUG to write lots of data/files

integer,parameter::nx=360,ny=288,nz=(385-1)/nproc + 1


integer,parameter::nz_tot=(nz-1)*nproc + 1
integer,parameter::nx2=3*nx/2,ny2=3*ny/2
integer,parameter::lh=nx/2+1,ld=2*lh,lh_big=nx2/2+1,ld_big=2*lh_big

  ! -- MMC : Bogus is used for debugging .
integer, parameter :: iBOGUS = -1234567890  !--NOT a new Apple product
real (rprec), parameter :: BOGUS = -1234567890._rprec

real(rprec),parameter::pi=3.1415926535897932384626433_rprec

real(kind=rprec),parameter::L_x=5.0_rprec*pi,L_y=4.0_rprec*pi
!real(rprec),parameter::z_i=1250._rprec, L_z=1500._rprec/nproc
real(rprec),parameter::z_i=1000._rprec, L_z=3000._rprec/nproc
! L_z is not nondimensionalized by z_i yet
! MM : For netural case => L_z=z_i

! set the aspect ratio of the box, already nondimensional
real(rprec),parameter::dz=L_z/z_i/(nz-1)
real(rprec),parameter::dx=L_x/nx,dy=L_y/ny

integer,parameter::nsteps=1 ! simulation steps + SCAL_init for 32^3 with dt_dim=0.25
integer,parameter::spectraCALC=0 ! 125000

!--Coriolis stuff; coriol=non-dim coriolis parameter,
! ug=horiz geostrophic vel, vg=transverse geostrophic vel
logical,parameter::coriolis_forcing=.true.

!--fpx3 way to get around sqrt issue below:
$define $ug_dim 1.0
$define $vg_dim 0.0
$define $_rprec _rprec
real (rprec),protected :: ug_dim = $ug_dim$_rprec
real (rprec),protected :: vg_dim = $vg_dim$_rprec

$perl: $u_star = sqrt( $ug_dim**2 + $vg_dim**2 )
real(rprec),protected::u_star=$u_star$_rprec
real(rprec),parameter::Pr=.4_rprec

!MM : I changed here dt_dim=0.04
real(rprec),parameter::dt_dim=0.05_rprec ! this is for unstable wt_s = 0.1, 64 cube case
real(rprec),protected:: dt !=dt_dim*u_star/z_i
! real(rprec),parameter::dt=dt_dim*u_star/z_i

real(rprec),protected::coriol !=1.3962634E-04*z_i/u_star !for latitude = 37.6N ? or 73N ?
! MM : I defined ug_global in Main.f90 variable in time :
 !real(rprec)::ug=ug_dim/u_star,vg=vg_dim/u_star
real(rprec),dimension(ld,ny,1:((nz-1)*nproc))::ug,vg


real(rprec),parameter::cp=1005.0_rprec       ! specific heat capacity at constant pressure (J/kg/K)
real(rprec),parameter::Le=2.5008E06          ! Latent heat of vaporization (J/kg)
real(rprec),parameter::Rv=461.5_rprec        ! Gas constant for water vapor (J/kg/K)
real(rprec),parameter::Rd=287.04_rprec        ! Gas constant for dry air (J/kg/K)
real(rprec),parameter::rho_d=1.2923_rprec      ! density of dry air (kg/m^3)
real(rprec),parameter::rho_w=1000.0_rprec     ! water density (kg/m^3)
real(rprec),parameter::pr_surf=100000.0_rprec    ! Surface pressure (Pa)
! MM
real(rprec),parameter::vonk=.4_rprec

! SKS p_count = output ones ; c_count = sampling
integer,parameter::c_count=20,p_count=10000 !p_count=800 => 5 mins for 64^3;dt=0.003
! integer,parameter::c_count=10,p_count=3600 !p_count=800 => 5 mins for 64^3;dt=0.003
! SKS
integer, parameter :: cs_count = 5  !--tsteps between dynamic Cs updates
logical,parameter::output=.true.
logical, parameter :: use_avgslice = .true.
! The parameter average_dim_num controls the number of dimensions averaged
! while outputting files from avgslice ..
! Possible values of average_dim_num: 2 => average over x,y,t ; 1 => average over y,t
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

! ubc: upper boundary condition: ubc=0 stress free lid, ubc=1 sponge
! damping method = 1: use Nieuwstadt's method, = 2: use Raleigh damping method
integer,parameter::ubc=1,damping_method=2

character (*), parameter :: lbc_mom = 'wall' !--'wall', 'stress free'

! prescribed inflow: constant or read from file
logical,parameter::inflow=.false.

real (rprec), parameter :: buff_end = 1._rprec    ! position of right end of buffer region, as a fraction of L_x
real (rprec), parameter :: buff_len = 0.125_rprec ! length of buffer region as a fraction of L_x
real (rprec), parameter :: face_avg = 1.0_rprec

! SKS
logical, parameter :: read_inflow_file = .false.
! logical, parameter :: read_inflow_file = .true.
! SKS

logical, parameter :: write_inflow_file = .false. !--records at position jx_s
integer, parameter :: jt_start_write = 15000

! forcing along top and bottom bdrys if inflow is true and force_top_bot is true,
! then the top & bottom velocities are forced to the inflow velocity
logical, parameter :: force_top_bot = .false.

logical, parameter :: use_mean_p_force = .false. ! .not.inflow
real (rprec), parameter :: mean_p_force = 1._rprec * z_i/L_z/nproc  !--usually just z_i/L_z

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
! SKS

integer,parameter::DYN_init=1, SCAL_init=1

! lbc=0: prescribed surface temperature, lbc=2 prescribed surface flux
! SKS
integer,parameter :: lbc=1, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0
! integer,parameter :: lbc=0, patch_flag=1, remote_flag=0,remote_homog_flag=0,remote_flux_homog_flag=0
! SKS
logical,parameter :: remote_to_patch_flag=.FALSE.! create a 2 patch T_s field using the remote-sensed data
! The corresponding remote_to_patch subroutine is in bottombc.f90
integer,parameter :: diurnal_forcing_flag=0, no_days=1
logical,parameter :: jan_diurnal_run=.false.,ug_diurnal_variation=.false.
logical,parameter :: GABLS_diurnal_test=.false.

   ! -- MMC : using previous simulations employing vel.out
logical,parameter :: initu=.false.,initsc=.false.,inilag=.true.,interp=.FALSE.
!logical,parameter :: initu=.true.,initsc=.true.,inilag=.true.,interp=.FALSE.
! initu   = true to read from a file; false to create with random noise
! initlag = true to initialize cs, FLM & FMM; false to read from vel.out

! Added a new parameter - passive_scalar for passive scalars with bldngs
logical,parameter :: passive_scalar=.false.,GABLS_test=.false.
logical,parameter :: test_phase=.FALSE., vec_map=.FALSE., smag_sc=.FALSE.
logical,parameter :: check_dt=.TRUE.
integer,parameter :: stencil_pts=4
logical,parameter :: coarse_grain_flag=.FALSE.
real(kind=rprec),parameter::g=9.81_rprec, inv_strength=0.0050_rprec,dTdz_top=0.0050_rprec
real(kind=rprec),parameter :: T_scale=400._rprec,T_init=300.00_rprec ! unstable
real(kind=rprec),parameter::inv_strength_q=-0.004_rprec,dqdz_top=-0.004_rprec
real(kind=rprec),parameter :: q_init=12.00_rprec 
real(kind=rprec) :: q_scale=100._rprec

namelist /scales_nml/ q_scale,ug_dim,vg_dim

!--------
contains
! ---------
  subroutine read_namelist
       integer::io
   open (983, file='input.nml')
  read (983, nml=scales_nml, iostat=io)
  close(983)
!  write(984,  nml=scales_nml)
!  write(*,  nml=scales_nml)

u_star = sqrt( ug_dim**2 + vg_dim**2 )
dt=dt_dim*u_star/z_i
coriol=1.3962634E-04*z_i/u_star
if(coord==0) then
write(*,  scales_nml)
write(*,*) u_star,dt,coriol
end if

end subroutine read_namelist

end module param
