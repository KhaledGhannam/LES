module scalars_module2
use types,only:rprec
use param 
use sim_param,only:u,v,w
use sgsmodule,only: Nu_t

!use main
 
implicit none

!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
!    average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile



contains

subroutine timestep_conditions(CFLu,CFLv,CFLw,visc_stab)
! This subroutine computes CFl and viscous stability and is called every wbase timesteps from main.f90
implicit none
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
! SKS
  real(kind=rprec),dimension(4,nproc)::max_var_tot_domain
! real(kind=rprec),dimension(nproc,4)::max_var_tot_domain
! SKS
$endif

real(kind=rprec) :: delta, u_res_max,v_res_max,w_res_max,uvw_res_max, nu_max
real(kind=rprec),dimension(1,4) :: max_vels
real(kind=rprec),intent(out) :: CFLu, CFLv, CFLw, visc_stab      ! u_res_max, nu_max, dt_upd
 
$if ($MPI)
  recvcounts = size (max_vels)
  displs = coord_of_rank * recvcounts 
  max_vels(1,1)=maxval(u(1:nx,1:ny,1:nz-1))
  max_vels(1,2)=maxval(v(1:nx,1:ny,1:nz-1))
  max_vels(1,3)=maxval(w(1:nx,1:ny,1:nz-1))
  max_vels(1,4)=maxval(Nu_t(1:nx,1:ny,1:nz-1))
  call mpi_gatherv (max_vels(1,1), size(max_vels), MPI_RPREC,                &
                    max_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  ! SKS
  u_res_max=sqrt(maxval(max_var_tot_domain(1,:))**2)
  v_res_max=sqrt(maxval(max_var_tot_domain(2,:))**2)
  w_res_max=sqrt(maxval(max_var_tot_domain(3,:))**2)
  nu_max=maxval(max_var_tot_domain(4,:))
$else
  u_res_max = sqrt(maxval(u(1:nx,1:ny,1:nz-1)**2+v(1:nx,1:ny,1:nz-1)**2+&
  w(1:nx,1:ny,1:nz-1)**2)) ! don't bother with interp here
  nu_max=maxval(Nu_t(1:nx,1:ny,1:nz-1))
$endif  

if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
    delta = min(dx, dy)
    CFLu = u_res_max*dt/dx
    CFLv = v_res_max*dt/dy
    CFLw = w_res_max*dt/dz
    visc_stab = dt*nu_max/(delta**2)

  !  dt_upd=0.08_rprec*delta/u_res_max
end if

end subroutine timestep_conditions

end module scalars_module2
