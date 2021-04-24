module atm_thermodynamics
use types,only:rprec
use sim_param
use param
! -----------------------
implicit none
! -----------------------------


$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension(ld,ny,$lbz:nz):: actual_T, actual_Tv, pr_atm,rel_hum, rel_hum_q,vapor_pr,sat_vapor_pr
real(kind=rprec),dimension(ld,ny,$lbz:nz):: sat_qmix, zlcl_all
real(kind=rprec),dimension(nx,ny):: zlcl_parcel,arg17
real(kind=rprec)::zlcl_parcel_ave, zlcl_parcel_ave2


! ------------
contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine thermo()
use test_filtermodule
use param
use sim_param

implicit none

integer:: nn,jx,jy,jz,k_global

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

do jz=1,nz-1
      k_global=jz+coord*(nz-1)
       
    actual_T(:,:,jz)=theta(:,:,jz) - (g/cp)*(k_global*dz*z_i)/T_scale
    actual_Tv(:,:,jz)=actual_T(:,:,jz)*(1+0.61_rprec*0.001_rprec*qmix(:,:,jz)*q_scale)
    sat_vapor_pr(:,:,jz)=611.2_rprec*exp(17.62_rprec*(actual_T(:,:,jz)*T_scale-273.2)/ &
                                (243.12+(actual_T(:,:,jz)*T_scale-273.2)))
    pr_atm(:,:,jz)=pr_surf*exp(-(k_global*dz*z_i)/(Rd*actual_Tv(:,:,jz)*T_scale/g))

    vapor_pr(:,:,jz)=(0.001_rprec*q_scale*qmix(:,:,jz)/(0.001_rprec*q_scale*qmix(:,:,jz)+0.622_rprec))* &
                           pr_atm(:,:,jz)
    sat_qmix(:,:,jz)=0.622_rprec*sat_vapor_pr(:,:,jz)/(pr_atm(:,:,jz)-sat_vapor_pr(:,:,jz))
    rel_hum(:,:,jz)=vapor_pr(:,:,jz)/sat_vapor_pr(:,:,jz)
    rel_hum_q(:,:,jz)=0.001_rprec*q_scale*qmix(:,:,jz)/sat_qmix(:,:,jz)
    zlcl_all(:,:,jz)= k_global*dz*z_i + (cp/g)*(actual_T(:,:,jz)*T_scale -55 - &
                             (1/(actual_T(:,:,jz)*T_scale-55) -dlog(rel_hum_q(:,:,jz))/2840)**(-1.0_rprec) )
end do



if ((jt_total .GE. 120000) .AND. (jt_total .LE. 144000) .AND.  modulo (jt_total, 40) == 0) then
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 1)) then

        zlcl_parcel(1:nx,1:ny)=1.5_rprec*dz*z_i + (cp/g)*(actual_T(1:nx,1:ny,1)*T_scale -55 - &
                             (1/(actual_T(1:nx,1:ny,1)*T_scale-55) - &
                             dlog(rel_hum_q(1:nx,1:ny,1))/2840)**(-1.0_rprec) )
        

      open(unit=8293,file=path//'output/zlcl_parcel.out',status="unknown",position="append")
                do jy=1,ny
                  write(8293,7187) (zlcl_parcel(jx,jy),jx=1,nx)
                end do
              close(8293)
end if
end if

if ((jt_total .GE. 1) .AND.  modulo (jt_total, 20) == 0) then
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 1)) then
        
        arg17(1:nx,1:ny)=1.5_rprec*dz*z_i + (cp/g)*(actual_T(1:nx,1:ny,1)*T_scale -55 - &
                             (1/(actual_T(1:nx,1:ny,1)*T_scale-55) - &
                             dlog(rel_hum_q(1:nx,1:nx,1))/2840)**(-1.0_rprec) )
                             
        zlcl_parcel_ave2=1.5_rprec*dz*z_i + (cp/g)*(sum(actual_T(1:nx,1:ny,1))/float(nx*ny)*T_scale -55 - &
                             (1/(sum(actual_T(1:nx,1:ny,1))/float(nx*ny)*T_scale-55) - &
                             dlog(sum(rel_hum_q(1:nx,1:nx,1))/float(nx*ny))/2840)**(-1.0_rprec) )
                             
        zlcl_parcel_ave=sum(arg17(1:nx,1:ny))/float(nx*ny)

             open (unit=8294,file=path//'output/zlcl_parcel_ave.out',status="unknown",position="append")
    write(8294,7187) zlcl_parcel_ave, zlcl_parcel_ave2
    close(8294) 
end if
end if

7187     format(1400(E14.5))

$if ($MPI)
    !--exchange ghost-node information
    !--send stuff up to ghost layer (0) (for z-derivs)
    !--nz should already be in sync with 1 level: done in project()
    call mpi_sendrecv (actual_T(1, 1, nz-1), ld*ny, MPI_RPREC, up, 727,  &
                       actual_T(1, 1, 0), ld*ny, MPI_RPREC, down, 727,   &
                       comm, status, ierr)

    call mpi_sendrecv (actual_Tv(1, 1, nz-1), ld*ny, MPI_RPREC, up, 728,  &
                       actual_Tv(1, 1, 0), ld*ny, MPI_RPREC, down, 728,   &
                       comm, status, ierr)

    call mpi_sendrecv (sat_vapor_pr(1, 1, nz-1), ld*ny, MPI_RPREC, up, 729,  &
                       sat_vapor_pr(1, 1, 0), ld*ny, MPI_RPREC, down, 729,   &
                       comm, status, ierr)
    call mpi_sendrecv (pr_atm(1, 1, nz-1), ld*ny, MPI_RPREC, up, 730,  &
                       pr_atm(1, 1, 0), ld*ny, MPI_RPREC, down, 730,   &
                       comm, status, ierr)
    call mpi_sendrecv (vapor_pr(1, 1, nz-1), ld*ny, MPI_RPREC, up, 731,  &
                       vapor_pr(1, 1, 0), ld*ny, MPI_RPREC, down, 731,   &
                       comm, status, ierr)
    call mpi_sendrecv (sat_qmix(1, 1, nz-1), ld*ny, MPI_RPREC, up, 732,  &
                       sat_qmix(1, 1, 0), ld*ny, MPI_RPREC, down, 732,   &
                       comm, status, ierr)

    call mpi_sendrecv (rel_hum(1, 1, nz-1), ld*ny, MPI_RPREC, up, 733,  &
                       rel_hum(1, 1, 0), ld*ny, MPI_RPREC, down, 733,   &
                       comm, status, ierr)
    call mpi_sendrecv (rel_hum_q(1, 1, nz-1), ld*ny, MPI_RPREC, up, 734,  &
                       rel_hum_q(1, 1, 0), ld*ny, MPI_RPREC, down, 734,   &
                       comm, status, ierr)
    call mpi_sendrecv (zlcl_all(1, 1, nz-1), ld*ny, MPI_RPREC, up, 735,  &
                       zlcl_all(1, 1, 0), ld*ny, MPI_RPREC, down, 735,   &
                       comm, status, ierr)

  $endif

  end subroutine thermo
end module atm_thermodynamics
