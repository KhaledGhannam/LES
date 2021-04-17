subroutine initial()
use types,only:rprec
use param
use sim_param,only: u,v,w,RHSx,RHSy,RHSz,theta,qmix
use sgsmodule,only: Cs_opt2,Cs_opt2_avg,F_LM,F_MM,F_QN,F_NN, &
                   G_LM,G_MM,G_QN,G_NN,Pr_t

use scalars_module,only: RHS_T,sgs_t3
use bottombc,only: psi_m,psi_m0,psi_h,psi_h0
use scalars_module_q,only: RHS_q,sgs_q3
use immersedbc,only:fx,fy,fz,u_des,v_des,w_des,n_bldg,bldg_pts

!----------------------------------------------
implicit none

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

character (200) :: fname, temp
real(kind=rprec),dimension(ld,ny,nz)::crap
real(kind=rprec)::z
integer::i,jz

Cs_opt2_avg=0._rprec
fx=0._rprec;fy=0._rprec;fz=0._rprec
u_des=0._rprec;v_des=0._rprec;w_des=0._rprec


if (S_FLAG) then
    fname = path // 'vel_sc_out/vel_sc.out'
else
    fname = path // 'vel_sc_out/vel.out'
end if

$if ($MPI)
  write (temp, '(".c",i0)') coord
  fname = trim (fname) // temp
$endif

open(11,file=fname,form='unformatted')

if(initu)then         ! if_initu_begins
  if(initsc) then     ! if_initsc_begins
  
    print *,'Reading initial velocity and temperature from file'
    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),theta(:,:,1:nz), &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz),RHSz(:, :, 1:nz),         &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m
      case (2:3)
        read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),   &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),        &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
                 
      case (4)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2 
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),  &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),       &
                    RHS_T(:,:,1:nz),sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM,  &
                    crap, crap
        end if
        
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    qmix(:,:,1:nz),RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), RHS_q(:,:,1:nz),sgs_t3(:,:,1), sgs_q3(:,:,1),psi_m,Cs_opt2
        else
       
         read (11)  u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),theta(:,:,1:nz),   &
                    qmix(:,:,1:nz),RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz),          &
                    RHS_T(:,:,1:nz),RHS_q(:,:,1:nz),sgs_t3(:,:,1),sgs_q3(:,:,1),&
                    psi_m,psi_m0,psi_h,psi_h0,Cs_opt2,&
                    F_LM,F_MM,F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t(:,:,1:nz)
        end if
      ! SKS
      case (6)
        read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz),   &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),        &
                  RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
      case (7)
        if (inilag) then  ! Not sure what all are needed when inilag = .true.
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2
        else
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),theta(:,:,1:nz), &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),      &
                    RHS_T(:,:,1:nz), sgs_t3(:,:,1), psi_m, Cs_opt2, F_LM, F_MM,&
                    F_QN, &
                    F_NN, &
                    G_LM, &
                    G_MM, &
                    G_QN, &
                    G_NN, &
                    Pr_t(:, :, 1:nz)
        end if
      ! SKS
      case default
        write (*,*) 'initial: invalid model number'
    end select
    
    do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
! SKS
     write(6,7789) coord,jz,z,(sum(u(:,:,jz))/float(nx*ny)),(sum(v(:,:,jz))/&
     float(nx*ny)),(sum(w(:,:,jz))/float(nx*ny)),&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,(sum(qmix(:,:,jz))/float(nx*ny))*q_scale

end do

7789 format('coord, jz, z, ubar, vbar, wbar,T_bar,q_bar:',(1x,I3,1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))


  else     ! if_init_sc continues (if false this follows read only velocity field)

    print *,'Reading initial velocity field from file'

    select case (model)
      case (1)
        read (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),       &
                  RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz)
      case (2:4)
        read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2
      case (5)
        if (inilag) then  !--not sure if Cs_opt2 should be there, just quickie
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2
        else
          read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2, F_LM, F_MM, F_QN, F_NN
        end if
      ! SKS
      case (6)
        read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                 RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                 Cs_opt2
      case (7)
        if (inilag) then  ! not sure what inilag does
          read (11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                    Cs_opt2
        else
          read(11) u(:, :, 1:nz),v(:, :, 1:nz),w(:, :, 1:nz),             &
                   RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),  &
                   Cs_opt2, F_LM, F_MM, F_QN, F_NN
        end if
      ! SKS
      case default
        write (*, *) 'initial: invalid model number'
    end select


   end if   ! if init_sc condition ends here (the else below is for initu being false)
  
else
  if (dns_bc) then
     print*, 'Creating initial velocity field with DNS BCs'
     call ic_dns()
  else
    print*, 'Creating initial fields'
    if (S_FLAG) then
         print *, 'Creating initial velocity & scalar fields'
         call ic_scal()
    else
       ! SKS
       ! If S_FLAG == 0 then ic() is called, else the above ones are called
       print*, 'Creating initial velocity field'
       call ic()
       ! SKS
    end if
  end if
end if



! bldg stuff
if (use_bldg) then
   open(1,file=path//'bldg.dat')
   read(1,*) n_bldg
   allocate(bldg_pts(5,n_bldg))
   do i=1,n_bldg
      read(1,*) bldg_pts(1:5,i)
      if(bldg_pts(5,i).ge.nz)then
         write(*,*)"lz should be less than nz"
         stop
      end if
   end do
   close(1)
end if

$if ($MPI)
  !--synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz 
  call mpi_sendrecv (u(1, 1, nz-1), ld*ny, MPI_RPREC, up, 1,  &
                     u(1, 1, 0), ld*ny, MPI_RPREC, down, 1,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, nz-1), ld*ny, MPI_RPREC, up, 2,  &
                     v(1, 1, 0), ld*ny, MPI_RPREC, down, 2,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, nz-1), ld*ny, MPI_RPREC, up, 3,  &
                     w(1, 1, 0), ld*ny, MPI_RPREC, down, 3,   &
                     comm, status, ierr)
  call mpi_sendrecv (u(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
                     u(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
                     comm, status, ierr)
  call mpi_sendrecv (v(1, 1, 1), ld*ny, MPI_RPREC, down, 5,  &
                     v(1, 1, nz), ld*ny, MPI_RPREC, up, 5,   &
                     comm, status, ierr)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, 6,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, 6,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, 7,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, 7,   &
                     comm, status, ierr)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, 8,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, 8,   &
                     comm, status, ierr)

  call mpi_sendrecv (qmix(1, 1, nz-1), ld*ny, MPI_RPREC, up, 9,  &
                     qmix(1, 1, 0), ld*ny, MPI_RPREC, down, 9,   &
                     comm, status, ierr)

  call mpi_sendrecv (qmix(1, 1, 1), ld*ny, MPI_RPREC, down, 10,  &
                     qmix(1, 1, nz), ld*ny, MPI_RPREC, up, 10,   &
                     comm, status, ierr)
$endif

if (USE_MPI .and. coord == 0) then
  !--set 0-level velocities to BOGUS
  u(:, :, $lbz) = BOGUS
  v(:, :, $lbz) = BOGUS
  w(:, :, $lbz) = BOGUS
  theta(:, :, $lbz) = BOGUS
  qmix(:, :, $lbz) = BOGUS
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ic_scal()

implicit none
real(kind=rprec),dimension(nz)::ubar,vbar,wbar
real(kind=rprec)::rms, noise, arg, arg2,theta_mean,q_mean
real(kind=rprec)::z,w_star,T_star,q_star,ran3
real(kind=rprec)::z_turb_limit,perturb_height_factor,z_inv,gam
integer::jx,jy,jz,seed,jz_abs
! SKS
real(kind=rprec),dimension(nz_tot-1)::utotinit,vtotinit,wtotinit,Ttotinit,qtotinit
real(kind=rprec),dimension(nz-1)::uinit,vinit,winit,Tinit,qinit
real(kind=rprec)::k_global

$if ($MPI)
  integer::recvcounts(nproc)
  integer::displs(nproc)
$endif
! SKS

       if (lbc .eq. 0) then
           theta_mean=T_init !T_s_min is dimensional while T_s is non-dimensional
           q_mean=q_init
           print *,'theta_mean = , q_mean= ',theta_mean, q_mean
           
       else
           theta_mean=T_init
           q_mean=q_init
         print *,'theta_mean = , q_mean= ',theta_mean, q_mean
       end if

       if (sum(wt_s(1:nx,1:ny))/float(nx*ny) .lt. 0._rprec) then
        
         perturb_height_factor=1.0_rprec
         z_inv=1.0_rprec*z_i
       else
        
         perturb_height_factor=1.0_rprec
         z_inv=1.0_rprec*z_i
    
       end if

       z_turb_limit=perturb_height_factor*z_i
      
      if (sum(wt_s(1:nx,1:ny))/float(nx*ny) .eq. 0.0_rprec) then

      w_star=(g/theta_mean*0.06_rprec*z_i)**(1._rprec/3._rprec)
      T_star=0.06_rprec/w_star
      q_star=T_star
      
      else
      w_star=sign((g/theta_mean*abs(sum(wt_s(1:nx,1:ny))/float(nx*ny))*z_i)** &
                   (1._rprec/3._rprec),sum(wt_s(1:nx,1:ny))/float(nx*ny))
      T_star=(sum(wt_s(1:nx,1:ny))/float(nx*ny))/w_star
      q_star=T_star
      end if


do jz=1,nz
         $if ($MPI)
            z = (coord*(nz-1) + jz - 0.5_rprec) * dz
         $else
            z = (real(jz)-0.5_rprec)*dz
         $endif
!c IC in equilibrium with rough surface (rough dominates in effective zo)
        arg2=z/(sum(zo(:,:))/float(nx*ny))
        arg=(1._rprec/vonk)*log(arg2)!-1./(2.*vonk*z_i*z_i)*z*z

    if (coriolis_forcing) then

         k_global = jz + coord*(nz-1)
       if (k_global > ((z_i/(L_z*nproc))*nproc*(nz-1))) then
                ubar(jz)=ug0/u_star
                vbar(jz)=vg0/u_star
       else
        gam=sqrt(u_star*coriol/(10*z_i))
        ubar(jz)=ug0*(1-exp(-gam*z*z_i)*cos(gam*z*z_i))/u_star
        vbar(jz)=ug0*exp(-gam*z*z_i)*sin(gam*z*z_i)/u_star
      
       end if
       
        wbar(jz)=0._rprec

    else
        ubar(jz)=arg
        vbar(jz)=0._rprec
        wbar(jz)=0._rprec
    end if

        if (z .gt. z_turb_limit) then
           ubar(jz)=ubar(jz-1)
           vbar(jz)=vbar(jz-1)
        end if
end do

  rms = 3._rprec
  do jz=1,nz
  $if ($MPI)
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
  $else
    jz_abs = jz
    z = (jz- 0.5_rprec) * dz * z_i
  $endif
  seed = -80 - jz_abs  !--trying to make consistent init for MPI
    do jy=1,ny
      do jx=1,nx
!c...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
!c...Taking std dev of vel as 1 at all heights

       if (z .LE. z_turb_limit) then
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         u(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+ubar(jz)
         
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         v(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+vbar(jz) 
         
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         w(jx,jy,jz)=noise*(1._rprec-z/z_i)*w_star/u_star+wbar(jz)
         
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         theta(jx,jy,jz)=(theta_mean+10._rprec*noise*(1-z/z_i)*T_star)/T_scale
         
         noise=rms/0.289_rprec*(ran3(seed)-0.5_rprec)
         qmix(jx,jy,jz)=(q_mean+10._rprec*noise*(1-z/z_i)*q_star)/q_scale 
         
      else
         u(jx,jy,jz)=ubar(jz)
         v(jx,jy,jz)=vbar(jz)
         w(jx,jy,jz)=wbar(jz)
 
         theta(jx,jy,jz)=(theta_mean+(z-z_inv)*inv_strength)/T_scale
         qmix(jx,jy,jz)=(q_mean+(z-z_inv)*inv_strength_q)/q_scale
             
      end if
        
       end do
     end do
  end do

  !...BC for W
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    w(1:nx, 1:ny, 1) = 0._rprec
  end if
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    w(1:nx, 1:ny, nz) = 0._rprec
  endif

  !...BC for U, V & T
  if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
    u(1:nx, 1:ny, nz) = u(1:nx, 1:ny, nz-1)
    v(1:nx, 1:ny, nz) = v(1:nx, 1:ny, nz-1)
    theta(1:nx, 1:ny, nz) = theta(1:nx, 1:ny, nz-1)+dTdz_top/T_scale*z_i*dz
    qmix(1:nx, 1:ny, nz) = qmix(1:nx, 1:ny, nz-1)+dqdz_top/q_scale*z_i*dz
   end if

!VK Display the mean vertical profiles of the initialized variables on the screen
!open(unit=44,file=path//'output/init_profiles.dat',status="unknown",position="append")

do jz=1,nz
     $if ($MPI)
       z = (coord*(nz-1) + jz - 0.5_rprec) * dz
     $else
       z = (jz - 0.5_rprec) * dz
     $endif
! SKS
     write(6,7781) coord,jz,z,(sum(u(:,:,jz))/float(nx*ny)),(sum(v(:,:,jz))/&
     float(nx*ny)),(sum(w(:,:,jz))/float(nx*ny)),&
     (sum(theta(:,:,jz))/float(nx*ny))*T_scale,(sum(qmix(:,:,jz))/float(nx*ny))*q_scale

end do

7781 format('coord, jz, z, ubar, vbar, wbar,T_bar,q_bar:',(1x,I3,1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))

! SKS
! This part written just to get the initial profiles 

open(unit=103,file=path//'output/initial_profiles.dat',status="unknown",position="append")

do jz=1,nz-1      
   uinit(jz) = (sum(u(:,:,jz))/float(nx*ny))*u_star
   vinit(jz) = (sum(v(:,:,jz))/float(nx*ny))*u_star
   winit(jz) = (sum(w(:,:,jz))/float(nx*ny))*u_star
   Tinit(jz) = (sum(theta(:,:,jz))/float(nx*ny))*T_scale
   qinit(jz) = (sum(qmix(:,:,jz))/float(nx*ny))*q_scale
end do
      
$if ($MPI) 
  recvcounts = size(uinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (uinit(1),size(uinit),MPI_RPREC, &
                    utotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(vinit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (vinit(1),size(vinit),MPI_RPREC, &
                    vtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(winit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (winit(1),size(winit),MPI_RPREC, &
                    wtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
  recvcounts = size(Tinit)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (Tinit(1),size(Tinit),MPI_RPREC, &
                    Ttotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)

  recvcounts = size(qinit)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (qinit(1),size(qinit),MPI_RPREC, &
                    qtotinit(1),recvcounts,displs,  &
                    MPI_RPREC,rank_of_coord(0),comm,ierr)
$else
  utotinit = uinit
  vtotinit = vinit
  wtotinit = winit
  Ttotinit = Tinit
  qtotinit = qinit
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  do jz=1,nz_tot-1
     write(103,8001) utotinit(jz),vtotinit(jz),wtotinit(jz),Ttotinit(jz),qtotinit(jz)
  end do
  write(103,*)
end if
8001  format(1400(E14.5))
! SKS    

end subroutine ic_scal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine initial
