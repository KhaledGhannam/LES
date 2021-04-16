module scalars_module
! HUMIDITY subroutines in place but not yet turned on !!
use types,only:rprec
use param 
use sim_param,only:u,v,w,theta,path,u_old,qmix
use bottombc ! Includes patches subroutine
use sgsmodule,only:Nu_t,Pr_t
implicit none
integer, parameter:: tag_counter = 273
logical, parameter:: SCALAR_DEBUG=.FALSE.
!!!!!!--------------------------------------------
! Part I of the scalar files - contains the basic subroutines
! Also look at scalars_module2.f90 for other subroutines !! 
! CONTAINS subroutines :
! theta_all_in_one - Performs the various steps for theta calculation
! humidity_all_in_one - Performs the various steps for humidity
! scalar_RHS_calc - computes the RHS side of the scalar evolution equation
! calcbeta - computes the buoyancy term for temperature
! step_scalar - time evolves the scalar
! obukhov - computes the obukhov similarity terms for use in scalars,wallstress and derivwall
! Authored by Vijayant Kumar
! Last modified - April 24, 2004
!!!!!!--------------------------------------------
 
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
!! KMG: added sgs_t1 and sgs_t2 for horizontal sgs heat fluxes Pi_1 and Pi_2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: beta_scal,Pr_
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dTdz,dTdx,dTdy
! Only ones needed for output..Might need to add x and y derivatives here in case they need to be outputted
! Right now they are in the "scalar"_all_in_one routines below !!
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS_Tf,RHS_T
real(kind=rprec), dimension(ld,ny,$lbz:nz):: sgs_t3,sgs_t1,sgs_t2 ! defines the surface sgs flux

real(kind=rprec),dimension(nx,ny)::L,wstar ! defines obukhov length and convective vel scale, w_star
real(kind=rprec),dimension(ld,ny)::T_s_filtered ! filtered T_s for calc of wT_s

! Now define local u_star in bottombc.f90 and calculate in wallstress.f90 and use that value everywhere else
integer, parameter:: obukhov_output=0 !Controls whether obukhov variables are outputted by scalar_slice
integer, parameter:: wt_s_vector_dim1=no_days*86400/300+1

real(kind=rprec),dimension(wt_s_vector_dim1,1) :: wt_s_vector
! Variables for heterogeneity analysis
! hetero_array_freqz = number of time steps equivalent to 20 seconds
integer,parameter:: hetero_array_freqz=int(20/dt_dim),hetero_count_out=p_count
integer,save::time_ind

!--MM 2015: Added Artificially Moving(Heating) Inversion:
integer, parameter :: nz_global = (nz-1) * nproc
integer ::k_global
real ::dt_avg
real(kind=rprec),dimension(ld,ny,$lbz:nz)::theta_old

contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine theta_all_in_one
use topbc,only:sponge
use test_filtermodule
implicit none
real(kind=rprec),dimension(nx,ny)::wt_s_current
real::dummy_t
integer::jz,ios,counter
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dTdx,dTdy
real(kind=rprec),dimension($lbz:nz)::sponge_theta

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        wt_s_current(1:nx,1:ny)=wt_s(1:nx,1:ny)
 end if

! SKS
  if (model == 1 .or. model == 2 .or. model == 3 .or. model == 4 .or. model == 5) then
     Pr_ = Pr
  else
     do jz = $lbz,nz ! ubc_jz
        Pr_(:,:,jz) = Pr_t(:,:,jz) ! Pr
     end do
  end if
! SKS

! Right now set the Prandtl num matrix equal to a constant Prandtl
! number as specified in param. could use the dynamic model ideal to compute Pr as well !!
! The plan-averaged dynamic Prandtl number model is already coded. just need to put it in !!

call filt_da(theta,dTdx,dTdy)
call ddz_uv (dTdz,theta)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   dTdz(:,:,Nz)=dTdz_top/T_scale*z_i ! Valid for temperature
end if

$if ($MPI)
! print *,'synchronizing in theta_all_in for coord = ',coord
! Need to synchronize w and dTdz across the processors for uvp node
! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+1,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+1,   &
                     comm, status, ierr)   
  call mpi_sendrecv (dTdz(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+2,  &
                     dTdz(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+2,   &
                     comm, status, ierr)   

! Also need to synchronize Nu_t across the processors for 
! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+3,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+3,   &
                     comm, status, ierr)   
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+4,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+4,   &
                     comm, status, ierr)
$endif

if (S_FLAG) then
   RHS_Tf=RHS_T


! Perform test filtering of T_s for calculation of surf fluxes
     if ((jt_total .eq. SCAL_init) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
       print *,'T_s b4 filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
       T_s_filtered(1:nx,1:ny)=T_s
       call test_filter(T_s_filtered,G_test)
       T_s=T_s_filtered(1:nx,1:ny)
       print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
     end if
     if ((jt_total .eq. nsteps-1) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
        print *,'T_s after filtering',sqrt(sum((T_s-sum(T_s)/float(nx*ny))**2))/float(nx*ny)
     end if

! call scalar_RHS_calc(T_s,z_os,RHS_T,sgs_t3,jt,psi_h,phi_h,Pr_,wt_s_current)
theta_old(:,:,:)=theta(:,:,:)

call scalar_RHS_calc(theta,dTdx,dTdy,dTdz,T_s,z_os,RHS_T,sgs_t1,sgs_t2,sgs_t3,wt_s_current)

! subroutine looks like this later
! subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)


if (ubc==1 .and. damping_method==2) then ! add the damping term to the scalar equation
   do jz=1,nz-1
   RHS_T(1:nx,1:ny,jz)=RHS_T(1:nx,1:ny,jz)-0.5_rprec*(sponge(jz)+sponge(jz+1))*&
                 (theta(1:nx,1:ny,jz)-sum(theta(1:nx,1:ny,jz))/(nx*ny))
   end do
end if

call calcbeta(theta,qmix) ! Calculates the buoyancy term which gets added to the vertical momentum equation
! call calcbeta(theta,beta_scal)

if (jt .eq. SCAL_INIT .and. (.not. initsc)) then
   RHS_Tf=RHS_T
end if

call step_scalar(theta,RHS_T,RHS_Tf)
end if



$if ($MPI)
  call mpi_sendrecv (theta(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+7,  &
                     theta(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+7,   &
                     comm, status, ierr)   
  call mpi_sendrecv (theta(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+8,  &
                     theta(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+8,   &
                     comm, status, ierr)
$endif

end subroutine theta_all_in_one

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine calcbeta (scalar,scalarv)
! subroutine calcbeta (scalar, beta_scal)
! This calculates the buoyancy term (beta_scal) to be added to the vertical momentum equation for temperature
! Authored by Vijayant Kumar
! Last updated April 14, 2004
implicit none
integer::i, j, k, jz_min
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz),intent(in)::scalar,scalarv
real(kind=rprec),dimension(0:nz)::scalar_bar
real(kind=rprec)::g_hat,above, below
!..Non-dimensionalize gravity
g_hat=g*(z_i/(u_star**2))
beta_scal=0._rprec



! Note Beta is stored on W nodes, but Theta is on UVP nodes
! We do not time-advance the ground nodes, so start at k=2
! VK: Inserted the averaging code inside this file itself rather than doing it in prof
do k=$lbz,nz
scalar_bar(k)=0.0    
   do j=1,ny
      do i=1,nx
        scalar_bar(k)=scalar_bar(k)+scalar(i,j,k)*(1+0.61_rprec*0.001_rprec*scalarv(i,j,k)*q_scale)
      end do
   end do
scalar_bar(k)=scalar_bar(k)/(nx*ny)
end do
! We do not time-advance the ground nodes, so start at k=2
! For the MPI case, this means that we start from jz=2 for
! coord=0 and jz=1 otherwise... enable by an if statment

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
 else
    jz_min = 1
 end if

do k=jz_min,Nz-1
      do j=1,Ny
             do i=1,nx
                above=(scalar(i,j,k)*(1+0.61_rprec*0.001_rprec*scalarv(i,j,k)*q_scale)-scalar_bar(k))/scalar_bar(k)
                below=(scalar(i,j,k-1)*(1+0.61_rprec*0.001_rprec*scalarv(i,j,k-1)*q_scale)-scalar_bar(k-1))/scalar_bar(k-1)
                beta_scal(i,j,k)=g_hat*(above + below)/2._rprec
             end do
      end do
end do
return
end subroutine calcbeta
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------

subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)
use test_filtermodule
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
integer::i,j,k,jz
integer::jz_min,ubc_jz
real::crap2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS,temp
real(kind=rprec),dimension(ld,ny,$lbz:nz):: scalar
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,sgs_vert,sgs_horx,sgs_hory
real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
 real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: RHS_m
real(kind=rprec),dimension(nx,ny):: ustar_local,S_Surf,surf_flux,z_os,wt_s2,surf_flux_current
real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux
real(kind=rprec),dimension (ptypes):: ustar
character (64) :: fname_hetero

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

         ustar_local(:,:)=ustar_avg(:,:)
          scalar_node_1=scalar(:,:,1)
          call test_filter(scalar_node_1,G_test) 
            if (lbc==1 .AND. scalar(1,1,1)<2) then
                     
                    wt_s2(1:nx,1:ny)=wt_s(1:nx,1:ny)
                    do i=1,Nx
                    do j=1,Ny
                    surf_flux(i,j)=wt_s2(i,j)/(T_scale*u_star)
                    end do
                    end do
  
            else if (lbc==0.and.scalar(1,1,1)<2) then
                    do i=1,Nx
                    do j=1,Ny 
                    surf_flux(i,j)=(S_Surf(i,j)-scalar_node_1(i,j))*vonk*ustar_local(i,j)&
                                   /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j) + psi_h0(i,j))
                   end do 
                   end do
             end if

!c....Now we have the lowest dsdz on the UVP nodes all others on w nodes
            do i=1,Nx
            do j=1,Ny
            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*dz/2._rprec)

            end do
            end do
end if !end for if (USE_MPI ...) block

 call dealias1(u,u_m)
 call dealias1(v,v_m)
 call dealias1(w,w_m)
 call dealias1(dsdx,dsdx_m)
 call dealias1(dsdy,dsdy_m)
 call dealias1(dsdz,dsdz_m)

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    jz_min = 2
 else
    jz_min = 1
 end if
 
    ubc_jz = nz-1

! Now compute the RHS term of the filtered scalar equation. 
! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

! xxxxx ------Comments valid for MPI case only ---------XXXX
! For MPI case, all the nodes have fluid nodes (1:nz-1) except for
! near-the wall processes (2:nz-1 for w nodes and 1:nz-1 for uvp nodes)
! and the top nodes (1:nz)
! The following loop executes from 2:nz-1 for the near-wall process
! and 1:nz-1 for the other processes. Also, the subsequent 2 loops
! take care of the first node for the near-wall process (coord = 0)
! and the topmost node for the top process (coord = nproc-1).
! Note also that we need to make ghost data available for dTdz and
! w for the topmost node (jz=n) within each process and therefore,
! this synchronization (MPI_SENDRECV()) has been done in the subroutine
! theta_all_in_one ()
! xxxxx --------- MPI Comment block ends ------------------XXXX

do k=jz_min,nz-1
  do j=1,Ny2
    do i=1,Nx2
    RHS_m(i,j,k)=u_m(i,j,k)*dsdx_m(i,j,k)+v_m(i,j,k)*dsdy_m(i,j,k)&
    +(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))/2._rprec
    end do
  end do
end do


 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny2
     do i=1,Nx2
      RHS_m(i,j,1)=u_m(i,j,1)*dsdx_m(i,j,1)+v_m(i,j,1)*dsdy_m(i,j,1)&
      +(0.5_rprec*w_m(i,j,2))*dsdz_m(i,j,2)
     end do
   end do
 end if
 
 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   do j=1,Ny2
     do i=1,Nx2
      RHS_m(i,j,Nz)=u_m(i,j,Nz)*dsdx_m(i,j,Nz)+v_m(i,j,Nz)*dsdy_m(i,j,Nz)
     end do
   end do
 end if


 call dealias2(RHS,RHS_m)

!MM -----------

!c...Now building the SGS part of the RHS.
! Here the sgs_term for scalars is built up using Nu_t from sgs_stag_W.f
! and dividing it by the turbulent Prandtl # specified in dimen.h
!c....Note: Since we bring the Convective term to RHS its sign changes.
!c....Below "Temp" is used for SGS flux; its divergence is added to RHS
!VK.. Nu_t is on w nodes everywhere except at z=dz/2.. while
!VK dsdx is on uvp nodes.. so, interpolate Nu_t as we want temp to
!VK be on uvp nodes
! All this valid only till Pr_ is a constant..
! This will need a makeover once Pr_ becomes dynamic as well...

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
       temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdx(i,j,1)
        sgs_horx(i,j,1)=-1.*temp(i,j,1)
     end do
   end do
end if

  do k=jz_min,ubc_jz
    do j=1,Ny
      do i=1,Nx
        temp(i,j,k)=(1./(0.5_rprec*(Pr_(i,j,k)+Pr_(i,j,k+1))))* &
                0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
         sgs_horx(i,j,k)=-1.*temp(i,j,k)
        ! temp(i,j,k)=(1./Pr_(i,j,k))*0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdx(i,j,k)
      end do
    end do
  end do


  call DDX (dtemp, temp)  

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,1) = (-1.*RHS(i,j,1) + dtemp(i,j,1))
       temp(i,j,1)=(1./Pr_(i,j,1))*Nu_t(i,j,1)*dsdy(i,j,1)
       sgs_hory(i,j,1)=-1.*temp(i,j,1)
     end do
   end do
end if

 do k=jz_min,ubc_jz
! Nu_t is on w nodes and dsdy is on uvp nodes for jz=2 to nz
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,k) = (-1.*RHS(i,j,k) + dtemp(i,j,k))
       temp(i,j,k)=(1./(0.5_rprec*(Pr_(i,j,k)+Pr_(i,j,k+1))))* &
               0.5_rprec*(Nu_t(i,j,k)+Nu_t(i,j,k+1))*dsdy(i,j,k)
       sgs_hory(i,j,k)=-1.*temp(i,j,k)
     end do
   end do
 end do


 call DDY (dtemp, temp)   
  
!c...Use MO flux at wall for the scalar sgs term !
!c Note that the total contribution to the scalar sgs term at
!c the first node comes from the surface flux computed above from
!c the specified heat flux, wt_s

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   do j=1,Ny
     do i=1,Nx
      RHS(i,j,1) = RHS(i,j,1) + dtemp(i,j,1)
      temp(i,j,1) = -1.*surf_flux(i,j)
      sgs_vert(i,j,1) = surf_flux(i,j)
     end do
   end do
 end if
! Note sgs_vert is -1*temp because sgs_vert is modeled as -Nu_t*dsdz/Pr
! while temp is the same term but w/o the minus sign due to the additional
! minus outside the scalar fluctuation flux term in RHS
! need to run this loop nz due to the nature of the differenetiation in ddz_w

 do k=jz_min,nz
   do j=1,Ny
     do i=1,Nx
       RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
       temp(i,j,k)=(1./Pr_(i,j,k))*Nu_t(i,j,k)*dsdz(i,j,k)
       sgs_vert(i,j,k)=-1.*temp(i,j,k)
     end do
   end do
 end do
!c...The SGS_z flux is on the W nodes, but DDZ_W will put it back on UVP nodes! 
!c Also note that sgs_vert(i,j,k) influences the computations in 
! OBUKHOV.f and is not involved in any computations in this routine.
! sgs_t3(i,j,1) (<w'theta'> is used for computing wt at the surface in OBUKHOV)

  call DDZ_w (dtemp, temp)
  do k=1,ubc_jz
    Do j=1,Ny
    do i=1,Nx
    RHS(i,j,k) = RHS(i,j,k) + dtemp(i,j,k)
    end do
    end do
  end do

end subroutine scalar_RHS_calc
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine step_scalar(scalar,RHS_pre,RHS_post)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar, RHS_pre, RHS_post
!real(kind=rprec)::wt_s_current
!cVK - This routine moves the scalar field (scalar in this case)
!cVK - forward in time using the scalar from previous time step
!cVK - and the RHS terms from the previous two time steps 
!cVK - using second order Adams-Bashforth scheme

! Note that in the last staments within this file, we set the value
! of scalar at the topmost node based on prescribed bc (inv_strength)
! and so, we could save some computation by only performing
! the scalar computation till Nz-1 global node...

do k=1,nz-1
      do j=1,ny
             do i=1,nx
                 scalar(i,j,k)= scalar(i,j,k)+dt*(1.5_rprec*RHS_pre(i,j,k)-0.5_rprec*RHS_post(i,j,k))
             end do
      end do
end do     

!VK Note that this subroutine was designed to be working with a set of scalars (incl.
!VK temperature and humidity and so, these boundary conditions as given below should
!VK be interpreted in the right context and not just for temperature
!VK For example, the first if block refers to a condition with humidity while the
!VK second and third statements are focussed to temperature

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
! if MPI - then clicks and is used for the process dealing wih the top nodes
! else in general is used for temp bc
! Note that L_z*nproc is the total vertical extent of the domain for the MPI and
! non-MPI cases ..(with MPI in place, we can not use L_z>z_i anymore)
     if ((L_z*nproc)>z_i) then ! for temperature and non-neutral case
         scalar(:,:,Nz)=scalar(:,:,Nz-1)+dTdz_top/T_scale*z_i*dz !ubc 
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
     else ! for everything else - neutral and passive scalars (may be modified depending on need)
         scalar(:,:,Nz)=scalar(:,:,Nz-1)
     end if
end if

end subroutine step_scalar

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine obukhov
use types,only:rprec
use sim_param,only:u,v,theta,qmix,path,w,txz,tyz
use scalars_module_q,only:sgs_q3
use bottombc !Includes patches subroutine, phi_m,psi_m,phi_h,psi_h,ustar_avg
use test_filtermodule
implicit none

integer:: jx,jy,i
real(kind=rprec), dimension(ld,ny):: wt_avg,wq_avg,theta_avg,u1,v1
real(kind=rprec), dimension(nx,ny):: x,x0,x0s,zeta,zeta0,zeta0s,u_avg
real(kind=rprec):: g_,wt_,wq_,ustar_,theta_,L_,zo_,wstar_avg
real(kind=rprec),save:: obuk_L,obuk_ustar,obuk_phi_m,obuk_phi_h,obuk_psi_m   
real(kind=rprec),save:: obuk_wt_sfc,obuk_psi_h,obuk_zo,obuk_wstar   


   if (jt .EQ. SCAL_init) then
    obuk_L=0._rprec
    obuk_ustar=0._rprec
    obuk_wstar=0._rprec
    obuk_phi_m=0._rprec
    obuk_phi_h=0._rprec
    obuk_psi_m=0._rprec
    obuk_psi_h=0._rprec
    obuk_zo=0._rprec
    psi_m=0._rprec
    psi_m0=0._rprec
    psi_h=0._rprec
    psi_h0=0._rprec
        if (.not. initsc) then
        sgs_t3(1:nx,1:ny,1)=wt_s(1:nx,1:ny)/(u_star*T_scale)
        end if
   end if

!  nondimensionalize g
   g_=(g*z_i)/(u_star**2)  ! end if

theta_avg=theta(:,:,1)*(1+0.61_rprec*0.001_rprec*qmix(:,:,1)*q_scale) 
wt_avg=sgs_t3(:,:,1) ! We need only the surface flux - defined by sgs
zo_=exp(sum(dlog(zo(1:nx,1:ny)))/float(nx*ny))

! averages over x-y plane @ z = 1
wt_=sum(sgs_t3(1:nx,1:ny,1))/float(nx*ny)
ustar_=sum(sqrt(u(1:nx,1:ny,1)**2+v(1:nx,1:ny,1)**2))/float(nx*ny)*vonK/&
(dlog(0.5_rprec*dz/zo_)-sum(psi_m(1:nx,1:ny))/float(nx*ny)+sum(psi_m0(1:nx,1:ny))/float(nx*ny))
theta_=sum(theta_avg(1:nx,1:ny))/float(nx*ny)


  u1=u(:,:,1)
  v1=v(:,:,1)
  call test_filter(u1,G_test_test)
  call test_filter(v1,G_test_test)

if (patch_flag==1) then  
  do jx=1,nx
    do jy=1,ny
      u_avg(jx,jy)=sqrt(u1(jx,jy)**2+v1(jx,jy)**2) 
      wt_avg(jx,jy)=sgs_t3(jx,jy,1)
      ustar_avg(jx,jy)=u_avg(jx,jy)*vonK/(dlog(0.5_rprec*dz/zo(jx,jy))-psi_m(jx,jy)+psi_m0(jx,jy))
      theta_avg(jx,jy)=theta_
    end do
  end do
end if
   

   do jx=1,nx
   do jy=1,ny
    L(jx,jy)=-ustar_avg(jx,jy)**3/(vonk*(g_/theta_avg(jx,jy))*wt_avg(jx,jy))
    zeta(jx,jy)=0.5_rprec*dz/L(jx,jy)
    zeta0(jx,jy)=zo(jx,jy)/L(jx,jy)
    zeta0s(jx,jy)=z_os(jx,jy)/L(jx,jy)
    
      if (zeta(jx,jy)<-20.0_rprec) then
              zeta(jx,jy)=-20.0_rprec
      end if

      if (zeta(jx,jy)>2.0_rprec) then
              zeta(jx,jy)=2.0_rprec
      end if

      if (zeta0(jx,jy)<-20.0_rprec) then
              zeta0(jx,jy)=-20.0_rprec
      end if

      if (zeta0(jx,jy)>2.0_rprec) then
              zeta0(jx,jy)=2.0_rprec
      end if
      
      if (zeta0s(jx,jy)<-20.0_rprec) then
              zeta0s(jx,jy)=-20.0_rprec
      end if

      if (zeta0s(jx,jy)>2.0_rprec) then
              zeta0s(jx,jy)=2.0_rprec
      end if


! for unstable conditions
      if ((L(jx,jy)<0._rprec) .and. (wt_avg(jx,jy) .ne. 0._rprec)) then
             x(jx,jy)=(1._rprec-16._rprec*zeta(jx,jy))**0.25_rprec

             psi_m(jx,jy)=2._rprec*dlog((1.+x(jx,jy))/2._rprec)+&
             dlog((1._rprec+x(jx,jy)**2)/2._rprec)-2._rprec*datan(x(jx,jy))+pi/2._rprec

             psi_h(jx,jy)=2._rprec*dlog((1._rprec+x(jx,jy)**2)/2._rprec)
     
             phi_m(jx,jy)=x(jx,jy)**(-1)
             phi_h(jx,jy)=x(jx,jy)**(-2)

             x0(jx,jy)=(1._rprec-16._rprec*zeta0(jx,jy))**.25_rprec
             x0s(jx,jy)=(1._rprec-16._rprec*zeta0s(jx,jy))**.25_rprec
             
             psi_m0(jx,jy)=2._rprec*dlog((1.+x0(jx,jy))/2._rprec)+&
             dlog((1._rprec+x0(jx,jy)**2)/2._rprec)-2._rprec*datan(x0(jx,jy))+pi/2._rprec
             
             psi_h0(jx,jy)=2._rprec*dlog((1._rprec+x0s(jx,jy)**2)/2._rprec)
! for stable condirions

           else if ((L(jx,jy)>0._rprec).and.(wt_avg(jx,jy).ne. 0._rprec)) then
                   
          phi_m(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          phi_h(jx,jy)=1._rprec+5.0_rprec*zeta(jx,jy)
          psi_m(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
          psi_h(jx,jy)=-1._rprec*5.0_rprec*zeta(jx,jy)
          psi_m0(jx,jy)=-1._rprec*5.0_rprec*zeta0(jx,jy)
          psi_h0(jx,jy)=-1._rprec*5.0_rprec*zeta0s(jx,jy)

      else
             psi_m(jx,jy)=0._rprec
             psi_m0(jx,jy)=0._rprec
             psi_h(jx,jy)=0._rprec
             psi_h0(jx,jy)=0._rprec
             phi_m(jx,jy)=1._rprec
             phi_h(jx,jy)=1._rprec

      end if ! (Loop5 ends)
  end do
  end do

  L_=-(ustar_**3)/(vonk*(g_/theta_)*wt_)
  wstar_avg=sign((g_/theta_*abs(wt_))**(1./3.),wt_)

! SKS
  if (mod(jt,100) == 0) then     ! 100 because wbase = 100

 
  write (6,7780) L_*z_i, ustar_*u_star, theta_*T_scale, wstar_avg*u_star, &
                (dz/2)/L_, wt_*u_star*T_scale
end if

7780 format ('L(m),ustar(m/s),thetav_1(K),wstar(m/s),z/L,wt_s(Km/s):',&
(6(1x,F15.7)))

!-------------------- OUTPUT ------------------------------
! Output the heat flux time series to a file to be used later
!    open (unit=47,file=path//'output/WT_sfc_tseries.out',status="unknown",position="append")
!    write(47,5168) (jt_total+1)*dt,wt_
 !   close(47)
!----------------------------------------------------------

  return
end subroutine obukhov 

end module scalars_module
