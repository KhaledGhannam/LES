module scalars_module_q
! HUMIDITY subroutines in place but not yet turned on !!
use types,only:rprec
use param 
use sim_param,only:u,v,w,qmix,path
use bottombc ! Includes patches subroutine
use sgsmodule,only:Nu_t,Pr_t
implicit none
integer, parameter:: tag_counter = 323
 
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
!! KMG: added sgs_t1 and sgs_t2 for horizontal sgs heat fluxes Pi_1 and Pi_2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: Pr_
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dqdz,dqdx,dqdy
! Only ones needed for output..Might need to add x and y derivatives here in case they need to be outputted
! Right now they are in the "scalar"_all_in_one routines below !!
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS_qf,RHS_q
real(kind=rprec), dimension(ld,ny,$lbz:nz):: sgs_q3,sgs_q1,sgs_q2 ! defines the surface sgs flux

real(kind=rprec),dimension(ld,ny)::q_s_filtered ! filtered T_s for calc of wT_s


real(kind=rprec),dimension(ld,ny,$lbz:nz)::qmix_old

contains
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine q_all_in_one
use topbc,only:sponge
use test_filtermodule
implicit none
real(kind=rprec),dimension(nx,ny)::wq_s_current
integer::jz,ios,counter
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::dqdx,dqdy
real(kind=rprec),dimension($lbz:nz)::sponge_q

 if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        wq_s_current(1:nx,1:ny)=wq_s(1:nx,1:ny)
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

call filt_da(qmix,dqdx,dqdy)
call ddz_uv (dqdz,qmix)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
   dqdz(:,:,Nz)=dqdz_top/q_scale*z_i ! Valid for temperature
end if

$if ($MPI)
! print *,'synchronizing in theta_all_in for coord = ',coord
! Need to synchronize w and dTdz across the processors for uvp node
! based computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (w(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+43,  &
                     w(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+43,   &
                     comm, status, ierr)   
  call mpi_sendrecv (dqdz(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+44,  &
                     dqdz(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+44,   &
                     comm, status, ierr)   

! Also need to synchronize Nu_t across the processors for 
! computation in scalar_RHS_calc (calculation of RHS_m)
  call mpi_sendrecv (Nu_t(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+60,  &
                     Nu_t(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+60,   &
                     comm, status, ierr)   
  call mpi_sendrecv (Nu_t(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+70,  &
                     Nu_t(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+70,   &
                     comm, status, ierr)
$endif

if (S_FLAG) then
   RHS_qf=RHS_q


! Perform test filtering of T_s for calculation of surf fluxes
     if ((jt_total .eq. SCAL_init) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
       print *,'q_s b4 filtering',sqrt(sum((q_s-sum(q_s)/float(nx*ny))**2))/float(nx*ny)
       q_s_filtered(1:nx,1:ny)=q_s
       call test_filter(q_s_filtered,G_test)
       q_s=q_s_filtered(1:nx,1:ny)
       print *,'q_s after filtering',sqrt(sum((q_s-sum(q_s)/float(nx*ny))**2))/float(nx*ny)
     end if
     if ((jt_total .eq. nsteps-1) .and. (lbc .eq. 0) .and. (remote_homog_flag .eq. 0)) then
        print *,'q_s after filtering',sqrt(sum((q_s-sum(q_s)/float(nx*ny))**2))/float(nx*ny)
     end if

! call scalar_RHS_calc(T_s,z_os,RHS_T,sgs_t3,jt,psi_h,phi_h,Pr_,wt_s_current)
qmix_old(:,:,:)=qmix(:,:,:)

call scalar_q_RHS_calc(qmix,dqdx,dqdy,dqdz,q_s,z_os,RHS_q,sgs_q1,sgs_q2,sgs_q3,wq_s_current)

! subroutine looks like this later
! subroutine scalar_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)


if (ubc==1 .and. damping_method==2) then ! add the damping term to the scalar equation
   do jz=1,nz-1
   RHS_q(1:nx,1:ny,jz)=RHS_q(1:nx,1:ny,jz)-0.5_rprec*(sponge(jz)+sponge(jz+1))*&
                 (qmix(1:nx,1:ny,jz)-sum(qmix(1:nx,1:ny,jz))/(nx*ny))
   end do
end if

if (jt .eq. SCAL_INIT .and. (.not. initsc)) then
   RHS_qf=RHS_q
end if

call step_scalar_q(qmix,RHS_q,RHS_qf)
end if



$if ($MPI)
  call mpi_sendrecv (qmix(1, 1, 1), ld*ny, MPI_RPREC, down, tag_counter+20,  &
                     qmix(1, 1, nz), ld*ny, MPI_RPREC, up, tag_counter+20,   &
                     comm, status, ierr)   
  call mpi_sendrecv (qmix(1, 1, nz-1), ld*ny, MPI_RPREC, up, tag_counter+21,  &
                     qmix(1, 1, 0), ld*ny, MPI_RPREC, down, tag_counter+21,   &
                     comm, status, ierr)
$endif

end subroutine q_all_in_one

!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------

subroutine scalar_q_RHS_calc(scalar,dsdx,dsdy,dsdz,S_Surf,z_os,RHS,sgs_horx,sgs_hory,sgs_vert,surf_flux_current)
use test_filtermodule
implicit none
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
integer::i,j,k,jz
integer::jz_min,ubc_jz
real:: crap2
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dsdx,dsdy,dsdz
real(kind=rprec),dimension(ld,ny,$lbz:nz):: RHS,temp
real(kind=rprec),dimension(ld,ny,$lbz:nz):: scalar
real(kind=rprec),dimension(ld,ny,$lbz:nz):: dtemp,sgs_vert,sgs_horx,sgs_hory
real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: u_m,v_m,w_m,dsdx_m,dsdy_m,dsdz_m
 real(kind=rprec),dimension(ld_big,ny2,$lbz:nz):: RHS_m
real(kind=rprec),dimension(nx,ny):: surf_flux_current,ustar_local,S_Surf,surf_flux,z_os,wq_s2
real(kind=rprec),dimension(ld,ny):: scalar_node_1 ! test filtered and used for computing surface flux

character (64) :: fname_hetero

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

         ustar_local(:,:)=ustar_avg(:,:)
          scalar_node_1=scalar(:,:,1)
          call test_filter(scalar_node_1,G_test) 
            if (lbc==1.and.scalar(1,1,1)<2) then
                     
                    wq_s2(1:nx,1:ny)=wq_s(1:nx,1:ny)
                    do i=1,Nx
                    do j=1,Ny
                    surf_flux(i,j)=wq_s2(i,j)/q_scale/u_star
                    end do
                    end do
  
            else if (lbc==0.and.scalar(1,1,1)<2) then
                    do i=1,Nx
                    do j=1,Ny 
                    surf_flux(i,j)=(S_Surf(i,j)-scalar_node_1(i,j))*vonk*ustar_local(i,j)&
                                   /(dlog(dz/(2._rprec*z_os(i,j)))-psi_h(i,j))
                   end do 
                   end do
             end if

!c....Now we have the lowest dsdz on the UVP nodes all others on w nodes
            do i=1,Nx
            do j=1,Ny
            dsdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar_local(i,j)*vonk*DZ/2._rprec)

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

end subroutine scalar_q_RHS_calc
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
!----------------------XXXXXXXXXXXVIJXXXXXXXXXXX-----------------------
subroutine step_scalar_q(scalar_q,RHS_pre_q,RHS_post_q)
implicit none
integer:: i,j,k
$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
real(kind=rprec),dimension(ld,ny,$lbz:nz)::scalar_q, RHS_pre_q, RHS_post_q

do k=1,nz-1
      do j=1,ny
             do i=1,nx
                 scalar_q(i,j,k)= scalar_q(i,j,k)+dt*(1.5_rprec*RHS_pre_q(i,j,k)-0.5_rprec*RHS_post_q(i,j,k))
             end do
      end do
end do     

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
     if ((L_z*nproc)>z_i) then ! for temperature and non-neutral case
         scalar_q(:,:,Nz)=scalar_q(:,:,Nz-1)+dqdz_top/q_scale*z_i*dz !ubc 
! inv_strength - refers to the slope of the inversion ~ 0.003K/Km for temperature
     else ! for everything else - neutral and passive scalars (may be modified depending on need)
         scalar_q(:,:,Nz)=scalar_q(:,:,Nz-1)
     end if
end if

end subroutine step_scalar_q


end module scalars_module_q
