module bottombc
use types,only:rprec
use param,only:nx,ny,ld,cp,Le,Rv,Rd,rho_d,rho_w

implicit none

integer,parameter::num_patch=1,ptypes=1
real(kind=rprec),parameter::zo1=0.10_rprec,zo2=0.10_rprec
integer,parameter::square=0
real(kind=rprec),parameter::zo_out=0.000010_rprec, zo_square=0.000010_rprec
! SKS
integer,parameter::Nx_sq=5,Ny_sq=3
real(kind=rprec),parameter::theta_s1=290._rprec,theta_s2=290._rprec,q_s1=12._rprec,q_s2=12._rprec

real(kind=rprec),dimension(nx,ny)::ustar_avg
! This is the common local ustar computed in obukhov and used in wallstress.f90
real(kind=rprec),dimension(nx,ny)::zo,T_s,q_s
real(kind=rprec),dimension(nx,ny)::z_os
real(kind=rprec),dimension(nx,ny)::zom,wt_s,wq_s
!TS add for non-neutral case
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m,psi_m0,phi_h,psi_h,psi_h0
!VK The obukhov similarity functions are computed using obukhov(scalars_module.f90) for non-neutral scenario
integer,dimension(nx,ny)::patch
integer,dimension(num_patch)::patchnum

contains
subroutine patches()
!VK This assigns momentum roughness, temperature and wetness for the different patches
!VK and fills the lookup tables patch and patchnum. This is called from routines 
!VK patch_or_remote.f90 depending whether to use remotely-sensed data or patches
use param
implicit none
integer::i, j,ix,jy, patchtype, begini, endi, type
!open(unit=77,file=path//'spatial_flux.dat',status='unknown')

! IMPORTANT: read fluxes in W/m^2

open(unit=9001,file='./wts.dat',status='unknown')
open(unit=9002,file='./wqs.dat',status='unknown')
open(unit=9003,file='./zom.dat',status='unknown')
     
      do j=1,ny
          read(9001,*) (wt_s(i,j),i=1,nx)
          read(9002,*) (wq_s(i,j),i=1,nx)
          read(9003,*) (zom(i,j),i=1,nx)
      end do
      close(9001)
      close(9002)
     close(9003)

 !         read(9001,*) wt_s
 !         read(9002,*) wq_s
 !         read(9003,*) zom

 !     close(9001)
 !     close(9002)
 !     close(9003)
 !if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !    print *, 'size of wt_s', size(wt_s,1), size(wt_s,2)
 !     print *, 'size of wq_s', size(wq_s,1), size(wq_s,2)
 !     print *, 'size of zom', size(zom,1), size(zom,2)
! end if

   do j=1,ny
   do i=1,nx
     zo(i,j)=zom(i,j)/z_i
     T_s(i,j)=theta_s1/T_scale
     q_s(i,j)=q_s1/q_scale
     z_os(i,j)=zom(i,j)/(z_i*7.39_rprec)
     wt_s(i,j)=wt_s(i,j)/(rho_d*cp)
     wq_s(i,j)=1000.0_rprec*wq_s(i,j)/(rho_d*Le)   ! multiply by 1000 for g/kg humidity
     
   end do       
   end do

!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

!open(unit=9004,file=path//'output/wts_hom_ef5_output.out',status="unknown",position="append")
!do jy=1,ny
!write(9004,7168) (wt_s(ix,jy),ix=1,nx)
!end do
!close(9004)

!open(unit=9005,file=path//'output/wqs_hom_ef5_output.out',status="unknown",position="append")
!do jy=1,ny
!write(9005,7168) (wq_s(ix,jy),ix=1,nx)
!end do
!close(9005)

!open(unit=9006,file=path//'output/zo_hom_ef5_output.out',status="unknown",position="append")
!do jy=1,ny
!write(9006,7168) (zo(ix,jy),ix=1,nx)
!end do
!close(9006)

!end if

7168     format(1400(E14.5))

end subroutine patches
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX
subroutine remote_to_patch(T_s_in,T_or_z)
implicit none
integer:: ii,jj,kk,T_or_z
real(kind=rprec),dimension(nx,ny):: T_s_in,dummy
real(kind=rprec):: sigma_multiplier,crap1,crap2,crap3
real(kind=rprec):: patch_cold, patch_hot

sigma_multiplier=1._rprec
dummy=0._rprec
if (T_or_z == 1) then
    crap1=sum(T_s_in)/float(nx*ny)
elseif (T_or_z == 2) then
    crap1=exp(sum(dlog(T_s_in))/float(nx*ny))
else
    print *,'Wrong choice of T_or_z in remote_to_patch().. STOPPING'
    stop
end if

crap2=sqrt(sum((T_s_in-crap1)**2)/float(nx*ny))

patch_hot=crap1+sigma_multiplier*crap2;
patch_cold=crap1-sigma_multiplier*crap2;

print *,'mean, std, patch_hot, patch_cold = ',crap1,crap2,patch_hot,patch_cold
!First do the patch business for temperature
if ((patch_hot .lt. 0._rprec) .OR. (patch_cold .lt. 0._rprec)) then
print *,'Hot & Cold patch calculation yields negative T_s'
print *,'Trying sigma_multiplier = 0.75'
     if (patch_cold < 0._rprec) then
        sigma_multiplier=0.75_rprec
        patch_cold=crap1-sigma_multiplier*crap2;
        patch_hot=crap1+sigma_multiplier*crap2;
        crap3=0.5_rprec*(patch_cold+patch_hot)
        print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
        if (patch_cold < 0._rprec) then
           print *,'sigma = 0.75 FAILED, Trying sigma_= 0.5'
           sigma_multiplier=0.5_rprec
           patch_cold=crap1-sigma_multiplier*crap2;
           patch_hot=crap1+sigma_multiplier*crap2;
           crap3=0.5_rprec*(patch_cold+patch_hot)
           print *,'NEW:mean, patch_hot, patch_cold = ',crap3,patch_hot,patch_cold
           if (patch_cold < 0._rprec) then
              print *,'sigma = 0.5 FAILED, STOPPING NOW...'
              print *,'This message is from the subroutine patch_to_remote in scalars_module2.f90.'
           end if
        end if
     end if
end if
! Assign to patches
if (T_or_z .eq. 1) then
  dummy(1:nx/2,:)=patch_hot
  dummy(nx/2+1:nx,:)=patch_cold
   print *,'2 patch T_s field generated: patch_hot,patch_cold = ',patch_hot,patch_cold
elseif (T_or_z .eq. 2) then
! Assign cold patch to the first patch and hot patch to the second one
! Note that this is exactly opposite to what we do for temperature as
! in REALITY, hot surfaces have low roughness and vice-versa
  dummy(1:nx/2,:)=patch_cold
  dummy(nx/2+1:nx,:)=patch_hot
   print *,'2 patch roughnesss field generated: patch_smooth,patch_rough = ',patch_cold,patch_hot
end if

T_s_in=dummy

end subroutine remote_to_patch
!!!xxxxxxxxxxxx--------VIJ-------------XXXXXXXXXXXXXXXXXXX

subroutine avgpatch(u_avg)
! computes the averaged value of a variable (at the wall) over a patch
! and assigns it to an nx X ny array

! sc: 
! note: this is inefficient: should calculate directly to u_avg --no maybe its good
use sim_param,only:u
implicit none
real(kind=rprec),dimension(nx,ny),intent(inout)::u_avg
integer::i,j,k
real(kind=rprec),dimension(ptypes)::temp

temp=0._rprec
do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      temp(patch(i,j))=temp(patch(i,j))+u(i,j,1)/real(patchnum(patch(i,j)))
   end if
end do
end do
end do

do j=1,ny
do i=1,nx
do k=1,ptypes
   if (patch(i,j).eq.k) then
      u_avg(i,j)=real(temp(k))
   end if
end do
end do
end do
end subroutine avgpatch

end module bottombc
