module bottombc
use types,only:rprec
use param,only:nx,ny,ld,cp,Le,Rv,Rd,rho_d,rho_w

implicit none


real(kind=rprec),parameter::zo1=0.10_rprec
real(kind=rprec),parameter::theta_s1=290._rprec,q_s1=12._rprec

real(kind=rprec),dimension(nx,ny)::ustar_avg
! This is the common local ustar computed in obukhov and used in wallstress.f90
real(kind=rprec),dimension(nx,ny)::zo,T_s,q_s
real(kind=rprec),dimension(nx,ny)::z_os
real(kind=rprec),dimension(nx,ny)::zom,wt_s,wq_s
!TS add for non-neutral case
real(kind=rprec),dimension(nx,ny)::phi_m,psi_m,psi_m0,phi_h,psi_h,psi_h0

!----------------------------------
contains
!----------------------------------

subroutine patches()

implicit none
integer::i, j,ix,jy


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


end module bottombc
