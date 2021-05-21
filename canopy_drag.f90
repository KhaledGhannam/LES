subroutine canopy_drag(lad,drag_coeff2, drag_u, drag_v, drag_w)

use types,only:rprec
use param
use sim_param,only:u,v,w,wind_xy,wind_z

implicit none

$if ($MPI)
  !--this dimensioning adds a ghost layer for finite differences
  !--its simpler to have all arrays dimensioned the same, even though
  !--some components do not need ghost layer
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(rprec),dimension(ld,ny,1:((nz-1)*nproc)), intent(in)::lad
real (rprec), dimension (ld, ny, $lbz:nz), intent (out) :: drag_coeff2, drag_u, drag_v, drag_w
integer::jx,jy,jz,Nzz,k_global,k
! ------------------------------------------------------------------
! --------------------------------------

 do k=1,nz-1
                   k_global=k+coord*(nz-1)

             if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord .eq. 0))) then
             wind_xy(:,:,1)=sqrt(u(:,:,1)**2+v(:,:,1)**2+(0.25_rprec*w(:,:,k+1))**2)
             wind_z(:,:,1)=0._rprec
             if (constant_Cd) then
                     drag_coeff2 (:,:,1)=0.4_rprec
             else
                     drag_coeff2(:,:,1)=min((wind_xy(:,:,1)/(0.22_rprec/u_star))**(-2._rprec/3._rprec),0.8_rprec)
             end if

             drag_u(:,:,1)=-drag_coeff2(:,:,1)*(0.5*(lad(:,:,1)+lad(:,:,2)))*wind_xy(:,:,1)*u(:,:,1)
             drag_v(:,:,1)=-drag_coeff2(:,:,1)*(0.5*(lad(:,:,1)+lad(:,:,2)))*wind_xy(:,:,1)*v(:,:,1)
             drag_w(:,:,1)=0._rprec

             else
             wind_xy(:,:,k)=sqrt(u(:,:,k)**2+v(:,:,k)**2+(0.5_rprec*(w(:,:,k)+w(:,:,k+1)))**2)
             wind_z(:,:,k)=sqrt((0.5*(u(:,:,k)+u(:,:,k-1)))**2+(0.5*(v(:,:,k)+v(:,:,k-1)))**2+w(:,:,k)**2)

             if (constant_Cd) then
                     drag_coeff2(:,:,k)=0.4
             else
                     drag_coeff2(:,:,k)=min((wind_xy(:,:,k)/(0.22_rprec/u_star))**(-2._rprec/3._rprec),0.8_rprec)
             end if

             drag_u(:,:,k)=-drag_coeff2(:,:,k)*(0.5*(lad(:,:,k_global)+lad(:,:,k_global+1)))*wind_xy(:,:,k)*u(:,:,k)
             drag_v(:,:,k)=-drag_coeff2(:,:,k)*(0.5*(lad(:,:,k_global)+lad(:,:,k_global+1)))*wind_xy(:,:,k)*v(:,:,k)
             drag_w(:,:,k)=-0.5*(drag_coeff2(:,:,k)+drag_coeff2(:,:,k-1))*lad(:,:,k_global)*wind_z(:,:,k)*w(:,:,k)
        end if
     end do



end subroutine canopy_drag

