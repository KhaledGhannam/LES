module output_slices

use types,only:rprec
use param
use sim_param
use scalars_module,only:sgs_t3
use scalars_module_q,only:sgs_q3
use bottombc
! -----------------------
implicit none
! -----------------------------


! -------------------------------------
contains
! --------------------------------------------

subroutine write_slices()
use types,only:rprec
use param
use sim_param
use scalars_module,only:sgs_t3
use scalars_module_q,only:sgs_q3
use bottombc


implicit none

integer::jx,jy,jz,Nzz,k_global,i,j,k
real(kind=rprec),dimension(nx,nz-1)::u_slice1_xz,w_slice1_xz,q_slice1_xz,theta_slice1_xz
real(kind=rprec),dimension(nx,nz-1)::u_slice2_xz,w_slice2_xz,q_slice2_xz,theta_slice2_xz
real(kind=rprec),dimension(nx,nz-1)::u_slice3_xz,w_slice3_xz,q_slice3_xz,theta_slice3_xz

real(kind=rprec),dimension(ny,nz-1)::u_slice1_yz,w_slice1_yz,q_slice1_yz,theta_slice1_yz
real(kind=rprec),dimension(ny,nz-1)::u_slice2_yz,w_slice2_yz,q_slice2_yz,theta_slice2_yz
real(kind=rprec),dimension(ny,nz-1)::u_slice3_yz,w_slice3_yz,q_slice3_yz,theta_slice3_yz

real(kind=rprec),dimension(nx,ny)::u_slice1_xy,q_slice1_xy, w_slice1_xy,theta_slice1_xy
real(kind=rprec),dimension(nx,ny)::u_slice2_xy,q_slice2_xy, w_slice2_xy,theta_slice2_xy
real(kind=rprec),dimension(nx,ny)::u_slice3_xy,q_slice3_xy, w_slice3_xy,theta_slice3_xy
real(kind=rprec),dimension(nx,ny)::u_slice4_xy,q_slice4_xy, w_slice4_xy,theta_slice4_xy

real(kind=rprec),dimension(:,:),allocatable::slice_out,slicey_out
real(kind=rprec),dimension(nx,ny)::arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8

do jz=1,nz-1

    ! if((jz .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord .eq. 0))) then
     !    arg1(1:nx,1:ny)=0.25_rprec*(w(1:nx,1:ny,1)-sum(w(1:nx,1:ny,1))/float(nx*ny) )
     !    arg2(1:nx,1:ny)=u(1:nx,1:ny,1)-sum(u(1:nx,1:ny,1))/float(nx*ny)
     !    arg3(1:nx,1:ny)=qmix(1:nx,1:ny,1)-sum(qmix(1:nx,1:ny,1))/float(nx*ny)
     !    arg4(1:nx,1:ny)=theta(1:nx,1:ny,1)-sum(theta(1:nx,1:ny,1))/float(nx*ny)

      !  else
        arg1(1:nx,1:ny)=w(1:nx,1:ny,jz)-sum(w(1:nx,1:ny,jz))/float(nx*ny)

        arg2(1:nx,1:ny)=u(1:nx,1:ny,jz)-sum(u(1:nx,1:ny,jz))/float(nx*ny)
        arg3(1:nx,1:ny)=qmix(1:nx,1:ny,jz)-sum(qmix(1:nx,1:ny,jz))/float(nx*ny)
        arg4(1:nx,1:ny)=theta(1:nx,1:ny,jz)-sum(theta(1:nx,1:ny,jz))/float(nx*ny)
        
     !   end if

     w_slice1_xz(1:nx,jz)=arg1(1:nx,90)
     u_slice1_xz(1:nx,jz)=arg2(1:nx,90)
     q_slice1_xz(1:nx,jz)=arg3(1:nx,90)
     theta_slice1_xz(1:nx,jz)=arg4(1:nx,90)
     
     w_slice2_xz(1:nx,jz)=arg1(1:nx,130)
     u_slice2_xz(1:nx,jz)=arg2(1:nx,130)
     q_slice2_xz(1:nx,jz)=arg3(1:nx,130)
     theta_slice2_xz(1:nx,jz)=arg4(1:nx,130)

     w_slice3_xz(1:nx,jz)=arg1(1:nx,200)
     u_slice3_xz(1:nx,jz)=arg2(1:nx,200)
     q_slice3_xz(1:nx,jz)=arg3(1:nx,200)
     theta_slice3_xz(1:nx,jz)=arg4(1:nx,200)


     w_slice1_yz(1:ny,jz)=arg1(55,1:ny)
     u_slice1_yz(1:ny,jz)=arg2(55,1:ny)
     q_slice1_yz(1:ny,jz)=arg3(55,1:ny)
     theta_slice1_yz(1:ny,jz)=arg4(55,1:ny)
     
     w_slice2_yz(1:ny,jz)=arg1(130,1:ny)
     u_slice2_yz(1:ny,jz)=arg2(130,1:ny)
     q_slice2_yz(1:ny,jz)=arg3(130,1:ny)
     theta_slice2_yz(1:ny,jz)=arg4(130,1:ny)

     w_slice3_yz(1:ny,jz)=arg1(300,1:ny)
     u_slice3_yz(1:ny,jz)=arg2(300,1:ny)
     q_slice3_yz(1:ny,jz)=arg3(300,1:ny)
     theta_slice3_yz(1:ny,jz)=arg4(300,1:ny)

     
     if (jz .eq. 1) then

          if (coord==1) then
              w_slice1_xy(1:nx,1:ny)=arg1(1:nx,1:ny)
              u_slice1_xy(1:nx,1:ny)=arg2(1:nx,1:ny)
              q_slice1_xy(1:nx,1:ny)=arg3(1:nx,1:ny)
              theta_slice1_xy(1:nx,1:ny)=arg4(1:nx,1:ny)

              open(unit=2293,file=path//'output/u_slice1_xy_coord1.out',status="unknown",position="append")
                do jy=1,ny
                  write(2293,7185) (u_slice1_xy(i,jy),i=1,nx)
                end do
              close(2293)


              open(unit=2294,file=path//'output/w_slice1_xy_coord1.out',status="unknown",position="append")
                do jy=1,ny
                  write(2294,7185) (w_slice1_xy(i,jy),i=1,nx)
                end do
              close(2294)


              open(unit=2295,file=path//'output/q_slice1_xy_coord1.out',status="unknown",position="append")
                do jy=1,ny
                  write(2295,7185) (q_slice1_xy(i,jy),i=1,nx)
                end do
              close(2295)


              open(unit=2296,file=path//'output/theta_slice1_xy_coord1.out',status="unknown",position="append")
                do jy=1,ny
                  write(2296,7185) (theta_slice1_xy(i,jy),i=1,nx)
                end do
              close(2296)

              
              else if (coord==8) then
              w_slice2_xy(1:nx,1:ny)=arg1(1:nx,1:ny)
              u_slice2_xy(1:nx,1:ny)=arg2(1:nx,1:ny)
              q_slice2_xy(1:nx,1:ny)=arg3(1:nx,1:ny)
              theta_slice2_xy(1:nx,1:ny)=arg4(1:nx,1:ny)

              open(unit=2623,file=path//'output/u_slice2_xy_coord8.out',status="unknown",position="append")
                do jy=1,ny
                  write(2623,7185) (u_slice2_xy(i,jy),i=1,nx)
                end do
              close(2623)


              open(unit=2624,file=path//'output/w_slice2_xy_coord8.out',status="unknown",position="append")
                do jy=1,ny
                  write(2624,7185) (w_slice2_xy(i,jy),i=1,nx)
                end do
              close(2624)


              open(unit=2625,file=path//'output/q_slice2_xy_coord8.out',status="unknown",position="append")
                do jy=1,ny
                  write(2625,7185) (q_slice2_xy(i,jy),i=1,nx)
                end do
              close(2625)


              open(unit=2626,file=path//'output/theta_slice2_xy_coord8.out',status="unknown",position="append")
                do jy=1,ny
                  write(2626,7185) (theta_slice2_xy(i,jy),i=1,nx)
                end do
              close(2626)

              
              else if (coord==45) then
              w_slice3_xy(1:nx,1:ny)=arg1(1:nx,1:ny)
              u_slice3_xy(1:nx,1:ny)=arg2(1:nx,1:ny)
              q_slice3_xy(1:nx,1:ny)=arg3(1:nx,1:ny)
              theta_slice3_xy(1:nx,1:ny)=arg4(1:nx,1:ny)

              open(unit=2823,file=path//'output/u_slice3_xy_coord45.out',status="unknown",position="append")
                do jy=1,ny
                  write(2823,7185) (u_slice3_xy(i,jy),i=1,nx)
                end do
              close(2823)


              open(unit=2824,file=path//'output/w_slice3_xy_coord45.out',status="unknown",position="append")
                do jy=1,ny
                  write(2824,7185) (w_slice3_xy(i,jy),i=1,nx)
                end do
              close(2824)


              open(unit=2825,file=path//'output/q_slice3_xy_coord45.out',status="unknown",position="append")
                do jy=1,ny
                  write(2825,7185) (q_slice3_xy(i,jy),i=1,nx)
                end do
              close(2825)


              open(unit=2826,file=path//'output/theta_slice3_xy_coord45.out',status="unknown",position="append")
                do jy=1,ny
                  write(2826,7185) (theta_slice3_xy(i,jy),i=1,nx)
                end do
              close(2826)
             
             
            else if (coord==120) then
              w_slice4_xy(1:nx,1:ny)=arg1(1:nx,1:ny)
              u_slice4_xy(1:nx,1:ny)=arg2(1:nx,1:ny)
              q_slice4_xy(1:nx,1:ny)=arg3(1:nx,1:ny)
              theta_slice4_xy(1:nx,1:ny)=arg4(1:nx,1:ny)

              open(unit=2827,file=path//'output/u_slice4_xy_coord120.out',status="unknown",position="append")
                do jy=1,ny
                  write(2827,7185) (u_slice4_xy(i,jy),i=1,nx)
                end do
              close(2827)


              open(unit=2828,file=path//'output/w_slice4_xy_coord120.out',status="unknown",position="append")
                do jy=1,ny
                  write(2828,7185) (w_slice4_xy(i,jy),i=1,nx)
                end do
              close(2828)


              open(unit=2829,file=path//'output/q_slice4_xy_coord120.out',status="unknown",position="append")
                do jy=1,ny
                  write(2829,7185) (q_slice4_xy(i,jy),i=1,nx)
                end do
              close(2829)


              open(unit=2830,file=path//'output/theta_slice4_xy_coord120.out',status="unknown",position="append")
                do jy=1,ny
                  write(2830,7185) (theta_slice4_xy(i,jy),i=1,nx)
                end do
              close(2830)
               
              end if
              
      end if
     end do

        allocate(slice_out(1:nx,1:(nz_tot-1)));
        call collocate_slices(u_slice1_xz,slice_out,8037,'u_slice1')
        call collocate_slices(w_slice1_xz,slice_out,8038,'w_slice1')
        call collocate_slices(q_slice1_xz,slice_out,8039,'q_slice1')
        call collocate_slices(theta_slice1_xz,slice_out,8040,'theta_slice1')
        
        call collocate_slices(u_slice2_xz,slice_out,8044,'u_slice2')
        call collocate_slices(w_slice2_xz,slice_out,8045,'w_slice2')
        call collocate_slices(q_slice2_xz,slice_out,8046,'q_slice2')
        call collocate_slices(theta_slice2_xz,slice_out,8047,'theta_slice2')
        
        call collocate_slices(u_slice3_xz,slice_out,8048,'u_slice3')
        call collocate_slices(w_slice3_xz,slice_out,8049,'w_slice3')
        call collocate_slices(q_slice3_xz,slice_out,8050,'q_slice3')
        call collocate_slices(theta_slice3_xz,slice_out,8051,'theta_slice3')
        
        deallocate(slice_out)
        
        allocate(slicey_out(1:ny,1:(nz_tot-1)));
        call collocate_slicesy(u_slice1_yz,slicey_out,8052,'u_slice1')
        call collocate_slicesy(w_slice1_yz,slicey_out,8053,'w_slice1')
        call collocate_slicesy(q_slice1_yz,slicey_out,8054,'q_slice1')
        call collocate_slicesy(theta_slice1_yz,slicey_out,8055,'theta_slice1')
        
        call collocate_slicesy(u_slice2_yz,slicey_out,8056,'u_slice2')
        call collocate_slicesy(w_slice2_yz,slicey_out,8057,'w_slice2')
        call collocate_slicesy(q_slice2_yz,slicey_out,8058,'q_slice2')
        call collocate_slicesy(theta_slice2_yz,slicey_out,8059,'theta_slice2')
        
        call collocate_slicesy(u_slice3_yz,slicey_out,8060,'u_slice3')
        call collocate_slicesy(w_slice3_yz,slicey_out,8061,'w_slice3')
        call collocate_slicesy(q_slice3_yz,slicey_out,8062,'q_slice3')
        call collocate_slicesy(theta_slice3_yz,slicey_out,8063,'theta_slice3')

        deallocate(slicey_out)
        
7185     format(1400(E14.5))

end subroutine write_slices





! --------------------------------------------------------
! -----------------------------------------------
subroutine collocate_slices(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,ind3,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(nx,nz-1)::avg_var_proc
real(kind=rprec),dimension(nx,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/xz_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1,1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
           do ind3=1,nz_tot-1
           
            write(file_ind,7307)(avg_var_tot_domain(ind1,ind3),ind1=1,nx)
           
           end do
        close(file_ind)
  end if
7307     format(1400(E14.5))

end subroutine collocate_slices
! -----------------------------------------------
! ----------------------------------------------

! --------------------------------------------------------
! -----------------------------------------------
subroutine collocate_slicesy(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,ind3,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(ny,nz-1)::avg_var_proc
real(kind=rprec),dimension(ny,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/yz_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1,1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
           do ind3=1,nz_tot-1
           
            write(file_ind,7308)(avg_var_tot_domain(ind1,ind3),ind1=1,ny)
           
           end do
        close(file_ind)
  end if
7308     format(1400(E14.5))

end subroutine collocate_slicesy
! -----------------------------------------------
! ----------------------------------------------


end module output_slices
