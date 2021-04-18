module io
use types,only:rprec
use param, only : ld, nx, ny, nz, write_inflow_file, path,  &
                  use_avgslice, USE_MPI, coord, rank, nproc,      &
                  average_dim_num, nz_tot,jt_total,p_count,dt,z_i,u_star

!use sim_param,only:L11t,L22t,L33t,Q11t,Q22t,Q33t
implicit none
private
public :: jt_total
public :: openfiles,output_loop,output_final,            &
          mean_u,mean_u2,mean_v,mean_v2,mean_w,mean_w2,  &
          inflow_read, inflow_write, avg_stats
public :: average_dim_select_flag, dim1_size, dim2_size,        &
          dim1_global, dim2_global, collocate_MPI_averages,     &
          compute_avg_var 

integer,parameter::num_hour_out=1
integer,parameter::base=50000,nwrite=base
! SKS
! Shouldnt be hard coded..base corresponds to the 
! time interval after which you want to write file
! So can be a factor of nsteps.
! SKS

logical,parameter:: io_spec=.false.,output_fields_3d_flag=.false.
integer,parameter::spec_write_freqz=600, fields_3d_write_freqz=6
integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! slice_inst sets the value of the y node for the chosen x-z inst slice
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process 
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical,parameter:: inst_slice_flag=.false.
integer,parameter:: num_vars=4 ! works currently only for u,v,w,T due to the size difference in Cs
integer,parameter:: slice_inst=(nz_tot-1)/2, inst_slice_freqz=5, inst_array_lim=200

logical, parameter :: cumulative_time = .true.
character (*), parameter :: fcumulative_time = path // 'total_time.dat'

integer, parameter :: n_avg_stats = 5000000 !--interval for updates in avg_stats
character (*), parameter :: end_hdr_avg = '# end header'


!! --------------------------------------------------------------------
integer,parameter :: average_dim_select_flag=1-(average_dim_num/2) 
! The variable average_dim_select_flag generates the following values based
! on the value of average_dim_num in param.f90 :- 
! a) average_dim_num = 2 : 
! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
! b) average_dim_num = 1 : 
! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
integer, parameter :: dim1_size=average_dim_select_flag*(nx-nz+1)+nz-1
integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
integer, parameter :: dim1_global=average_dim_select_flag*(nx-nz_tot+1)+nz_tot-1
integer, parameter :: dim2_global=average_dim_select_flag*(nz_tot-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

!!!!  time_spec>0 output time series spectrum (need additional calcu.)
integer,parameter::time_spec=0
integer::n_obs, jt_total_init
integer,allocatable::obs_pt(:,:)

!!!!  io_mean=.true. output small domain time-averaged velocity
logical,parameter::io_mean=.false.
integer,parameter::jx_pls=1,jx_ple=nx,width=0
integer,parameter::jy_pls=ny/2-width,jy_ple=ny/2+width
real(kind=rprec),dimension(jx_pls:jx_ple,jy_pls:jy_ple,nz)::&
     mean_u,mean_v,mean_w,mean_u2,mean_v2,mean_w2

!!!!  io_lambda2
logical,parameter::io_lambda2=.false.
real(kind=rprec),dimension(nx,ny,nz)::lam2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine openfiles()
use sim_param,only:path
implicit none

!--to hold file names
character (64) :: temp
character (64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
                  fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

integer::i

logical :: exst

if (cumulative_time) then

  inquire (file=fcumulative_time, exist=exst)
  if (exst) then
    open (1, file=fcumulative_time)
    read (1, *) jt_total
    jt_total_init=jt_total 
    close (1)
  else  !--assume this is the first run on cumulative time
    write (*, *) 'file ', fcumulative_time, ' not found'
    write (*, *) 'assuming jt_total = 0'
    jt_total = 0
    jt_total_init=jt_total 
  end if

end if

!if ((.not. USE_MPI) .or. (USE_MPI .and. rank == 0)) then
!  open(13,file=path//'output/check_ke.out',status="unknown",position="append")
  ! SKS
!  open(14,file=path//'output/Cs2byPr_t.out',status="unknown",position="append")
!  open(84,file=path//'output/spec_u.out',status="unknown",position="append")
!  open(85,file=path//'output/spec_v.out',status="unknown",position="append")
!  open(86,file=path//'output/spec_w.out',status="unknown",position="append")
!  open(87,file=path//'output/spec_T.out',status="unknown",position="append")
!  open(88,file=path//'output/spec_freq.out',status="unknown",position="append")
!  open(111,file=path//'output/tkeVsz.dat',status="unknown",position="append")
!  ! SKS
!end if

!if(time_spec.gt.0)then
!  open(15,file=path//'output/velspec.out',form='unformatted',position='append')
!  if(jt_total.eq.0)rewind(15)
!endif

!if(io_mean)then
!  open(51,file=path//'output/mean_u.out',form='unformatted',position='append')
!  if(jt_total.eq.0)then
!    rewind(51)
!    write(51)jx_pls,jx_ple,jy_pls,jy_ple
!  endif
!endif

!fCS1plan = path // 'output/CS1plan.out'
!fCS2plan = path // 'output/CS2plan.out'
!fCS4plan = path // 'output/CS4plan.out'
!fVISCplan = path // 'output/VISCplan.out'
!fDISSplan = path // 'output/DISSplan.out'
!fCS1Vplan = path // 'output/CS1Vplan.out'
!fCS2Vplan = path // 'output/CS2Vplan.out'
!fCS4Vplan = path // 'output/CS4Vplan.out'

!$if ($MPI)
!  !--append coordinate identifiers
!  write (temp, '(".c",i0)') coord
!  fCS1plan = trim (fCS1plan) // temp
!  fCS2plan = trim (fCS2plan) // temp
!  fCS4plan = trim (fCS4plan) // temp
!  fVISCplan = trim (fVISCplan) // temp
!  fDISSplan = trim (fDISSplan) // temp
!  fCS1Vplan = trim (fCS1Vplan) // temp
!  fCS2Vplan = trim (fCS2Vplan) // temp
!  fCS4Vplan = trim (fCS4Vplan) // temp
!$endif

!if(time_spec.gt.0)then
!open(1,file=path//'obs.pt')
!read(1,*)n_obs
!allocate(obs_pt(1:2,n_obs))
!do i=1,n_obs
!read(1,*)obs_pt(1:2,i)
!enddo
!close(1)
!endif

end subroutine openfiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop
use param,only:output,dt,c_count,S_FLAG,SCAL_init,jt,jan_diurnal_run
use sim_param,only:path,u,v,w,dudz,dudx,p,&
     RHSx,RHSy,RHSz,theta, txx, txy, txz, tyy, tyz, tzz
use sgsmodule,only:Cs_opt2
use scalars_module,only:sgs_t3
implicit none
real(kind=rprec),dimension(nz)::u_ndim
character(len=20)::req

! SKS
character (100) :: fname, temp  ! With 64 the code was giving an error !
! character (64) :: fname, temp ! sort of segmentation fault i guess.
! SKS

integer::jx,jy,jz
integer:: fields_3d_write_freqz_temp

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif

real(kind=rprec),dimension(ld,$lbz:nz,num_vars*inst_array_lim),save::inst_slice_array
integer,save:: inst_slice_counter

!jt_total=jt_total+1  !--moved into main
!call calculate_mean

if ((use_avgslice) .and. (mod(jt,c_count)==0)) then
       if ((S_FLAG) .and. (jt.GE.SCAL_init)) then ! Output scalar variables
         call MM_XYZ_Out()
       elseif (.not. S_FLAG) then
         !call avgslice()
       !  call MM_budget_slice()
       end if
end if

if (output) then
  if ((mod(jt_total,base)==0).and.(jt_total.ge.1)) then
    if (S_FLAG) then
   !    write (fname, '(a,i6.6,a)') path // 'output/vel_sc', jt_total, '.out'
    else
       write (fname, '(a,i6.6,a)') path // 'output/vel', jt_total, '.out'
    end if
    $if ($MPI)
       write (temp, '(".c",i0)') coord
       fname = trim (fname) // temp
    $endif

    open(1,file=fname,form='unformatted')
    call checkpoint (1)
    close(1)
 !   if (io_mean) call io_mean_out
  end if
end if

 if (S_FLAG) then ! If Block 1
      if ((inst_slice_flag) .AND. mod(jt_total, inst_slice_freqz)==0) then !If Block 2
        if (jt .eq. inst_slice_freqz) inst_slice_counter=1
        if (jt .eq. inst_slice_freqz) print *,'inst_slice_counter = ',inst_slice_counter
            
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
                print *,'inst_slice_counter = ',inst_slice_counter
            end if

       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+1) = u(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+2) = v(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+3) = w(:,slice_inst,:)
       inst_slice_array(:,:,(inst_slice_counter-1)*num_vars+4) = theta(:,slice_inst,:)

         if (mod(inst_slice_counter,inst_array_lim) .eq. 0) then !If Block 3 begins
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'INSIDE:inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
          write(fname,'(A,i6.6,A)')path//'output/fields_3d/inst_slices_uvwT_till_',jt_total,'.bin'
          $if ($MPI)
            write (temp, '(".c",i0)') coord
            fname = trim (fname) // temp
          $endif
           open(1,file=fname,form='unformatted')
           write(1) real(inst_slice_array)
           close(1)
           inst_slice_array=0._rprec; inst_slice_counter=0;
         end if ! If Block 3 ends

       inst_slice_counter = inst_slice_counter+1 !increment the slice counter by 1
            if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then !If Block 4
                print *,'inst_slice_counter = ',inst_slice_counter
            end if !If Block 4 ends
     end if ! If Block 2 ends
       
       fields_3d_write_freqz_temp=50

    if ((output_fields_3d_flag) .and. mod(jt_total,fields_3d_write_freqz_temp)==0) then !If Block 5 begins

    write(fname,'(A,i6.6,A)')path//'output/fields_3d/o_uvwT_',jt_total,'.bin'
    $if ($MPI)
      write (temp, '(".c",i0)') coord
      fname = trim (fname) // temp
    $endif
    open(1,file=fname,form='unformatted')
    write(1) real(u),real(v),real(w),real(theta); close(1)
    
    end if ! If Block 5 ends
 end if ! If Block 1 ends

!  if ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then 
!   if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
!  end if

!  if (time_spec.gt.0) call timeseries_spec

end subroutine output_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_avg_var(avg_var,raw_var,output_dim)
use param,only:nx,ny,nz,c_count,p_count
integer :: output_dim
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var
real(kind=rprec),dimension(1:nx,1:ny,1:nz-1):: raw_var
real(kind=rprec):: fr
     
fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
select case (output_dim)
  case(1) !average over y and t
      avg_var=avg_var+fr*sum(raw_var(1:nx,1:ny,1:nz-1),dim=2)/ny
  case(2) ! average over x,y and t
   !   avg_var(:,1)=avg_var(:,1)+fr*sum(sum(raw_var(1:nx,1:ny,1:nz-1),dim=1),dim=2)/(nx*ny)
end select
end subroutine compute_avg_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

  avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,                &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
      if (average_dim_num .eq. 1) then
         do ind2=1,nz_tot-1
          write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
         end do
      else if (average_dim_num .eq. 2) then
         write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
      end if
         close(file_ind)
  end if

5168     format(1400(E14.5))
end subroutine collocate_MPI_averages

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_N(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(dim1_size,dim2_size)::avg_var_proc
real(kind=rprec),dimension(dim1_global,dim2_global)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (avg_var_proc(1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1, 1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
        if (average_dim_num .eq. 1) then
           do ind2=1,nz_tot-1
            write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,ind2),ind1=1,nx)
           end do
        else if (average_dim_num .eq. 2) then
           write(file_ind,5168) jt*dt,(avg_var_tot_domain(ind1,1),ind1=1,nz_tot-1)
        end if
        close(file_ind)
  end if
5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!xxxxxxxxxx----------VIJ----------XXXXXXXXXXXXXXXXXXXXX
!!!xxxxxxxxxx----------MM----------XXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--This Subroutine Give the output for whole domain (MM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MM_XYZ_Out
use sim_param
use param,only:dz,p_count,c_count,jt,jt_total,T_scale,q_scale,u_star,z_i
use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
use scalars_module,only: L,wstar,dTdz,sgs_t1,sgs_t2,sgs_t3,dTdx,dTdy
use scalars_module_q,only: dqdz,sgs_q2,sgs_q1,sgs_q3,dqdx,dqdy
use atm_thermodynamics
implicit none
integer::i,j,k


! ----------------  thermodynamics ---------------
real(kind=rprec),dimension(nx,ny,nz-1),save::aactual_T,aactual_Tv,apr_atm,arel_hum,arel_hum_q
real(kind=rprec),dimension(nx,ny,nz-1),save::avapor_pr,asat_vapor_pr,asat_qmix,azlcl_all
! -----------------------------
! First, second, and third individual moments
real(kind=rprec),dimension(nx,ny,nz-1),save::ap,au,av,aw,atheta,aq,p2,u2,v2,w2,theta2,q2,u3,v3,w3,q3,T3
!real(kind=rprec),dimension(nx,ny,nz-1),save::apf,pf2
! Derivatives
real(kind=rprec),dimension(nx,ny,nz-1),save::adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,adTdx,adTdy,adTdz
real(kind=rprec),dimension(nx,ny,nz-1),save::adqdx,adqdy,adqdz,adpdx,adpdy,adpdz
! 3rd moments and fluxes
real(kind=rprec),dimension(nx,ny,nz-1),save::tke,awe,aue,ave,auw,avw,awt,awq,aut,avt,auq,avq,auv
real(kind=rprec),dimension(nx,ny,nz-1),save::awuu,awvv,awqq,awtt,awwu,awwv,awwt,awwq
!real(kind=rprec),dimension(nx,ny,nz-1),save::pi1_t,pi2_t,pi3_t,pi1_q,pi2_q,pi3_q
!SGS
real(kind=rprec),dimension(nx,ny,nz-1),save::atxx,atxy,atxz,atyy,atyz,atzz
!real(kind=rprec),dimension(nx,ny,nz-1),save::atpx,atpy,atpz,aqpx,aqpy,aqpz
real(kind=rprec),dimension(nx,ny,nz-1),save::aup,avp,awp,atp,aqp,atq,aupx,aupy,aupz,avpx,avpy,avpz
real(kind=rprec),dimension(nx,ny,nz-1),save::awpx,awpy,awpz,at11ux,at12uy,at13uz
real(kind=rprec),dimension(nx,ny,nz-1),save::at21vx,at22vy,at23vz,at31wx,at32wy,at33wz
real(kind=rprec),dimension(nx,ny,nz-1),save::apux,apuy,apuz,apvx,apvy,apvz,apwx,apwy,apwz,aptx,apty,aptz,apqx,apqy,apqz
real(kind=rprec),dimension(nx,ny,nz-1),save::at11s11,at12s12,at13s13,at22s22,at23s23,at33s33
!real(kind=rprec),dimension(nx,ny,nz-1),save::api1_tx,api2_ty,api3_tz,api1_qx,api2_qy,api3_qz
real(kind=rprec),dimension(nx,ny,nz-1),save::as11,as12,as13,as22,as23,as33
!real(kind=rprec),dimension(nx,ny,nz-1),save::asgst1,asgst2,asgst3,asgsq1,asgsq2,asgsq3
real(kind=rprec),dimension(nx,ny,nz-1),save::asgst3,asgsq3
real(kind=rprec),dimension(nx,ny,nz-1),save::aut11,aut12,aut13,avt21,avt22,avt23,awt31,awt32,awt33
!real(kind=rprec),dimension(nx,ny,nz-1),save::awt1,awt2,awt3,awt4,awq1,awq2,awq3,awq4,auw1,auw2,auw3,auw4
!real(kind=rprec),dimension(nx,ny,nz-1),save::avw1,avw2,avw3,avw4
real(kind=rprec),allocatable,dimension(:,:,:)::avg_out
real(kind=rprec)::fr
real(kind=rprec),dimension(nx,ny)::arg1,arg2,arg3,arg4,arg4f,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13
real(kind=rprec),dimension(nx,ny)::ux,uy,uz,vx,vy,vz,wx,wy,wz,px,py,pz,tx,ty,tz,qx,qy,qz
real(kind=rprec),dimension(nx,ny)::S11,S22,S33,S12,S13,S23
real(rprec),parameter::delbar=2._rprec
real(rprec)::denomsig=0.5_rprec/((delbar**(2.0_rprec/3.0_rprec))-1._rprec)
real(rprec)::denomsigq=0.5_rprec/(((2._rprec*delbar)**(2.0_rprec/3.0_rprec))-1._rprec)

character (len=256) :: local_filename
!

!if(.not.allocated(aactual_T)) allocate(aactual_T(nx,ny,nz-1))

if(jt_total .LE. c_count) then
         au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;
         !apf=0._rprec;atheta=0._rprec;
         atheta=0._rprec;
        aq=0._rprec; p2=0._rprec; u2=0._rprec;v2=0._rprec;
       ! pf2=0._rprec;
        w2=0._rprec;theta2=0._rprec; q2=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
        q3=0._rprec; T3=0._rprec;adudx=0._rprec;adudy=0._rprec;adudz=0._rprec;
        advdx=0._rprec;advdy=0._rprec;advdz=0._rprec;adwdx=0._rprec;adwdy=0._rprec;adwdz=0._rprec;
        adTdx=0._rprec;adTdy=0._rprec;adTdz=0._rprec;adqdx=0._rprec;adqdy=0._rprec;adqdz=0._rprec;
         adpdx=0._rprec;adpdy=0._rprec;adpdz=0._rprec;tke=0._rprec;awe=0._rprec;aue=0._rprec;ave=0._rprec;
         auw=0._rprec;avw=0._rprec;awt=0._rprec;awq=0._rprec;aut=0._rprec;avt=0._rprec;auq=0._rprec;
         avq=0._rprec;auv=0._rprec;awuu=0._rprec;awvv=0._rprec;awqq=0._rprec;awtt=0._rprec;awwu=0._rprec;
         awwv=0._rprec;awwt=0._rprec;awwq=0._rprec;
         !pi1_t=0._rprec;pi2_t=0._rprec;pi3_t=0._rprec;
         !pi1_q=0._rprec;pi2_q=0._rprec;pi3_q=0._rprec;
         atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
        atzz=0._rprec;aup=0._rprec;avp=0._rprec;awp=0._rprec;atp=0._rprec;aqp=0._rprec;atq=0._rprec;
        aupx=0._rprec;aupy=0._rprec;aupz=0._rprec;avpx=0._rprec;avpy=0._rprec;avpz=0._rprec;
        awpx=0._rprec;awpy=0._rprec;awpz=0._rprec;
        !atpx=0._rprec;atpy=0._rprec;atpz=0._rprec;
        !aqpx=0._rprec;aqpy=0._rprec;aqpz=0._rprec;
        at11ux=0._rprec;at12uy=0._rprec;at13uz=0._rprec;
        at21vx=0._rprec;at22vy=0._rprec;at23vz=0._rprec;at31wx=0._rprec;at32wy=0._rprec;at33wz=0._rprec;
        apux=0._rprec;apuy=0._rprec;apuz=0._rprec;apvx=0._rprec;apvy=0._rprec;apvz=0._rprec;
        apwx=0._rprec;apwy=0._rprec;apwz=0._rprec;
        aptx=0._rprec;apty=0._rprec;aptz=0._rprec;
        apqx=0._rprec;apqy=0._rprec;apqz=0._rprec;
        at11s11=0._rprec;at12s12=0._rprec;at13s13=0._rprec;at22s22=0._rprec;at23s23=0._rprec;at33s33=0._rprec;
        !api1_tx=0._rprec;api2_ty=0._rprec;api3_tz=0._rprec;api1_qx=0._rprec;api2_qy=0._rprec;api3_qz=0._rprec;
        as11=0._rprec;as12=0._rprec;as13=0._rprec;as22=0._rprec;as23=0._rprec;as33=0._rprec;
        aut11=0._rprec;aut12=0._rprec;aut13=0._rprec;avt21=0._rprec;avt22=0._rprec;avt23=0._rprec;
        awt31=0._rprec;awt32=0._rprec;awt33=0._rprec;
        !asgst1=0._rprec;asgst2=0._rprec;asgst3=0._rprec;asgsq1=0._rprec;asgsq2=0._rprec;asgsq3=0._rprec;
        asgst3=0._rprec;asgsq3=0._rprec;
        arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg4f=0._rprec;arg5=0._rprec;arg6=0._rprec;
        arg7=0._rprec;arg8=0._rprec;arg9=0._rprec;arg10=0._rprec;arg11=0._rprec;
        arg12=0._rprec;ux=0._rprec;uy=0._rprec;uz=0._rprec;vx=0._rprec;vy=0._rprec;
        !awt1=0._rprec;awt2=0._rprec;awt3=0._rprec;awt4=0._rprec;awq1=0._rprec;awq2=0._rprec;
        !awq3=0._rprec;awq4=0._rprec;auw1=0._rprec;auw2=0._rprec;auw3=0._rprec;auw4=0._rprec;
        !avw1=0._rprec;avw2=0._rprec;avw3=0._rprec;avw4=0._rprec;
        aactual_T=0._rprec;aactual_Tv=0._rprec;apr_atm=0._rprec;arel_hum=0._rprec;arel_hum_q=0._rprec;
        avapor_pr=0._rprec;asat_vapor_pr=0._rprec;asat_qmix=0._rprec;azlcl_all=0._rprec;
end if



fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
do k=1,Nz-1
!do i=1,Nx
!do j=1,Ny



       if((k .eq. 1) .AND. ((.not. USE_MPI) .or. (USE_MPI .and. coord .eq. 0))) then
         arg1(1:nx,1:ny)=u(1:nx,1:ny,k)
         arg2(1:nx,1:ny)=v(1:nx,1:ny,k)
         arg3(1:nx,1:ny)=0.25_rprec*w(1:nx,1:ny,k+1)
         arg4(1:nx,1:ny)=p(1:nx,1:ny,k)-0.5_rprec*(arg1(1:nx,1:ny)*arg1(1:nx,1:ny)+&
                           arg2(1:nx,1:ny)*arg2(1:nx,1:ny)+arg3(1:nx,1:ny)*arg3(1:nx,1:ny))
                           
        ! arg4f(1:nx,1:ny)= (p(1:nx,1:ny,k)- 0.5_rprec*(arg1*arg1+arg2*arg2+arg3*arg3)-&
        !         (1._rprec/3._rprec)*(denomsig)*(L11t(1:nx,1:ny,k)+L22t(1:nx,1:ny,k)+L33t(1:nx,1:ny,k)))
                 
                 
         arg5(1:nx,1:ny)=(txz(1:nx,1:ny,k)+txz(1:nx,1:ny,k+1))/2.0_rprec
         arg6(1:nx,1:ny)=(tyz(1:nx,1:ny,k)+tyz(1:nx,1:ny,k+1))/2.0_rprec
         arg7(1:nx,1:ny)=theta(1:nx,1:ny,k)
         arg8(1:nx,1:ny)=qmix(1:nx,1:ny,k)

         arg9(1:nx,1:ny)=(u(1:nx,1:ny,k)-sum(u(1:nx,1:ny,k))/(nx*ny))
         arg13(1:nx,1:ny)=(v(1:nx,1:ny,k)-sum(v(1:nx,1:ny,k))/(nx*ny))
         arg10(1:nx,1:ny)=0.25_rprec*(w(1:nx,1:ny,k+1)-sum(w(1:nx,1:ny,k+1))/(nx*ny))
         arg11(1:nx,1:ny)=(theta(1:nx,1:ny,k)-sum(theta(1:nx,1:ny,k))/(nx*ny))
         arg12(1:nx,1:ny)=(qmix(1:nx,1:ny,k)-sum(qmix(1:nx,1:ny,k))/(nx*ny))

         ux(1:nx,1:ny)=dudx(1:nx,1:ny,k)
         uy(1:nx,1:ny)=dudy(1:nx,1:ny,k)
         uz(1:nx,1:ny)=dudz(1:nx,1:ny,k)

         vx(1:nx,1:ny)=dvdx(1:nx,1:ny,k)
         vy(1:nx,1:ny)=dvdy(1:nx,1:ny,k)
         vz(1:nx,1:ny)=dvdz(1:nx,1:ny,k)

         wx(1:nx,1:ny)=((dwdx(1:nx,1:ny,k+1))/2.0_rprec)
         wy(1:nx,1:ny)=((dwdy(1:nx,1:ny,k+1))/2.0_rprec)
         wz(1:nx,1:ny)=dwdz(1:nx,1:ny,k)

         px(1:nx,1:ny)=dpdx(1:nx,1:ny,k) - arg1(1:nx,1:ny)*ux(1:nx,1:ny) - &
                        arg2(1:nx,1:ny)*vx(1:nx,1:ny) - arg3(1:nx,1:ny)*wx(1:nx,1:ny)
                        
         py(1:nx,1:ny)=dpdy(1:nx,1:ny,k) - arg1(1:nx,1:ny)*uy(1:nx,1:ny) - &
                        arg2(1:nx,1:ny)*vy(1:nx,1:ny) - arg3(1:nx,1:ny)*wy(1:nx,1:ny)
                        
         pz(1:nx,1:ny)=dpdz(1:nx,1:ny,k) - arg1(1:nx,1:ny)*uz(1:nx,1:ny) - &
                         arg2(1:nx,1:ny)*vz(1:nx,1:ny) - arg3(1:nx,1:ny)*wz(1:nx,1:ny)

         tx(1:nx,1:ny)=dTdx(1:nx,1:ny,k)
         ty(1:nx,1:ny)=dTdy(1:nx,1:ny,k)
         tz(1:nx,1:ny)=dTdz(1:nx,1:ny,k)

         qx(1:nx,1:ny)=dqdx(1:nx,1:ny,k)
         qy(1:nx,1:ny)=dqdy(1:nx,1:ny,k)
         qz(1:nx,1:ny)=dqdz(1:nx,1:ny,k)

         else

         arg1(1:nx,1:ny)=(u(1:nx,1:ny,k)+u(1:nx,1:ny,k-1))/2.0_rprec
         arg2(1:nx,1:ny)=(v(1:nx,1:ny,k)+v(1:nx,1:ny,k-1))/2.0_rprec
         arg3(1:nx,1:ny)=w(1:nx,1:ny,k)
         arg4(1:nx,1:ny)=((p(1:nx,1:ny,k)+p(1:nx,1:ny,k-1))/2.0_rprec - &
                 0.5_rprec*(arg1(1:nx,1:ny)*arg1(1:nx,1:ny)+ &
                     arg2(1:nx,1:ny)*arg2(1:nx,1:ny)+arg3(1:nx,1:ny)*arg3(1:nx,1:ny)))
    !     arg4f=((p(1:nx,1:ny,k)+p(1:nx,1:ny,k-1))/2.0_rprec - &
     !             0.5_rprec*(arg1*arg1+arg2*arg2+arg3*arg3)- &
      !           (1._rprec/3._rprec)*(denomsig)* &
      !           (L11t(1:nx,1:ny,k)+L22t(1:nx,1:ny,k)+L33t(1:nx,1:ny,k)))
         arg5(1:nx,1:ny)=txz(1:nx,1:ny,k)
         arg6(1:nx,1:ny)=tyz(1:nx,1:ny,k)
         arg7(1:nx,1:ny)=(theta(1:nx,1:ny,k)+theta(1:nx,1:ny,k-1))/2.0_rprec
         arg8(1:nx,1:ny)=(qmix(1:nx,1:ny,k)+qmix(1:nx,1:ny,k-1))/2.0_rprec

         arg9(1:nx,1:ny)=0.5_rprec*(u(1:nx,1:ny,k)-sum(u(1:nx,1:ny,k))/(nx*ny)+&
                            u(1:nx,1:ny,k-1)-sum(u(1:nx,1:ny,k-1))/(nx*ny))
                            
         arg13(1:nx,1:ny)=0.5_rprec*(v(1:nx,1:ny,k)-sum(v(1:nx,1:ny,k))/(nx*ny)+&
                            v(1:nx,1:ny,k-1)-sum(v(1:nx,1:ny,k-1))/(nx*ny))
                            
         arg10(1:nx,1:ny)=(w(1:nx,1:ny,k)-sum(w(1:nx,1:ny,k))/(nx*ny))
         
         arg11(1:nx,1:ny)=0.5_rprec*(theta(1:nx,1:ny,k)-sum(theta(1:nx,1:ny,k))/(nx*ny)+&
                          theta(1:nx,1:ny,k-1)- sum(theta(1:nx,1:ny,k-1))/(nx*ny))
                          
         arg12(1:nx,1:ny)=0.5_rprec*(qmix(1:nx,1:ny,k)-sum(qmix(1:nx,1:ny,k))/(nx*ny)+&
                          qmix(1:nx,1:ny,k-1)- sum(qmix(1:nx,1:ny,k-1))/(nx*ny))


         ux(1:nx,1:ny)=((dudx(1:nx,1:ny,k)+dudx(1:nx,1:ny,k-1))/2._rprec)
         uy(1:nx,1:ny)=((dudy(1:nx,1:ny,k)+dudy(1:nx,1:ny,k-1))/2._rprec)
         uz(1:nx,1:ny)=dudz(1:nx,1:ny,k)

         vx(1:nx,1:ny)=((dvdx(1:nx,1:ny,k)+dvdx(1:nx,1:ny,k-1))/2._rprec)
         vy(1:nx,1:ny)=((dvdy(1:nx,1:ny,k)+dvdy(1:nx,1:ny,k-1))/2._rprec)
         vz(1:nx,1:ny)=dvdz(1:nx,1:ny,k)

         wx(1:nx,1:ny)=dwdx(1:nx,1:ny,k)
         wy(1:nx,1:ny)=dwdy(1:nx,1:ny,k)
         wz(1:nx,1:ny)=((dwdz(1:nx,1:ny,k)+dwdz(1:nx,1:ny,k-1))/2._rprec)

         px(1:nx,1:ny)=((dpdx(1:nx,1:ny,k)+dpdx(1:nx,1:ny,k-1))/2._rprec) - &
                         arg1(1:nx,1:ny)*ux(1:nx,1:ny) - arg2(1:nx,1:ny)*vx(1:nx,1:ny) -&
                          arg3(1:nx,1:ny)*wx(1:nx,1:ny)
                          
         py(1:nx,1:ny)=((dpdy(1:nx,1:ny,k)+dpdy(1:nx,1:ny,k-1))/2._rprec) - &
                         arg1(1:nx,1:ny)*uy(1:nx,1:ny) - arg2(1:nx,1:ny)*vy(1:nx,1:ny) -&
                          arg3(1:nx,1:ny)*wy(1:nx,1:ny)
                          
         pz(1:nx,1:ny)=dpdz(1:nx,1:ny,k) - arg1(1:nx,1:ny)*uz(1:nx,1:ny) - &
                        arg2(1:nx,1:ny)*vz(1:nx,1:ny) - arg3(1:nx,1:ny)*wz(1:nx,1:ny)

         tx(1:nx,1:ny)=((dTdx(1:nx,1:ny,k)+dTdx(1:nx,1:ny,k-1))/2._rprec)
         ty(1:nx,1:ny)=((dTdy(1:nx,1:ny,k)+dTdy(1:nx,1:ny,k-1))/2._rprec)
         tz(1:nx,1:ny)=dTdz(1:nx,1:ny,k)

         qx(1:nx,1:ny)=((dqdx(1:nx,1:ny,k)+dqdx(1:nx,1:ny,k-1))/2._rprec)
         qy(1:nx,1:ny)=((dqdy(1:nx,1:ny,k)+dqdy(1:nx,1:ny,k-1))/2._rprec)
         qz(1:nx,1:ny)=dqdz(1:nx,1:ny,k)

       end if


! -----   Thermodynamics -----------

aactual_T(1:nx,1:ny,k)=aactual_T(1:nx,1:ny,k)+fr*actual_T(1:nx,1:ny,k)*T_scale
aactual_Tv(1:nx,1:ny,k)=aactual_Tv(1:nx,1:ny,k)+fr*actual_Tv(1:nx,1:ny,k)*T_scale
apr_atm(1:nx,1:ny,k)=apr_atm(1:nx,1:ny,k)+fr*pr_atm(1:nx,1:ny,k)
arel_hum(1:nx,1:ny,k)=arel_hum(1:nx,1:ny,k)+fr*rel_hum(1:nx,1:ny,k)
arel_hum_q(1:nx,1:ny,k)=arel_hum_q(1:nx,1:ny,k)+fr*rel_hum_q(1:nx,1:ny,k)

avapor_pr(1:nx,1:ny,k)=avapor_pr(1:nx,1:ny,k)+fr*vapor_pr(1:nx,1:ny,k)
asat_vapor_pr(1:nx,1:ny,k)=asat_vapor_pr(1:nx,1:ny,k)+fr*sat_vapor_pr(1:nx,1:ny,k)
asat_qmix(1:nx,1:ny,k)=asat_qmix(1:nx,1:ny,k)+fr*sat_qmix(1:nx,1:ny,k)
azlcl_all(1:nx,1:ny,k)=azlcl_all(1:nx,1:ny,k)+fr*zlcl_all(1:nx,1:ny,k)

! First moments

       au(1:nx,1:ny,k)=au(1:nx,1:ny,k)+fr*arg1(1:nx,1:ny)
       av(1:nx,1:ny,k)=av(1:nx,1:ny,k)+fr*arg2(1:nx,1:ny)
       aw(1:nx,1:ny,k)=aw(1:nx,1:ny,k)+fr*arg3(1:nx,1:ny)
       ap(1:nx,1:ny,k)=ap(1:nx,1:ny,k)+fr*arg4(1:nx,1:ny)
  !     apf(1:nx,1:ny,k)=apf(1:nx,1:ny,k)+fr*arg4f
       atheta(1:nx,1:ny,k)=atheta(1:nx,1:ny,k)+fr*arg7(1:nx,1:ny)
       aq(1:nx,1:ny,k)=aq(1:nx,1:ny,k)+fr*arg8(1:nx,1:ny)

! Variances and skewnesses

       u2(1:nx,1:ny,k)=u2(1:nx,1:ny,k)+fr*arg1(1:nx,1:ny)*arg1(1:nx,1:ny)
       v2(1:nx,1:ny,k)=v2(1:nx,1:ny,k)+fr*arg2(1:nx,1:ny)*arg2(1:nx,1:ny)
       w2(1:nx,1:ny,k)=w2(1:nx,1:ny,k)+fr*arg3(1:nx,1:ny)*arg3(1:nx,1:ny)
       p2(1:nx,1:ny,k)=p2(1:nx,1:ny,k)+fr*arg4(1:nx,1:ny)*arg4(1:nx,1:ny)
   !    pf2(1:nx,1:ny,k)=pf2(1:nx,1:ny,k)+fr*arg4f*arg4f
       theta2(1:nx,1:ny,k)=theta2(1:nx,1:ny,k)+fr*arg7(1:nx,1:ny)*arg7(1:nx,1:ny)
       q2(1:nx,1:ny,k)=q2(1:nx,1:ny,k)+fr*arg8(1:nx,1:ny)*arg8(1:nx,1:ny)
       tke(1:nx,1:ny,k)=tke(1:nx,1:ny,k)+fr*(0.5_rprec*(arg1(1:nx,1:ny)*arg1(1:nx,1:ny)+&
                          arg2(1:nx,1:ny)*arg2(1:nx,1:ny)+arg3(1:nx,1:ny)*arg3(1:nx,1:ny)))

       u3(1:nx,1:ny,k)=u3(1:nx,1:ny,k)+fr*arg1(:,:)*arg1(:,:)*arg1(:,:)
       v3(1:nx,1:ny,k)=v3(1:nx,1:ny,k)+fr*arg2(:,:)*arg2(:,:)*arg2(:,:)
       w3(1:nx,1:ny,k)=w3(1:nx,1:ny,k)+fr*arg3(:,:)*arg3(:,:)*arg3(:,:)
       T3(1:nx,1:ny,k)=T3(1:nx,1:ny,k)+fr*arg7(:,:)*arg7(:,:)*arg7(:,:)
       q3(1:nx,1:ny,k)=q3(1:nx,1:ny,k)+fr*arg8(:,:)*arg8(:,:)*arg8(:,:)

! Fluxes and third-order moments

      awe(1:nx,1:ny,k)=awe(1:nx,1:ny,k)+fr*arg3(:,:)*(0.5_rprec*(arg1(:,:)*arg1(:,:)+&
                         arg2(:,:)*arg2(:,:)+arg3(:,:)*arg3(:,:)))
      ave(1:nx,1:ny,k)=ave(1:nx,1:ny,k)+fr*arg2(:,:)*(0.5_rprec*(arg1(:,:)*arg1(:,:)+&  
                        arg2(:,:)*arg2(:,:)+arg3(:,:)*arg3(:,:)))
      aue(1:nx,1:ny,k)=aue(1:nx,1:ny,k)+fr*arg1(:,:)*(0.5_rprec*(arg1(:,:)*arg1(:,:)+&
                        arg2(:,:)*arg2(:,:)+arg3(:,:)*arg3(:,:)))
      auw(1:nx,1:ny,k)=auw(1:nx,1:ny,k)+fr*arg1(:,:)*arg3(:,:)
      avw(1:nx,1:ny,k)=avw(1:nx,1:ny,k)+fr*arg2(:,:)*arg3(:,:)
      awt(1:nx,1:ny,k)=awt(1:nx,1:ny,k)+fr*arg3(:,:)*arg7(:,:)
      awq(1:nx,1:ny,k)=awq(1:nx,1:ny,k)+fr*arg3(:,:)*arg8(:,:)
      auv(1:nx,1:ny,k)=auv(1:nx,1:ny,k)+fr*arg1(:,:)*arg2(:,:)
      aut(1:nx,1:ny,k)=aut(1:nx,1:ny,k)+fr*arg1(:,:)*arg7(:,:)
      avt(1:nx,1:ny,k)=avt(1:nx,1:ny,k)+fr*arg2(:,:)*arg7(:,:)
      auq(1:nx,1:ny,k)=auq(1:nx,1:ny,k)+fr*arg1(:,:)*arg8(:,:)
      avq(1:nx,1:ny,k)=avq(1:nx,1:ny,k)+fr*arg2(:,:)*arg8(:,:)

      awuu(1:nx,1:ny,k)=awuu(1:nx,1:ny,k)+fr*arg3(:,:)*arg1(:,:)*arg1(:,:)
      awvv(1:nx,1:ny,k)=awvv(1:nx,1:ny,k)+fr*arg3(:,:)*arg2(:,:)*arg2(:,:)
      awtt(1:nx,1:ny,k)=awtt(1:nx,1:ny,k)+fr*arg3(:,:)*arg7(:,:)*arg7(:,:)
      awqq(1:nx,1:ny,k)=awqq(1:nx,1:ny,k)+fr*arg3(:,:)*arg8(:,:)*arg8(:,:)
      awwu(1:nx,1:ny,k)=awwu(1:nx,1:ny,k)+fr*arg3(:,:)*arg3(:,:)*arg1(:,:)
      awwv(1:nx,1:ny,k)=awwv(1:nx,1:ny,k)+fr*arg3(:,:)*arg3(:,:)*arg2(:,:)
      awwt(1:nx,1:ny,k)=awwt(1:nx,1:ny,k)+fr*arg3(:,:)*arg3(:,:)*arg7(:,:)
      awwq(1:nx,1:ny,k)=awwq(1:nx,1:ny,k)+fr*arg3(:,:)*arg3(:,:)*arg8(:,:)

      !pi1_t(1:nx,1:ny,k)=pi1_t(1:nx,1:ny,k)+fr*(sgs_t1(1:nx,1:ny,k))*arg7
      !pi2_t(1:nx,1:ny,k)=pi2_t(1:nx,1:ny,k)+fr*(sgs_t2(1:nx,1:ny,k))*arg7
      !pi3_t(1:nx,1:ny,k)=pi3_t(1:nx,1:ny,k)+fr*(sgs_t3(1:nx,1:ny,k))*arg7


      !pi1_q(1:nx,1:ny,k)=pi1_q(1:nx,1:ny,k)+fr*(sgs_q1(1:nx,1:ny,k))*arg8
      !pi2_q(1:nx,1:ny,k)=pi2_q(1:nx,1:ny,k)+fr*(sgs_q2(1:nx,1:ny,k))*arg8
      !pi3_q(1:nx,1:ny,k)=pi3_q(1:nx,1:ny,k)+fr*(sgs_q3(1:nx,1:ny,k))*arg8

      !asgst1(1:nx,1:ny,k)=asgst1(1:nx,1:ny,k)+fr*(sgs_t1(1:nx,1:ny,k))
      !asgst2(1:nx,1:ny,k)=asgst2(1:nx,1:ny,k)+fr*(sgs_t2(1:nx,1:ny,k))
      asgst3(1:nx,1:ny,k)=asgst3(1:nx,1:ny,k)+fr*(sgs_t3(1:nx,1:ny,k))

      !asgsq1(1:nx,1:ny,k)=asgsq1(1:nx,1:ny,k)+fr*(sgs_q1(1:nx,1:ny,k))
      !asgsq2(1:nx,1:ny,k)=asgsq2(1:nx,1:ny,k)+fr*(sgs_q2(1:nx,1:ny,k))
      asgsq3(1:nx,1:ny,k)=asgsq3(1:nx,1:ny,k)+fr*(sgs_q3(1:nx,1:ny,k))


! -----------------Derivatives -----------

     adudx(1:nx,1:ny,k)=adudx(1:nx,1:ny,k)+fr*ux(:,:)
     adudy(1:nx,1:ny,k)=adudy(1:nx,1:ny,k)+fr*uy(:,:)
     adudz(1:nx,1:ny,k)=adudz(1:nx,1:ny,k)+fr*uz(:,:)

     advdx(1:nx,1:ny,k)=advdx(1:nx,1:ny,k)+fr*vx(:,:)
     advdy(1:nx,1:ny,k)=advdy(1:nx,1:ny,k)+fr*vy(:,:)
     advdz(1:nx,1:ny,k)=advdz(1:nx,1:ny,k)+fr*vz(:,:)

     adwdx(1:nx,1:ny,k)=adwdx(1:nx,1:ny,k)+fr*wx(:,:)
     adwdy(1:nx,1:ny,k)=adwdy(1:nx,1:ny,k)+fr*wy(:,:)
     adwdz(1:nx,1:ny,k)=adwdz(1:nx,1:ny,k)+fr*wz(:,:)

     adTdx(1:nx,1:ny,k)=adTdx(1:nx,1:ny,k)+fr*tx(:,:)
     adTdy(1:nx,1:ny,k)=adTdy(1:nx,1:ny,k)+fr*ty(:,:)
     adTdz(1:nx,1:ny,k)=adTdz(1:nx,1:ny,k)+fr*tz(:,:)

     adqdx(1:nx,1:ny,k)=adqdx(1:nx,1:ny,k)+fr*qx(:,:)
     adqdy(1:nx,1:ny,k)=adqdy(1:nx,1:ny,k)+fr*qy(:,:)
     adqdz(1:nx,1:ny,k)=adqdz(1:nx,1:ny,k)+fr*qz(:,:)

     adpdx(1:nx,1:ny,k)=adpdx(1:nx,1:ny,k)+fr*px(:,:)
     adpdy(1:nx,1:ny,k)=adpdy(1:nx,1:ny,k)+fr*py(:,:)
     adpdz(1:nx,1:ny,k)=adpdz(1:nx,1:ny,k)+fr*pz(:,:)


!------------ SGS -----------------------------

      atxx(1:nx,1:ny,k)=atxx(1:nx,1:ny,k)+fr*txx(1:nx,1:ny,k)
      atxy(1:nx,1:ny,k)=atxy(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)
      atxz(1:nx,1:ny,k)=atxz(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)
      atyz(1:nx,1:ny,k)=atyz(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)
      atyy(1:nx,1:ny,k)=atyy(1:nx,1:ny,k)+fr*tyy(1:nx,1:ny,k)
      atzz(1:nx,1:ny,k)=atzz(1:nx,1:ny,k)+fr*tzz(1:nx,1:ny,k)

! ------- pressure transport ---------

      aup(1:nx,1:ny,k)=aup(1:nx,1:ny,k)+fr*arg1(:,:)*arg4(:,:)
      avp(1:nx,1:ny,k)=avp(1:nx,1:ny,k)+fr*arg2(:,:)*arg4(:,:)
      awp(1:nx,1:ny,k)=awp(1:nx,1:ny,k)+fr*arg3(:,:)*arg4(:,:)
      atp(1:nx,1:ny,k)=atp(1:nx,1:ny,k)+fr*arg7(:,:)*arg4(:,:)
      aqp(1:nx,1:ny,k)=aqp(1:nx,1:ny,k)+fr*arg8(:,:)*arg4(:,:)
      atq(1:nx,1:ny,k)=atq(1:nx,1:ny,k)+fr*arg7(:,:)*arg8(:,:)
      aupx(1:nx,1:ny,k)=aupx(1:nx,1:ny,k)+fr*arg1(:,:)*px(:,:)
      aupy(1:nx,1:ny,k)=aupy(1:nx,1:ny,k)+fr*arg1(:,:)*py(:,:)
      aupz(1:nx,1:ny,k)=aupz(1:nx,1:ny,k)+fr*arg1(:,:)*pz(:,:)
      avpx(1:nx,1:ny,k)=avpx(1:nx,1:ny,k)+fr*arg2(:,:)*px(:,:)
      avpy(1:nx,1:ny,k)=avpy(1:nx,1:ny,k)+fr*arg2(:,:)*py(:,:)
      avpz(1:nx,1:ny,k)=avpz(1:nx,1:ny,k)+fr*arg2(:,:)*pz(:,:)
      awpx(1:nx,1:ny,k)=awpx(1:nx,1:ny,k)+fr*arg3(:,:)*px(:,:)
      awpy(1:nx,1:ny,k)=awpy(1:nx,1:ny,k)+fr*arg3(:,:)*py(:,:)
      awpz(1:nx,1:ny,k)=awpz(1:nx,1:ny,k)+fr*arg3(:,:)*pz(:,:)

 !-------------------------------------------
      at11ux(1:nx,1:ny,k)=at11ux(1:nx,1:ny,k)+fr*txx(1:nx,1:ny,k)*dudx(1:nx,1:ny,k)
      at12uy(1:nx,1:ny,k)=at12uy(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)*dudy(1:nx,1:ny,k)
      at13uz(1:nx,1:ny,k)=at13uz(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)*dudz(1:nx,1:ny,k)

      at21vx(1:nx,1:ny,k)=at21vx(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)*dvdx(1:nx,1:ny,k)
      at22vy(1:nx,1:ny,k)=at22vy(1:nx,1:ny,k)+fr*tyy(1:nx,1:ny,k)*dvdy(1:nx,1:ny,k)
      at23vz(1:nx,1:ny,k)=at23vz(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)*dvdz(1:nx,1:ny,k)

      at31wx(1:nx,1:ny,k)=at31wx(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)*dwdx(1:nx,1:ny,k)
      at32wy(1:nx,1:ny,k)=at32wy(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)*dwdy(1:nx,1:ny,k)
      at33wz(1:nx,1:ny,k)=at33wz(1:nx,1:ny,k)+fr*tzz(1:nx,1:ny,k)*dwdz(1:nx,1:ny,k)

      apux(1:nx,1:ny,k)=apux(1:nx,1:ny,k)+fr*arg4*ux(:,:)
      apuy(1:nx,1:ny,k)=apuy(1:nx,1:ny,k)+fr*arg4*uy(:,:)
      apuz(1:nx,1:ny,k)=apuz(1:nx,1:ny,k)+fr*arg4*uz(:,:)

      apvx(1:nx,1:ny,k)=apvx(1:nx,1:ny,k)+fr*arg4(:,:)*vx(:,:)
      apvy(1:nx,1:ny,k)=apvy(1:nx,1:ny,k)+fr*arg4(:,:)*vy(:,:)
      apvz(1:nx,1:ny,k)=apvz(1:nx,1:ny,k)+fr*arg4(:,:)*vz(:,:)

      apwx(1:nx,1:ny,k)=apwx(1:nx,1:ny,k)+fr*arg4(:,:)*wx(:,:)
      apwy(1:nx,1:ny,k)=apwy(1:nx,1:ny,k)+fr*arg4(:,:)*wy(:,:)
      apwz(1:nx,1:ny,k)=apwz(1:nx,1:ny,k)+fr*arg4(:,:)*wz(:,:)

      aptx(1:nx,1:ny,k)=aptx(1:nx,1:ny,k)+fr*arg4(:,:)*tx(:,:)
      apty(1:nx,1:ny,k)=apty(1:nx,1:ny,k)+fr*arg4(:,:)*ty(:,:)
      aptz(1:nx,1:ny,k)=aptz(1:nx,1:ny,k)+fr*arg4(:,:)*tz(:,:)

      !atpx(1:nx,1:ny,k)=atpx(1:nx,1:ny,k)+fr*arg7*px
      !atpy(1:nx,1:ny,k)=atpy(1:nx,1:ny,k)+fr*arg7*py
      !atpz(1:nx,1:ny,k)=atpz(1:nx,1:ny,k)+fr*arg7*pz

      !aqpx(1:nx,1:ny,k)=atpx(1:nx,1:ny,k)+fr*arg8*px
      !aqpy(1:nx,1:ny,k)=atpy(1:nx,1:ny,k)+fr*arg8*py
      !aqpz(1:nx,1:ny,k)=atpz(1:nx,1:ny,k)+fr*arg8*pz

      apqx(1:nx,1:ny,k)=apqx(1:nx,1:ny,k)+fr*arg4(:,:)*qx(:,:)
      apqy(1:nx,1:ny,k)=apqy(1:nx,1:ny,k)+fr*arg4(:,:)*qy(:,:)
      apqz(1:nx,1:ny,k)=apqz(1:nx,1:ny,k)+fr*arg4(:,:)*qz(:,:)


! --------------------- Strains -----------------
        S11(:,:)=ux(:,:)
        S12(:,:)=0.5_rprec*(uy(:,:)+vx(:,:))
        S13(:,:)=0.5_rprec*(uz(:,:)+wx(:,:))
        S22(:,:)=vy(:,:)
        S23(:,:)=0.5_rprec*(vz(:,:)+wy(:,:))
        S33(:,:)=wz(:,:)

        as11(1:nx,1:ny,k)=as11(1:nx,1:ny,k)+fr*S11(:,:)
        as12(1:nx,1:ny,k)=as12(1:nx,1:ny,k)+fr*S12(:,:)
        as13(1:nx,1:ny,k)=as13(1:nx,1:ny,k)+fr*S13(:,:)
        as22(1:nx,1:ny,k)=as22(1:nx,1:ny,k)+fr*S22(:,:)
        as23(1:nx,1:ny,k)=as23(1:nx,1:ny,k)+fr*S23(:,:)
        as33(1:nx,1:ny,k)=as33(1:nx,1:ny,k)+fr*S33(:,:)

        at11s11(1:nx,1:ny,k)=at11s11(1:nx,1:ny,k)+fr*txx(1:nx,1:ny,k)*S11(:,:)
        at12s12(1:nx,1:ny,k)=at12s12(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)*S12(:,:)
        at13s13(1:nx,1:ny,k)=at13s13(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)*S13(:,:)
        at22s22(1:nx,1:ny,k)=at22s22(1:nx,1:ny,k)+fr*tyy(1:nx,1:ny,k)*S22(:,:)
        at23s23(1:nx,1:ny,k)=at23s23(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)*S23(:,:)
        at33s33(1:nx,1:ny,k)=at33s33(1:nx,1:ny,k)+fr*tzz(1:nx,1:ny,k)*S33(:,:)

        aut11(1:nx,1:ny,k)=aut11(1:nx,1:ny,k)+fr*txx(1:nx,1:ny,k)*arg1(:,:)
        aut12(1:nx,1:ny,k)=aut12(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)*arg1(:,:)
        aut13(1:nx,1:ny,k)=aut13(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)*arg1(:,:)

        avt21(1:nx,1:ny,k)=avt21(1:nx,1:ny,k)+fr*txy(1:nx,1:ny,k)*arg2(:,:)
        avt22(1:nx,1:ny,k)=avt22(1:nx,1:ny,k)+fr*tyy(1:nx,1:ny,k)*arg2(:,:)
        avt23(1:nx,1:ny,k)=avt23(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)*arg2(:,:)

        awt31(1:nx,1:ny,k)=awt31(1:nx,1:ny,k)+fr*txz(1:nx,1:ny,k)*arg3(:,:)
        awt32(1:nx,1:ny,k)=awt32(1:nx,1:ny,k)+fr*tyz(1:nx,1:ny,k)*arg3(:,:)
        awt33(1:nx,1:ny,k)=awt33(1:nx,1:ny,k)+fr*tzz(1:nx,1:ny,k)*arg3(:,:)

! -------- sgs flux scalar gradient interactions ----

        !api1_tx(1:nx,1:ny,k)=api1_tx(1:nx,1:ny,k)+fr*(sgs_t1(1:nx,1:ny,k))*tx
        !api2_ty(1:nx,1:ny,k)=api2_ty(1:nx,1:ny,k)+fr*(sgs_t2(1:nx,1:ny,k))*ty
        !api3_tz(1:nx,1:ny,k)=api3_tz(1:nx,1:ny,k)+fr*(sgs_t3(1:nx,1:ny,k))*tz

        !api1_qx(1:nx,1:ny,k)=api1_qx(1:nx,1:ny,k)+fr*(sgs_q1(1:nx,1:ny,k))*qx
        !api2_qy(1:nx,1:ny,k)=api2_qy(1:nx,1:ny,k)+fr*(sgs_q2(1:nx,1:ny,k))*qy
        !api3_qz(1:nx,1:ny,k)=api3_qz(1:nx,1:ny,k)+fr*(sgs_q3(1:nx,1:ny,k))*qz



      !  if ((arg10 > 0.0_rprec) .AND. (arg9 > 0.0_rprec)) then
      !          auw1(1:nx,1:ny,k)=auw1(1:nx,1:ny,k)+fr*arg9*arg10
      !          auw2(1:nx,1:ny,k)=auw2(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw3(1:nx,1:ny,k)=auw3(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw4(1:nx,1:ny,k)=auw4(1:nx,1:ny,k)+fr*0.0_rprec
      !  else if ((arg10 > 0.0_rprec) .AND. (arg9 < 0.0_rprec)) then

      !          auw1(1:nx,1:ny,k)=auw1(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw2(1:nx,1:ny,k)=auw2(1:nx,1:ny,k)+fr*arg9*arg10
      !          auw3(1:nx,1:ny,k)=auw3(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw4(1:nx,1:ny,k)=auw4(1:nx,1:ny,k)+fr*0.0_rprec

      !   else if ((arg10 < 0.0_rprec) .AND. (arg9 < 0.0_rprec)) then

      !          auw1(1:nx,1:ny,k)=auw1(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw2(1:nx,1:ny,k)=auw2(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw3(1:nx,1:ny,k)=auw3(1:nx,1:ny,k)+fr*arg9*arg10
      !          auw4(1:nx,1:ny,k)=auw4(1:nx,1:ny,k)+fr*0.0_rprec

      !  else
      !          auw1(1:nx,1:ny,k)=auw1(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw2(1:nx,1:ny,k)=auw2(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw3(1:nx,1:ny,k)=auw3(1:nx,1:ny,k)+fr*0.0_rprec
      !          auw4(1:nx,1:ny,k)=auw4(1:nx,1:ny,k)+fr*arg9*arg10
      !  end if


      ! if ((arg10 > 0.0_rprec) .AND. (arg13 > 0.0_rprec)) then
      !          avw1(1:nx,1:ny,k)=avw1(1:nx,1:ny,k)+fr*arg13*arg10
      !          avw2(1:nx,1:ny,k)=avw2(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw3(1:nx,1:ny,k)=avw3(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw4(1:nx,1:ny,k)=avw4(1:nx,1:ny,k)+fr*0.0_rprec
                
      !  else if ((arg10 > 0.0_rprec) .AND. (arg13 < 0.0_rprec)) then

      !          avw1(1:nx,1:ny,k)=avw1(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw2(1:nx,1:ny,k)=avw2(1:nx,1:ny,k)+fr*arg13*arg10
      !          avw3(1:nx,1:ny,k)=avw3(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw4(1:nx,1:ny,k)=avw4(1:nx,1:ny,k)+fr*0.0_rprec

      !   else if ((arg10 < 0.0_rprec) .AND. (arg13 < 0.0_rprec)) then

      !          avw1(1:nx,1:ny,k)=avw1(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw2(1:nx,1:ny,k)=avw2(1:nx,1:ny,k)+fr*0.0_rprec
      !          avw3(1:nx,1:ny,k)=avw3(1:nx,1:ny,k)+fr*arg13*arg10
      !          avw4(1:nx,1:ny,k)=avw4(1:nx,1:ny,k)+fr*0.0_rprec

      !  else
      !          avw1(1:nx,1:ny,k)=avw1(1:nx,1:ny,k)+fr*0.0_rprec
       !         avw2(1:nx,1:ny,k)=avw2(1:nx,1:ny,k)+fr*0.0_rprec
       !         avw3(1:nx,1:ny,k)=avw3(1:nx,1:ny,k)+fr*0.0_rprec
       !         avw4(1:nx,1:ny,k)=avw4(1:nx,1:ny,k)+fr*arg13*arg10
       ! end if



        !if ((arg10 > 0.0_rprec) .AND. (arg11 > 0.0_rprec)) then
        !        awt1(1:nx,1:ny,k)=awt1(1:nx,1:ny,k)+fr*arg11*arg10
        !        awt2(1:nx,1:ny,k)=awt2(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt3(1:nx,1:ny,k)=awt3(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt4(1:nx,1:ny,k)=awt4(1:nx,1:ny,k)+fr*0.0_rprec

        !else if ((arg10 > 0.0_rprec) .AND. (arg11 < 0.0_rprec)) then

        !        awt1(1:nx,1:ny,k)=awt1(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt2(1:nx,1:ny,k)=awt2(1:nx,1:ny,k)+fr*arg11*arg10
        !        awt3(1:nx,1:ny,k)=awt3(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt4(1:nx,1:ny,k)=awt4(1:nx,1:ny,k)+fr*0.0_rprec

        ! else if ((arg10 < 0.0_rprec) .AND. (arg11 < 0.0_rprec)) then

        !        awt1(1:nx,1:ny,k)=awt1(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt2(1:nx,1:ny,k)=awt2(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt3(1:nx,1:ny,k)=awt3(1:nx,1:ny,k)+fr*arg11*arg10
        !        awt4(1:nx,1:ny,k)=awt4(1:nx,1:ny,k)+fr*0.0_rprec

        !else
        !        awt1(1:nx,1:ny,k)=awt1(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt2(1:nx,1:ny,k)=awt2(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt3(1:nx,1:ny,k)=awt3(1:nx,1:ny,k)+fr*0.0_rprec
        !        awt4(1:nx,1:ny,k)=awt4(1:nx,1:ny,k)+fr*arg11*arg10
        !end if



        !if ((arg10 > 0.0_rprec) .AND. (arg12 > 0.0_rprec)) then
        !        awq1(1:nx,1:ny,k)=awq1(1:nx,1:ny,k)+fr*arg12*arg10
        !        awq2(1:nx,1:ny,k)=awq2(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq3(1:nx,1:ny,k)=awq3(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq4(1:nx,1:ny,k)=awq4(1:nx,1:ny,k)+fr*0.0_rprec
        !else if ((arg10 > 0.0_rprec) .AND. (arg12 < 0.0_rprec)) then

        !        awq1(1:nx,1:ny,k)=awq1(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq2(1:nx,1:ny,k)=awq2(1:nx,1:ny,k)+fr*arg12*arg10
        !        awq3(1:nx,1:ny,k)=awq3(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq4(1:nx,1:ny,k)=awq4(1:nx,1:ny,k)+fr*0.0_rprec

        ! else if ((arg10 < 0.0_rprec) .AND. (arg12 < 0.0_rprec)) then

        !        awq1(1:nx,1:ny,k)=awq1(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq2(1:nx,1:ny,k)=awq2(1:nx,1:ny,k)+fr*0.0_rprec
        !        awq3(1:nx,1:ny,k)=awq3(1:nx,1:ny,k)+fr*arg12*arg10
        !        awq4(1:nx,1:ny,k)=awq4(1:nx,1:ny,k)+fr*0.0_rprec

        !else
         !       awq1(1:nx,1:ny,k)=awq1(1:nx,1:ny,k)+fr*0.0_rprec
         !       awq2(1:nx,1:ny,k)=awq2(1:nx,1:ny,k)+fr*0.0_rprec
         !       awq3(1:nx,1:ny,k)=awq3(1:nx,1:ny,k)+fr*0.0_rprec
         !       awq4(1:nx,1:ny,k)=awq4(1:nx,1:ny,k)+fr*arg12*arg10
        !end if

end do
!end do
!end do
!---------------

if (mod(jt,p_count)==0) then
        allocate(avg_out(1:nx,1:ny,1:(nz_tot-1)));

        call collocate_MPI_averages_SHH2(au,avg_out,120,'u')
        call collocate_MPI_averages_SHH2(av,avg_out,121,'v')
        call collocate_MPI_averages_SHH2(aw,avg_out,122,'w')
        call collocate_MPI_averages_SHH2(atheta,avg_out,125,'theta')
        call collocate_MPI_averages_SHH2(aq,avg_out,126,'q')
        call collocate_MPI_averages_SHH2(u2,avg_out,129,'u2')
        call collocate_MPI_averages_SHH2(v2,avg_out,130,'v2')
        call collocate_MPI_averages_SHH2(w2,avg_out,131,'w2')
        call collocate_MPI_averages_SHH2(theta2,avg_out,132,'T2')
        call collocate_MPI_averages_SHH2(q2,avg_out,133,'q2')
        call collocate_MPI_averages_SHH2(auw,avg_out,161,'uw')
        call collocate_MPI_averages_SHH2(avw,avg_out,162,'vw')
        call collocate_MPI_averages_SHH2(awt,avg_out,163,'wt')
        call collocate_MPI_averages_SHH2(awq,avg_out,164,'wq')
        call collocate_MPI_averages_SHH2(atxz,avg_out,186,'txz')
        call collocate_MPI_averages_SHH2(atyz,avg_out,188,'tyz')
        call collocate_MPI_averages_SHH2(asgst3,avg_out,258,'sgst3')
        call collocate_MPI_averages_SHH2(asgsq3,avg_out,261,'sgsq3')
        
        deallocate(avg_out)
        
        !!VK Zero out the outputted averages !!

          au=0._rprec;av=0._rprec;aw=0._rprec;
          atheta=0._rprec;
        aq=0._rprec; u2=0._rprec;v2=0._rprec;
        w2=0._rprec;theta2=0._rprec; q2=0._rprec;
         auw=0._rprec;avw=0._rprec;awt=0._rprec;awq=0._rprec;
         atxz=0._rprec; atyz=0._rprec;
        asgst3=0._rprec;asgsq3=0._rprec;
        arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;
        arg7=0._rprec;arg8=0._rprec;

end if

          
if (mod(jt,3*p_count)==0) then
        allocate(avg_out(1:nx,1:ny,1:(nz_tot-1))); 
        
        call collocate_MPI_averages_SHH2(ap,avg_out,123,'p')
        call collocate_MPI_averages_SHH2(p2,avg_out,127,'p2')
        call collocate_MPI_averages_SHH2(u3,avg_out,134,'u3')
        call collocate_MPI_averages_SHH2(v3,avg_out,135,'v3')
        call collocate_MPI_averages_SHH2(w3,avg_out,136,'w3')
        call collocate_MPI_averages_SHH2(q3,avg_out,137,'q3')
        call collocate_MPI_averages_SHH2(T3,avg_out,138,'T3')
        call collocate_MPI_averages_SHH2(adudx,avg_out,139,'dudx')
        call collocate_MPI_averages_SHH2(adudy,avg_out,140,'dudy')
        call collocate_MPI_averages_SHH2(adudz,avg_out,141,'dudz')
        call collocate_MPI_averages_SHH2(advdx,avg_out,142,'dvdx')
        call collocate_MPI_averages_SHH2(advdy,avg_out,143,'dvdy')
        call collocate_MPI_averages_SHH2(advdz,avg_out,144,'dvdz')
        call collocate_MPI_averages_SHH2(adwdx,avg_out,145,'dwdx')
        call collocate_MPI_averages_SHH2(adwdy,avg_out,146,'dwdy')
        call collocate_MPI_averages_SHH2(adwdz,avg_out,147,'dwdz')
        call collocate_MPI_averages_SHH2(adTdx,avg_out,148,'dTdx')
        call collocate_MPI_averages_SHH2(adTdy,avg_out,149,'dTdy')
        call collocate_MPI_averages_SHH2(adTdz,avg_out,150,'dTdz')
        call collocate_MPI_averages_SHH2(adqdx,avg_out,151,'dqdx')
        call collocate_MPI_averages_SHH2(adqdy,avg_out,152,'dqdy')
        call collocate_MPI_averages_SHH2(adqdz,avg_out,153,'dqdz')
        call collocate_MPI_averages_SHH2(adpdx,avg_out,154,'dpdx')
        call collocate_MPI_averages_SHH2(adpdy,avg_out,155,'dpdy')
        call collocate_MPI_averages_SHH2(adpdz,avg_out,156,'dpdz')
        call collocate_MPI_averages_SHH2(tke,avg_out,157,'tke')
        call collocate_MPI_averages_SHH2(awe,avg_out,158,'we')
        call collocate_MPI_averages_SHH2(aue,avg_out,159,'ue')
        call collocate_MPI_averages_SHH2(ave,avg_out,160,'ve')

        call collocate_MPI_averages_SHH2(aut,avg_out,165,'ut')
        call collocate_MPI_averages_SHH2(avt,avg_out,166,'vt')
        call collocate_MPI_averages_SHH2(auq,avg_out,167,'uq')
        call collocate_MPI_averages_SHH2(avq,avg_out,168,'vq')
        call collocate_MPI_averages_SHH2(auv,avg_out,169,'uv')
        call collocate_MPI_averages_SHH2(awuu,avg_out,170,'wuu')
        call collocate_MPI_averages_SHH2(awvv,avg_out,171,'wvv')
        call collocate_MPI_averages_SHH2(awtt,avg_out,172,'wtt')
        call collocate_MPI_averages_SHH2(awqq,avg_out,173,'wqq')
        call collocate_MPI_averages_SHH2(awwu,avg_out,174,'wwu')
        call collocate_MPI_averages_SHH2(awwv,avg_out,175,'wwv')
        call collocate_MPI_averages_SHH2(awwt,avg_out,176,'wwt')
        call collocate_MPI_averages_SHH2(awwq,avg_out,177,'wwq')
        !call collocate_MPI_averages_SHH2(pi1_t,avg_out,178,'pi1t')
        !call collocate_MPI_averages_SHH2(pi2_t,avg_out,179,'pi2t')
        !call collocate_MPI_averages_SHH2(pi3_t,avg_out,180,'pi3t')
        !call collocate_MPI_averages_SHH2(pi1_q,avg_out,181,'pi1q')
        !call collocate_MPI_averages_SHH2(pi2_q,avg_out,182,'pi2q')
        !call collocate_MPI_averages_SHH2(pi3_q,avg_out,183,'pi3q')
        call collocate_MPI_averages_SHH2(atxx,avg_out,184,'txx')
        call collocate_MPI_averages_SHH2(atxy,avg_out,185,'txy')
        call collocate_MPI_averages_SHH2(atyy,avg_out,187,'tyy')
        call collocate_MPI_averages_SHH2(atzz,avg_out,189,'tzz')
        call collocate_MPI_averages_SHH2(aup,avg_out,190,'up')
        call collocate_MPI_averages_SHH2(avp,avg_out,191,'vp')
        call collocate_MPI_averages_SHH2(awp,avg_out,192,'wp')
        call collocate_MPI_averages_SHH2(atp,avg_out,193,'tp')
        call collocate_MPI_averages_SHH2(aqp,avg_out,194,'qp')
        call collocate_MPI_averages_SHH2(atq,avg_out,195,'tq')
        call collocate_MPI_averages_SHH2(aupx,avg_out,196,'upx')
        call collocate_MPI_averages_SHH2(aupy,avg_out,197,'upy')
        call collocate_MPI_averages_SHH2(aupz,avg_out,198,'upz')
        call collocate_MPI_averages_SHH2(avpx,avg_out,199,'vpx')
        call collocate_MPI_averages_SHH2(avpy,avg_out,200,'vpy')
        call collocate_MPI_averages_SHH2(avpz,avg_out,201,'vpz')
        call collocate_MPI_averages_SHH2(awpx,avg_out,202,'wpx')
        call collocate_MPI_averages_SHH2(awpy,avg_out,203,'wpy')
        call collocate_MPI_averages_SHH2(awpz,avg_out,204,'wpz')
        call collocate_MPI_averages_SHH2(at11ux,avg_out,205,'t11ux')
        call collocate_MPI_averages_SHH2(at12uy,avg_out,206,'t12uy')
        call collocate_MPI_averages_SHH2(at13uz,avg_out,207,'t13uz')
        call collocate_MPI_averages_SHH2(at21vx,avg_out,208,'t21vx')
        call collocate_MPI_averages_SHH2(at22vy,avg_out,209,'t22vy')
        call collocate_MPI_averages_SHH2(at23vz,avg_out,210,'t23vz')
        call collocate_MPI_averages_SHH2(at31wx,avg_out,211,'t31wx')
        call collocate_MPI_averages_SHH2(at32wy,avg_out,212,'t32wy')
        call collocate_MPI_averages_SHH2(at33wz,avg_out,213,'t33wz')
        call collocate_MPI_averages_SHH2(apux,avg_out,214,'pux')
        call collocate_MPI_averages_SHH2(apuy,avg_out,215,'puy')
        call collocate_MPI_averages_SHH2(apuz,avg_out,216,'puz')
        call collocate_MPI_averages_SHH2(apvx,avg_out,217,'pvx')
        call collocate_MPI_averages_SHH2(apvy,avg_out,218,'pvy')
        call collocate_MPI_averages_SHH2(apvz,avg_out,219,'pvz')
        call collocate_MPI_averages_SHH2(apwx,avg_out,220,'pwx')
        call collocate_MPI_averages_SHH2(apwy,avg_out,221,'pwy')
        call collocate_MPI_averages_SHH2(apwz,avg_out,222,'pwz')
        call collocate_MPI_averages_SHH2(aptx,avg_out,223,'ptx')
        call collocate_MPI_averages_SHH2(apty,avg_out,224,'pty')
        call collocate_MPI_averages_SHH2(aptz,avg_out,225,'ptz')
        call collocate_MPI_averages_SHH2(apqx,avg_out,226,'pqx')
        call collocate_MPI_averages_SHH2(apqy,avg_out,227,'pqy')
        call collocate_MPI_averages_SHH2(apqz,avg_out,228,'pqz')
        call collocate_MPI_averages_SHH2(at11s11,avg_out,229,'t11s11')
        call collocate_MPI_averages_SHH2(at12s12,avg_out,230,'t12s12')
        call collocate_MPI_averages_SHH2(at13s13,avg_out,231,'t13s13')
        call collocate_MPI_averages_SHH2(at22s22,avg_out,232,'t22s22')
        call collocate_MPI_averages_SHH2(at23s23,avg_out,233,'t23s23')
        call collocate_MPI_averages_SHH2(at33s33,avg_out,234,'t33s33')
        !call collocate_MPI_averages_SHH2(api1_tx,avg_out,235,'pi1_tx')
        !call collocate_MPI_averages_SHH2(api2_ty,avg_out,236,'pi2_ty')
        !call collocate_MPI_averages_SHH2(api3_tz,avg_out,237,'pi3_tz')
        !call collocate_MPI_averages_SHH2(api1_qx,avg_out,238,'pi1_qx')
        !call collocate_MPI_averages_SHH2(api2_qy,avg_out,239,'pi2_qy')
        !call collocate_MPI_averages_SHH2(api3_qz,avg_out,240,'pi3_qz')
        call collocate_MPI_averages_SHH2(as11,avg_out,241,'s11')
        call collocate_MPI_averages_SHH2(as12,avg_out,242,'s12')
        call collocate_MPI_averages_SHH2(as13,avg_out,243,'s13')
        call collocate_MPI_averages_SHH2(as22,avg_out,244,'s22')
        call collocate_MPI_averages_SHH2(as23,avg_out,245,'s23')
        call collocate_MPI_averages_SHH2(as33,avg_out,246,'s33')
        call collocate_MPI_averages_SHH2(aut11,avg_out,247,'ut11')
        call collocate_MPI_averages_SHH2(aut12,avg_out,248,'ut12')
        call collocate_MPI_averages_SHH2(aut13,avg_out,249,'ut13')
        call collocate_MPI_averages_SHH2(avt21,avg_out,250,'vt21')
        call collocate_MPI_averages_SHH2(avt22,avg_out,251,'vt22')
        call collocate_MPI_averages_SHH2(avt23,avg_out,252,'vt23')
        call collocate_MPI_averages_SHH2(awt31,avg_out,253,'wt31')
        call collocate_MPI_averages_SHH2(awt32,avg_out,254,'wt32')
        call collocate_MPI_averages_SHH2(awt33,avg_out,255,'wt33')
        !call collocate_MPI_averages_SHH2(asgst1,avg_out,256,'sgst1')
        !call collocate_MPI_averages_SHH2(asgst2,avg_out,257,'sgst2')

        !call collocate_MPI_averages_SHH2(auw1,avg_out,262,'uw1')
        !call collocate_MPI_averages_SHH2(auw2,avg_out,263,'uw2')
        !call collocate_MPI_averages_SHH2(auw3,avg_out,264,'uw3')
        !call collocate_MPI_averages_SHH2(auw4,avg_out,265,'uw4')
        !call collocate_MPI_averages_SHH2(awt1,avg_out,266,'wt1')
        !call collocate_MPI_averages_SHH2(awt2,avg_out,267,'wt2')
        !call collocate_MPI_averages_SHH2(awt3,avg_out,268,'wt3')
        !call collocate_MPI_averages_SHH2(awt4,avg_out,269,'wt4')
        !call collocate_MPI_averages_SHH2(awq1,avg_out,270,'wq1')
        !call collocate_MPI_averages_SHH2(awq2,avg_out,271,'wq2')
        !call collocate_MPI_averages_SHH2(awq3,avg_out,272,'wq3')
        !call collocate_MPI_averages_SHH2(awq4,avg_out,273,'wq4')
        !call collocate_MPI_averages_SHH2(atpx,avg_out,274,'tpx')
        !call collocate_MPI_averages_SHH2(atpy,avg_out,275,'tpy')
        !call collocate_MPI_averages_SHH2(atpz,avg_out,276,'tpz')
        call collocate_MPI_averages_SHH2(aactual_T,avg_out,277,'actual_T')
        call collocate_MPI_averages_SHH2(aactual_Tv,avg_out,278,'actual_Tv')
        call collocate_MPI_averages_SHH2(apr_atm,avg_out,279,'pr_atm')
        call collocate_MPI_averages_SHH2(arel_hum,avg_out,280,'rel_hum')
        call collocate_MPI_averages_SHH2(arel_hum_q,avg_out,281,'rel_hum_q')
        call collocate_MPI_averages_SHH2(avapor_pr,avg_out,282,'vapor_pr')
        call collocate_MPI_averages_SHH2(asat_vapor_pr,avg_out,283,'sat_vapor_pr')
        call collocate_MPI_averages_SHH2(asat_qmix,avg_out,284,'sat_qmix')
        call collocate_MPI_averages_SHH2(azlcl_all,avg_out,285,'zlcl_all')
        !call collocate_MPI_averages_SHH2(avw1,avg_out,286,'vw1')
        !call collocate_MPI_averages_SHH2(avw2,avg_out,287,'vw2')
        !call collocate_MPI_averages_SHH2(avw3,avg_out,288,'vw3')
        !call collocate_MPI_averages_SHH2(avw4,avg_out,289,'vw4')
        !call collocate_MPI_averages_SHH2(aqpx,avg_out,290,'qpx')
        !call collocate_MPI_averages_SHH2(aqpy,avg_out,291,'qpy')
        !call collocate_MPI_averages_SHH2(aqpz,avg_out,292,'qpz')

     
     
     deallocate(avg_out)
!
!!VK Zero out the outputted averages !!

ap=0._rprec;
        p2=0._rprec; 
       ! pf2=0._rprec;
        u3=0._rprec;v3=0._rprec;w3=0._rprec;
        q3=0._rprec; T3=0._rprec;adudx=0._rprec;adudy=0._rprec;adudz=0._rprec;
        advdx=0._rprec;advdy=0._rprec;advdz=0._rprec;adwdx=0._rprec;adwdy=0._rprec;adwdz=0._rprec;
        adTdx=0._rprec;adTdy=0._rprec;adTdz=0._rprec;adqdx=0._rprec;adqdy=0._rprec;adqdz=0._rprec;
         adpdx=0._rprec;adpdy=0._rprec;adpdz=0._rprec;tke=0._rprec;awe=0._rprec;aue=0._rprec;ave=0._rprec;
         aut=0._rprec;avt=0._rprec;auq=0._rprec;
         avq=0._rprec;auv=0._rprec;awuu=0._rprec;awvv=0._rprec;awqq=0._rprec;awtt=0._rprec;awwu=0._rprec;
         awwv=0._rprec;awwt=0._rprec;awwq=0._rprec;
         !pi1_t=0._rprec;pi2_t=0._rprec;pi3_t=0._rprec;
         !pi1_q=0._rprec;pi2_q=0._rprec;pi3_q=0._rprec;
         atxx=0._rprec;atyy=0._rprec;
        atzz=0._rprec;aup=0._rprec;avp=0._rprec;awp=0._rprec;atp=0._rprec;aqp=0._rprec;atq=0._rprec;
        aupx=0._rprec;aupy=0._rprec;aupz=0._rprec;avpx=0._rprec;avpy=0._rprec;avpz=0._rprec;
        awpx=0._rprec;awpy=0._rprec;awpz=0._rprec;
        !atpx=0._rprec;atpy=0._rprec;atpz=0._rprec;
        !aqpx=0._rprec;aqpy=0._rprec;aqpz=0._rprec;
        at11ux=0._rprec;at12uy=0._rprec;at13uz=0._rprec;
        at21vx=0._rprec;at22vy=0._rprec;at23vz=0._rprec;at31wx=0._rprec;at32wy=0._rprec;at33wz=0._rprec;
        apux=0._rprec;apuy=0._rprec;apuz=0._rprec;apvx=0._rprec;apvy=0._rprec;apvz=0._rprec;
        apwx=0._rprec;apwy=0._rprec;apwz=0._rprec;
        aptx=0._rprec;apty=0._rprec;aptz=0._rprec;
        apqx=0._rprec;apqy=0._rprec;apqz=0._rprec;
        at11s11=0._rprec;at12s12=0._rprec;at13s13=0._rprec;at22s22=0._rprec;at23s23=0._rprec;at33s33=0._rprec;
        !api1_tx=0._rprec;api2_ty=0._rprec;api3_tz=0._rprec;api1_qx=0._rprec;api2_qy=0._rprec;api3_qz=0._rprec;
        as11=0._rprec;as12=0._rprec;as13=0._rprec;as22=0._rprec;as23=0._rprec;as33=0._rprec;
        aut11=0._rprec;aut12=0._rprec;aut13=0._rprec;avt21=0._rprec;avt22=0._rprec;avt23=0._rprec;
        awt31=0._rprec;awt32=0._rprec;awt33=0._rprec;
        !asgst1=0._rprec;asgst2=0._rprec;asgst3=0._rprec;asgsq1=0._rprec;asgsq2=0._rprec;asgsq3=0._rprec;
        
        arg1=0._rprec;arg2=0._rprec;arg3=0._rprec;arg4=0._rprec;arg4f=0._rprec;arg5=0._rprec;arg6=0._rprec;
        arg7=0._rprec;arg8=0._rprec;arg9=0._rprec;arg10=0._rprec;arg11=0._rprec;
        arg12=0._rprec;ux=0._rprec;uy=0._rprec;uz=0._rprec;vx=0._rprec;vy=0._rprec;
        !awt1=0._rprec;awt2=0._rprec;awt3=0._rprec;awt4=0._rprec;awq1=0._rprec;awq2=0._rprec;
        !awq3=0._rprec;awq4=0._rprec;auw1=0._rprec;auw2=0._rprec;auw3=0._rprec;auw4=0._rprec;
        !avw1=0._rprec;avw2=0._rprec;avw3=0._rprec;avw4=0._rprec;
        aactual_T=0._rprec;aactual_Tv=0._rprec;apr_atm=0._rprec;arel_hum=0._rprec;arel_hum_q=0._rprec;
        avapor_pr=0._rprec;asat_vapor_pr=0._rprec;asat_qmix=0._rprec;azlcl_all=0._rprec;

end if
5168     format(1400(E15.6))

end subroutine MM_XYZ_Out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--The following subroutine does the collocation of the MPI arrays for
! SHHOutput Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_SHH(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,ind3,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(nx,ny,nz-1)::avg_var_proc
real(kind=rprec),dimension(nx,ny,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (avg_var_proc(1,1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1,1,1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
           do ind3=1,nz_tot-1
           do ind2=1,ny
            write(file_ind,5168)(avg_var_tot_domain(ind1,ind2,ind3),ind1=1,nx)
           end do
           end do
        close(file_ind)
  end if


5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_SHH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--The following subroutine does the collocation of the MPI arrays for
! SHHOutput Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collocate_MPI_averages_SHH2(avg_var_proc,avg_var_tot_domain,file_ind,filename_str)
!subroutine collocate_MPI_averages(avg_var_proc,avg_var_tot_domain,file_ind)
use param
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer :: ind1,ind2,ind3,file_ind
character (*),intent(in) :: filename_str
character (len=256) :: local_filename
real(kind=rprec),dimension(nx,ny,nz-1)::avg_var_proc
real(kind=rprec),dimension(nx,ny,nz_tot-1)::avg_var_tot_domain

local_filename=path//'output/aver_'//trim(filename_str)//'.out'

avg_var_tot_domain=0._rprec
$if ($MPI)
  recvcounts = size (avg_var_proc)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (avg_var_proc(1,1,1), size (avg_var_proc), MPI_RPREC,       &
                    avg_var_tot_domain(1,1,1), recvcounts, displs, MPI_RPREC, &
                    rank_of_coord(0), comm, ierr)
$else
  avg_var_tot_domain=avg_var_proc
$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
        open(file_ind,file=trim(local_filename),status="unknown",position="append")
           do ind3=1,250
           do ind2=1,ny
            write(file_ind,5168)(avg_var_tot_domain(ind1,ind2,ind3),ind1=1,nx)
           end do
           end do
        close(file_ind)
  end if


5168     format(1400(E14.5))
end subroutine collocate_MPI_averages_SHH2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--assumes lun is open and positioned correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (lun)

use param,only:nz,S_FLAG
use sim_param,only:u,v,w,RHSx,RHSy,RHSz,theta,qmix
use sgsmodule,only:Cs_opt2,F_LM,F_MM,F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t
use scalars_module,only:RHS_T,sgs_t3
use bottombc,only:psi_m,psi_m0,psi_h,psi_h0
use scalars_module_q,only:RHS_q,sgs_q3

implicit none
integer,intent(in)::lun

if (S_FLAG) then ! WITH SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),theta(:,:,1:nz),   &
               qmix(:,:,1:nz),RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz),          &
               RHS_T(:,:,1:nz),RHS_q(:,:,1:nz),sgs_t3(:,:,1),sgs_q3(:,:,1),&
               psi_m,psi_m0,psi_h,psi_h0,Cs_opt2,&
               F_LM,F_MM,F_QN,F_NN,G_LM,G_MM,G_QN,G_NN,Pr_t(:,:,1:nz)
else ! No SCALARS
   write (lun) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),          &
               RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz), &
               Cs_opt2,F_LM,F_MM,F_QN,F_NN
end if

end subroutine checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_final(jt, lun_opt)
implicit none

integer,intent(in)::jt
integer, intent (in), optional :: lun_opt  !--if present, write to unit lun
integer, parameter :: lun_default = 11
integer::jx,jy,jz
integer :: lun
logical :: opn

if (present (lun_opt)) then
  lun = lun_opt
else
  lun = lun_default
end if

inquire (unit=lun, opened=opn)

if (.not. opn) then
  write (*, *) 'output_final: lun=', lun, ' is not open'
  stop
end if

rewind (lun)

call checkpoint (lun)

close (lun)

if ((cumulative_time) .and. (lun == lun_default)) then
  !--only do this for true final output, not intermediate recording
  open (1, file=fcumulative_time)
  write (1, *) jt_total
  close (1)
end if

end subroutine output_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_lambda2_out
use sim_param,only:path
implicit none
character(len=24)::fname
call lambda2()
write(fname,'(A13,i6.6,A4)')path//'output/lam-',jt_total,'.out'
open(1,file=fname,form='unformatted')
write(1)nx,ny,nz
write(1)real(lam2)
close(1)
end subroutine io_lambda2_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lambda2()
use types,only:rprec
use sim_param,only:u,v,w,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use param,only:dx,dy,dz
implicit none
real(kind=rprec)::S11,S22,S33,S12,O12,S13,O13,S23,O23,&
     ux,uy,uz,vx,vy,vz,wx,wy,wz
integer::jx,jy,jz
! following used for eispack call...
integer::neis,nmeis,matzeis,ierreis,iv1eis(3)
double precision,dimension(3,3)::aeis,zeis
double precision,dimension(3)::wreis,wieis,fv1eis
double precision::ave

! assignments for eispack call
neis=3
nmeis=3
matzeis=0
ierreis=0
lam2=0._rprec

! at level z=dz/2.  (Where dudz and dvdz are on UVP nodes)
jz=1
do jy=1,ny
do jx=1,nx              
   ux=dudx(jx,jy,1)  ! uvp-node
   uy=dudy(jx,jy,1)  ! uvp-node
   uz=dudz(jx,jy,1)  ! uvp-node
   vx=dvdx(jx,jy,1)  ! uvp-node
   vy=dvdy(jx,jy,1)  ! uvp-node
   vz=dvdz(jx,jy,1)  ! uvp-node 
! special case
   wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  ! uvp-node
   wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  ! uvp-node
   wz=dwdz(jx,jy,1)  ! uvp-node
   S11=ux          ! uvp-node
   S12=0.5_rprec*(uy+vx) ! uvp-node
! taken care of with wall stress routine
   S13=0.5_rprec*(uz+wx) ! uvp
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! uvp-node
! taken care of with wall stress routine 
   S23=0.5_rprec*(vz+wy) ! uvp
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! uvp-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
  write (*, *) 'rg temporarily removed, sorry'; stop
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
! calculate derivatives/strain on w-nodes
do jz=2,nz-1  
do jy=1,ny
do jx=1,nx              
   ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))  ! w-node
   uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))  ! w-node
   uz=dudz(jx,jy,jz)  ! w-node
   vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))  ! w-node
   vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))  ! w-node
   vz=dvdz(jx,jy,jz)  ! w-node
   wx=dwdx(jx,jy,jz)  ! w-node
   wy=dwdy(jx,jy,jz)  ! w-node
   wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))  ! w-node
   S11=ux          ! w-node
   S12=0.5_rprec*(uy+vx) ! w-node
   S13=0.5_rprec*(uz+wx) ! w-node
   O12=0.5_rprec*(uy-vx) ! w-node
   O13=0.5_rprec*(uz-wx) ! w-node
   S22=vy          ! w-node
   S23=0.5_rprec*(vz+wy) ! w-node
   O23=0.5_rprec*(vz-wy) ! w-node
   S33=wz          ! w-node
   aeis(1,1)=s11*s11+s12*s12+s13*s13-O12*O12-O13*O13
   aeis(1,2)=s11*s12+s12*s22+s13*s23-O13*O23
   aeis(1,3)=s11*s13+s12*s23+s13*s33+O12*O23
   aeis(2,2)=s12*s12+s22*s22+s23*s23-O12*O12-O23*O23
   aeis(2,3)=s12*s13+s22*s23+s23*s33-O12*O13
   aeis(3,3)=s13*s13+s23*s23+s33*s33-O13*O13-O23*O23
   aeis(2,1)=aeis(1,2)
   aeis(3,1)=aeis(1,3)
   aeis(3,2)=aeis(2,3)
   if(wreis(1).ge.wreis(2).and.wreis(1).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(1),kind=rprec)
   elseif(wreis(2).ge.wreis(1).and.wreis(2).le.wreis(3)) then
      lam2(jx,jy,jz)=real(wreis(2),kind=rprec)
   else
      lam2(jx,jy,jz)=real(wreis(3),kind=rprec)
   endif
end do
end do
end do

print*,'minmax',minval(lam2),maxval(lam2)
end subroutine lambda2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine io_mean_out
implicit none
write(51)real(mean_u),real(mean_u2),real(mean_v),real(mean_v2),&
     real(mean_w),real(mean_w2)
!--the nz/4*3 stuff has to go
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz/4*3)=0._rprec
end subroutine io_mean_out

subroutine calculate_mean
use sim_param,only:u,v,w
use sgsmodule,only:Cs_opt2,Cs_opt2_avg
implicit none
Cs_opt2_avg(:,:,:)=Cs_opt2_avg(:,:,:)+Cs_opt2(:,:,:)/nwrite
!TS
mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_u2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     u(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_v2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     v(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)/nwrite
mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)=&
     mean_w2(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)+&
     w(jx_pls:jx_ple,jy_pls:jy_ple,1:nz)**2/nwrite
end subroutine calculate_mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timeseries_spec
use sim_param,only:u,v,w,theta
implicit none
integer::jx,jy,jz,i
if(mod(jt_total,time_spec)==0.and.jt_total.gt.2000)then
jx=NX/8
jy=NY/2+1
jz=NZ/2
endif
end subroutine timeseries_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine post_spec(jt_local)
use sim_param,only:path,u,v,w,theta
use param
use fft
implicit none
real(kind=rprec),dimension(nx/2,nz)::spectra_u,spectra_v,spectra_w,&
     spectra_theta
real(kind=rprec),dimension(4,nx/2,nz-1)::spectra_uvwT
real(kind=rprec),dimension(4,nx/2,nz_tot-1)::spectra_uvwT_tot
integer,intent(in)::jt_local
integer::k,jz,z
character(len=64)::fname1,fname2,fname3,fname4
$if ($MPI)
  $define $lbz 0
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$else
  $define $lbz 1
$endif

write(fname1,'(a,a)') path//'output/spec_x','.dat'
open(82,file=fname1,form='formatted')
do jz=1,nz-1
   z=(jz-0.5_rprec)*dz*z_i
   write(82,*) (real(kx(k,1)/z_i*z),k=1,nx/2)
   call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
   call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
   call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
   call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))
enddo
   close(82)
$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
                    spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
                    MPI_RPREC, rank_of_coord(0), comm, ierr)
$else
  spectra_uvwT_tot=spectra_uvwT
$endif

if((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
   write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
   open(83,file=fname1,form='unformatted')
   write(83) real(spectra_uvwT_tot(:,1:nx/2,:))
   close(83)
end if

end subroutine post_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spectrum(u, spec)
use fft
implicit none      
real(kind=rprec),dimension(ld,ny),intent(in)::u
real(kind=rprec),dimension(nx/2),intent(out)::spec  !--assumes Nyquist is 0

integer::jy,jz,k
real(kind=rprec),dimension(nx)::vel_r,vel_c

integer*8, save :: plan
logical, save :: init = .false.

if (.not. init) then
  call rfftw_f77_create_plan(plan,nx,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
  init = .true.
end if

! initialize
spec(:)=0._rprec
do jy=1,ny
   vel_r(:)= u(1:nx,jy)/real(nx,kind=rprec)
! check this normaliztion-part of forward; call the fft
   call rfftw_f77_one(plan,vel_r,vel_c)
! compute magnitudes the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
   spec(1)=spec(1)+0.5*vel_c(1)*vel_c(1)
   do k=2,nx/2
      spec(k)=spec(k)+vel_c(k)*vel_c(k)+vel_c(nx+2-k)*vel_c(nx+2-k)
   end do

   !--assume Nyquist is 0
   !spec(nx/2+1)=spec(nx/2+1)+vel_c(nx/2+1)*vel_c(nx/2+1)
end do
spec(:)=spec(:)/real(Ny,kind=rprec) ! for average over Ny
end subroutine spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avg_stats ()
use param
use sim_param, only : u, v, w, txz
use fft, only : kx
implicit none

!--choose naming convention that does not conflict with qpost
character (*), parameter :: fubar_avg = 'output/ubar-avg_stats.dat'
character (*), parameter :: fupr2bar_avg = 'output/upr2bar-avg_stats.dat'
character (*), parameter :: fstressbar_avg = 'output/stressbar-avg_stats.dat'
character (*), parameter :: fEozbar_avg = 'output/Eozbar-avg_stats.dat'

integer, parameter :: hdr_len = 256
logical, parameter :: DEBUG = .false.
character (hdr_len) :: Eozbar_hdr

$if ($MPI)
  $define $lbz 0
$else
  $define $lbz 1
$endif
$if ($MPI)
  integer :: recvcounts(nproc)
  integer :: displs(nproc)
$endif
integer, save :: n_ubar_avg
integer, save :: n_upr2bar_avg
integer, save :: n_stressbar_avg
integer, save :: n_Eozbar_avg
integer :: jz

logical, save :: init = .false.

real (rprec) :: z
real (rprec) :: zu(1, nz_tot-1)
real (rprec) :: kz_z(2, nx/2)
real (rprec), save :: ubar_avg(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_avg(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_avg(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_avg(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
!--tot is a temp for current stats at nz_tot size
real (rprec), save :: ubar_tot(1, nz_tot-1)      !--1 is <u>
real (rprec), save :: upr2bar_tot(3, nz_tot-1)   !--<u'^2>, <v'^2>, <w'^2>
real (rprec), save :: stressbar_tot(3, nz_tot-1) !--1 is <u'w'>, 2 is <txz>, 3 is <u'w'> + <txz>
real (rprec), save :: Eozbar_tot(1, nx/2, nz_tot-1)  !--E11(k1,z)/z
real (rprec) :: upr(nx, ny), vpr(nx, ny), wpr(nx, ny)
real (rprec) :: ubar(nz-1), vbar(nz-1), wbar(nz-1)
real (rprec) :: upr2bar(3, nz-1)
real (rprec) :: stressbar(3, nz-1)
real (rprec) :: Eozbar(nx/2, nz-1)
!---------------------------------------------------------------------

!--check whether or not to actually do anything
!--motivation for doing this way is that it cleans up interface in main
if (modulo (jt, n_avg_stats) /= 0) goto 001  !--do nothing, exit cleanly

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then  
  if (.not. init) then  !--initialization

    call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
    call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
    call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg) 
    do jz = 1, nz-2
      call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                     leaveopn='yes')
    end do
    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)

    init = .true.

  end if
end if

!--calculation of current stats
do jz = $lbz, nz-1

  ubar(jz) = sum (u(1:nx, 1:ny, jz)) / (nx * ny)
  vbar(jz) = sum (v(1:nx, 1:ny, jz)) / (nx * ny)
  wbar(jz) = sum (w(1:nx, 1:ny, jz)) / (nx * ny)

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
       (jz == 1) ) then
    upr = 0._rprec
    vpr = 0._rprec
    wpr = 0._rprec
  else
    !--see qpost for u/w-node interpolation
    !--convention will be to put up, vp, wp on w-nodes
    upr = 0.5_rprec * (u(1:nx, 1:ny, jz) - ubar(jz) +  &
                       u(1:nx, 1:ny, jz-1) - ubar(jz-1))
    vpr = 0.5_rprec * (v(1:nx, 1:ny, jz) - vbar(jz) +  &
                       v(1:nx, 1:ny, jz-1) - vbar(jz-1))
    wpr = w(1:nx, 1:ny, jz) - wbar(jz)
  end if
 
  upr2bar(1, jz) = sum (upr**2) / (nx * ny)
  upr2bar(2, jz) = sum (vpr**2) / (nx * ny)
  upr2bar(3, jz) = sum (wpr**2) / (nx * ny)

  stressbar(1, jz) = sum (upr * wpr) / (nx * ny) 
  stressbar(2, jz) = sum (txz(1:nx, 1:ny, jz)) / (nx * ny)
  stressbar(3, jz) = sum (stressbar(1:2, jz))

  !--energy spectra
  call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
  z = (jz - 0.5_rprec) * dz
  Eozbar(:, jz) = Eozbar(:, jz) / z

end do

!--collect current stats into nz_tot sized arrays
$if ($MPI)

  if (DEBUG) then
    write (*, *) coord, ': ubar(1) = ', ubar(1)
  end if

  recvcounts = size (ubar)
  displs = coord_of_rank * recvcounts 
  call mpi_gatherv (ubar(1), size (ubar), MPI_RPREC,                &
                    ubar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (upr2bar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (upr2bar(1, 1), size (upr2bar), MPI_RPREC,          &
                    upr2bar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)  

  recvcounts = size (stressbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (stressbar(1, 1), size (stressbar), MPI_RPREC,        &
                    stressbar_tot(1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

  recvcounts = size (Eozbar)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (Eozbar(1, 1), size (Eozbar), MPI_RPREC,              &
                    Eozbar_tot(1, 1, 1), recvcounts, displs, MPI_RPREC,  &
                    rank_of_coord(0), comm, ierr)

$else

  ubar_tot(1, :) = ubar
  upr2bar_tot = upr2bar
  stressbar_tot = stressbar
  Eozbar_tot(1, :, :) = Eozbar

$endif

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !--calculation of cumulative average stats
  ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
  n_ubar_avg = n_ubar_avg + 1

  upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
                (n_upr2bar_avg + 1)
  n_upr2bar_avg = n_upr2bar_avg + 1

  stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
                  (n_stressbar_avg + 1)
  n_stressbar_avg = n_stressbar_avg + 1

  Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
  n_Eozbar_avg = n_Eozbar_avg + 1

  !--prepare list of z-coordinates
  forall (jz=1:nz_tot-1) zu(1, jz) = (jz - 0.5_rprec) * dz
  !--prepare  header, optional

  !--write out to file
  call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
  call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
  call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)

  !--this is a bit awkward: maybe just have another routine to do it right
  Eozbar_hdr = 'zone' !--this is for tecplot... 
  kz_z(1, :) = kx(1:nx/2, 1) * zu(1, 1)
  kz_z(2, :) = zu(1, 1)
  call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, 1),  &
                  hdr=Eozbar_hdr) 

  do jz = 2, nz_tot - 1
    kz_z(1, :) = kx(1:nx/2, 1) * zu(1, jz)
    kz_z(2, :) = zu(1, jz)
    call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z, Eozbar_avg(:, :, jz),  &
                    hdr=Eozbar_hdr, position='append') 
  end do

end if

001 continue  !--exit cleanly

end subroutine avg_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_ccol  !--num. coord columns: x, y, etc.
real (rprec), intent (out) :: a_avg(:, :)
integer, intent (out) :: n_avg
character (*), optional, intent (in) :: leaveopn
character (128) :: buff
logical :: exst, opn
integer :: j
real (rprec) :: z(n_ccol)

!---------------------------------------------------------------------
inquire (file=file_avg, exist=exst, opened=opn)

if (exst) then

  if (.not. opn) then
    open (1, file=file_avg)
    read (1, '(a)') buff

    if (buff(1:1) == '#') then
      read (buff(2:), *) n_avg
    else
      write (*, *) 'avg_stats: error'
      write (*, *) trim (file_avg), ' does not have expected format on line 1'
      stop  !--need to replace all stops with nice mpi exits
    end if
  end if

  !--skip data header lines here
  do
    read (1, '(a)') buff
    if (trim (buff) == trim (end_hdr_avg)) exit
  end do

  do j = 1, size (a_avg, 2)
    read (1, *) z, a_avg(:, j)  !--z is just placeholder here
  end do

  if (present (leaveopn)) then
    if (leaveopn /= 'yes') close (1)  !--case sensitive here
  else
    close (1)
  end if

else

  n_avg = 0
  a_avg = 0._rprec

end if

end subroutine init_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
implicit none

character (*), intent (in) :: file_avg
integer, intent (in) :: n_avg
real (rprec), intent (in) :: x(:, :)  !--coord columns for x, y, etc
real (rprec), intent (in) :: a_avg(:, :)

character (*), optional, intent (in) :: hdr
character (*), optional, intent (in) :: position

character (64) :: r_fmt, fmt
character (32) :: posn

integer :: j

!---------------------------------------------------------------------

!--check sizes compatible
if (size (x, 2) /= size (a_avg, 2)) then
  write (*, *) 'write_avg: error with sizes of x, a_avg'
  stop
end if

if (present (position)) then
  posn = position
else
  posn = 'rewind'
end if

open (1, file=file_avg, position=posn)

if (trim (posn) /= 'append') then  !--case sensitive
  write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
end if

if (present (hdr)) then
  !--write data header, if present
  write (1, '(a)') trim (hdr)
end if

!--write something to indicate end of header, always do this
write (1, '(a)') end_hdr_avg

!--define output format
write (r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                           '.', precision (1._rprec)
write (fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                         '(1x,', trim (r_fmt), '))'

!--write to file
do j = 1, size (a_avg, 2)
  write (1, fmt) x(:, j), a_avg(:, j)
end do

close (1)

end subroutine write_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_write ()
use param, only : jt_total, model, jt_start_write, buff_end,  &
                  read_inflow_file, write_inflow_file
use sgsmodule, only : F_MM, F_LM, F_QN, F_NN
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_write'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: field_file = 'output/inflow.vel.out'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80
integer, parameter :: field_lun = 81

logical, parameter :: DEBUG = .false.

character (64) :: fname

integer, save :: rec = 0
integer :: nrec
integer :: iolen
integer :: iend, iend_w

logical, save :: initialized = .false.
logical :: opn, exst

!---------------------------------------------------------------------

!--option check
if ( read_inflow_file .and. write_inflow_file ) then
  write (*, *) sub // ': cannot have read_inflow_file and write_inflow_file'
  stop
end if

!--check consistency with inflow_cond
iend = floor (buff_end * nx + 1._rprec)
iend_w = modulo (iend - 1, nx) + 1

if (.not. initialized) then

  inquire ( unit=lun, exist=exst, opened=opn )
  if ( .not. exst ) then
    write (*, *) sub // ': lun = ', lun, ' does not exist'
    stop
  end if
  if (opn) then
    write (*, *) sub // ': lun = ', lun, ' is already open'
    stop
  end if

  if ( USE_MPI ) then
      write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
  else
      write ( fname, '(a)' ) inflow_file
  end if
  
  inquire ( file=fname, exist=exst, opened=opn )
  if (exst .and. opn) then
    write (*, *) sub // ': file = ', trim (fname), ' is already open'
    stop
  end if
  
  !--figure out the record length
  if ( model.eq.4 ) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if (model.eq.5) then
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    inquire (iolength=iolen) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if

  !--always add to existing inflow_file
  !--inflow_file records always start at 1
  if ( exst ) then
      !--figure out the number of records already in file
      call len_da_file (fname, iolen, nrec)
      write (*, *) sub // ': #records in ' // trim (fname) // '= ', nrec
      rec = nrec
  else
      rec = 0
  end if
  
  !--using direct-access file to allow implementation of 'inflow recycling'
  !  more easily
  !--may want to put in some simple checks on ny, nz
  open (unit=lun, file=fname, access='direct', action='write',  &
        recl=iolen)

  initialized = .true.

end if

if (jt_total == jt_start_write) then  !--write entire flow field out
  inquire (unit=field_lun, exist=exst, opened=opn)
  if (exst .and. .not. opn) then
    open (unit=field_lun, file=field_file, form='unformatted')
    call output_final (jt_total, field_lun)
  else
    write (*, *) sub // ': problem opening ' // field_file
    stop
  end if
end if

if (jt_total >= jt_start_write) then
  rec = rec + 1
  if ( model.eq.4 ) then
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:)
  else if ( model.eq.5) then 
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), F_MM(1,:,:), F_LM(1,:,:), F_QN(1,:,:), F_NN(1,:,:)
  else
    write (unit=lun, rec=rec) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:)
  end if
  if ( DEBUG ) write (*, *) sub // ': wrote record ', rec
end if

end subroutine inflow_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inflow_read ()
use param, only : model, ny, nz, pi, nsteps, jt_total, buff_end
use sgsmodule, only : FMM_hold, FLM_hold, FQN_hold, FNN_hold
use sim_param, only : u, v, w
implicit none

character (*), parameter :: sub = 'inflow_read'
character (*), parameter :: inflow_file = 'output/inflow_BC.out'
character (*), parameter :: debug_file = 'inflow_read_debug.dat'
character (*), parameter :: MPI_suffix = '.c'

integer, parameter :: lun = 80  !--inflow_write for now
integer, parameter :: lun_DEBUG = 88

integer, parameter :: l_blend = 300  !--length of blending zone (recycling)
                                     !--should correspond to integral scale
                                     !--this is number of t-steps
logical, parameter :: recycle = .false.

logical, parameter :: DEBUG = .false.

character (32) :: fmt
character (64) :: fname

!--check for consistency with sim_param here
!--could define a fortran integer lbz in sim_param, and make it visible
!  here, however, this may have complications elsewhere where the name lbz
!  is used.
$if ( $MPI )
    $define $lbz 0
$else
    $define $lbz 1
$endif

integer :: jy, jz
integer :: iend, iend_w
integer :: i
integer :: iolen
integer, save :: rec
integer, save :: nrec
integer :: recp

logical, save :: init_DEBUG = .false.
logical, save :: initialized = .false.
logical :: exst, opn

real (rprec) :: wgt

real (rprec) :: u_tmp(ny, $lbz:nz), v_tmp(ny, $lbz:nz), w_tmp(ny, $lbz:nz)

!---------------------------------------------------------------------

iend = floor ( buff_end * nx + 1.0_rprec )
iend_w = modulo ( iend - 1, nx ) + 1

if ( .not. initialized ) then

    inquire ( unit=lun, exist=exst, opened=opn )
    if ( .not. exst ) then
        write (*, *) sub // ': lun = ', lun, ' does not exist'
        stop
    end if
    if ( opn ) then
        write (*, *) sub // ': lun = ', lun, ' is already open'
        stop
    end if

    if ( USE_MPI ) then
        write ( fname, '(a,a,i0)' ) trim (inflow_file), MPI_suffix, coord
    else
        write ( fname, '(a)' ) inflow_file
    end if
    
    inquire ( file=fname, exist=exst, opened=opn )
    if ( exst ) then
        if ( opn ) then
            write (*, *) sub // ': file = ', fname, ' is already open'
            stop
        end if
    else
        write (*, *) sub // ': file = ', fname, ' does not exist'
        stop
    end if

    !--can only reach this point if exst and .not. opn
  
    !--figure out the record length
    if ( model.eq.4 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
    else if ( model.eq.5 ) then
        inquire ( iolength=iolen ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FMM_hold, FQN_hold, FNN_hold 
    else
        inquire ( iolength=iolen ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :)
    endif
 
    !--figure out the number of records
    call len_da_file ( fname, iolen, nrec )

    write (*, *) sub // ': number of records = ', nrec

    if ( recycle ) then
        !--check minimum length requirement
        !  checks that there are some points that will be non-blended
        
        if ( 2 * (l_blend - 1) > nrec ) then
            write (*, *) sub // ': ', fname, 'is too short to recycle'
            stop
        end if
    end if

    open ( unit=lun, file=fname, access='direct', action='read',  &
           recl=iolen )

    !--file always starts a record 1, but in continued runs, we may need to
    !  access a record that is not 1 to begin with
    !--actually, with wrap-around of records, it means the reading can start
    !  at any point in the file and be OK
    !--intended use: jt_total = 1 here at start of set of runs reading
    !  from the inflow_file, so the first record read will be record 1
    rec = jt_total - 1

    initialized = .true.

end if
rec = rec + 1
if ( recycle ) then
    rec = modulo ( rec - 1, nrec - l_blend + 1 ) + 1
else
    rec = modulo ( rec - 1, nrec ) + 1
end if

if ( model.eq.4 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold
else if ( model.eq.5 ) then
    read ( unit=lun, rec=rec ) u(iend_w,:,:), v(iend_w,:,:), w(iend_w,:,:), FMM_hold, FLM_hold, FQN_hold, FNN_hold
else
    read ( unit=lun, rec=rec ) u(iend_w, :, :), v(iend_w, :, :), w(iend_w, :, :) 
endif
if ( DEBUG ) write (*, *) sub // ' : read record ', rec
    
if ( recycle ) then
    if ( rec < l_blend ) then
        recp = nrec - l_blend + 1 + rec
        wgt = 0.5_rprec * ( 1.0_rprec -                              &
                            cos ( pi * real (rec, rprec) / l_blend ) )
            !--wgt = 0+ when rec = 1
            !  wgt = 1- when rec = l_blend
        read ( unit=lun, rec=recp ) u_tmp, v_tmp, w_tmp
        u(iend_w, :, :) = wgt * u(iend_w, :, :) + (1.0_rprec - wgt) * u_tmp
        v(iend_w, :, :) = wgt * v(iend_w, :, :) + (1.0_rprec - wgt) * v_tmp
        w(iend_w, :, :) = wgt * w(iend_w, :, :) + (1.0_rprec - wgt) * w_tmp
    end if
end if

if ( DEBUG ) then  !--write out slices as an ascii time series
    if ( .not. init_DEBUG ) then
        inquire ( unit=lun_DEBUG, exist=exst, opened=opn )
        if ( exst .and. (.not. opn) ) then
            if ( USE_MPI ) then
                open ( unit=lun_DEBUG, file=debug_file // MPI_suffix )
            else
                open ( unit=lun_DEBUG, file=debug_file )
            end if
        
            write ( lun_DEBUG, '(a)' ) 'variables = "y" "z" "t" "u" "v" "w"'
            write ( lun_DEBUG, '(3(a,i0))' ) 'zone, f=point, i= ', ny,  &
                                             ', j= ', nz,               &
                                             ', k= ', nsteps
        else
            write (*, *) sub // ': problem opening debug file'
            stop
        end if
        init_DEBUG = .true.
    end if

    fmt = '(3(1x,i0),3(1x,es12.5))'
    do jz = 1, nz
        do jy = 1, ny
            write ( lun_DEBUG, fmt ) jy, jz, jt_total, u(iend_w, jy, jz),  &
                                     v(iend_w, jy, jz), w(iend_w, jy, jz)
        end do
    end do
end if

end subroutine inflow_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--finds number of records on existing direct-access unformatted file
!--taken from Clive Page's comp.lang.fortran posting (12/16/2003), 
!  under the topic counting number of records in a Fortran direct file
!--minor changes/renaming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine len_da_file(fname, lenrec, length)
implicit none
character (*), intent(in) :: fname  ! name of existing direct-access file
integer, intent(in)       :: lenrec ! record length (O/S dependent units)
integer, intent(out) :: length      ! number of records.
!
character (1) :: cdummy
integer :: lunit, nlo, nhi, mid, kode
logical :: exists, open
!
! find a free unit on which to open the file
!
do lunit = 99, 1, -1
  !--units to skip (compiler dependent)
  select case (lunit)
    case (5:6)
      !--do nothing
    case default
      inquire(unit=lunit, exist=exists, opened=open)
      if(exists .and. .not. open) exit
  end select
end do
open(unit=lunit, file=fname, access="direct", recl=lenrec, iostat=kode)
if(kode /= 0) then
  print *, 'error in len_da_file: ', trim(fname), ' does not exist'
  return
end if
!
! expansion phase
!
mid = 1
do
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode /= 0) exit
  mid = 2 * mid
end do
!
! length is between mid/2 and mid, do binary search to refine
!
nlo = mid/2
nhi = mid
do while(nhi - nlo > 1)
  mid = (nlo + nhi) / 2
  read(lunit, rec=mid, iostat=kode) cdummy
  if(kode == 0) then
     nlo = mid
  else
     nhi = mid
  end if
end do
length = nlo
close(unit=lunit)
return
end subroutine len_da_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module io
