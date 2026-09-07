!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2024 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!!
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!!
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!!
!! Authors:
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================

!> Module containing the subroutines to perform probes data collection
!! during postprocessing
module mod_post_probes

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi, one_4pi, ascii_real

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning  , info, printout, dust_time, t_realtime, new_file_unit

use mod_geometry, only: &
  t_geo_component, destroy_elements

use mod_parse, only: &
  t_parse, &
  getstr, getrealarray, &
  countoption

use mod_hdf5_io, only: &
  h5loc, &
  open_hdf5_file, &
  close_hdf5_file , &
  open_hdf5_group, &
  close_hdf5_group, &
  read_hdf5

use mod_stringtools, only: &
  LowCase

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_wake, only: &
  t_wake

use mod_actuatordisk, only: &
  t_actdisk

use mod_post_load, only: &
  load_refs, load_res, load_wake_post

use mod_tecplot_out, only: &
  tec_out_probes

use mod_dat_out, only: &
  dat_out_probes_header

use mod_wind, only: &
  variable_wind
  
use mod_surfpan, only: &
  t_surfpan

implicit none

public :: post_probes

private

character(len=max_char_len), parameter :: this_mod_name = 'mod_post_probes'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_probes( sbprms , basename , data_basename , an_name , ia , &
                      out_frmt , components_names , all_comp , &
                      an_start , an_end , an_step )
  type(t_parse), pointer                                   :: sbprms
  character(len=*) , intent(in)                            :: basename
  character(len=*) , intent(in)                            :: data_basename
  character(len=*) , intent(in)                            :: an_name
  integer          , intent(in)                            :: ia
  character(len=*) , intent(in)                            :: out_frmt
  character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
  logical , intent(inout)                                  :: all_comp
  integer , intent(in)                                     :: an_start , an_end , an_step

  type(t_geo_component), allocatable                       :: comps(:), comps_old(:)
  integer                                                  :: nstep
  real(wp), allocatable                                    :: probe_vars(:,:,:)
  real(wp), allocatable                                    :: time(:)
  character(len=max_char_len), allocatable                 :: probe_var_names(:)
  character(len=max_char_len), allocatable                 :: probe_loc_names(:)
  character(len=max_char_len)                              :: in_type , str_a , filename_in , var_name , vars_str
  integer                                                  :: n_probes , n_vars 
  real(wp), allocatable                                    :: rr_probes(:,:)
  logical                                                  :: probe_vel , probe_p , probe_vort, probe_cp
  integer                                                  :: fid_out , i_var , nprint
  integer                                                  :: ie , ip , ic , it , ires, nelems_comp
  integer(h5loc)                                           :: floc , ploc
  character(len=max_char_len)                              :: filename
  real(wp), allocatable                                    :: points(:,:), points_old(:,:)
  real(wp), allocatable                                    :: points_virtual(:,:), points_virtual_old(:,:)
  integer                                                  :: nelem, nelem_virtual

  real(wp)                                                 :: u_inf(3)
  real(wp)                                                 :: P_inf , rho
  real(wp)                                                 :: vel_probe(3) = 0.0_wp , vort_probe(3) = 0.0_wp
  real(wp)                                                 :: v(3) = 0.0_wp, phi = 0.0_wp
  real(wp), allocatable , target                           :: sol(:)
  real(wp)                                                 :: pres_probe, pot_probe = 0.0_wp, pot_probe_old = 0.0_wp
  real(wp)                                                 :: t, t_old

  real(wp), allocatable                                    :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                    :: vort(:), cp(:), vort_old(:), cp_old(:)

  type(t_wake)                                             :: wake, wake_old
  type(t_elem_p), allocatable                              :: wake_elems(:), wake_elems_old(:)
  integer                                                  :: ivar, ierr
  
  logical                                                  :: ex
  real(wp)                                                 :: dummy

  character(len=max_char_len), parameter                   :: this_sub_name = 'post_probes'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' probes'//nl
  call printout(trim(msg))

  ! Select all the components
  !   components_names is allocated in load_components_postpro()
  !   and deallocated in dust_post at the end of each analysis
  if ( allocated(components_names) ) then
    call warning(trim(this_sub_name), trim(this_mod_name), &
        'All the components are used. <Components> input &
        &is ignored, and deallocated.' )
    deallocate(components_names)
  end if
  all_comp = .true.

  ! Read probe coordinates: point_list or from_file
  in_type =  getstr(sbprms,'input_type')
  select case ( trim(in_type) )
    case('point_list')
      n_probes = countoption(sbprms,'point')
      allocate(rr_probes(3,n_probes))
      do ip = 1 , n_probes
        rr_probes(:,ip) = getrealarray(sbprms,'point',3)
      end do
    case('from_file')
      filename_in = getstr(sbprms,'file')   ! N probes and then their coordinates
      fid_out = 21
      open(unit=fid_out,file=trim(filename_in))
      read(fid_out,*) n_probes
      allocate(rr_probes(3,n_probes))
      do ip = 1 , n_probes
        read(fid_out,*) rr_probes(:,ip)
      end do
      close(fid_out)
    case default
      write(str_a,*) ia
      call error('dust_post','','Unknown InputType: '//trim(in_type)//&
                  ' for analysis n.'//trim(str_a)//'.'//nl//&
                  'It must either "point_list" or "from_file".')
  end select

  ! Read variables to save : velocity | pressure | vorticity
  ! TODO: add Cp
  probe_vel = .false. ; probe_p = .false. ; probe_vort = .false. ; probe_cp = .false.
  n_vars = countoption(sbprms,'variable')

  if ( n_vars .eq. 0 ) then ! default: velocity | pressure | vorticity
    probe_vel = .true. ; probe_p = .true. ; probe_vort = .true.
  else
    do i_var = 1 , n_vars
      var_name = getstr(sbprms,'variable') ; call LowCase(var_name)
      select case(trim(var_name))
        case ( 'velocity' ) ; probe_vel = .true.
        case ( 'pressure' ) ; probe_p   = .true.
        !case ( 'vorticity') ; probe_vort= .true.
        case ( 'cp'       ) ; probe_cp = .true.
        case ( 'all') ; probe_vel = .true. ; probe_p   = .true. ; !probe_vort= .true.
        case default
        write(str_a,*) ia
        call error('dust_post','','Unknown Variable: '//trim(var_name)//&
                    ' for analysis n.'//trim(str_a)//'.'//nl//&
                    'Choose "velocity", "pressure" .') !, "vorticity".')
      end select
    end do
  end if
  
  if (probe_cp) probe_p = .true. ! if pressure coefficient is requested, compute pressure and adimens.
  
  ! Find the number of fields to be plotted ( vec{vel} , p , vec{omega})
  nprint = 0
  if(probe_vel)  nprint = nprint+3
  if(probe_p)    nprint = nprint+1
  if(probe_cp)   nprint = nprint+1

  ! Write .dat file header
  allocate(probe_var_names(nprint))
  allocate(probe_loc_names(n_probes))
  probe_var_names = ''
  i_var = 1
  if(probe_vel) then
    probe_var_names(i_var)     = 'ux'
    probe_var_names(i_var + 1) = 'uy'
    probe_var_names(i_var + 2) = 'uz'
    i_var = i_var + 3
  endif
  if(probe_p) then
    probe_var_names(i_var) = 'p'
    i_var = i_var + 1
  endif 
  if(probe_cp) then
    probe_var_names(i_var) = 'cp'
  end if
  !if(probe_vort) then
  !  probe_var_names(i_var)     = 'omx'
  !  probe_var_names(i_var + 1) = 'omy'
  !  probe_var_names(i_var + 2) = 'omz'
  !  i_var = i_var + 3
  !endif
  vars_str = ''
  do i_var = 1,size(probe_var_names)
    vars_str = trim(vars_str)//'  '//trim(probe_var_names(i_var))
  enddo

  do ip = 1,n_probes
    write(probe_loc_names(ip),'(A,F10.5,A,F10.5,A,F10.5)') &
        'x=',rr_probes(1,ip), &
        'y=',rr_probes(2,ip), &
        'z=',rr_probes(3,ip)
  enddo

  !Allocate where the solution will be stored
  nstep = (an_end-an_start)/an_step + 1
  allocate(probe_vars(nprint, nstep, n_probes))
  allocate(time(nstep))


  ! load the geo components just once
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  !> TODO: here get the run id
  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual,  floc, &
                                components_names,  all_comp)
  call close_hdf5_file(floc)

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  if (probe_p) then ! geo for previous timestep (should be the same, but to avoid issues with allocatables...)
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !> TODO: here get the run id
    call load_components_postpro(comps_old, points_old, points_virtual_old, nelem, nelem_virtual, floc, &
                                components_names,  all_comp)  
    call close_hdf5_file(floc)
    call prepare_geometry_postpro(comps_old)    
  end if
  
  ! Allocate and point to sol
  ! TODO why?
  allocate(sol(nelem)) ; sol = 0.0_wp
  ip = 0
  do ic = 1 , size(comps)
    do ie = 1 , size(comps(ic)%el)
      ip = ip + 1
      comps(ic)%el(ie)%mag => sol(ip)
    end do
  end do

  !time history
  ires = 0
  !it_old = -1
  do it =an_start, an_end, an_step
    ires = ires+1
    ! Open the result file ----------------------
    if (probe_p) then ! load previous result for unsteady contribution
      ! check with wake allocation between steps
      !if (it .eq. it_old+1) then ! do not reload and update, use data from previous iteration
      !  comps_old = comps
      !  points_old = points
      !  wake_old = wake
      !  wake_elems_old = wake_elems
      !  vort_old = vort
      !  cp_old = cp
      !  t_old = t
      !else 
      
      ! load and update previous timestep
        write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it-1,'.h5'
        inquire(file=filename, exist=ex)
        if (.not. ex) then
          call error(this_mod_name,this_sub_name,'Output data '//filename//&
                    &' from the previous step is not available, cannot compute&
                    & unsteady contribution to pressure.')
        end if
        call open_hdf5_file(trim(filename),floc)
        ! Load the references and move the points ---
        call load_refs(floc,refs_R,refs_off)
        
        call update_points_postpro(comps_old, points_old, points_virtual_old, refs_R, refs_off, &
                                filen = trim(filename) )
        ! Load the results --------------------------
        call load_res(floc, comps_old, vort_old, cp_old, t_old)

        call load_wake_post(floc, wake_old, wake_elems_old)
        call close_hdf5_file(floc)
        
        ! Re-compute u_vort from loaded old wake
        do ic = 1,size(comps_old)
          nelems_comp = comps_old(ic)%nelems
        !$omp parallel do private(ie) 
          do ie = 1, nelems_comp
            select type(el =>comps_old(ic)%el(ie))
              type is(t_surfpan)
                call el%get_vort_vel(wake_old%vort_p)
            end select
          end do
        !$omp end parallel do
        end do
      !end if
    end if    
    
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
    call open_hdf5_file(trim(filename),floc)
    
    ! Load u_inf --------------------------------
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(u_inf,'u_inf',ploc)
    call read_hdf5(P_inf,'P_inf',ploc)
    call read_hdf5(rho,'rho_inf',ploc)
    call close_hdf5_group(ploc)

    sim_param%u_inf = u_inf
    sim_param%P_inf = P_inf
    sim_param%rho_inf = rho

    ! Load the references and move the points ---
    call load_refs(floc,refs_R,refs_off)
    call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                              filen = trim(filename) )
    ! Load the results --------------------------
    call load_res(floc, comps, vort, cp, t)

    call load_wake_post(floc, wake, wake_elems)
    call close_hdf5_file(floc)
    
    ! Re-compute u_vort from loaded wake
    do ic = 1,size(comps)
      nelems_comp = comps(ic)%nelems
      !$omp parallel do private(ie) 
      do ie = 1, nelems_comp
        select type(el =>comps(ic)%el(ie))
          type is(t_surfpan)
            call el%get_vort_vel(wake%vort_p)
        end select
      end do
      !$omp end parallel do
    end do

    time(ires) = t
    
    ! Compute velocity --------------------------
    do ip = 1 , n_probes ! probes

      vel_probe = 0.0_wp ; pres_probe = 0.0_wp ; vort_probe = 0.0_wp
      pot_probe = 0.0_wp ; pot_probe_old = 0.0_wp
      i_var = 1
      if ( probe_vel .or. probe_p ) then

      ! compute velocity
      do ic = 1,size(comps)
!$omp parallel do private(ie, v) reduction(+:vel_probe)
        do ie = 1 , size( comps(ic)%el )
          call comps(ic)%el(ie)%compute_vel( rr_probes(:,ip) , v )
          vel_probe = vel_probe + v*one_4pi

        end do
!$omp end parallel do
      end do

!$omp parallel do private( ie, v) reduction(+:vel_probe)
      do ie = 1, size(wake_elems)
        call wake_elems(ie)%p%compute_vel( rr_probes(:,ip) , v )
        vel_probe = vel_probe + v*one_4pi
      enddo
!$omp end parallel do

      vel_probe = vel_probe + variable_wind(rr_probes(:,ip), t)
    end if

    if(probe_vel) then
      probe_vars(i_var:i_var+2, ires, ip) = vel_probe
      i_var = i_var+3
    endif

    ! compute pressure
    if ( probe_p .or. probe_cp ) then
      ! Bernoulli equation
      ! rho * dphi/dt + P + 0.5*rho*V^2 = P_infty + 0.5*rho*V_infty^2
  
      ! compute current phi
      do ic = 1,size(comps)
!$omp parallel do private(ie, phi) reduction(+:pot_probe)
        do ie = 1 , size( comps(ic)%el )
          call comps(ic)%el(ie)%compute_pot(phi, dummy, rr_probes(:,ip), 1, 2)
          pot_probe = pot_probe + phi
        end do
!$omp end parallel do
      end do
      
! TODO check parallel with select type
!$omp parallel do private( ie, phi) reduction(+:pot_probe)
      do ie = 1, size(wake_elems)
        select type(el => wake_elems(ie)%p)
          class is (c_pot_elem)
            call el%compute_pot(phi, dummy, rr_probes(:,ip), 1, 2)
            pot_probe = pot_probe + phi
        end select  
      enddo
!$omp end parallel do
    
      ! compute previous phi
      do ic = 1,size(comps_old)
!$omp parallel do private(ie, phi) reduction(+:pot_probe_old)
        do ie = 1 , size( comps_old(ic)%el )
          call comps_old(ic)%el(ie)%compute_pot(phi, dummy, rr_probes(:,ip), 1, 2)
          pot_probe_old = pot_probe_old + phi
        end do
!$omp end parallel do
      end do
!$omp parallel do private( ie, phi) reduction(+:pot_probe_old)
      do ie = 1, size(wake_elems_old)
        select type(el => wake_elems_old(ie)%p)
          class is (c_pot_elem)
            call el%compute_pot(phi, dummy, rr_probes(:,ip), 1, 2)
            pot_probe_old = pot_probe_old + phi
        end select  
      enddo
!$omp end parallel do

      pres_probe = P_inf + 0.5_wp*rho*norm2(u_inf)**2 - 0.5_wp*rho*norm2(vel_probe)**2 -&
                  & rho*(pot_probe-pot_probe_old)/(t-t_old)

      if (probe_p) then ! output pressure coefficient
        probe_vars(i_var, ires, ip) = pres_probe
        i_var = i_var + 1
      endif 
      
      if (probe_cp) then 
        if (norm2(u_inf) .gt. 1e-6_wp) then
          probe_vars(i_var, ires, ip) = (pres_probe - P_inf)/(0.5_wp*rho*norm2(u_inf)**2)
        else ! todo add u_ref
          call error(this_mod_name,this_sub_name,'Pressure coefficient requested, but u_inf =0')
        end if
      end if
    end if

    ! compute vorticity
    !if ( probe_vort ) then
    !  vort_probe = 0.0_wp
    !  probe_vars(i_var:i_var+2, ires, ip) = vort_probe
    !  i_var = i_var+3
    !end if

  end do  ! probes
  
    !it_old = it
    
    ! TODO check if wake needs to be destroyed/deallocated between iterations
  end do

  !Output the results in the correct format
  select case (trim(out_frmt))

    case ('dat')
      call new_file_unit(fid_out, ierr)
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.dat'
      open(unit=fid_out,file=trim(filename))
      call dat_out_probes_header( fid_out , rr_probes , vars_str, nstep )

      do ires = 1, size(time)
        write(fid_out,'('//ascii_real//')',advance='no') time(ires)
        do ip = 1, n_probes
          do ivar = 1, size(probe_vars,1)
            write(fid_out,'('//ascii_real//')',advance='no') probe_vars(ivar,ires,ip)
            write(fid_out,'(A)',advance='no') ' '
          enddo
        enddo
        write(fid_out,*) ' '
      enddo
      close(fid_out)

    case ('tecplot')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
      call tec_out_probes(filename, time, probe_vars, probe_var_names, probe_loc_names)
    case default
      call error('dust_post','','Unknown format '//trim(out_frmt)//&
                ' for probe output')
    end select

  do ic = 1 , size(comps)
    do ie = 1 , size(comps(ic)%el)
      ip = ip + 1
      comps(ic)%el(ie)%mag => null()
    end do
  end do
  deallocate(sol)

  deallocate(probe_vars,time)
  deallocate(probe_var_names,probe_loc_names)
  deallocate(rr_probes)

  call destroy_elements(comps)
  deallocate(comps, points,components_names)
  
  if (allocated(comps_old)) then 
    call destroy_elements(comps_old)
    deallocate(comps_old)
  end if
  if (allocated(points_old)) deallocate(points_old)

  write(msg,'(A,I0,A)') nl//'++++++++++ Probes done'//nl
  call printout(trim(msg))

end subroutine post_probes

! ----------------------------------------------------------------------

end module mod_post_probes
