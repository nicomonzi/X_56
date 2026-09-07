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
!!          Andrea Colli
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================

!> This is the main file of the DUST solver

program dust

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_sim_param, only: &
  t_sim_param, sim_param, create_param_main, init_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_basename, &
  check_file_exists

use mod_geometry, only: &
  t_geo, &
  create_geometry, update_geometry, &
  t_tedge,  destroy_geometry, destroy_elements

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_doublet, only: &
  initialize_doublet

use mod_surfpan, only: &
  t_surfpan , initialize_surfpan

use mod_vortlatt, only: &
  t_vortlatt 

use mod_stripe, only: &
  t_stripe

use mod_nlvl_solver, only: &
  solve_nlvl_correction

use mod_liftlin, only: &
  update_liftlin, t_liftlin_p, &
  build_ll_kernel, &
  solve_liftlin

use mod_actuatordisk, only: &
  update_actdisk

use mod_virtual, only: &
  c_elem_virtual, t_elem_virtual_p

use mod_vortline, only: &
  initialize_vortline

use mod_vortpart, only: &
  initialize_vortpart

use mod_c81, only: &
  t_aero_tab

use mod_linsys_vars, only: &
  t_linsys

use mod_linsys, only: &
  initialize_linsys, assemble_linsys, solve_linsys, destroy_linsys, &
  factorize_kutta, dump_linsys 

use mod_pressure_equation, only: &
  dump_linsys_pres, press_normvel_der, initialize_pressure_sys, &
  assemble_pressure_sys, solve_pressure_sys

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

use mod_wake, only: &
  t_wake, initialize_wake, update_wake, &
  prepare_wake, load_wake, complete_wake, destroy_wake

use mod_vtk_out, only: &
  vtk_out_bin

use mod_tecplot_out, only: &
  tec_out_sol_bin

use mod_hdf5_io, only: &
  h5loc, initialize_hdf5, destroy_hdf5, new_hdf5_file, open_hdf5_file, &
  close_hdf5_file, new_hdf5_group, open_hdf5_group, close_hdf5_group, &
  write_hdf5, write_hdf5_attr, read_hdf5, read_hdf5_al, append_hdf5

use mod_dust_io, only: &
  save_status, load_solution

use mod_viscosity, only: &
  viscosity_effects

use mod_octree, only: &
  initialize_octree, destroy_octree, sort_particles, t_octree, &
  apply_multipole_panels

use mod_math, only: & 
  cross, dot, vec2mat, invmat_banded, linear_interp

#if USE_PRECICE
  use mod_precice, only: &
    t_precice
#endif

use mod_wind, only: &
  variable_wind

implicit none

!> Run-id
integer                           :: run_id(10)

!> Input
!> Main parameters parser
type(t_parse)                     :: prms
character(len=*), parameter       :: input_file_name_def = 'dust.in'
character(len=max_char_len)       :: input_file_name
character(len=max_char_len)       :: target_file
character(len=extended_char_len)  :: message

!> Time parameters
real(wp)                          :: time
integer                           :: it, nstep, nout
real(wp)                          :: t_last_out, t_last_debug_out
real(wp)                          :: time_no_out, time_no_out_debug
logical                           :: time_2_out, time_2_debug_out
logical                           :: output_start
real(wp)                          :: dt_debug_out
logical                           :: already_solv_restart

!> Main variables
!> All the implicit elements, sorted first static then moving
type(t_impl_elem_p), allocatable  :: elems(:)
!> All the explicit elements
type(t_expl_elem_p), allocatable  :: elems_expl(:)
!> Only the lifting line elements
type(t_liftlin_p), allocatable    :: elems_ll(:)
!> Only the actuator disk elements
type(t_expl_elem_p), allocatable  :: elems_ad(:)
!> All the elements (panels+ll)
type(t_pot_elem_p), allocatable   :: elems_tot(:)
!> All the non corrected elements (panels+vl+ll+a-vl_nl)
type(t_pot_elem_p), allocatable   :: elems_non_corr(:)
!> All the corrected vortex lattice elements 
type(t_pot_elem_p), allocatable   :: elems_corr(:) 
!> Elem virtual 
type(t_elem_virtual_p), allocatable :: elems_virtual(:)

!> Geometry
type(t_geo)                       :: geo
!> Trailing edge
type(t_tedge)                     :: te
!> Airfoil table data
type(t_aero_tab), allocatable     :: airfoil_data(:)
!> Linear system
type(t_linsys)                    :: linsys
!> Wake
type(t_wake)                      :: wake
!> Timing vars
real(t_realtime)                  :: t1 , t0, t00, t11, t22
!> I/O prefixes
character(len=max_char_len)       :: frmt, frmt_vl
character(len=max_char_len)       :: basename_debug

real(wp), allocatable             :: res_old(:)
real(wp), allocatable             :: surf_vel_SurfPan_old(:,:)
real(wp), allocatable             :: nor_SurfPan_old(:,:)
real(wp), allocatable             :: al_kernel(:,:), al_v(:)
real(wp)                          :: theta_cen(3), R_cen(3, 3) 

!> Kutta correction 
real(wp), allocatable             :: delta_pres_te(:), delta_pres_te_perturbed(:) 
real(wp), allocatable             :: delta_pres_te_iter(:), delta_pres_te_iter_old(:)
real(wp), allocatable             :: mag_pert(:,:) 
real(wp), allocatable             :: delta_mag_te(:), delta_mag_te_old(:) 
real(wp), allocatable             :: delta_mag_te_iter(:)
real(wp), allocatable             :: delta_mag_te_lower_old(:), delta_mag_te_upper_old(:)
real(wp), allocatable             :: mag_te_tmp_lower(:), mag_te_tmp_upper(:)
real(wp), allocatable             :: rhs_tmp(:), res_tmp(:) 
real(wp), allocatable             :: A_tmp(:,:)
real(wp), allocatable             :: jacobi(:,:)
integer                           :: it_pan, n_pan_te

!> VL viscous correction
integer                           :: i_el, j_el, k_el, i_c, i_s, i, sel, i_p, i_c2, i_s2, ih
integer                           :: it_vl, it_stall
real(wp)                          :: tol, diff, max_diff 
real(wp)                          :: d_cd(3), vel(3), v(3), a_v, area_stripe, dforce_stripe(3)
real(wp)                          :: nor(3), tang_cen(3), u_v, q_inf

!> relaxation 
real(wp), allocatable             :: residual_vl(:), residual_vl_old(:), residual_vl_delta(:), gamma_tmp(:)
real(wp)                          :: rel_aitken, den_aitken

!> octree parameters
type(t_octree)                    :: octree

#if USE_PRECICE
  !> PreCICE data 
  type(t_precice)                 :: precice
  integer                         :: bool
  logical                         :: precice_convergence
  integer                         :: j
  real(wp)                        :: sum_force(3)
  integer                         :: it_precice = 0
#endif


call printout(nl//'>>>>>> DUST beginning >>>>>>'//nl)

t00 = dust_time()

call get_run_id(run_id)

!> Modules initialization 
call initialize_hdf5()

!> Input reading 
if(command_argument_count().gt.0) then
  call get_command_argument(1, value = input_file_name)
else
  input_file_name = input_file_name_def
endif

call printout(nl//'Reading input parameters from file "'//&
                trim(input_file_name)//'"'//nl)

call create_param_main(prms)

!> Get the parameters and print them out
call printout(nl//'====== Input parameters: ======')
call check_file_exists(input_file_name,'dust main')
call prms%read_options(input_file_name, printout_val=.true.)

!> Initialize all the parameters
nout = 0  !> Reset the numbering for output files
output_start = getlogical(prms, 'output_start')
call init_sim_param(sim_param, prms, nout, output_start)

!> Remaining parameters
if(countoption(prms, 'dt_debug_out') .lt. 1) then
  dt_debug_out = sim_param%dt_out
else
  dt_debug_out = getreal(prms, 'dt_debug_out')
endif

if(countoption(prms, 'basename_debug') .lt. 1) then
  basename_debug = sim_param%basename
else
  basename_debug = getstr(prms,'basename_debug')
endif

sim_param%basename_debug = basename_debug

#if USE_PRECICE
!> Create PreCICE dust participant
! Do it here because it needs to read precice_config path from dust.in
call precice%initialize()
#endif


!> Parameters Initializations 
call initialize_doublet()
call initialize_vortline()
call initialize_surfpan()

!> Check that tend .gt. tinit
if ( sim_param%tend .le. sim_param%t0 ) then
  write(message,*) 'The end time of the simulation',sim_param%tend,'is&
                  & lower or equal to the start time',sim_param%t0,'.&
                  & Remembver that when restarting without resetting the&
                  & time the start time is taken from the restart result!'
  call error('dust','',message)
end if


!> Printout the parameters
if (sim_param%debug_level .ge. 3) then
  write(message,*) 'Initial time tstart: ', sim_param%t0;             call printout(message)
  write(message,*) 'Final time tend:     ', sim_param%tend;           call printout(message)
  write(message,*) 'Time step dt:        ', sim_param%dt;             call printout(message)
  write(message,*) 'Output interval:     ', sim_param%dt_out;         call printout(message)
  write(message,*) 'Output first step:   ', output_start;             call printout(message)
  write(message,*) 'Debug level:         ', sim_param%debug_level;    call printout(message)
  write(message,*) 'Free stream velocity:', sim_param%u_inf;          call printout(message)
  write(message,*) 'Maximum wake panels: ', sim_param%n_wake_panels;  call printout(message)
  write(message,*) 'Results basename:    ', trim(sim_param%basename); call printout(message)
  write(message,*) 'Debug basename:      ', trim(basename_debug);     call printout(message)
endif


!> Check that the basenames are valid 
call check_basename(trim(sim_param%basename),'dust main')
if(sim_param%debug_level.ge.10)  call check_basename(trim(basename_debug),'dust main')

!> Simulation parameters 
nstep = sim_param%n_timesteps
allocate(sim_param%time_vec(sim_param%n_timesteps))
sim_param%time_vec = (/ ( sim_param%t0 + &
          real(i-1,wp)*sim_param%dt, i=1,sim_param%n_timesteps ) /)

!> Geometry creation 
call printout(nl//'====== Geometry Creation ======')

t0 = dust_time()
target_file = trim(sim_param%basename)//'_geo.h5'

call create_geometry(sim_param%GeometryFile, sim_param%ReferenceFile, &
                    input_file_name, geo, te, elems, elems_expl, elems_ad, &
                    elems_ll, elems_virtual, elems_corr, elems_non_corr, elems_tot, airfoil_data, target_file, run_id)

t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Created geometry in: ' , t1 - t0,' s.'
  call printout(message)
endif

if(sim_param%debug_level .ge. 15) &
  call debug_printout_geometry_minimal(elems, geo, basename_debug, 0)
if(sim_param%debug_level .ge. 15) &
  call debug_ll_printout_geometry(elems_ll, geo, basename_debug, 0)

!> TODO: check whether to move these calls before, and precisely what they do
if(sim_param%debug_level .ge. 7) call ignoredParameters(prms)
call finalizeParameters(prms)

!> Initialize PreCICE mesh and fields 
#if USE_PRECICE
  te%t_hinged = te%t
  call precice%initialize_mesh( geo )
  call precice%initialize_fields()
  call precicef_initialize()
  call precice%update_elems(geo, elems_tot, elems_virtual, te ) 
#endif

!> Initialization 
if(sim_param%use_fmm) then
  call printout(nl//'====== Initializing Octree ======')
  t0 = dust_time()
  call initialize_octree(sim_param%BoxLength, sim_param%NBox, &
                        sim_param%OctreeOrigin, sim_param%NOctreeLevels, &
                        sim_param%MinOctreePart, sim_param%MultipoleDegree, &
                        sim_param%octree_vortex_rad, octree)
  t1 = dust_time()
  if(sim_param%debug_level .ge. 1) then
    write(message,'(A,F9.3,A)') 'Initialized octree in: ' , t1 - t0,' s.'
    call printout(message)
  endif
endif

call printout(nl//'====== Initializing Wake ======')

call initialize_wake(wake, geo, te, sim_param%n_wake_panels, &
      sim_param%n_wake_panels, sim_param%n_wake_particles)

call printout(nl//'====== Initializing Linear System ======')
t0 = dust_time()
call initialize_linsys(linsys, geo, elems, elems_expl, wake ) 

t1 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Initialized linear system in: ' , t1 - t0,' s.'
  call printout(message)
endif

!> Restart 
if (sim_param%restart_from_file) then
  call load_solution(sim_param%restart_file, geo%components, geo%refs)
  call load_wake(sim_param%restart_file, wake, elems_tot)
  already_solv_restart = .true.
else ! Set to zero the intensity of all the singularities

  do i_el = 1 ,size(elems)      ! implicit elements (vr, sp)
      elems(i_el)%p%mag = 0.0_wp
  end do
  
  do i_el = 1 , size(elems_expl) ! explicit elements (ll, ad)
    elems_expl(i_el)%p%mag = 0.0_wp
  end do
  if (size(elems_ll) .gt. 0) then
    do i_el = 1, size(elems_ll)
      elems_ll(i_el)%p%Gamma_old = 0.0_wp
    end do
  endif
  already_solv_restart = .false.
endif

!> Keep the first-step potential solution path unchanged: geometry update
!> is handled inside the normal timestep flow.

t22 = dust_time()
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') nl//'------ Completed all preliminary operations &
                              &in: ' , t22 - t00,' s.'
  call printout(message)
endif

!> Build kernel for regularisation (averaging) process for LL 
if (( size(elems_ll) .gt. 0 )) then
  allocate(al_v(size(elems_ll))) 
  al_v = 0.0_wp  
  call build_ll_kernel(elems_ll, al_v, al_kernel)
  deallocate(al_v)
end if

!> Initialize coupling 
#if USE_PRECICE
  call precicef_is_coupling_ongoing(precice%is_ongoing)  
  write(*,*) ' is coupling ongoing: ', precice%is_ongoing
  precice_convergence = .true.
  call precicef_get_max_time_step_size(precice%dt_precice)
  precice%dt_precice = min(precice%dt_precice, sim_param%dt)
  write(*,*) ' dt_precice         : ', precice%dt_precice
  
  !> Before entering the time cycle we need to actually initialize the
  ! position of the elements, so we have to query mbdyn
  !> Read data from structural solver
  
  
  do i = 1, size(precice%fields)
    if ( trim(precice%fields(i)%fio) .eq. 'read' ) then
      if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
        call precicef_read_data(  precice%mesh%mesh_name, &
                                  precice%fields(i)%fname, &
                                  precice%mesh%nnodes  , &
                                  precice%mesh%node_ids, &
                                  precice%dt_precice, &
                                  precice%fields(i)%fdata(1,:) )
      elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
          call precicef_read_data(  precice%mesh%mesh_name, &
                                    precice%fields(i)%fname, &
                                    precice%mesh%nnodes  , &
                                    precice%mesh%node_ids, &
                                    precice%dt_precice, &
                                    precice%fields(i)%fdata )
      endif
    end if
  end do

  !> Update dust geometry ( elems and first wake panels )
  call precice%update_elems( geo, elems_tot, elems_virtual, te )

  !> Update geo_data()
  do i_el = 1, size(elems_tot)
    call elems_tot(i_el)%p%calc_geo_data( &
                          geo%points(:,elems_tot(i_el)%p%i_ver) )
  end do

  !> Update geo_data()
  do i_el = 1, size(elems_virtual)
    call elems_virtual(i_el)%p%calc_geo_data_virtual( &
              geo%points_virtual(:,elems_virtual(i_el)%p%i_ver) )
  end do

  !> Update near-field wake
  call precice%update_near_field_wake( geo, wake, te )
  
  !> Store the second row of the wake
  wake%old_second_row = wake%pan_w_points(:,:,2)
#endif


!=========================== Time Cycle ==============================
!> General overview:
!> - build and solve systems
!> - compute loads
!> - save data
!> - update and prepare for next step

call printout(nl//'////////// Performing Computations //////////')
time = sim_param%t0
sim_param%time_old = sim_param%t0 + 1
t_last_out = time
t_last_debug_out = time
time_no_out = 0.0_wp
time_no_out_debug = 0.0_wp

if (sim_param%kutta_correction) then
  n_pan_te = size(wake%pan_gen_elems_id,2) 
  allocate(delta_pres_te(n_pan_te)); delta_pres_te = 0.0_wp
  allocate(delta_pres_te_perturbed(n_pan_te)); delta_pres_te_perturbed = 0.0_wp
  allocate(delta_pres_te_iter(n_pan_te)); delta_pres_te_iter = 0.0_wp
  allocate(delta_pres_te_iter_old(n_pan_te)); delta_pres_te_iter_old = 0.0_wp

  allocate(delta_mag_te(n_pan_te)); delta_mag_te = 0.0_wp
  allocate(delta_mag_te_old(n_pan_te)); delta_mag_te_old = 0.0_wp
  allocate(mag_te_tmp_upper(n_pan_te)); mag_te_tmp_upper = 0.0_wp
  allocate(mag_te_tmp_lower(n_pan_te)); mag_te_tmp_lower = 0.0_wp
  allocate(delta_mag_te_upper_old(n_pan_te)); delta_mag_te_upper_old = 0.0_wp
  allocate(delta_mag_te_lower_old(n_pan_te)); delta_mag_te_lower_old = 0.0_wp
  allocate(delta_mag_te_iter(n_pan_te)); delta_mag_te_iter = 0.0_wp

  allocate(mag_pert(n_pan_te,n_pan_te)); mag_pert = 1.0_wp  
  allocate(rhs_tmp(size(linsys%b))); rhs_tmp = 0.0_wp 
  allocate(res_tmp(size(linsys%b))); res_tmp = 0.0_wp
  allocate(A_tmp(size(linsys%A,1),size(linsys%A,1))); A_tmp = 0.0_wp
  allocate(jacobi(n_pan_te,n_pan_te)); jacobi = 0.0_wp
endif 

allocate(surf_vel_SurfPan_old(geo%nSurfpan,3)) ; surf_vel_SurfPan_old = 0.0_wp
allocate(     nor_SurfPan_old(geo%nSurfpan,3)) ;      nor_SurfPan_old = 0.0_wp
allocate(res_old(size(elems)))
res_old = 0.0_wp

!===========> Start time cycle 
t11 = dust_time()
it = 0
#if USE_PRECICE
it = 1
  do while ( ( it .lt. nstep ) .and. ( precice%is_ongoing .eq. 1 ) ) 
#else
  do while ( ( it .lt. nstep ) )
    it = it + 1
#endif

    sim_param%time_old = sim_param%time

    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,I5,A,I5,A,F9.4)') nl//'--> Step ',it,' of ', &
                                        nstep, ' simulation time: ', time
      call printout(message)
      t22 = dust_time()
      write(message,'(A,F9.3,A)') 'Elapsed wall time: ', t22 - t00
      call printout(message)
    endif

#if USE_PRECICE
    if ( precice_convergence ) then
#endif

    call init_timestep(time)

#if USE_PRECICE
      it_precice = it_precice + 1
      precice_convergence = .false.
    end if
#endif

    if (sim_param%time_varying_u_inf) then
      !Interpolation of Uinf from file on DUST time
      do i = 1, 3
        call linear_interp(sim_param%u_inf_comps(i,:), sim_param%u_inf_time, &
                          time, sim_param%u_inf(i)) 
      end do 
    end if

    if ( mod( it-1, sim_param%ndt_update_wake ) .eq. 0 ) then
      call prepare_wake(wake, elems_tot, octree)
    end if
    
    if(already_solv_restart) then
    ! prepare_wake zeros uvort for ALL elems_tot, so it has to be computed again for a restart.
    ! TODO check unwanted behaviour of prepare_wake and limit its access to elems_tot 
      do i_el = 1 , size(elems)
        select type( el => elems(i_el)%p ) ; type is (t_surfpan)
          call el%get_vort_vel(wake%vort_p)
        end select
      end do
    end if
    
    call update_liftlin(elems_ll,linsys)
    call update_actdisk(elems_ad,linsys)

    !> Debug geometry printing
    if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
              call debug_printout_geometry(elems, geo, basename_debug, it)
    if((sim_param%debug_level .ge. 16).and.time_2_debug_out)&
              call debug_ll_printout_geometry(elems_ll, geo, basename_debug, it)

#if USE_PRECICE

    call precicef_requires_writing_checkpoint(bool )
    if ( bool .eq. 1 ) then ! Save old state
      !> Save old state: forces and moments
      do j = 1, size(precice%fields)
        if ( trim(precice%fields(j)%fio) .eq. 'write' ) then
          precice%fields(j)%cdata = precice%fields(j)%fdata
        end if
      end do
    end if

    !> Read data from structural solver
    do i = 1, size(precice%fields)
      if ( trim(precice%fields(i)%fio) .eq. 'read' ) then
        if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
          call precicef_read_data(  precice%mesh%mesh_name, &
                                    precice%fields(i)%fname, &
                                    precice%mesh%nnodes  , &
                                    precice%mesh%node_ids, &
                                    precice%dt_precice, &
                                    precice%fields(i)%fdata(1,:) )
        elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
            call precicef_read_data(  precice%mesh%mesh_name, &
                                      precice%fields(i)%fname, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%dt_precice, &
                                      precice%fields(i)%fdata )
        endif
      end if
    end do

    !> Update dust geometry ( elems and first wake panels )
    call precice%update_elems( geo, elems_tot, elems_virtual, te )

    !> Update geo_data()
    do i_el = 1, size(elems_tot)
      call elems_tot(i_el)%p%calc_geo_data( &
                            geo%points(:,elems_tot(i_el)%p%i_ver) )
    end do

    do i_el = 1, size(elems_virtual)
      call elems_virtual(i_el)%p%calc_geo_data_virtual( &
                geo%points_virtual(:,elems_virtual(i_el)%p%i_ver) )
    end do

    !> Update near-field wake
    call precice%update_near_field_wake( geo, wake, te )

    !> Update dt--> mbdyn should take care of the dt and send it to precice (TODO)
#else
#endif
  do i_el = 1 , geo%nSurfPan
      select type ( el => elems(i_el)%p ) ; class is ( t_surfpan )
#if USE_PRECICE
      if (geo%components(el%comp_id)%coupling) then
        call el%create_local_velocity_stencil(el%R_cen)
        !> chtls stencil
        call el%create_chtls_stencil(el%R_cen)
      else 
        call el%create_local_velocity_stencil( &    
                geo%refs(geo%components(el%comp_id)%ref_id)%R_g)
        !> chtls stencil
        call el%create_chtls_stencil( &             
                geo%refs(geo%components(el%comp_id)%ref_id)%R_g)
      endif 
#else 
      call el%create_local_velocity_stencil( &    
              geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
      !> chtls stencil
      call el%create_chtls_stencil( &             
              geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
#endif 
      end select
  end do
    !> Calculate the normal velocity derivative for the pressure equation
    call press_normvel_der(geo, elems, surf_vel_SurfPan_old)
  
  !> If the simulation is restarted the solution part can be skipped 
  ! since it's loaded from the restart file
  ! TODO: not working for coupled simulations, so far
#if USE_PRECICE
#else
  if(.not. already_solv_restart) then
#endif    
    !>-------------- Assemble the system ------
    t0 = dust_time()
    sel = size(elems) ! total number of elements
    call assemble_linsys(linsys, geo, elems, elems_expl, wake)
    !call assemble_pressure_sys(linsys, geo, elems, wake)
    t1 = dust_time()
    
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Assembled linear system in: ' , t1 - t0,' s.'
      call printout(message)
    endif


    !debug output of the system
    if ((sim_param%debug_level .ge. 50).and.time_2_debug_out) then
      write(frmt,'(I4.4)') it
      call dump_linsys(linsys, trim(basename_debug)//'A_'//trim(frmt)//'.dat', &
                              trim(basename_debug)//'b_'//trim(frmt)//'.dat' )
      call dump_linsys_pres(linsys, &
                              trim(basename_debug)//'Apres_'//trim(frmt)//'.dat', &
                              trim(basename_debug)//'bpres_'//trim(frmt)//'.dat')
    endif

    if (sim_param%kutta_correction .and. geo%nSurfPan .gt. 0) then  
      call factorize_kutta(linsys) 
    end if 

    !------ Solve the system ------
    t0 = dust_time()
    if (linsys%rank .gt. 0) then
      call solve_linsys(linsys)
    endif
    t1 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)')  'Solved linear system in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    !> compute dGamma_dt for unsteady contribution
!$omp parallel do private(i_el)
    do i_el = 1 , sel
      elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
    end do
!$omp end parallel do


!debug print of the results
if (sim_param%debug_level .ge. 20 .and. time_2_debug_out) &
                      call debug_printout_result(linsys, basename_debug, it)

  !------ Update the explicit part ------  % v-----implicit elems: p,v
  if ( size(elems_ll) .gt. 0 ) then
    if ( trim(sim_param%llSolver) .eq. 'GammaMethod' ) then ! Gamma-method
      if (sim_param%time .gt. sim_param%time_old)  then
        do i_el = 1, size(elems_ll)
          elems_ll(i_el)%p%Gamma_old_old = elems_ll(i_el)%p%Gamma_old
          elems_ll(i_el)%p%Gamma_old = elems_ll(i_el)%p%mag
        enddo
      endif
      call solve_liftlin(elems_ll, elems_tot, elems , elems_ad , &
              (/ wake%pan_p, wake%rin_p/), wake%vort_p, airfoil_data, it)
  else
    call error('dust','dust',' Wrong string for LLsolver. &
          &This parameter should have been set equal to "GammaMethod" (default) &
          &in init_sim_param() routine. Something went wrong. Stop')
    end if
  end if


#if USE_PRECICE
  !$omp parallel do private(i_el)
    do i_el = 1, sel
      ! ifort bugs workaround:
      ! apparently it is not possible to call polymorphic methods inside
      ! select cases for intel, need to call these for all elements and for the
      ! vortex lattices it is going to be a dummy empty function call
      if (geo%components(elems(i_el)%p%comp_id)%coupling) then  
        !> calculate the pressure using the relative orientation matrix
        call elems(i_el)%p%compute_pres(elems(i_el)%p%R_cen)  ! update surf_vel field too
      else !> non coupled component 
          call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
              geo%refs(geo%components(elems(i_el)%p%comp_id)%ref_id)%R_g)
      endif  
      call elems(i_el)%p%compute_dforce()      
    end do
  !$omp end parallel do
#else
  
  !$omp parallel do private(i_el)
    do i_el = 1 , sel
      call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
              geo%refs(geo%components(elems(i_el)%p%comp_id)%ref_id)%R_g)
      call elems(i_el)%p%compute_dforce()  
    end do
  !$omp end parallel do
#endif

  ! ifort bugs workaround:
  ! since even if the following calls looks thread safe, they mess up with
  ! ifort and parallel runs, so the cycle is executed another time just for the
  ! vortex lattices
if ( geo%nVortLatt .gt. 0) then
  !$omp parallel do private(i_el)
    do i_el = 1 , sel      
      select type(el => elems(i_el)%p)        
        class is(t_vortlatt)    
          ! compute vel at 1/4 chord (some approx, see the comments in the fcn)      
          call el%get_vel_ctr_pt( elems_tot, (/ wake%pan_p, wake%rin_p/), wake%vort_p)
      end select
    end do 
  !$omp end parallel do

  do i_el = 1 , sel      
    select type(el => elems(i_el)%p)        
      class is(t_vortlatt)         
        !> compute dforce using AVL formula 
        call el%compute_dforce_jukowski(.true.) 
        !> update the pressure field, p = df.n / area
        el%pres = sum(el%dforce * el%nor)/el%area
      end select
  end do
end if 
! Unsteady Kutta Condition
  ! 1. get the delta p at the trailing edge for each panel element  
  ! 2. compute the jacobi matrix deltagamma = (1 + diag(-beta))*gamma 
  !    a. calculate new rhs using linear distribution at the trailing edge
  !    b. solve the linear system to get delta mu = A^-1 * (rhs - T_L*deltagamma) 
  !    c. compute the new p 
  !    update jacobian matrix = (p_tilde - p)/(- gamma*beta) 
  ! 3. start newton iteration to get the new gamma 
  
  if (geo%nSurfPan .gt. 0 .and. sim_param%kutta_correction .and. (it .gt. sim_param%kutta_startstep - 1)) then  
    !> save the temporary values
    rhs_tmp = linsys%b 
    A_tmp = linsys%A 
    res_tmp = linsys%res
    
    if  (it .gt. sim_param%kutta_startstep) then

      linsys%b = rhs_tmp - matmul(linsys%TR, delta_mag_te_old)

      !> Solve the factorized system
      if (linsys%rank .gt. 0) then ! A matrix includes the steady kutta condition
        linsys%skip = .true.
        call solve_linsys(linsys)
      end if
      
      !$omp parallel do private(k_el)
      do k_el = 1 , sel      
        elems(k_el)%p%didou_dt = (linsys%res(k_el) - res_old(k_el)) / sim_param%dt
      enddo 
      !$omp end parallel do

      !> compute delta pressure at the trailing edge 
      
#if USE_PRECICE
      !$omp parallel do private(j_el)
        do j_el = 1, n_pan_te
          if (associated(wake%pan_gen_elems(2,j_el)%p)) then 
            if (geo%components(te%e(1, j_el)%p%comp_id)%coupling) then  
              !> calculate the pressure using the relative orientation matrix
              !> upper side of the panel
              call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres(te%e(1, j_el)%p%R_cen) 
              !> lower side of the panel  
              call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres(te%e(2, j_el)%p%R_cen) 
            else !> non coupled component
              !> upper side of the panel 
              call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
              !> lower side of the panel
              call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
            endif 
          endif 
        end do 
      !$omp end parallel do      
#else
      !$omp parallel do private(j_el)
        do j_el = 1, n_pan_te
          if (associated(wake%pan_gen_elems(2,j_el)%p)) then 
            !> upper side of the panel
            call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                  geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
            !> lower side of the panel
            call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres( &     
                  geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
          endif 
        end do
      !$omp end parallel do 
#endif    
    endif 

    ! check i_el for multiple components   
    do i_el = 1, n_pan_te  
      if (associated(wake%pan_gen_elems(2,i_el)%p)) then !> surface panels 
        delta_mag_te(i_el) = elems(wake%pan_gen_elems_id(1,i_el))%p%mag - &
                              elems(wake%pan_gen_elems_id(2,i_el))%p%mag
        delta_pres_te(i_el) = elems(wake%pan_gen_elems_id(1,i_el))%p%pres - &
                              elems(wake%pan_gen_elems_id(2,i_el))%p%pres
        mag_te_tmp_upper(i_el) = elems(wake%pan_gen_elems_id(1,i_el))%p%mag
        mag_te_tmp_lower(i_el) = elems(wake%pan_gen_elems_id(2,i_el))%p%mag
      endif 
    enddo 
#if USE_PRECICE 
    !update jacobian only at the first iteration
    if ((it .eq. sim_param%kutta_startstep .or. mod(it, sim_param%kutta_update_jacobian) .eq. 0) & 
        .and. it_precice .eq. 1 ) then 
#else
    if (it .eq. sim_param%kutta_startstep .or. mod(it, sim_param%kutta_update_jacobian) .eq. 0) then 
#endif 
      !> initialize the perturbation matrix: for each columns we have the 
      !> intensities from the steady kutta condition for each trailing edge panel element 
      
      do i_el = 1, n_pan_te  ! for all panel elements 
        if (associated(wake%pan_gen_elems(2,i_el)%p)) then 
          mag_pert(i_el, :) = delta_mag_te
        endif
      enddo 

      !> Starting the perturbation on the diagonal terms of mag_pert 

      do i_el = 1, n_pan_te  ! for all panel elements
        if(associated(wake%pan_gen_elems(2,i_el)%p)) then 
          !> perturbation of the circulation
          mag_pert(i_el, i_el) = (1.0_wp - sim_param%kutta_beta)*(delta_mag_te(i_el)) 

        !> rhs perturbed 
        linsys%b = rhs_tmp - matmul(linsys%TL, mag_pert(i_el, :)) - matmul(linsys%TR, delta_mag_te_old)      
        
        !> Solve the factorized system: A_free_wake is factorized only once and then preserved 
        !> for A_free_wake and rhs_perturbed  
        if (linsys%rank .gt. 0) then ! A matrix includes the steady kutta condition
            linsys%skip = .true.
            linsys%A = linsys%A_wake_free
            call solve_linsys(linsys) 
        end if
        !> update unsteady term 
        !$omp parallel do private(k_el)
        do k_el = 1 , sel      
          elems(k_el)%p%didou_dt = (linsys%res(k_el) - res_old(k_el)) / sim_param%dt
        enddo 
        !$omp end parallel do 
        !> compute perturbed pressure at trailing edge 

#if USE_PRECICE
        !$omp parallel do private(j_el)
          do j_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
              if (geo%components(elems(wake%pan_gen_elems_id(1,j_el))%p%comp_id)%coupling) then  
                !> calculate the pressure using the relative orientation matrix
                !> upper side of the panel
                call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres(te%e(1, j_el)%p%R_cen) 
                !> lower side of the panel  
                call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres(te%e(2, j_el)%p%R_cen) 
              else !> non coupled component
                !> upper side of the panel 
                call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                      geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
                !> lower side of the panel
                call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres( &     
                      geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
              endif   
            endif    
          enddo
        !$omp end parallel do
#else
        !$omp parallel do private(j_el)
          do j_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
              !> upper side of the panel
              call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
              !> lower side of the panel
              call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
            endif 
          enddo
        !$omp end parallel do
#endif   
        !> compute the perturbed pressure difference
        do j_el = 1, n_pan_te
          if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
            delta_pres_te_perturbed(j_el) = ((elems(wake%pan_gen_elems_id(1,j_el))%p%pres - & 
                                              elems(wake%pan_gen_elems_id(2,j_el))%p%pres))
            jacobi(j_el, i_el) = (delta_pres_te_perturbed(j_el) - delta_pres_te(j_el))/ & 
                                  (mag_pert(i_el, i_el)-delta_mag_te(i_el)) 
          endif 
        enddo
        !> restore the initial values to the unperturbed ones 
        linsys%b = rhs_tmp  
        linsys%res = res_tmp
        linsys%A = A_tmp 
      endif
      enddo
      !> inverse of jacobian matrix (only the non zero part is inverted, the rest is zero) 
      call invmat_banded(jacobi, size(jacobi,1))
    endif !> update jacobi matrix

    !> Newton Raphson iteration to get the new circulation/pressure difference
    if (it .gt. sim_param%kutta_startstep) then 
      !> Newton-Raphson iteration to get the new circulation 
      it_pan = 1 
      tol = sim_param%kutta_tol + 1.0e-6_wp
      delta_pres_te_iter = delta_pres_te
      delta_mag_te_iter = delta_mag_te  

      do while (tol .gt. sim_param%kutta_tol .and. it_pan .lt. sim_param%kutta_maxiter)
!        
        delta_mag_te_iter = delta_mag_te_iter - matmul(jacobi, delta_pres_te_iter)
        linsys%b = rhs_tmp - matmul(linsys%TL, delta_mag_te_iter) - matmul(linsys%TR, delta_mag_te_old)     

        !> Solve the factorized system
        if (linsys%rank .gt. 0) then
          linsys%A = linsys%A_wake_free
          call solve_linsys(linsys) 
        end if

        !> update unsteady term 
        !$omp parallel do private(k_el)
        do k_el = 1 , sel      
          elems(k_el)%p%didou_dt = (linsys%res(k_el) - res_old(k_el)) / sim_param%dt
        enddo 
        !$omp end parallel do 

        !> compute perturbed pressure at trailing edge
        
#if USE_PRECICE
        !$omp parallel do private(j_el)
          do j_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
              if (geo%components(elems(wake%pan_gen_elems_id(1,j_el))%p%comp_id)%coupling) then  
                !> calculate the pressure using the relative orientation matrix
                !> upper side of the panel
                call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres(te%e(1, j_el)%p%R_cen) 
                !> lower side of the panel  
                call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres(te%e(2, j_el)%p%R_cen) 
              else !> non coupled component
                !> upper side of the panel 
                call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                      geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
                !> lower side of the panel
                call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres( &     
                      geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
              endif    
            end if  
          enddo
        !$omp end parallel do
#else   
        !$omp parallel do private(j_el)
          do j_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
              !> upper side of the panel
              call elems(wake%pan_gen_elems_id(1,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(1, j_el)%p%comp_id)%ref_id)%R_g)
              !> lower side of the panel
              call elems(wake%pan_gen_elems_id(2,j_el))%p%compute_pres( &     
                    geo%refs(geo%components(te%e(2, j_el)%p%comp_id)%ref_id)%R_g)
            endif 
          end do
        !$omp end parallel do
#endif
        !$omp parallel do private(j_el)
        do j_el = 1, n_pan_te
          if(associated(wake%pan_gen_elems(2,j_el)%p)) then 
            delta_pres_te_iter(j_el) = (elems(wake%pan_gen_elems_id(1,j_el))%p%pres - &
                                        elems(wake%pan_gen_elems_id(2,j_el))%p%pres) 
          endif 
        enddo
        !$omp end parallel do

        tol = maxval(abs(delta_pres_te_iter))
        it_pan = it_pan + 1

        if(it_pan .eq. sim_param%kutta_maxiter) then
          call warning('dust','dust','max iteration reached for kutta condition:&
                      &increase kutta_maxiter!') 
          write(message,*) 'Last iteration error: ', tol
          call printout(message)
        endif
        
      enddo !> while 
      
    endif !>  sim_param%kutta_startstep

!> update pressure field 
#if USE_PRECICE
  !$omp parallel do private(i_el, theta_cen, R_cen)
    do i_el = 1, sel
      ! ifort bugs workaround:
      ! apparently it is not possible to call polymorphic methods inside
      ! select cases for intel, need to call these for all elements and for the
      ! vortex lattices it is going to be a dummy empty function call
      if (geo%components(elems(i_el)%p%comp_id)%coupling) then  
        !> calculate the pressure using the relative orientation matrix
        call elems(i_el)%p%compute_pres(elems(i_el)%p%R_cen)  ! update surf_vel field too
      else !> non coupled component 
          call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
              geo%refs(geo%components(elems(i_el)%p%comp_id)%ref_id)%R_g)
      endif  
      call elems(i_el)%p%compute_dforce()      
    end do
  !$omp end parallel do
#else
  
  !$omp parallel do private(i_el)
    do i_el = 1 , sel
      call elems(i_el)%p%compute_pres( &     ! update surf_vel field too
              geo%refs(geo%components(elems(i_el)%p%comp_id)%ref_id)%R_g)
      call elems(i_el)%p%compute_dforce()  
    end do
  !$omp end parallel do
#endif
    linsys%skip = .false.
    linsys%A = A_tmp !> restore the A matrix for the next time iteration
  endif !> sim_param%kutta 

!> Zeroing pressure and force for lateral hinge panels (i.e. blending region) to avoid unphysical forces
!  Crude hardcoded workaround, checking angle between panel normal and hinge axis for >10 deg and then
!  checking if the panel is actually close to the hinge (hardcoded 2*max edge) by looking at the distance
!  in the hinge axis direction. 
!  Todo: fix it by adding a reference to the actual blending region 
!  Warning: might give issues with large dihedrals/deflections near the hinge or for unconventional geometries
!!$omp parallel do private(i_el, ih)
!  do i_el = 1, sel 
!    if (size(geo%components(elems(i_el)%p%comp_id)%hinge) .gt. 0) then 
!      do ih = 1, size(geo%components(elems(i_el)%p%comp_id)%hinge)
!        if(abs(dot_product(elems(i_el)%p%nor, geo%components(elems(i_el)%p%comp_id)%hinge(ih)%ref%h(:,1))) .gt. sin(10.0_wp*pi/180.0_wp) .and. &
!          min(abs(dot_product(elems(i_el)%p%cen-geo%components(elems(i_el)%p%comp_id)%hinge(ih)%ref%rr0,geo%components(elems(i_el)%p%comp_id)%hinge(ih)%ref%h(:,1))),&
!              abs(dot_product(elems(i_el)%p%cen-geo%components(elems(i_el)%p%comp_id)%hinge(ih)%ref%rr1,geo%components(elems(i_el)%p%comp_id)%hinge(ih)%ref%h(:,2))))&
!              .lt. 2.0_wp*maxval(elems(i_el)%p%edge_len))then                        
!                elems(i_el)%p%dforce = 0.0_wp 
!                elems(i_el)%p%pres = 0.0_wp
!        endif  
!      enddo   
!    endif 
!  end do 
!!$omp end parallel do 

  call solve_nlvl_correction(geo, wake, elems, elems_non_corr, elems_corr, linsys, airfoil_data, &
                             res_old, sel, it, time_2_debug_out, basename_debug)

  ! Explicit elements:
  ! - liftlin: _pres and _dforce computed in solve_liftin()
  ! - actdisk: avg delta_pressure and force computed here,
  !            to include thier effects in postpro (e.g. integral loads)
!$omp parallel do private(i_el)
  do i_el = 1 , size(elems_ad)
    call elems_ad(i_el)%p%compute_pres( &     ! update surf_vel field too
          geo%refs( geo%components(elems_ad(i_el)%p%comp_id)%ref_id )%R_g)
    call elems_ad(i_el)%p%compute_dforce()
  end do
!$omp end parallel do

#if USE_PRECICE
    sum_force = 0.0_wp
    do i = 1, size(elems_tot)
      sum_force = sum_force + elems_tot(i)%p%dforce
    end do

    !> Update force and moments to be passed to the structural solver
    call precice%update_force(geo, elems_tot)

    call precicef_is_coupling_ongoing(precice%is_ongoing)
    if (precice%is_ongoing .eq. 1) then
      call precicef_advance( precice%dt_precice )
    end if

    !> Write force and moments to structural solver
    do i = 1, size(precice%fields)
      if (trim(precice%fields(i)%fio) .eq. 'write') then
        if (trim(precice%fields(i)%ftype) .eq. 'scalar') then
          call precicef_write_data(   precice%mesh%mesh_name, &
                                      precice%fields(i)%fname, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%fields(i)%fdata(1,:) )
        elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
          call precicef_write_data(   precice%mesh%mesh_name, &
                                      precice%fields(i)%fname, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%fields(i)%fdata )
        endif
      end if
    end do

    !===========> precice convergence check
    ! if not converging skip everything and goto end of time cycle to repeat step
    ! if converged continue to the update section
    call precicef_requires_reading_checkpoint(bool)
    if ( bool .eq. 1 ) then 
      !> timestep not converged
      !> Reload checkpoint state
      do j = 1, size(precice%fields)
        if ( trim(precice%fields(j)%fio) .eq. 'write' ) then
          precice%fields(j)%fdata = precice%fields(j)%cdata
        end if
      end do
      it_precice = it_precice + 1
    else ! else contains everything down to the end of the time cycle (l. 1310)
      !> timestep converged
      precice_convergence = .true.
#else
#endif


    !> Print the results 
    if(time_2_out)  then
      nout = nout+1
      call save_status(geo, wake, nout, time, run_id)
    endif
  
  !> if the simulation has been restarted it should jump here because
  ! the solution comes from the restart file 
#if USE_PRECICE
#else
  endif ! already solved because restarted
#endif
  already_solv_restart = .false.

    !> Viscous Effects and Flow Separations 
    ! some computation of surface quantities and vorticity to be released.
    ! Free vortices will be introduced in prepare_wake(), some lines below
    if(sim_param%use_ve) call viscosity_effects( geo , elems , te )

    !> Wake update comes next, but it's for the next step, so if the 
    ! simulation is coupled we need to query mbdyn again for the updated
    ! positions
#if USE_PRECICE
    !> Read data from structural solver
    do i = 1, size(precice%fields)
      if ( trim(precice%fields(i)%fio) .eq. 'read' ) then
        if ( trim(precice%fields(i)%ftype) .eq. 'scalar' ) then
          call precicef_read_data(  precice%mesh%mesh_name, &
                                    precice%fields(i)%fname, &
                                    precice%mesh%nnodes  , &
                                    precice%mesh%node_ids, &
                                    precice%dt_precice, &
                                    precice%fields(i)%fdata(1,:) )
        elseif ( trim(precice%fields(i)%ftype) .eq. 'vector' ) then
            call precicef_read_data(  precice%mesh%mesh_name, &
                                      precice%fields(i)%fname, &
                                      precice%mesh%nnodes  , &
                                      precice%mesh%node_ids, &
                                      precice%dt_precice, &
                                      precice%fields(i)%fdata )
        endif
      end if
    end do

    !> Update dust geometry
    call precice%update_elems( geo, elems_tot,  elems_virtual, te)

    !> Update dt--> mbdyn should take care of the dt and send it to precice (TODO)
#endif


    !> Treat the wake: this needs to be done after output
    !  in practice the update is for the next iteration;
    !  this means that in the following routines the wake points are already
    !  updated to the next step, so we can immediatelt
    !  add the new particles if needed
    
    t0 = dust_time()
    if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
      
      !> update the trailing edge to go into the wake (average)
      if (sim_param%kutta_correction .and. it .gt. sim_param%kutta_startstep - 1) then
        !$omp parallel do private(i_el)
          do i_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,i_el)%p)) then 
              elems(wake%pan_gen_elems_id(1,i_el))%p%mag  = (elems(wake%pan_gen_elems_id(1,i_el))%p%mag + &
                                                              delta_mag_te_upper_old(i_el))/2.0_wp
              elems(wake%pan_gen_elems_id(2,i_el))%p%mag  = (elems(wake%pan_gen_elems_id(2,i_el))%p%mag + &
                                                              delta_mag_te_lower_old(i_el))/2.0_wp
            else !> todo 
            endif
          enddo
        !$omp end parallel do
      endif

      call update_wake(wake, geo, elems_tot, octree)
      
      !> restore the trailing edge to its original intensity 
      if (sim_param%kutta_correction .and. it .gt. sim_param%kutta_startstep - 1) then
        !$omp parallel do private(i_el)
          do i_el = 1, n_pan_te
            if(associated(wake%pan_gen_elems(2,i_el)%p)) then 
              elems(wake%pan_gen_elems_id(1,i_el))%p%mag  = (elems(wake%pan_gen_elems_id(1,i_el))%p%mag*2.0_wp - &
                                                              delta_mag_te_upper_old(i_el))
              elems(wake%pan_gen_elems_id(2,i_el))%p%mag  = (elems(wake%pan_gen_elems_id(2,i_el))%p%mag*2.0_wp - &
                                                              delta_mag_te_lower_old(i_el))
            else !> todo 
            endif 
          enddo
        !$omp end parallel do
      endif
    
    end if
    t1 = dust_time()

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Updated wake in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    !> Pressure integral equation 
    !> save old velocity on the surfpan (before updating the geom, few lines below)
    do i_el = 1 , geo%nSurfPan
      select type ( el => elems(i_el)%p ) ; class is ( t_surfpan )
        surf_vel_SurfPan_old( geo%idSurfPanG2L(i_el) , :) = el%ub   ! el%surf_vel
            nor_SurfPan_old( geo%idSurfPanG2L(i_el) , :) = el%nor   ! el%surf_vel
      end select
    end do

#if USE_PRECICE
    !> Update geo_data()
    do i_el = 1, size(elems_tot)
      call elems_tot(i_el)%p%calc_geo_data( &
                            geo%points(:,elems_tot(i_el)%p%i_ver) )
    end do
    do i_el = 1, size(elems_virtual)
      call elems_virtual(i_el)%p%calc_geo_data_virtual( &
                            geo%points_virtual(:,elems_virtual(i_el)%p%i_ver) )
    end do

    !> Update near-field wake
    call precice%update_near_field_wake( geo, wake, te )
#endif

    t0 = dust_time()
    !> update geometry 
    if(it .lt. nstep) then
      time = min(sim_param%tend, sim_param%time_vec(it+1))
      !> Update geometry
      call update_geometry(geo, te, time, .false., .true.)
      if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
        call complete_wake(wake, geo, elems_tot, elems_virtual, te, octree, it)
      end if
    endif
    t1 = dust_time() 

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Updated geometry in: ' , t1 - t0,' s.'
      call printout(message)
    endif
    !> Save old solution (at the previous dt) of the linear system
    res_old = linsys%res

    
    if (geo%nSurfPan .gt. 0 .and. sim_param%kutta_correction .and. & 
        (it .gt. sim_param%kutta_startstep - 1)) then  
      delta_mag_te_old = delta_mag_te
      delta_mag_te_upper_old = mag_te_tmp_upper
      delta_mag_te_lower_old = mag_te_tmp_lower
    endif 


    !> Update nor_old (moved from geo/mod_geo.f90/update_geometry(), l.2220 approx
!$omp parallel do private(i_el)
    do i_el = 1 , size(elems_tot)
      elems_tot(i_el)%p%nor_old = elems_tot(i_el)%p%nor
    end do
!$omp end parallel do

    !> Save alpha_old for dynamic stall 
    do i_c = 1, size(geo%components)
      if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
          trim(geo%components(i_c)%aero_correction) .eq. 'true' .and. &
          sim_param%vl_dynstall) then 
!$omp parallel do private(i_s)
          do i_s = 1, size(geo%components(i_c)%stripe)
            geo%components(i_c)%stripe(i_s)%alpha_old = geo%components(i_c)%stripe(i_s)%alpha 
          end do 
!$omp end parallel do  
      end if
    end do 

#if USE_PRECICE
      ! *** to do *** dirty implementation. it-update moved into #if USE_PRECICE
      ! statement, to avoid double time counter update, when the code is not coupled
      ! with external softwares through PRECICE.
      ! #if .not. USE_PRECICE, it-update occurs at the beginning of the time loop,
      ! approximately at l.615.
      ! Is there a reason why it-update should occur in two different places of the
      ! time loop (at the begin w/o precice, at the end w/ precice)?
      
      !> Precice iters are done, so we can store the old points
      wake%update_old_second_row = .true.
      
      !> Update n. time step
      it = it + 1
      !> reset it_precice counter 
      it_precice = 0
    endif ! End of the if statement that check whether the timestep
          ! has converged or not (l. 1115)
#endif

  enddo !> while do 
call printout(nl//'\\\\\\\\\\  Computations Finished \\\\\\\\\\')
deallocate(res_old)
!> End Time Cycle 

!> Finalize PreCICE 
#if USE_PRECICE
  call precicef_finalize()
#endif

!> Cleanup 
call destroy_wake(wake)
call destroy_octree(octree)
call destroy_linsys(linsys)
call destroy_elements(geo)
!call destroy_geometry(geo, elems_tot)
call destroy_hdf5()

t22 = dust_time()
!> Debug
if(sim_param%debug_level .ge. 1) then
  write(message,'(A,F9.3,A)') 'Completed all computations in ',t22-t00,' s'
  call printout(message)
endif
call printout(nl//'<<<<<< DUST end <<<<<<'//nl)

contains

!> Functions (maybe put everything in a dedicated module) 
subroutine get_run_id (run_id)
  integer, intent(out) :: run_id(10)
  real(wp)             :: randr
  integer              :: maxi, randi

  !> First 8 values are the date and time
  call date_and_time(VALUES=run_id(1:8))

  !> Last 3 values are 2 random integers
  maxi = huge(maxi)
  call random_number(randr)
  randi = int(randr*real(maxi,wp))
  run_id(9) = randi

  call random_number(randr)
  randi = int(randr*real(maxi,wp))
  run_id(10) = randi
end subroutine



!> Perform preliminary procedures each timestep, mainly chech if it is time
!! to perform output or not
subroutine init_timestep(t)
  real(wp), intent(in) :: t

  sim_param%time = t

  !if (real(t-t_last_out) .ge. real(sim_param%dt_out)) then
  if (real(time_no_out) .ge. real(sim_param%dt_out)) then
    time_2_out = .true.
    t_last_out = t
    time_no_out = 0.0_wp
  else
    time_2_out = .false.
  endif

  !if (real(t-t_last_debug_out) .ge. real(dt_debug_out)) then
  if (real(time_no_out_debug) .ge. real(dt_debug_out)) then
    time_2_debug_out = .true.
    t_last_debug_out = t
    time_no_out_debug = 0.0_wp
  else
    time_2_debug_out = .false.
  endif

  !If it is the last timestep output the solution, unless dt_out is set
  !longer than the whole execution, declaring implicitly that no output is
  !required.
  if((it .eq. nstep) .and. (sim_param%dt_out .le. sim_param%tend)) then
    time_2_out = .true.
    time_2_debug_out = .true.
  endif

  !if requested, print the output also at the first step (t0)
  if ((t.eq.sim_param%t0) .and. output_start) then
    t_last_out = t
    t_last_debug_out = t
    time_2_out = .true.
    time_2_debug_out = .true.
  endif

  time_no_out = time_no_out + sim_param%dt
  time_no_out_debug = time_no_out_debug + sim_param%dt

end subroutine init_timestep

subroutine debug_printout_result(linsys, basename, it)
  type(t_linsys),   intent(in) :: linsys
  character(len=*), intent(in) :: basename
  integer,          intent(in) :: it
  real(wp), allocatable        :: res(:,:)
  character(len=max_char_len)  :: sit

  allocate(res(1,linsys%rank+linsys%n_expl))
  res(1,:) = (/linsys%res,linsys%res_expl(:,1)/)
  write(sit,'(I4.4)') it
  call write_basic(res,trim(basename)//'_result_'//trim(sit)//'.dat')
  deallocate(res)

end subroutine debug_printout_result

subroutine debug_printout_geometry(elems, geo, basename, it)
  type(t_impl_elem_p),   intent(in) :: elems(:)
  type(t_geo),      intent(in)      :: geo
  character(len=*), intent(in)      :: basename
  integer,          intent(in)      :: it
  real(wp), allocatable             :: norm(:,:), cent(:,:), velb(:,:)
  integer, allocatable              :: el(:,:), conn(:,:)
  character(len=max_char_len)       :: sit
  integer                           :: ie, iv
  !> surf_vel and vel_phi
  real(wp), allocatable             :: surf_vel(:,:), vel_phi(:,:)
  integer                           :: i_e
  !> chtls
  real(wp)                          :: f(5)
  integer                           :: n_neigh

  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  !> surf_vel and vel_phi
  allocate( surf_vel(3,size(elems)), vel_phi(3,size(elems)) )
  surf_vel = -666.6_wp ; vel_phi = -666.6_wp

  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    velb(:,ie) = elems(ie)%p%ub
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
    do iv=1,elems(ie)%p%n_ver
      if(associated(elems(ie)%p%neigh(iv)%p)) then
        conn(iv, ie) = elems(ie)%p%neigh(iv)%p%id
      else
        conn(iv, ie) = 0
      endif
    enddo

    !> surf_vel and vel_phi for surfpan only
    select type( el => elems(ie)%p )
      class is(t_surfpan)
        surf_vel(:,ie) = el%surf_vel   ! elems(ie)%p%surf_vel
        !>  2. CHTLS method ------
        vel_phi(:,ie) = 0.0_wp
        f = 0.0_wp ; n_neigh = 0
        do i_e = 1 , el%n_ver
          if ( associated(el%neigh(i_e)%p) ) then
            n_neigh = n_neigh + 1
            f(n_neigh) = - ( el%neigh(i_e)%p%mag - el%mag )
          end if
        end do
        f(n_neigh+1) = sum(el%nor * (-variable_wind(el%cen, sim_param%time) - el%uvort + el%ub) )
        vel_phi(:,ie) = matmul( el%chtls_stencil , f(1:n_neigh+1) )
    end select
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(velb,       trim(basename)//'_mesh_velb_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_mesh_elems_'  //trim(sit)//'.dat')
  call write_basic(conn,       trim(basename)//'_mesh_conn_'   //trim(sit)//'.dat')
  deallocate(norm, cent, el, conn, velb)

  ! surf_vel and vel_phi
  call write_basic(surf_vel,   trim(basename)//'_mesh_surfvel_'  //trim(sit)//'.dat')
  call write_basic(vel_phi ,   trim(basename)//'_mesh_velphi_'   //trim(sit)//'.dat')
  deallocate(surf_vel,vel_phi)

end subroutine debug_printout_geometry

!------------------------------------------------------------------------------

subroutine debug_printout_geometry_minimal(elems,geo,basename, it)
  type(t_impl_elem_p),   intent(in) :: elems(:)
  type(t_geo),      intent(in)      :: geo
  character(len=*), intent(in)      :: basename
  integer,          intent(in)      :: it
  real(wp), allocatable             :: norm(:,:), cent(:,:)
  integer, allocatable              :: el(:,:)
  character(len=max_char_len)       :: sit
  integer                           :: ie
  integer(h5loc)                    :: h5fid

  allocate(norm(3,size(elems)), cent(3,size(elems)))
  allocate(el(4,size(elems))); el = 0

  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_mesh_elems_'  //trim(sit)//'.dat')

  call new_hdf5_file(trim(basename)//'_geo_'  //trim(sit)//'.h5',h5fid)
  call write_hdf5(geo%points,'points',h5fid)
  call close_hdf5_file(h5fid)

  deallocate(norm, cent, el)

end subroutine debug_printout_geometry_minimal

! ---------------------------------------------------------------------

subroutine debug_ll_printout_geometry(elems, geo, basename, it)
  type(t_liftlin_p),   intent(in)  :: elems(:)
  type(t_geo),      intent(in)     :: geo
  character(len=*), intent(in)     :: basename
  integer,          intent(in)     :: it
  real(wp), allocatable            :: norm(:,:), cent(:,:), velb(:,:)
  integer, allocatable             :: el(:,:), conn(:,:)
  character(len=max_char_len)      :: sit
  integer                          :: ie, iv
  
  allocate(norm(3,size(elems)), cent(3,size(elems)), velb(3,size(elems)))
  allocate(el(4,size(elems))); el = 0
  allocate(conn(4,size(elems))); conn = -666;
  
  do ie=1,size(elems)
    norm(:,ie) = elems(ie)%p%nor
    cent(:,ie) = elems(ie)%p%cen
    velb(:,ie) = elems(ie)%p%ub
    el(1:elems(ie)%p%n_ver,ie) = elems(ie)%p%i_ver
    do iv=1,elems(ie)%p%n_ver
      if(associated(elems(ie)%p%neigh(iv)%p)) then
        conn(iv, ie) = elems(ie)%p%neigh(iv)%p%id
      else
        conn(iv, ie) = 0
      endif
    enddo
  enddo

  write(sit,'(I4.4)') it
  call write_basic(geo%points, trim(basename)//'_ll_mesh_points_'//trim(sit)//'.dat')
  call write_basic(norm,       trim(basename)//'_ll_mesh_norm_'  //trim(sit)//'.dat')
  call write_basic(velb,       trim(basename)//'_ll_mesh_velb_'  //trim(sit)//'.dat')
  call write_basic(cent,       trim(basename)//'_ll_mesh_cent_'  //trim(sit)//'.dat')
  call write_basic(el,         trim(basename)//'_ll_mesh_elems_'  //trim(sit)//'.dat')
  call write_basic(conn,       trim(basename)//'_ll_mesh_conn_'   //trim(sit)//'.dat')
  deallocate(norm, cent, el, conn, velb)

end subroutine debug_ll_printout_geometry

end program dust




