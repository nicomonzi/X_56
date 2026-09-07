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
!!          Federico Gentile
!!          Matteo Dall'Ora       
!!          Alessandro Cocco
!!=========================================================================

!> This is the main file of the DUST solver

program dust_test

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_sim_param , only: &
  init_sim_param_particle, sim_param, t_sim_param, create_param_test_particle

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_basename, &
  check_file_exists

use mod_aeroel, only: &
  c_elem, c_vort_elem, & 
  t_elem_p, t_vort_elem_p

use mod_vortpart, only: &
  initialize_vortpart

use mod_wake_test, only: &
  t_wake, update_wake, initialize_wake, &
  prepare_wake, complete_wake

use mod_vtk_out, only: &
  vtk_out_bin

use mod_tecplot_out, only: &
  tec_out_sol_bin

use mod_hdf5_io, only: &
  h5loc, initialize_hdf5, destroy_hdf5, new_hdf5_file, open_hdf5_file, &
  close_hdf5_file, new_hdf5_group, open_hdf5_group, close_hdf5_group, &
  write_hdf5, write_hdf5_attr, read_hdf5, read_hdf5_al, append_hdf5

use mod_dust_io_test, only: &
  save_status

use mod_octree_test, only: &
  initialize_octree, destroy_octree, sort_particles, t_octree

use mod_math, only: & 
  cross, dot, vec2mat, invmat_banded

use mod_wind, only: &
  variable_wind

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

use mod_basic_io, only: &
  read_real_array_from_file 

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
logical                           :: time_2_out, time_2_debug_out, output_start 
real(wp)                          :: dt_debug_out
logical                           :: already_solv_restart 

!> Wake
type(t_wake)                      :: wake
!> Timing vars
real(t_realtime)                  :: t1 , t0, t00, t11, t22
!> I/O prefixes
character(len=max_char_len)       :: frmt, frmt_vl
character(len=max_char_len)       :: basename_debug

real(wp), allocatable             :: al_kernel(:,:), al_v(:), particlesMat(:,:)
real(wp)                          :: theta_cen(3), R_cen(3, 3) 
integer                           :: i

!> octree parameters
type(t_octree)                    :: octree

dt_debug_out = sim_param%dt_out
call printout(nl//'>>>>>> DUST_TEST beginning >>>>>>'//nl)

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

call create_param_test_particle(prms)

!> Get the parameters and print them out
call printout(nl//'====== Input parameters: ======')
call check_file_exists(input_file_name,'dust_test')
call prms%read_options(input_file_name, printout_val=.true.)

!> Initialize all the parameters
nout = 0  !> Reset the numbering for output files
output_start = getlogical(prms, 'output_start')
call init_sim_param_particle(sim_param, prms, nout, output_start) 

call read_real_array_from_file (8 , sim_param%particles_file, particlesMat)
sim_param%part_n0 = size(particlesMat,1)

allocate(sim_param%part_pos0(sim_param%part_n0, 3))
allocate(sim_param%part_vort0_dir(sim_param%part_n0, 3))
allocate(sim_param%part_vort0_mag(sim_param%part_n0))
allocate(sim_param%part_vol(sim_param%part_n0))

sim_param%part_pos0      = particlesMat(:,1:3)
sim_param%part_vort0_dir = particlesMat(:,4:6)
sim_param%part_vort0_mag = particlesMat(:,7)
sim_param%part_vol       = particlesMat(:,8)

!> Check that tend .gt. tinit
if ( sim_param%tend .le. sim_param%t0 ) then
  write(message,*) 'The end time of the simulation',sim_param%tend,'is&
                  & lower or equal to the start time',sim_param%t0,'.&
                  & Remembver that when restarting without resetting the&
                  & time the start time is taken from the restart result!'
  call error('dust_test','',message)
end if


!> Check that the basenames are valid 
call check_basename(trim(sim_param%basename),'dust test main')

!> Simulation parameters 
nstep = sim_param%n_timesteps
allocate(sim_param%time_vec(sim_param%n_timesteps))
sim_param%time_vec = (/ ( sim_param%t0 + &
          real(i-1,wp)*sim_param%dt, i=1,sim_param%n_timesteps ) /)

!> Initialization 
if(sim_param%use_fmm) then
  call printout(nl//'====== Initializing Octree ======')
  call initialize_octree(sim_param%BoxLength, sim_param%NBox, &
                        sim_param%OctreeOrigin, sim_param%NOctreeLevels, &
                        sim_param%MinOctreePart, sim_param%MultipoleDegree, &
                        sim_param%octree_vortex_rad, octree)
endif

call printout(nl//'====== Initializing Wake ======')

call initialize_wake(wake)

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

!===========> Start time cycle 
t11 = dust_time()
it = 0


do while ( ( it .lt. nstep ) )
  it = it + 1
  sim_param%time_old = sim_param%time
  write(message,'(A,I5,A,I5,A,F9.4)') nl//'--> Step ',it,' of ', &
                                        nstep, ' simulation time: ', time
  call printout(message)
  t22 = dust_time()
  write(message,'(A,F9.3,A)') 'Elapsed wall time: ', t22 - t00
  call printout(message)

  call init_timestep(time)

  if ( mod( it-1, sim_param%ndt_update_wake ) .eq. 0 ) then
      call prepare_wake(wake, octree)
  end if


!> Print the results 
    if(time_2_out)  then
      nout = nout+1
      call save_status(wake, nout, time, run_id)
    endif

    !> Treat the wake: this needs to be done after output
    !  in practice the update is for the next iteration;
    !  this means that in the following routines the wake points are already
    !  updated to the next step, so we can immediatelt
    !  add the new particles if needed
    
    t0 = dust_time()
    if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then     
      call update_wake(wake, octree)   
    end if
    t1 = dust_time()

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Updated wake in: ' , t1 - t0,' s.'
      call printout(message)
    endif

    t0 = dust_time()
    if(it .lt. nstep) then
      time = min(sim_param%tend, sim_param%time_vec(it+1))
      if ( mod( it, sim_param%ndt_update_wake ) .eq. 0 ) then
        call complete_wake(wake, octree) 
      end if
    endif
    t1 = dust_time() 

    !> debug message
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Completed wake in: ' , t1 - t0,' s.'
      call printout(message)
    endif
enddo !> while do 


call printout(nl//'\\\\\\\\\\  Computations Finished \\\\\\\\\\')
!> End Time Cycle 

!> Cleanup 
call destroy_octree(octree)
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

  time_no_out = time_no_out + sim_param%dt
  time_no_out_debug = time_no_out_debug + sim_param%dt

end subroutine init_timestep

! ---------------------------------------------------------------------

end program dust_test




