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


!> Module to handle the input/output of the solution of DUST in hdf5 format,
!! meant principally to exchange data with the post processor and restart
module mod_dust_io_test

use mod_param, only: &
  wp, max_char_len, nl

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

use mod_aeroel, only: &
  c_elem, c_vort_elem, &
  t_elem_p, t_vort_elem_p

use mod_hdf5_io, only: &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  write_hdf5_attr, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5


use mod_wake_test, only: &
  t_wake

use mod_stringtools, only: &
  stricmp

use mod_version, only: &
  git_sha1, version

use mod_parse, only: &
  t_parse, getstr, getint, getintarray, getreal, getrealarray, getlogical, countoption, &
  getsuboption,&
  finalizeparameters, t_link, check_opt_consistency


!----------------------------------------------------------------------

implicit none

public :: save_status

private


character(len=*), parameter :: this_mod_name = 'mod_dust_io_test'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine save_status(wake, it, time, run_id)
  type(t_wake), intent(in)          :: wake
  integer, intent(in)               :: it
  real(wp), intent(in)              :: time
  integer, intent(in)               :: run_id(10)

  integer(h5loc)                    :: floc, gloc1, gloc2, gloc3, gloc4, ploc
  character(len=max_char_len)       :: sit
  integer                           :: iref, nref
  character(len=max_char_len)       :: ref_name
  integer                           :: ie, ne, ip
  real(wp), allocatable             :: turbvisc(:), v_rad(:)
  real(wp), allocatable             :: points_w(:,:,:), cent(:,:,:) , vel_w(:,:,:)
  real(wp), allocatable             :: vort_v(:,:)
  integer,  allocatable             :: conn_pe(:)
  real(wp), allocatable             :: ori(:,:)


  !> create the output file
  write(sit,'(I4.4)') it
  call new_hdf5_file(trim(sim_param%basename)//'_res_'//trim(sit)//'.h5', &
                      floc)
  call write_hdf5_attr(git_sha1, 'git_sha1', floc)
  call write_hdf5_attr(version, 'version', floc)
  call sim_param%save_param(floc)
  call write_hdf5(time,'time',floc)

  call new_hdf5_group(floc, 'Parameters', ploc)
  call write_hdf5(sim_param%u_inf,'u_inf', ploc)
  call write_hdf5(sim_param%P_inf,'P_inf', ploc)
  call write_hdf5(sim_param%rho_inf,'rho_inf', ploc)
  call write_hdf5(sim_param%mu_inf,'mu_inf', ploc)
  call write_hdf5(sim_param%a_inf,'a_inf', ploc)  
  call close_hdf5_group(ploc)

  ! 1) %%%% Wake:
  ! just print the whole points, the solution and the starting row
  ! connectivity to build connectivity after
  !=== Particles ===
  call new_hdf5_group(floc, 'ParticleWake', gloc1)
  allocate(points_w(3,wake%n_prt,1))
  allocate(vort_v(3,wake%n_prt))
  allocate(turbvisc(wake%n_prt))
  allocate(vel_w(3,wake%n_prt,1))
  allocate(v_rad(wake%n_prt))

  do ip = 1, wake%n_prt
    points_w(:,ip,1) = wake%part_p(ip)%p%cen
    vel_w(:,ip,1) = wake%part_p(ip)%p%vel
    vort_v(:,ip) = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag
    turbvisc(ip) = wake%part_p(ip)%p%turbvisc
    v_rad(ip) = wake%part_p(ip)%p%r_Vortex
  enddo
  !write(*,*) 'points_w = ', points_w
  call write_hdf5(points_w(:,:,1),'WakePoints',gloc1)
  call write_hdf5(   vel_w(:,:,1),'WakeVels'  ,gloc1)
  call write_hdf5( turbvisc,'turbvisc'  ,gloc1)
  call write_hdf5( v_rad,'VortexRad'  ,gloc1)
  call write_hdf5(vort_v,'WakeVort',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, vort_v, vel_w, turbvisc, v_rad)

  ! 3) %%%% References
!  call new_hdf5_group(floc, 'References', gloc1)
!  nref = size(geo%refs)
!  call write_hdf5(nref, 'NReferences', gloc1)
!  do iref = 0, nref-1
!    write(ref_name,'(A,I3.3)')'Ref',iref
!    call new_hdf5_group(gloc1, trim(ref_name), gloc2)
!
!    call write_hdf5(geo%refs(iref)%tag, 'Tag', gloc2)
!    call write_hdf5(geo%refs(iref)%of_g, 'Offset',gloc2)
!    call write_hdf5(geo%refs(iref)%R_g, 'R',gloc2)
!    call write_hdf5(geo%refs(iref)%f_g, 'Vel',gloc2)
!    call write_hdf5(geo%refs(iref)%G_g, 'RotVel',gloc2)
!    call write_hdf5(geo%refs(iref)%relative_pol_pos, 'RelativePolPos',gloc2)
!    call write_hdf5(geo%refs(iref)%relative_rot_pos, 'RelativeRotPos',gloc2)
!
!    call close_hdf5_group(gloc2)
!  enddo
!  call close_hdf5_group(gloc1)
    
  call close_hdf5_file(floc)
  
end subroutine save_status

!----------------------------------------------------------------------

end module mod_dust_io_test
