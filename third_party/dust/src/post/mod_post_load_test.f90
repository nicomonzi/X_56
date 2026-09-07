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

!> Module containing the tools to load the solution for postprocessing
module mod_post_load_test

use mod_param, only: &
  wp, nl, max_char_len, pi

use mod_handling, only: &
  error, internal_error, warning, info, printout

use mod_hdf5_io, only: &
  h5loc, &
  open_hdf5_file, &
  close_hdf5_file, &
  open_hdf5_group, &
  close_hdf5_group, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5, &
  get_dset_dimensions_hdf5

use mod_wake_test, only: &
  t_wake


implicit none

public :: load_wake_post, load_wake_viz


private

character(len=*), parameter :: this_mod_name = 'mod_post_load_test'
character(len=max_char_len) :: msg

contains

!----------------------------------------------------------------------

subroutine load_wake_viz(floc, vppoints,  vpvort, &
                          vpvort_v, v_rad, vpturbvisc)
  integer(h5loc), intent(in)                   :: floc
  real(wp), allocatable, intent(out)           :: vppoints(:,:)
  real(wp), allocatable, intent(out)           :: vpvort(:)
  real(wp), allocatable, intent(out)           :: vpvort_v(:,:)
  real(wp), allocatable, intent(out)           :: v_rad(:)

  real(wp), allocatable, intent(out), optional :: vpturbvisc(:)
  
  integer(h5loc)                               :: gloc
  logical                                      :: got_dset
  real(wp), allocatable                        :: wvort_read(:,:)
  integer                                      :: ip


  got_dset = check_dset_hdf5('ParticleWake',floc)
  if(got_dset) then

    call open_hdf5_group(floc,'ParticleWake',gloc)
    call read_hdf5_al(vppoints,'WakePoints',gloc)
    call read_hdf5_al(wvort_read,'WakeVort',gloc)
    call read_hdf5_al(v_rad,'VortexRad',gloc)
    if(present(vpturbvisc)) call read_hdf5_al(vpturbvisc,'turbvisc',gloc)
    allocate(vpvort(size(wvort_read,2)))
    allocate(vpvort_v(3,size(wvort_read,2))) !vorticity vector 
    do ip = 1,size(vpvort)
      vpvort(ip) = norm2(wvort_read(:,ip))
      vpvort_v(:,ip) = wvort_read(:,ip)/vpvort(ip)
      !write(*,*) 'vorticity vector', vpvort_v(:,ip) !, vpvort(ip)
    enddo
    
    !deallocate(wvort_read)
    call move_alloc(wvort_read, vpvort_v)
    call close_hdf5_group(gloc)
  endif

end subroutine load_wake_viz
!----------------------------------------------------------------------

!> Load the wake for the postprocessing.
!! do everything
subroutine load_wake_post(floc, wake)
  integer(h5loc), intent(in)                :: floc
  type(t_wake), target, intent(out)         :: wake

  integer(h5loc)                            :: gloc
  real(wp), allocatable                     :: vppoints(:,:), vpvort(:,:)
  real(wp), allocatable                     :: v_rad(:)
  integer                                   :: p1 , p2
  integer                                   :: ip , iw, id, ir, iconn, i


  !=== Particles + Line vortex ===
  call open_hdf5_group(floc, 'ParticleWake', gloc)
  call read_hdf5_al(vppoints,'WakePoints',gloc)
  call read_hdf5_al(vpvort,'WakeVort',gloc)
  call read_hdf5_al(v_rad,'VortexRad',gloc)
  call close_hdf5_group(gloc)

  wake%n_prt = size(vpvort,2)
  wake%nmax_prt = size(vpvort,2)

  allocate(wake%wake_parts(wake%nmax_prt))
  allocate(wake%prt_ivort(wake%nmax_prt))
  allocate(wake%part_p(wake%n_prt))
  if(wake%n_prt .gt. 0) then
    allocate(wake%vort_p(wake%n_prt))
  endif

  do ip = 1,wake%n_prt
    wake%wake_parts(ip)%cen = vppoints(:,ip)
    wake%wake_parts(ip)%mag => wake%prt_ivort(ip)
    wake%wake_parts(ip)%mag = norm2(vpvort(:,ip))
    wake%wake_parts(ip)%r_Vortex = v_rad(ip)
    if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
      wake%wake_parts(ip)%dir = vpvort(:,ip)/wake%wake_parts(ip)%mag
    else
      wake%wake_parts(ip)%dir = vpvort(:,ip)
    endif
    wake%wake_parts(ip)%free = .false.
    wake%part_p(ip)%p => wake%wake_parts(ip)
    wake%vort_p(ip)%p => wake%wake_parts(ip)
  enddo

  deallocate(vppoints, vpvort, v_rad)

end subroutine

!----------------------------------------------------------------------

end module mod_post_load_test
