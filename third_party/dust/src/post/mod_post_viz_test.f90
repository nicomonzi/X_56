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

!> Module containing the subroutines to perform visualizations during
!! postprocessing
module mod_post_viz_test

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, &
  countoption

use mod_hdf5_io, only: &
  h5loc, &
  open_hdf5_file, &
  close_hdf5_file, & 
  open_hdf5_group, &
  close_hdf5_group, &
  read_hdf5

use mod_stringtools, only: &
  LowCase, isInList

use mod_tecplot_out, only: &
  tec_out_viz

use mod_vtk_out_test, only: &
  vtk_out_viz


use mod_vtk_utils, only: &
  t_output_var, add_output_var, copy_output_vars, clear_output_vars

use mod_wind, only: &
  variable_wind

use mod_post_load_test, only: &
  load_wake_viz

implicit none

public :: post_viz

private

character(len=*), parameter :: this_mod_name='mod_post_viz_test'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_viz( sbprms , basename , data_basename , an_name , ia , &
                      out_frmt , an_start , an_end , an_step, average, an_avg)
                      
  type(t_parse), pointer                                    :: sbprms
  character(len=*) , intent(in)                             :: basename
  character(len=*) , intent(in)                             :: data_basename
  character(len=*) , intent(in)                             :: an_name
  integer          , intent(in)                             :: ia
  character(len=*) , intent(in)                             :: out_frmt
  integer , intent(in)                                      :: an_start , an_end , an_step, an_avg
  logical, intent(in)                                       :: average

  character(len=max_char_len)                               :: filename, filename_virtual, filename_in
  integer(h5loc)                                            :: floc , ploc
  logical                                                   :: out_vort, out_vort_vec, out_vel, out_cp, out_press
  logical                                                   :: out_wake, out_surfvel, out_vrad
  logical                                                   :: out_dforce, out_dmom
  logical                                                   :: out_turbvisc
  integer                                                   :: n_var , i_var
  character(len=max_char_len), allocatable                  :: var_names(:)
  real(wp), allocatable                                     :: points(:,:), points_exp(:,:)
  real(wp), allocatable                                     :: vppoints(:,:), vpvort(:)
  real(wp), allocatable                                     :: vpvort_v(:,:), vpturbvisc(:), v_rad(:)

  real(wp), allocatable                                     :: points_ave(:,:)

  real(wp)                                                  :: u_inf(3)
  real(wp)                                                  :: P_inf , rho

  real(wp), allocatable                                     :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                     :: vort(:), cp(:), vel(:), press(:), surfvel(:,:)

  type(t_output_var), allocatable                           :: out_vars(:), ave_out_vars(:)
  type(t_output_var), allocatable                           :: out_vars_w(:), out_vars_vp(:)
  type(t_output_var), allocatable                           :: out_vars_virtual(:)
  integer                                                   :: nprint , nprint_w, nelem_out, nelem_vp

  integer                                                   :: it, ires
  real(wp)                                                  :: t
  character(len=*), parameter                               :: this_sub_name='post_viz'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' visualization'//nl
  call printout(trim(msg))

  ! Print the wake or not
  out_wake = getlogical(sbprms,'wake')

  !Check which variables to analyse
  out_vort = .false.; 
  out_vort_vec = .false.;
  out_turbvisc = .false.;
  out_vel = .false.; 
  out_vrad = .false.; 

  n_var = countoption(sbprms, 'variable')
  
  allocate(var_names(n_var))
  do i_var = 1, n_var
    var_names(i_var) = getstr(sbprms, 'variable') ; call LowCase(var_names(i_var))
  enddo
  out_vort     = isInList('vorticity',          var_names) ! Always lower case string in the code !
  out_vort_vec = isInList('vorticity_vector',   var_names) ! Always lower case string in the code !
  out_vel      = isInList('velocity' ,          var_names)
  out_turbvisc = isInList('turbulent_viscosity',var_names)
  out_vrad     = isInList('vortex_rad',         var_names) 

  nprint = 0; nprint_w = 0
  if(out_vort)      nprint = nprint+1
  if(out_vort_vec)  nprint = nprint+1
  if(out_vel)       nprint = nprint+1  
  if(out_turbvisc)  nprint = nprint+1
  if(out_vrad)      nprint = nprint+1
  if(out_wake) then
    allocate(out_vars_vp(nprint))
  endif

  if(out_wake .and. average) call error(this_sub_name, this_mod_name, &
  'Cannot output an averaged wake visualization. Remove the wake or avoid &
  &averaging')

  ! Time loop
  ires = 0
  do it = an_start, an_end, an_step

    ires = ires+1

    ! Open the file
    write(filename_in,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
    call open_hdf5_file(trim(filename_in),floc)

    ! Load free-stream parameters
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(u_inf,'u_inf',ploc)
    call read_hdf5(P_inf,'P_inf',ploc)
    call read_hdf5(rho,'rho_inf',ploc)
    call close_hdf5_group(ploc)
    
    if(.not. average) then
      ! Output filename
      write(filename,'(A,I4.4)') trim(basename)//'_'//trim(an_name)//'-',it

      if(out_turbvisc) then
        call load_wake_viz(floc, vppoints, vpvort, &
                          vpvort_v, v_rad, vpturbvisc)
      else
        call load_wake_viz(floc, vppoints, vpvort, &
                            vpvort_v, v_rad)
      endif
      nelem_vp = size(vppoints,2)
  

      i_var = 1
      if(out_vort) then
        call add_output_var(out_vars_vp(i_var), vpvort, &
                'Singularity_Intensity',.false.)
        i_var = i_var +1
      endif
      if(out_turbvisc) then
        call add_output_var(out_vars_vp(i_var), vpturbvisc, &
                'Turbulent_Viscosity',.false.)
        i_var = i_var +1
      endif
      if(out_vort_vec) then
        call add_output_var(out_vars_vp(i_var), vpvort_v, &
                'Vorticity',.false.)
        i_var = i_var +1
      endif
      if(out_vrad) then
        call add_output_var(out_vars_vp(i_var), v_rad, &
                'VortexRad',.false.)
        i_var = i_var +1
      endif 
      
      !Output the results (with wake)
      select case (trim(out_frmt))
        case ('vtk')
          filename = trim(filename)//'.vtu'
          call  vtk_out_viz(filename, &
                              vp_rr=vppoints, vp_vars=out_vars_vp)                     
        case default
          call error('dust_post','','Unknown format '//trim(out_frmt)//&
                    ' for visualization output')
      end select

      call clear_output_vars(out_vars_vp)

    endif !not average

    call close_hdf5_file(floc)

    if (allocated(vort   ) ) deallocate(vort )
    if (allocated(vel    ) ) deallocate(vel  )
    if (allocated(v_rad  ) ) deallocate(v_rad)
    if (allocated(vpturbvisc)) deallocate(vpturbvisc)


  end do ! Time loop

  deallocate(out_vars_vp)

  write(msg,'(A,I0,A)') nl//'++++++++++ Visualization done'//nl
  call printout(trim(msg))

end subroutine post_viz

! ----------------------------------------------------------------------

end module mod_post_viz_test
