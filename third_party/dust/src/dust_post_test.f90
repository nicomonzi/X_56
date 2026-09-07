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

!> This is the main file of the DUST postprocessor
program dust_post_test

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_sim_param, only: &
  create_param_post, sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit, &
  check_basename, check_file_exists

use mod_aeroel, only: &
  c_elem, c_vort_elem, &
  t_elem_p, t_vort_elem_p
  
use mod_vortpart, only: &
  initialize_vortpart

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, getreal, getint, &
  getrealarray, getintarray , &
  ignoredParameters, finalizeParameters, &
  countoption, getsuboption, t_link, check_opt_consistency, &
  print_parse_debug

use mod_hdf5_io, only: &
  initialize_hdf5, destroy_hdf5, &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5

use mod_stringtools, only: &
  LowCase, isInList, stricmp

use mod_tecplot_out, only: &
  tec_out_viz, tec_out_probes, tec_out_box, tec_out_loads

use mod_vtk_out_test, only: &
  vtk_out_viz , vtr_write

use mod_dat_out, only: &
  dat_out_probes_header, &
  dat_out_loads_header

use mod_math, only: &
  cross

use mod_post_viz_test, only: &
  post_viz

use mod_post_probes_test, only: &
  post_probes

use mod_post_flowfield_test, only: &
  post_flowfield

implicit none

!Input
character(len=*), parameter               :: input_file_name_def = 'dust_post_test.in'
character(len=max_char_len)               :: input_file_name

!Geometry parameters
type(t_parse)                             :: prms
type(t_parse), pointer                    :: sbprms , bxprms
integer                                   :: n_analyses, ia             
character(len=max_char_len)               :: basename, data_basename
character(len=max_char_len)               :: an_name, an_type
integer                                   :: an_start, an_end, an_step, an_avg
character(len=max_char_len), allocatable  :: var_names(:)
character(len=max_char_len)               :: lowstr
character(len=max_char_len)               :: out_frmt
logical                                   :: average

call printout(nl//'>>>>>> DUST POSTPROCESSOR beginning >>>>>>'//nl)
call initialize_hdf5()

!------ Input reading ------

if(command_argument_count().gt.0) then
  call get_command_argument(1,value=input_file_name)
else
  input_file_name = input_file_name_def
endif

call printout(nl//'Reading input parameters from file "'//&
                trim(input_file_name)//'"'//nl)

call create_param_post(prms, sbprms , bxprms)

call check_file_exists(input_file_name, 'dust postprocessor')
call prms%read_options(input_file_name, printout_val=.false.)

basename = getstr(prms,'basename')
data_basename = getstr(prms,'data_basename')

sim_param%RankineRad = getreal(prms, 'rankine_rad')
sim_param%VortexRad = getreal(prms, 'vortex_rad')
sim_param%CutoffRad  = getreal(prms, 'cutoff_rad')


n_analyses = countoption(prms,'analysis')

!> Check that the basenames are valid
call check_basename(trim(basename),'dust postprocessor')
call check_basename(trim(data_basename),'dust postprocessor')

!> Cycle on all the analyses
do ia = 1, n_analyses

  !> Get general parameter for the analysis
  call getsuboption(prms,'Analysis',sbprms)
  an_type  = getstr(sbprms,'type')
  call LowCase(an_type)
  
  an_name  = getstr(sbprms,'name')
  out_frmt = getstr(sbprms,'format')
  call LowCase(out_frmt)
  
  an_start = getint(sbprms,'start_res')
  an_end   = getint(sbprms,'end_res')
  an_step  = getint(sbprms,'step_res')
  average  = getlogical(sbprms, 'average')
  
  if (countoption(sbprms, 'avg_res') .eq. 0) then
    an_avg = an_start
  else
    an_avg  = getint(sbprms,'avg_res')
  endif
  
  !> Fork the different kind of analyses
  select case( trim(an_type) )

    !> Visualizations
    case('viz')
      call post_viz( sbprms , basename , data_basename , an_name , ia , &
                    out_frmt , an_start , an_end , an_step, average, an_avg )
    case('probes')
      call post_probes( sbprms , basename , data_basename , an_name , ia , &
                            out_frmt, an_start , an_end , an_step )
    !> Flow Field
    case('flow_field')
      call post_flowfield( sbprms , basename , data_basename , an_name , ia , &
                          out_frmt , an_start , an_end , an_step, average)
                          
    case default

      call error('dust_post','','Unknown type of analysis: '//trim(an_type))
  end select

  if(allocated(var_names)) deallocate(var_names)

  sbprms => null()

enddo !analysis

call destroy_hdf5()
call printout(nl//'<<<<<< DUST POSTPROCESSOR end       <<<<<<'//nl)

!----------------------------------------------------------------------

end program dust_post_test
