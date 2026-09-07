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
module mod_post_viz

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

use mod_geometry, only: &
  t_geo, t_geo_component, destroy_elements

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro, prepare_geometry_postpro, &
  expand_actdisk_postpro

use mod_tecplot_out, only: &
  tec_out_viz

use mod_vtk_out, only: &
  vtk_out_viz

use mod_post_load, only: &
  load_refs, load_res, load_wake_viz , check_if_components_exist

use mod_vtk_utils, only: &
  t_output_var, add_output_var, copy_output_vars, clear_output_vars

use mod_wind, only: &
  variable_wind

implicit none

public :: post_viz

private

character(len=*), parameter :: this_mod_name='mod_post_viz'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_viz( sbprms , basename , data_basename , an_name , ia , &
                      out_frmt , components_names , all_comp , &
                      an_start , an_end , an_step, average, an_avg)
  type(t_parse), pointer                                    :: sbprms
  character(len=*) , intent(in)                             :: basename
  character(len=*) , intent(in)                             :: data_basename
  character(len=*) , intent(in)                             :: an_name
  integer          , intent(in)                             :: ia
  character(len=*) , intent(in)                             :: out_frmt
  character(len=max_char_len), allocatable , intent(inout)  :: components_names(:)
  logical , intent(in)                                      :: all_comp
  integer , intent(in)                                      :: an_start , an_end , an_step, an_avg
  logical, intent(in)                                       :: average

  type(t_geo_component), allocatable                        :: comps(:)
  character(len=max_char_len)                               :: filename, filename_virtual, filename_in
  integer(h5loc)                                            :: floc , ploc
  logical                                                   :: out_vort, out_vort_vec, out_vel, out_cp, out_press
  logical                                                   :: out_wake, out_surfvel, out_vrad, out_pid
  logical                                                   :: out_dforce, out_dmom, out_sectional
  logical                                                   :: out_turbvisc
  logical                                                   :: separate_wake
  integer                                                   :: n_var , i_var
  character(len=max_char_len), allocatable                  :: var_names(:)
  real(wp), allocatable                                     :: points(:,:), points_exp(:,:) , wpoints(:,:)
  real(wp), allocatable                                     :: points_virtual(:,:), points_virtual_exp(:,:)
  real(wp), allocatable                                     :: vppoints(:,:), vpvort(:)
  real(wp), allocatable                                     :: vpvort_v(:,:), vpturbvisc(:), v_rad(:), vp_parentid(:)
  integer , allocatable                                     :: elems(:,:) , elems_virtual(:,:), welems(:,:)
  integer                                                   :: nelem , nelem_virtual, nelem_w, nelem_vp

  real(wp), allocatable                                     :: points_ave(:,:)

  real(wp)                                                  :: u_inf(3)
  real(wp)                                                  :: P_inf , rho

  real(wp), allocatable                                     :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                     :: vort(:), cp(:), vel(:), press(:), surfvel(:,:)
  real(wp), allocatable                                     :: force(:,:), moment(:,:), sectional(:,:)
  real(wp), allocatable                                     :: wvort(:)

  type(t_output_var), allocatable                           :: out_vars(:), ave_out_vars(:)
  type(t_output_var), allocatable                           :: out_vars_w(:), out_vars_vp(:)
  type(t_output_var), allocatable                           :: out_vars_virtual(:)
  integer                                                   :: nprint , nprint_w, nelem_out

  integer                                                   :: it, ires
  real(wp)                                                  :: t
  character(len=*), parameter                               :: this_sub_name='post_viz'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' visualization'//nl
  call printout(trim(msg))

  ! Print the wake or not
  out_wake = getlogical(sbprms,'wake')
  separate_wake = getlogical(sbprms,'separate_wake')

  !Check which variables to analyse 
  out_vort      = .false.
  out_vort_vec  = .false.
  out_vel       = .false. 
  out_surfvel   = .false.
  out_press     = .false. 
  out_cp        = .false.
  out_turbvisc  = .false.
  out_vrad      = .false.
  out_dforce    = .false.
  out_sectional = .false. 

  out_dmom      = .false.  
  out_pid       = .false.
  n_var = countoption(sbprms, 'variable')
  
  allocate(var_names(n_var))
  do i_var = 1, n_var
    var_names(i_var) = getstr(sbprms, 'variable') ; call LowCase(var_names(i_var))
  enddo

  out_vort     = isInList('vorticity',          var_names) 
  out_vort_vec = isInList('vorticity_vector',   var_names) 
  out_vel      = isInList('velocity' ,          var_names)
  out_surfvel  = isInList('surface_velocity',   var_names)
  out_press    = isInList('pressure' ,          var_names)
  out_cp       = isInList('cp'       ,          var_names)
  out_turbvisc = isInList('turbulent_viscosity',var_names)
  out_vrad     = isInList('vortex_rad',         var_names)
  out_dforce   = isInList('force',              var_names) 
  out_sectional= IsInList('sectional',          var_names)  
  out_dmom     = isInList('moment',             var_names)
  out_pid      = isInList('parent_id',          var_names)

  nprint = 0; nprint_w = 0
  if(out_vort)      nprint = nprint+1
  if(out_vort_vec)  nprint = nprint+1
  if(out_cp)        nprint = nprint+1
  if(out_surfvel)   nprint = nprint+1
  if(out_vel)       nprint = nprint+1  
  if(out_press)     nprint = nprint+1  
  if(out_turbvisc)  nprint = nprint+1
  if(out_vrad)      nprint = nprint+1
  if(out_dforce)    nprint = nprint+1
  if(out_sectional) nprint = nprint+1
  if(out_dmom)      nprint = nprint+1 
  if(out_pid)       nprint = nprint+1
  
  allocate(out_vars(nprint))
  allocate(out_vars_virtual(0))

  if(average) allocate(ave_out_vars(nprint))
  
  if(out_wake) then
    allocate(out_vars_w(nprint))
    allocate(out_vars_vp(nprint))
  endif


  ! Load the components (just once)
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)

  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
                                components_names, all_comp)

  call close_hdf5_file(floc)

  if(out_wake .and. average) call error(this_sub_name, this_mod_name, &
  'Cannot output an averaged wake visualization. Remove the wake or avoid &
  &averaging')

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

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

    !> Load the references
    call load_refs(floc,refs_R,refs_off)

    !> Move the points
    call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                filen = trim(filename_in) )
    !> Expand the actuator disks
    call expand_actdisk_postpro(comps, points, points_virtual, points_exp, points_virtual_exp, elems, elems_virtual)

    if(average .and. it .eq. an_avg) then
      ! Save the points of this iteration for the average visualization
        allocate(points_ave(size(points_exp,1),size(points_exp,2)))
        points_ave = points
    endif

    ! Load the results 
    ! TODO: better checks on optional variables
    if(out_sectional .and. out_surfvel) then
      call load_res(floc, comps, vort, press, t, force, moment, sectional, surfvel)
    else if(out_sectional) then
      call load_res(floc, comps, vort, press, t, force, moment, sec_force=sectional)
    else if(out_surfvel) then
      call load_res(floc, comps, vort, press, t, force, moment, surfvel=surfvel)
    else
      call load_res(floc, comps, vort, press, t, force, moment)
    endif
    
    !Prepare the variable for output
    nelem_out = size(vort)
    ! TODO: compute/or read pressure and velocity field. Now set equal to zero
    allocate(  vel(size( vort,1)) ) ; vel = 0.0_wp
    
    if(out_cp) then
      ! pressure coefficient
      allocate(   cp(size(press,1)) ) ;  cp = 0.0_wp
      if (norm2(u_inf) .gt. 1e-10_wp) then
        cp = (press - P_inf)/(0.5_wp*rho*norm2(u_inf)**2)
      else
        call error(this_sub_name, this_mod_name, &
              'Cp was requested, but u_inf is zero.')
      endif
    endif

    i_var = 1
    if(out_vort) then
      call add_output_var(out_vars(i_var), vort, 'Singularity_Intensity', &
                          .false.)
      i_var = i_var +1
    endif
    if(out_vrad) then
      call add_output_var(out_vars(i_var), vort, 'Vortex_Rad', &
                          .false.)
      i_var = i_var +1
    endif
    if(out_vort_vec) then
      call add_output_var(out_vars(i_var), vort, 'Vorticity_Vector', &
                          .false.)
      i_var = i_var +1
    endif    
    if(out_cp) then
      call add_output_var(out_vars(i_var), cp, 'Cp',.false.)
      i_var = i_var +1
    endif
    if(out_surfvel) then
      call add_output_var(out_vars(i_var), surfvel, 'Surface_Velocity',.false.)
      i_var = i_var +1
    endif
    if(out_vel) then
      call add_output_var(out_vars(i_var), vel, 'Velocity',.false.)
      i_var = i_var +1
    endif
    if(out_press) then
      call add_output_var(out_vars(i_var), press, 'Pressure',.false.)
      i_var = i_var +1
    endif
    if(out_dforce) then 
      call add_output_var(out_vars(i_var), force, 'Force',.false.)
      i_var = i_var +1
    endif
    if(out_sectional) then 
      call add_output_var(out_vars(i_var), sectional, 'Sectional',.false.)
      i_var = i_var +1
    endif
    if(out_dmom) then 
      call add_output_var(out_vars(i_var), moment, 'Moment',.false.)
      i_var = i_var +1
    endif

    
    if(average) then
      if( ires .eq. 1) then
        call copy_output_vars(out_vars, ave_out_vars, .true.)
      endif
      do i_var = 1,size(ave_out_vars)
      ave_out_vars(i_var)%var = ave_out_vars(i_var)%var * &
                                (real(ires-1,wp)/real(ires,wp)) + &
                                out_vars(i_var)%var/real(ires,wp)
      enddo
    endif

    if(.not. average) then
      ! Output filename
      write(filename,'(A,I4.4)') trim(basename)//'_'//trim(an_name)//'-',it
      write(filename_virtual,'(A,I4.4)') trim(basename)//'_virtual_'//trim(an_name)//'-',it
      if (out_wake) then

        if(out_turbvisc) then
          call load_wake_viz(floc, wpoints, welems, wvort, vppoints, vpvort, &
                            vpvort_v, v_rad, vp_parentid, vpturbvisc)
        else
          call load_wake_viz(floc, wpoints, welems, wvort, vppoints, vpvort, &
                              vpvort_v, v_rad, vp_parentid)
        endif
        nelem_w = size(welems,2)
        nelem_vp = size(vppoints,2)

        i_var = 1
        if(out_vort) then
          call add_output_var(out_vars_w(i_var), wvort, &
                  'Singularity_Intensity',.false.)
          call add_output_var(out_vars_vp(i_var), vpvort, &
                  'Singularity_Intensity',.false.)
          i_var = i_var +1
        endif
        if(out_cp) then
          call add_output_var(out_vars_w(i_var), cp, &
                  'Cp',.true.)
          call add_output_var(out_vars_vp(i_var), cp, &
                  'Cp',.true.)
          i_var = i_var +1
        endif
        if(out_surfvel) then
          call add_output_var(out_vars_w(i_var), surfvel, &
                  'Surface_Velocity',.true.)
          call add_output_var(out_vars_vp(i_var), surfvel, &
                  'Surface_Velocity',.true.)
          i_var = i_var +1
        endif
        if(out_vel) then
          call add_output_var(out_vars_w(i_var), vel, &
                  'Velocity',.true.)
          call add_output_var(out_vars_vp(i_var), vel, &
                  'Velocity',.true.)
          i_var = i_var +1
        endif
        if(out_press) then
          call add_output_var(out_vars_w(i_var), press, &
                  'Pressure',.true.)
          call add_output_var(out_vars_vp(i_var), press, &
                  'Pressure',.true.)
          i_var = i_var +1
        endif
        if(out_turbvisc) then
          call add_output_var(out_vars_w(i_var), vpturbvisc, &
                  'Turbulent_Viscosity',.true.)
          call add_output_var(out_vars_vp(i_var), vpturbvisc, &
                  'Turbulent_Viscosity',.false.)
          call add_output_var(out_vars(i_var), vpturbvisc, &
                  'Turbulent_Viscosity',.true.)
          i_var = i_var +1
        endif
        if(out_vort_vec) then
          call add_output_var(out_vars_vp(i_var), vpvort_v, &
                  'Vorticity',.false.)
          call add_output_var(out_vars_w(i_var), vpvort_v, &
                  'Vorticity',.true.)

          i_var = i_var +1
        endif
        if(out_vrad) then
          call add_output_var(out_vars_vp(i_var), v_rad, &
                  'VortexRad',.false.)
          call add_output_var(out_vars_w(i_var), v_rad, &
                  'VortexRad',.true.)
          !call add_output_var(out_vars(i_var), v_rad, &
          !        'VortexRad',.true.)
          i_var = i_var +1
        endif
        if (out_dforce) then 
          call add_output_var(out_vars_w(i_var), force, &
                              'Force',.true.)
          call add_output_var(out_vars_vp(i_var), force, &
                              'Force',.true.)
          i_var = i_var +1
        endif
        if (out_dmom) then 
          call add_output_var(out_vars_w(i_var), moment, &
                              'Moment', .true.)
          call add_output_var(out_vars_vp(i_var), moment, &
                              'Moment',.true.)
          i_var = i_var +1
        endif 
        if (out_sectional) then 
          call add_output_var(out_vars_w(i_var), sectional, &
                              'Sectional', .true.)
          call add_output_var(out_vars_vp(i_var), sectional, &
                              'Sectional',.true.)
          i_var = i_var +1
        endif
        if(out_pid) then
          call add_output_var(out_vars_vp(i_var), vp_parentid, &
                  'ParentId',.false.)
          call add_output_var(out_vars_w(i_var), vp_parentid, &
                  'ParentId',.true.)
          call add_output_var(out_vars(i_var), vp_parentid, &
                  'ParentId',.true.)
          i_var = i_var +1
        endif
        
        !Output the results (with wake)
        select case (trim(out_frmt))
          case ('tecplot')
            filename = trim(filename)//'.plt'
            filename_virtual = trim(filename_virtual)//'.plt'
            call  tec_out_viz(filename, t, &
                      points_exp, elems, out_vars, &
                      w_rr=wpoints, w_ee=welems, w_vars=out_vars_w, &
                      vp_rr=vppoints, vp_vars=out_vars_vp)
            call  tec_out_viz(filename_virtual, t, &
                      points_virtual_exp, elems_virtual, out_vars_virtual)
          case ('vtk')
            filename = trim(filename)//'.vtu'
            filename_virtual = trim(filename_virtual)//'.vtu'
            call  vtk_out_viz(filename, &
                        points_exp, elems, out_vars, &
                      w_rr=wpoints, w_ee=welems, w_vars=out_vars_w, &
                      vp_rr=vppoints, vp_vars=out_vars_vp, &
                      separate_wake = separate_wake)
            call  vtk_out_viz(filename_virtual, &
                      points_virtual_exp, elems_virtual, out_vars_virtual)
          case default
            call error('dust_post','','Unknown format '//trim(out_frmt)//&
                      ' for visualization output')
        end select

        deallocate (wpoints, welems,  wvort)
        call clear_output_vars(out_vars_w)
        call clear_output_vars(out_vars_vp)

      else

        !Output the results (without wake)
        select case (trim(out_frmt))
          case ('tecplot')
            filename = trim(filename)//'.plt'
            filename_virtual = trim(filename_virtual)//'.plt'
            call  tec_out_viz(filename, t, &
                          points_exp, elems, out_vars)
            call  tec_out_viz(filename_virtual, t, &
                          points_virtual_exp, elems_virtual, out_vars_virtual)
          case ('vtk')
            filename = trim(filename)//'.vtu'
            filename_virtual = trim(filename_virtual)//'.vtu'
            call  vtk_out_viz(filename, &
                          points_exp, elems, out_vars)
            call  vtk_out_viz(filename_virtual, &
                          points_virtual_exp, elems_virtual, out_vars_virtual)                          
          case default
            call error('dust_post','','Unknown format '//trim(out_frmt)//&
                        ' for visualization output')
        end select

      endif !output wake

    endif !not average

    call close_hdf5_file(floc)

    deallocate(refs_R, refs_off)
    !deallocate(print_var_names, print_vars)
    call clear_output_vars(out_vars)

    if (allocated(vort   ) ) deallocate(vort )
    if (allocated(press  ) ) deallocate(press)
    if (allocated(surfvel) ) deallocate(surfvel)
    if (allocated(vel    ) ) deallocate(vel  )
    if (allocated(cp     ) ) deallocate(cp   )
    if (allocated(v_rad  ) ) deallocate(v_rad)
    if (allocated(vpturbvisc)) deallocate(vpturbvisc)
    if (allocated(force  ) ) deallocate(force)
    if (allocated(moment)) deallocate(moment)

  end do ! Time loop

  !Print the average
  if(average) then
    write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_ave'
    !Output the results (without wake)
    select case (trim(out_frmt))
      case ('tecplot')
        filename = trim(filename)//'.plt'
        call  tec_out_viz(filename, t, &
                    points_ave, elems, ave_out_vars)
      case ('vtk')
        filename = trim(filename)//'.vtu'
        call  vtk_out_viz(filename, &
                    points_ave, elems, ave_out_vars)
      case default
        call error('dust_post','','Unknown format '//trim(out_frmt)//&
                    ' for visualization output')
    end select
    call clear_output_vars(ave_out_vars)
    deallocate(ave_out_vars)
    deallocate(points_ave)
  endif

  deallocate(points, points_virtual)
  call destroy_elements(comps)
  deallocate(comps)
  deallocate(out_vars)
  if(out_wake) deallocate(out_vars_w, out_vars_vp)


  write(msg,'(A,I0,A)') nl//'++++++++++ Visualization done'//nl
  call printout(trim(msg))

end subroutine post_viz

! ----------------------------------------------------------------------

end module mod_post_viz
