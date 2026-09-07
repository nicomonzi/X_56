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

!> Module containing the subroutines to perform output of specific data
!! required for aeroacoustics computations
module mod_post_aa

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi, ascii_real

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
  read_hdf5, &
  read_hdf5_attr

use mod_stringtools, only: &
  LowCase, isInList

use mod_geometry, only: &
  t_geo, t_geo_component, destroy_elements

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro, prepare_geometry_postpro, &
  expand_actdisk_postpro

use mod_post_load, only: &
  load_refs, load_res, load_wake_viz , check_if_components_exist

use mod_dat_out, only: &
  dat_out_aa_header, dat_out_aa

use mod_surfpan, only: &
  t_surfpan

implicit none

public :: post_aeroacoustics

private

character(len=max_char_len), parameter :: this_mod_name='mod_post_aa'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_aeroacoustics( sbprms, basename, data_basename, an_name, ia, &
                      out_frmt, components_names, all_comp, an_start, an_end, &
                      an_step, average )
  type(t_parse), pointer                                   :: sbprms
  character(len=max_char_len) , intent(in)                 :: basename
  character(len=max_char_len) , intent(in)                 :: data_basename
  character(len=max_char_len) , intent(in)                 :: an_name
  integer          , intent(in)                            :: ia
  character(len=max_char_len) , intent(in)                 :: out_frmt
  character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
  logical , intent(in)                                     :: all_comp
  integer , intent(in)                                     :: an_start , an_end , an_step
  logical, intent(in)                                      :: average

  type(t_geo_component), allocatable  :: comps(:)
  character(len=max_char_len)         :: filename
  integer(h5loc)                      :: floc , ploc
  real(wp), allocatable               :: points(:,:), points_virtual(:,:)
  integer                             :: nelem, nelem_virtual, n_time
  real(wp)                            :: time
  real(wp)                            :: p_inf, a_inf, rho_inf, mu_inf
  real(wp)                            :: u_inf(3)
  integer                             :: fid_out, fid_time
  integer                             :: i_comp, ie, ierr
  real(wp), allocatable               :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable               :: refs_G(:,:,:), refs_f(:,:)
  real(wp), allocatable               :: vort(:), press(:), surfvel(:,:)
  integer                             :: it, ires, imult, num_el
  real(wp)                            :: t
  integer                             :: iv_u, iv_l
  real(wp)                            :: cen_u(3), cen_l(3)
  real(wp)                            :: nor_u(3), nor_l(3)
  real(wp)                            :: vel_u(3), vel_l(3)
  real(wp)                            :: area_u, area_l
  logical                             :: mult, isopen
  character(len=max_char_len)         :: comp_root, last_mult, compname
  character(len=max_char_len), parameter  :: this_sub_name='post_aeroacoustics'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' aeroacoustics'//nl
  call printout(trim(msg))

  ! Load the components (just once)
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)

  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
                                components_names, all_comp)

  call close_hdf5_file(floc)


  if( average ) call error(this_sub_name, this_mod_name, &
    'Cannot output an averaged aeroacoustics analysis')
  if(trim(out_frmt) .ne. 'dat') call error(this_sub_name, this_mod_name, &
    'aeroacoustics analysis available only in dat ascii file format')

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  ! Output time filename
  write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_time.dat'
  call new_file_unit(fid_time, ierr)
  open(unit=fid_time,file=trim(filename))

  ! Time loop
  ires = 0
  n_time = (an_end-an_start)/an_step + 1 ! int general eger division 
  do it = an_start, an_end, an_step
    ires = ires+1

    ! Open the file
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
    call open_hdf5_file(trim(filename),floc)

    ! Load free-stream parameters
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(time,'time',floc)
    call read_hdf5(rho_inf, 'rho_inf', ploc)
    call read_hdf5(mu_inf, 'mu_inf', ploc)
    call read_hdf5(u_inf, 'u_inf', ploc)
    call read_hdf5(a_inf, 'a_inf', ploc)
    call read_hdf5(p_inf, 'P_inf', ploc)   

    call close_hdf5_group(ploc)

    write(fid_time, '('//ascii_real//')') time

    ! Load the references
    call load_refs(floc, refs_R, refs_off, refs_G, refs_f)

    ! Move the points
    call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, refs_G, refs_f, &
                                filen = trim(filename) )

    !Load the results
    call load_res(floc, comps, vort, press, t, surfvel=surfvel)
    ! num_el is the total number of elements (lines of the dat file, useful for parsing)
    num_el = 0
    do i_comp = 1,size(comps)   
      num_el = num_el + size(comps(i_comp)%el)
    enddo 

    !cycle on components
    do i_comp = 1,size(comps)
    associate(comp => comps(i_comp))

        ! check if the file unit is still open from a previous file
        last_mult = trim(compname)

        !open the file
        if (i_comp .eq. 1) then
          ! Output filename
          write(filename,'(A,I0,A)') trim(basename)//'_'//trim(an_name)//&
          '-',it,'.dat'
          call new_file_unit(fid_out, ierr)
          open(unit=fid_out,file=trim(filename))
          !write the header here
          call dat_out_aa_header( fid_out, num_el, time, p_inf, rho_inf, a_inf, mu_inf, u_inf)

        endif

      !cycle on elements
      do ie = 1,size(comp%el)

        cen_u = 0.0_wp; cen_l = 0.0_wp
        nor_u = 0.0_wp; nor_l = 0.0_wp
        vel_u = 0.0_wp; vel_l = 0.0_wp
        area_u = 0.0_wp; area_l = 0.0_wp

        call get_vl_virtual_ids(comp, ie, iv_u, iv_l)

        if (iv_u .gt. 0) then
          cen_u = comp%el_virtual(iv_u)%cen
          nor_u = comp%el_virtual(iv_u)%nor
          vel_u = comp%el_virtual(iv_u)%ub
          area_u = comp%el_virtual(iv_u)%area
        end if

        if (iv_l .gt. 0) then
          cen_l = comp%el_virtual(iv_l)%cen
          nor_l = comp%el_virtual(iv_l)%nor
          vel_l = comp%el_virtual(iv_l)%ub
          area_l = comp%el_virtual(iv_l)%area
        end if

        call dat_out_aa( fid_out, comp%el(ie)%cen, comp%el(ie)%nor, &
                comp%el(ie)%ub, comp%el(ie)%area, comp%el(ie)%pres, rho_inf, &
                cen_u, cen_l, nor_u, nor_l, area_u, area_l, vel_u, vel_l)
      enddo
    end associate
    enddo

    if (.not. mult .or. i_comp .eq. size(comps)) then
      !close the file
      close(fid_out)
    endif

    deallocate(refs_R, refs_off, refs_G, refs_f)

    if (allocated(vort ) ) deallocate(vort )
    if (allocated(press) ) deallocate(press)
    if (allocated(press) ) deallocate(surfvel)

    call close_hdf5_file(floc)

  end do ! Time loop

  close(fid_time)

  deallocate(points)
  call destroy_elements(comps)
  deallocate(comps)

  write(msg,'(A,I0,A)') nl//'++++++++++ Aeroacoustics done'//nl
  call printout(trim(msg))

end subroutine post_aeroacoustics

! ----------------------------------------------------------------------

subroutine get_vl_virtual_ids(comp, ie_vl, i_up, i_low)
  type(t_geo_component), intent(in) :: comp
  integer, intent(in)               :: ie_vl
  integer, intent(out)              :: i_up, i_low

  integer                           :: n_span, n_chord_vl
  integer                           :: n_chord_virtual_tot, n_chord_virtual_half
  integer                           :: i_span, i_chord, i_chord_virtual
  real(wp)                          :: eta

  i_up = 0
  i_low = 0

  if (trim(comp%comp_el_type) .ne. 'v') return

  n_span = comp%parametric_nelems_span
  n_chord_vl = comp%parametric_nelems_chor
  n_chord_virtual_tot = comp%parametric_nelems_chor_virtual

  if (n_span .le. 0 .or. n_chord_vl .le. 0 .or. n_chord_virtual_tot .le. 1) return
  if (mod(n_chord_virtual_tot,2) .ne. 0) return

  n_chord_virtual_half = n_chord_virtual_tot / 2
  if (n_chord_virtual_half .le. 0) return

  i_span = (ie_vl-1) / n_chord_vl + 1
  i_chord = ie_vl - (i_span-1) * n_chord_vl

  if (i_span .lt. 1 .or. i_span .gt. n_span) return

  if (n_chord_vl .eq. 1) then
    i_chord_virtual = 1
  else
    eta = real(i_chord-1,wp) / real(n_chord_vl-1,wp)
    i_chord_virtual = 1 + nint(eta * real(n_chord_virtual_half-1,wp))
  end if

  i_chord_virtual = max(1, min(n_chord_virtual_half, i_chord_virtual))

  ! Virtual panels are ordered as lower TE->LE then upper LE->TE.
  i_low = (i_span-1) * n_chord_virtual_tot + (n_chord_virtual_half - i_chord_virtual + 1)
  i_up = (i_span-1) * n_chord_virtual_tot + (n_chord_virtual_half + i_chord_virtual)

  if (i_low .lt. 1 .or. i_low .gt. size(comp%el_virtual)) i_low = 0
  if (i_up .lt. 1 .or. i_up .gt. size(comp%el_virtual)) i_up = 0

end subroutine get_vl_virtual_ids

! ----------------------------------------------------------------------

function is_multiple(comp_name, name_root, imult) result(ismult)
  character(len=max_char_len), intent(in)  :: comp_name
  character(len=max_char_len), intent(out) :: name_root
  integer, intent(out)                     :: imult
  logical                                  :: ismult

  character(len=max_char_len), parameter   :: nums='0123456789'
  character(len=max_char_len), parameter   :: underscore='_'
  integer                                  :: strlen, check

  ismult = .false.
  strlen = len_trim(comp_name)

  if (strlen .lt. 5) then
    return
  endif

  check = verify(comp_name(strlen-3:strlen-2), underscore)
  check = check + verify(comp_name(strlen-1:strlen),nums)

  if (check .eq. 0) then
    ismult = .true.
    name_root = comp_name(1:strlen-4)
    read(comp_name(strlen-1:strlen),*) imult
    imult = imult-1
  endif
end function


end module mod_post_aa
