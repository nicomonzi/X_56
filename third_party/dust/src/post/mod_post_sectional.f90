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
!!          Alessandro Cocco
!!=========================================================================

!> Module containing the subroutines to perform sectional loads
!! analysis during postprocessing
module mod_post_sectional

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, internal_error, warning, info, printout, dust_time, t_realtime, &
  new_file_unit

use mod_parse, only: &
  t_parse, &
  getrealarray, getreal, getint, getlogical, &
  getsuboption,  countoption

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

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro

use mod_geometry, only: &
  t_geo, t_geo_component, destroy_elements

use mod_post_load, only: &
  load_refs, load_res, load_ll, load_vl

use mod_dat_out, only: &
  dat_out_sectional, dat_out_sectional_ll,  dat_out_sectional_vl

use mod_tecplot_out, only: &
  tec_out_sectional

use mod_math, only: &
  cross, vec2mat

implicit none

! type rof box sectional loads
type t_box_secloads
  integer                :: nelems
  integer , allocatable  :: elems(:)
  real(wp), allocatable  :: fracs(:)
  real(wp), allocatable  :: cen(:,:)
end type t_box_secloads

public :: post_sectional

private

character(len=*), parameter :: this_mod_name = 'mod_post_sectional'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------
!TODO (to be tested):
! - parametric components aligned with y-axis
! - tol_y_cen is hard-coded. write tol_y_cen as an input
! - preCICE coupled components 
subroutine post_sectional (sbprms, bxprms, basename, data_basename, an_name, &
                          ia, out_frmt, components_names, all_comp, &
                          an_start, an_end, an_step, average )
type(t_parse), pointer                                  :: sbprms
type(t_parse), pointer                                  :: bxprms
character(len=*) , intent(in)                           :: basename
character(len=*) , intent(in)                           :: data_basename
character(len=*) , intent(in)                           :: an_name
integer          , intent(in)                           :: ia
character(len=*) , intent(in)                           :: out_frmt
character(len=max_char_len), allocatable, intent(inout) :: components_names(:)
logical , intent(in)                                    :: all_comp
integer , intent(in)                                    :: an_start , an_end , an_step
logical , intent(in)                                    :: average

type(t_geo_component), allocatable                      :: comps(:)
character(len=max_char_len)                             :: cname
integer(h5loc)                                          :: floc, gloc, cloc
real(wp), allocatable                                   :: refs_R(:,:,:), refs_off(:,:)
real(wp)                                                :: R_cen(3,3)
real(wp), allocatable                                   :: vort(:), cp(:)
real(wp), allocatable                                   :: ll_data(:,:,:), ll_data_ave(:,:,:)
real(wp), allocatable                                   :: vl_data(:,:,:), vl_data_ave(:,:,:)
real(wp), allocatable                                   :: points(:,:), points_virtual(:,:)
integer                                                 :: nelem, nelem_virtual
integer                                                 :: n_comp , n_comp_tot  
integer                                                 :: i_comp , id_comp , ax_coor , ref_id
character(len=max_char_len), allocatable                :: all_components_names(:)
character(len=max_char_len), allocatable                :: components_names_tmp(:)
integer                                                 :: nelem_span , nelem_chor , n_sect , n_time

real(wp), allocatable                                   :: time(:)
real(wp)                                                :: t
integer                                                 :: ires

!> sectional loads: paramteric components -------------------
real(wp)                                                :: axis_dir(3) , axis_nod(3)
real(wp) , allocatable                                  :: sec_loads(:,:,:)
real(wp) , allocatable                                  :: sec_loads_ave(:,:,:)
integer                                                 :: is                                 
integer, parameter                                      :: n_loads = 4   ! F and moment around an axis
integer, parameter                                      :: n_ll_data = 9, n_vl_data = 12
real(wp) , allocatable                                  :: ref_mat(:,:) , off_mat(:,:)
real(wp) , allocatable                                  :: y_cen(:) , y_span(:), chord(:)
real(wp) , parameter                                    :: tol_y_cen = 1.0e-3_wp
real(wp) , allocatable                                  :: r_axis(:,:) , r_axis_bas(:,:)
real(wp)                                                :: F_bas(3) , F_bas1(3)
real(wp)                                                :: M_bas(3)
integer                                                 :: ie , ic , it , i1
logical                                                 :: print_ll, print_vl

!> sectional loads: box -------------------------------------
character(len=max_char_len)                             :: comp_input
integer                                                 :: n_box
logical                                                 :: reshape_box
real(wp)                                                :: nCoordMaxBox
real(wp) , allocatable                                  :: box_coord(:,:) , box_coord_tmp(:,:)
real(wp)                                                :: vVec(3) , bVec(3) , wVec(3)  
real(wp)                                                :: baseLen(2) , heigLen(2) , spanLen
real(wp)                                                :: nVec(3) , ref_node(3) 
real(wp)                                                :: nCoordMax , nCoordMin
real(wp)                                                :: nCoord , dnCoord
real(wp) , allocatable                                  :: nCoordSec(:)
real(wp) , allocatable                                  :: nCoordCen(:,:)
real(wp)                                                :: normalLateralFaces(3,4)
real(wp)                                                :: refPoiLateralFaces(3,4)
real(wp)                                                :: b1Vec(3) , b2Vec(3)
type(t_box_secloads) , allocatable                      :: box_secloads(:)
real(wp)                                                :: distance(4)
real(wp) , parameter                                    :: nCoord_relOff = 0.05_wp
real(wp)                                                :: ac(3)
! "extra-param". TODO: check if they are called before #### box
real(wp)                                                :: nCoordVert(4)   ! TRI o QUAD elements
integer                                                 :: secVert(4)   ! TRI o QUAD elements
real(wp)                                                :: interSectPoints(3,2)
real(wp)                                                :: interSectAreas(2) , interSectCen(3,2)
integer                                                 :: nInterSect , index2
real(wp)                                                :: node1(3) , node2(3) , box_secloads_cen(3)
integer                                                 :: iSec1 , iSec2 , sec1_nVer , sec2_nVer
integer                                                 :: nNodeInt(2,2)
integer                                                 :: nver , iv, chord_start, chord_end, ii

character(len=max_char_len)                             :: filename

character(len=*), parameter :: this_sub_name = 'post_sectional'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' sectional loads'//nl
  call printout(trim(msg))


  ! Some warnings and errors -------------------------------
  ! WARNING: sectional loads are computed for one components only
  !  at time. If no component is specified --> return
  n_comp = countoption(sbprms,'component')
  if ( n_comp .le. 0 ) then
    call warning(this_mod_name, this_sub_name, 'No component specified for &
                &sectional_loads analysis. Skipped analysis.')
    return
  else if ( n_comp .ge. 2 ) then
    call warning(this_mod_name, this_sub_name, &
        'More than one component specified &
      &for ''sectional_loads'' analysis: just the first one is considered. &
      &Please run another ''sectional_analysis'' if you need it on more than &
      &component.')
  end if

  !> check if components_names is allocated (it should alwasy be allocated)
  if ( .not. allocated(components_names) ) then
    call internal_error(this_mod_name,this_sub_name, &
                        'components_name not allocated.')
  end if

  !------------ Check if the component really exists ---------- !
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp_tot,'NComponents',gloc)

  allocate(all_components_names(n_comp_tot))
  do i_comp = 1 , n_comp_tot
    write(cname,'(A,I3.3)') 'Comp',i_comp
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(all_components_names(i_comp),'CompName',cloc)    
    call close_hdf5_group(cloc)
  end do
  call close_hdf5_group(gloc)

  id_comp = -333
  do i_comp = 1 , n_comp_tot
      if ( trim(components_names(1)) .eq. trim(all_components_names(i_comp)) ) then
        id_comp = i_comp
      end if
  end do
  if ( id_comp .eq. -333 ) then ! component not found
    call printout('  Component not found. The available components are: ')
    do i_comp = 1 , n_comp_tot
      write(msg,'(A,I0,A)') '  comp. ' , i_comp , &
                          ':  '//trim(all_components_names(i_comp))
      call printout(trim(msg))
    end do
    call error(this_sub_name,this_mod_name,'No valid component found.')
  end if

  write(msg,*) '  Analysing component : ' , trim(components_names(1))//nl
  call printout(trim(msg))
  
  !> load only the first component just once 
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)

  allocate(components_names_tmp(1))
  components_names_tmp(1) = trim(components_names(1))

  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
                              components_names_tmp,  all_comp)
  call close_hdf5_file(floc)

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  ! Read the axis for the computation of the sectional loads
  axis_dir = getrealarray(sbprms,'axis_dir',3)
  axis_nod = getrealarray(sbprms,'axis_nod',3)

  ! if a box is defined ---> use the box
  comp_input = trim(comps(1)%comp_input )
  n_box = countoption(sbprms,'box_sect')
  if ( n_box .eq. 1 ) comp_input = 'cgns'
  if (n_box .gt. 1) call error(this_sub_name, this_mod_name, 'More than one&
  & box defined for sectional loads')

  !> decide on printing lifting lines data
  print_ll = getlogical(sbprms,'lifting_line_data')
  if (print_ll .and. trim(comps(1)%comp_el_type) .ne. 'l')  then
    print_ll = .false.
    call warning(this_sub_name, this_mod_name, 'Cannot output lifting &
          &line data for a non lifting line component, output of lifting &
          &line data skipped')
  endif

  !> decide on printing corrected vortex lattice data
  print_vl = getlogical(sbprms,'vortex_lattice_data')
  
  if (print_vl .and. (trim(comps(1)%comp_el_type) .ne. 'v' .or.  &
    trim(comps(1)%aero_correction) .ne. 'true'))  then
    print_vl = .false.
    call warning(this_sub_name, this_mod_name, 'Cannot output corrected &
          &vortex lattice data for a non vortex lattice component, output of vortex &
          &lattice data skipped')
  endif


  select case( trim(comp_input) )
  
  case ( 'parametric', 'pointwise' )

    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
    ax_coor = 2   ! 2. the parametric elements is defined
                  !    along the y-axis

    nelem_span = comps(id_comp)%parametric_nelems_span
    nelem_chor = comps(id_comp)%parametric_nelems_chor
    n_sect = nelem_span

    do i1 = 1, size(comps(id_comp)%loc_points, 2) 
      comps(id_comp)%loc_points(:, i1) = matmul(comps(id_comp)%coupling_node_rot,comps(id_comp)%loc_points(:,i1))
    end do
    
    ! ######################################################################
    ! Find the coordinates of the reference points in the local reference frame

    ! Find the coordinate along the axis of the sections -----------  --------
    ! the <ax_coor> coordinate y_cen of each section is built with the first
    ! panel in chord. Then the distance between the <ax_coor> of the centre
    ! of the other panels and y_cen is checked
    !TODO: find y-coord is general enough for all the 'parametric' components
    allocate(y_cen(n_sect),y_span(n_sect),chord(n_sect))
    ie = 0
    chord_start = 0
    chord_end = 0
    ii = 0
    do is = 1 , n_sect
      
      ie = ie + 1
      ii = ii + 1
      chord_start = (ii-1)*(nelem_chor + 1) + 1
      chord_end = ii*(nelem_chor + 1) 
      y_cen(is) = &
      sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver)) / &
                  real(comps(id_comp)%el(ie)%n_ver,wp)
      y_span(is) = &
      abs ( comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver(1) )&
      - comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver(2) ) )
      !> local chord as projection of the profile on x-z plane
      chord(is) = sqrt((minval(comps(id_comp)%loc_points(ax_coor - 1,chord_start:chord_end)) - & 
                        maxval(comps(id_comp)%loc_points(ax_coor - 1,chord_start:chord_end)))**2  + &
                        (minval(comps(id_comp)%loc_points(ax_coor + 1,chord_start:chord_end)) - &  
                        maxval(comps(id_comp)%loc_points(ax_coor + 1,chord_start:chord_end)))**2)
      do ic = 2 , nelem_chor ! check
        ie = ie + 1
        if (abs( y_cen(is) - &
          sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver))&
            / real(comps(id_comp)%el(ie)%n_ver,wp) ) &
            .gt. tol_y_cen ) then
          call error(this_sub_name,this_mod_name,'Wrong section definition &
          &for component '//trim(comps(id_comp)%comp_name)//' for sectional&
          & loads analysis.')
        end if
      end do
    end do

    ! Find the coordinate of the reference points on the axis --------------
    !  ( with coord. y_cen )
    if ( abs(axis_dir(2)) .lt. 1e-6_wp ) then
      call error('dust_post','','Wrong definition of the axis in&
            & sectional_loads analysis: abs(axis_dir(2)) .lt. 1e-6.&
            & STOP')
    end if

    allocate(r_axis(3,n_sect), r_axis_bas(3,n_sect))
    do is = 1 , n_sect
      r_axis(:,is) = axis_nod + &
              ( y_cen(is) - axis_nod(ax_coor) )/axis_dir(ax_coor) * axis_dir
    end do
    ! only initialisation here:
    ! - r_axis: coordinates in the local ref.frame
    ! - r_axis_bas: coordinates in the base ref.frame
    r_axis_bas = r_axis

    ! Find the ref_id or the reference frame where loads are projected -----
    ref_id = comps(id_comp)%ref_id
    write(msg,'(A,I0,A)') '   Employing the local reference frame for moment &
      &and force projection.'//nl//'   Reference tag: '&
      //trim(comps(id_comp)%ref_tag)//'  (ID: ',ref_id,')'
    call printout(trim(msg))

    ! allocate tmp array to store the results --------------
    if (an_end .lt. an_start) then 
      call error('dust_post', 'Wrong definition in time vector')
    endif
    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( sec_loads(n_time,n_sect,n_loads) )
    sec_loads = -333.0_wp ! initialisation for DEBUG
    allocate( time(n_time) ) ; time = -333.0_wp
    allocate( ref_mat(n_time,9) , off_mat(n_time,3) )

    if(print_ll) then  
      allocate(ll_data(n_time,n_sect,n_ll_data))
    elseif (print_vl) then  
        allocate(vl_data(n_time,n_sect,n_vl_data))
    endif 
  
    ires = 0
    do it = an_start, an_end, an_step
      ires = ires + 1

      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)

      ! Load the references and move the points ---
      call load_refs(floc, refs_R, refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                  filen = trim(filename) )
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)
      
      if(print_ll) then  
        call load_ll(floc, comps, ll_data(ires,:,:))
      elseif(print_vl) then  
        call load_vl(floc, comps, vl_data(ires,:,:))
      endif 

      call close_hdf5_file(floc)

      ! from local coordinates to coordinates in the base ref.frame
      do is = 1 , n_sect
        r_axis_bas(:,is) = matmul( refs_R(:,:,ref_id) , r_axis(:,is) ) + &
                          refs_off(:,ref_id)
      end do

      ! compute sectional loads. Loop over the panels of each section
      ie = 0
      do is = 1 , n_sect ! loop over sections
        ! force
        F_bas = 0.0_wp ; M_bas = 0.0_wp
        do ic = 1 , nelem_chor ! loop over chord
          ie = ie + 1
          F_bas1 = comps(id_comp)%el(ie)%dforce
          F_bas  = F_bas + F_bas1
#if USE_PRECICE
          if (comps(id_comp)%coupling) then 
            call vec2mat(comps(id_comp)%el(ie)%ori, R_cen)
            R_cen = matmul(comps(id_comp)%coupling_node_rot, R_cen)
            if (trim(comps(id_comp)%comp_el_type) .eq. 'l') then 
              ac = sum ( comps(id_comp)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
              M_bas  = M_bas + cross( matmul(transpose(R_cen), ac), F_bas1) &
                            + comps(id_comp)%el(ie)%dmom 
            else 
              M_bas  = M_bas + cross( matmul(transpose(R_cen), &
                              comps(id_comp)%el(ie)%cen) , F_bas1 )      &
                              + comps(id_comp)%el(ie)%dmom 
            endif
          else 
            if (trim(comps(id_comp)%comp_el_type) .eq. 'l') then 
              ac = sum ( comps(id_comp)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
              M_bas  = M_bas + cross(ac - r_axis_bas(:,is) , F_bas1 )      &
                        + comps(id_comp)%el(ie)%dmom 
            else
              M_bas  = M_bas + cross( comps(id_comp)%el(ie)%cen &
                        - r_axis_bas(:,is) , F_bas1 )      &
                        + comps(id_comp)%el(ie)%dmom 
            endif 
          endif
#else 
          if (trim(comps(id_comp)%comp_el_type) .eq. 'l') then 
            ac = sum ( comps(id_comp)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
            M_bas  = M_bas + cross(ac - r_axis_bas(:,is) , F_bas1 ) &
                      + comps(id_comp)%el(ie)%dmom 
          else
            M_bas  = M_bas + cross( comps(id_comp)%el(ie)%cen &
                      - r_axis_bas(:,is) , F_bas1 )           &
                      + comps(id_comp)%el(ie)%dmom 
          endif 
#endif 
        end do ! loop over chord

        ! From global to local coordinates of forces and moments
#if USE_PRECICE
        if (comps(id_comp)%coupling) then 
          call vec2mat(comps(id_comp)%el(ie)%ori, R_cen)
          R_cen = matmul(comps(id_comp)%coupling_node_rot, R_cen)
          sec_loads(ires,is,1:3) = matmul( &
                transpose( R_cen), F_bas ) / y_span(is)
          ! moment ( only the component around the <ax_coord> )
          M_bas = matmul( transpose(R_cen) , M_bas)
          sec_loads(ires,is,4) = M_bas(ax_coor) / y_span(is)
        else
          sec_loads(ires,is,1:3) = matmul( &
            transpose( refs_R(:,:, ref_id) ) , F_bas ) / y_span(is)
          ! moment ( only the component around the <ax_coord> )
          M_bas = matmul( &
                transpose( refs_R(:,:, ref_id) ) , M_bas )
          sec_loads(ires,is,4) = M_bas(ax_coor) / y_span(is)
        endif  
#else 
        sec_loads(ires,is,1:3) = matmul( &
            transpose( refs_R(:,:, ref_id) ) , F_bas ) / y_span(is)
        ! moment ( only the component around the <ax_coord> )
        M_bas = matmul( &
                transpose( refs_R(:,:, ref_id) ) , M_bas )
        sec_loads(ires,is,4) = M_bas(ax_coor) / y_span(is)
#endif 
      end do ! loop over sections

      ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
      off_mat(ires,:) = refs_off(:,ref_id)
      time(ires) = t

    end do

    if(average) then
      allocate(sec_loads_ave(1,size(sec_loads,2),size(sec_loads,3)))
      sec_loads_ave(1,:,:) = sum(sec_loads, 1)/real(size(sec_loads,1),wp)
      if (print_ll) then
        allocate(ll_data_ave(1,size(ll_data,2),size(ll_data,3)))
        ll_data_ave(1,:,:) = sum(ll_data, 1)/real(size(ll_data,1),wp)
      elseif (print_vl) then
        allocate(vl_data_ave(1,size(vl_data,2),size(vl_data,3)))
        vl_data_ave(1,:,:) = sum(vl_data, 1)/real(size(vl_data,1),wp)
      endif 

      select case(trim(out_frmt))
        case('dat')
            write(filename,'(A)') trim(basename)//'_'//trim(an_name)
            call dat_out_sectional ( filename, components_names(1), y_cen, &
                                    y_span, chord, time(1:1), sec_loads_ave, ref_mat, &
                                    off_mat, average )
          if(print_ll) then 
            call dat_out_sectional_ll(filename, components_names(1), &
                                    y_cen, y_span, chord, time(1:1), ll_data_ave, average)
          elseif(print_vl) then
            call dat_out_sectional_vl(filename, components_names(1), &
                                    y_cen, y_span, chord, time(1:1), vl_data_ave, average)
          endif 
          
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_ave.plt'
          if (print_ll) then
            call tec_out_sectional (filename, time(1:1), sec_loads_ave, y_cen, &
                                    y_span, chord, ll_data_ave)
          elseif (print_vl) then
            call tec_out_sectional (filename, time(1:1), sec_loads_ave, y_cen, &
                                    y_span, chord, vl_data_ave)
          else 
            call tec_out_sectional (filename, time(1:1), sec_loads_ave, y_cen, &
                                                                      y_span, chord )
          endif
        end select
      deallocate(sec_loads_ave)
    else
      select case(trim(out_frmt))
        case('dat')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)
          call dat_out_sectional ( filename, components_names(1), y_cen, &
                                  y_span, chord, time, sec_loads, &
                                  ref_mat, off_mat, average )
          if(print_ll) then 
            call dat_out_sectional_ll (filename, components_names(1),&
                                        y_cen, y_span, chord, time, ll_data, average )
          elseif (print_vl) then 
            call dat_out_sectional_vl (filename, components_names(1),&
                                        y_cen, y_span, chord, time, vl_data, average )
          endif 
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
          if (print_ll) then
            call tec_out_sectional (filename, time, sec_loads, y_cen, y_span, chord, &
                                    ll_data)
          elseif (print_vl) then 
            call tec_out_sectional (filename, time, sec_loads, y_cen, y_span, chord, &
                                    vl_data)
          else
            call tec_out_sectional ( filename, time, sec_loads, y_cen, y_span, chord)
          endif
      end select
    endif

    deallocate(r_axis,r_axis_bas)
    deallocate(y_cen, y_span)
    deallocate(sec_loads, time, ref_mat, off_mat)



  case ( 'cgns' )

    !TODO: move to a subroutine

    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded

    ! Read input for box and build the 8 nodes of the box -----------------------------
    call getsuboption(sbprms,'BoxSect',bxprms)
    allocate(box_coord(8,3),box_coord_tmp(8,3))
    ref_node = getrealarray(bxprms,'ref_node',3)
    vVec = getrealarray(bxprms,'face_vec',3)
    baseLen = getrealarray(bxprms,'face_bas',2)
    heigLen = getrealarray(bxprms,'face_hei',2)
    bVec = getrealarray(bxprms,'span_vec',3)
    spanLen = getreal(bxprms,'span_len')
    n_sect = getint(bxprms,'num_sect')
    reshape_box = getlogical(bxprms,'reshape_box')
    nVec = getrealarray(sbprms,'axis_mom',3)

    ! Normalise
    vVec = vVec / norm2(vVec)
    bVec = bVec / norm2(bVec)
    nVec = nVec / norm2(nVec)

    wVec = cross(vVec,nVec) ; wVec = wVec / norm2(wVec)

    box_coord(1,:) = ref_node
    box_coord(2,:) = box_coord(1,:) + baseLen(1) * vVec
    box_coord(3,:) = box_coord(2,:) + heigLen(1) * wVec
    box_coord(4,:) = box_coord(1,:) + heigLen(1) * wVec

    box_coord(5,:) = box_coord(1,:) + spanLen * bVec
    box_coord(6,:) = box_coord(5,:) + baseLen(2) * vVec
    box_coord(7,:) = box_coord(6,:) + heigLen(2) * wVec
    box_coord(8,:) = box_coord(5,:) + heigLen(2) * wVec


    ! Find if the whole component is contained in the box and/or reshape
    !  the box  ( in the direction perpendicolar to the nVec = axis_mom )
    ! nCoord: distance of a point from the first reference plane
    ie = 1
    nCoordMax = sum ( ( comps(id_comp)%loc_points( :, ie ) &
     - box_coord(1,:) ) * nVec )
    nCoordMin = nCoordMax !0.0_wp
    do ie = 2 , size( comps(id_comp)%loc_points , 2 ) ! size(comps(id_comp)%el)
      nCoord = sum ( ( comps(id_comp)%loc_points( :, ie ) &
        - box_coord(1,:) ) * nVec )
      if ( nCoord .gt. nCoordMax ) nCoordMax = nCoord
      if ( nCoord .lt. nCoordMin ) nCoordMin = nCoord
    end do

    ! Add small offset to the min and max nCoord ---------------------
    dnCoord = ( nCoordMax - nCoordMin ) / real(n_sect,wp)
    nCoordMin = nCoordMin - dnCoord * nCoord_relOff
    nCoordMax = nCoordMax + dnCoord * nCoord_relOff
    dnCoord = ( nCoordMax - nCoordMin ) / real(n_sect,wp)

    !check. TODO: add the offset to min,max() in the following lines
    nCoordMaxBox = sum(bVec*nVec) * spanLen
    nCoordMin = max(nCoordMin,0.0_wp)
    nCoordMax = min(nCoordMax,nCoordMaxBox)

    ! Reshape the box, if it is "too large" --------------------------
    if ( reshape_box ) then

!     ! Update box corners: update only if the box is too large
!     allocate(box_coord_tmp(8,3))
      box_coord_tmp = 0.0_wp
      box_coord_tmp(1,:) = box_coord(1,:) + ( box_coord(5,:) - &
                           box_coord(1,:) ) * nCoordMin / spanLen
      box_coord_tmp(2,:) = box_coord(2,:) + ( box_coord(6,:) - &
                           box_coord(2,:) ) * nCoordMin / spanLen
      box_coord_tmp(3,:) = box_coord(3,:) + ( box_coord(7,:) - &
                           box_coord(3,:) ) * nCoordMin / spanLen
      box_coord_tmp(4,:) = box_coord(4,:) + ( box_coord(8,:) - &
                           box_coord(4,:) ) * nCoordMin / spanLen
      box_coord_tmp(5,:) = box_coord(1,:) + ( box_coord(5,:) - &
                           box_coord(1,:) ) * nCoordMax / spanLen
      box_coord_tmp(6,:) = box_coord(2,:) + ( box_coord(6,:) - &
                           box_coord(2,:) ) * nCoordMax / spanLen
      box_coord_tmp(7,:) = box_coord(3,:) + ( box_coord(7,:) - &
                           box_coord(3,:) ) * nCoordMax / spanLen
      box_coord_tmp(8,:) = box_coord(4,:) + ( box_coord(8,:) - &
                           box_coord(4,:) ) * nCoordMax / spanLen

      box_coord = box_coord_tmp

    end if ! end reshape

    ! Find the n-coord of the sections and the reference points at the centre
    ! of the sections
    ! obs: so far, the nCoord has been defined w.r.t. the original first face
    !      of the box. From now on, nCoordSec, nCoordCen are defined with
    !      respect to the first section, that can be moved
    allocate(nCoordSec(n_sect+1))
    allocate(nCoordCen(3,n_sect)) ; nCoordCen = 0.0_wp

    nCoordSec(1) = 0.0_wp ! + nCoordMin
    do is = 1 , n_sect
      nCoordSec(is+1) = (nCoordMax-nCoordMin) * real(is,wp) / real(n_sect,wp)
      nCoordCen(2,is) = (nCoordMax-nCoordMin) * ( real(is,wp) - 0.5_wp)  &
                                                            / real(n_sect,wp)
    end do

    !TODO: use generalised version of this formula.
    ! this is valid only if nVec = yVec
    ! TODO: ask if this offset is required: + ref_node(2)
    allocate(y_cen( n_sect)) ; y_cen = nCoordCen(2,:)
    allocate(y_span(n_sect)) ; y_span= nCoordSec(2:) - nCoordSec(1:n_Sect)

    allocate(r_axis(3,n_sect),r_axis_bas(3,n_sect))
    do is = 1 , n_sect
      r_axis(:,is) = axis_nod + &
              ( y_cen(is) - axis_nod(2) )/axis_dir(2) * axis_dir
    end do

    ! only initialisation here:
    ! - r_axis: coordinates in the local ref.frame
    ! - r_axis_bas: coordinates in the base ref.frame
    r_axis_bas = r_axis

    ! Define the unit normal vectors to the 4 planar lateral faces of the
    ! boxes and the reference point, used as the origin for measuring
    ! distance from the plane
    refPoiLateralFaces(:,1) = box_coord(1,:)
    refPoiLateralFaces(:,2) = box_coord(1,:)
    refPoiLateralFaces(:,3) = box_coord(3,:)
    refPoiLateralFaces(:,4) = box_coord(4,:)

    b1Vec = box_coord(6,:) - box_coord(2,:) ; b1Vec = b1Vec / norm2(b1Vec)
    b2Vec = box_coord(8,:) - box_coord(4,:) ; b2Vec = b2Vec / norm2(b2Vec)
    normalLateralFaces(:,1) =  cross( bVec,wVec) / norm2(cross( bVec,wVec))
    normalLateralFaces(:,2) =  cross( vVec,bVec) / norm2(cross( vVec,bVec))
    normalLateralFaces(:,3) = -cross(b1Vec,wVec) / norm2(cross(b1Vec,wVec))
    normalLateralFaces(:,4) =  cross(b2Vec,vVec) / norm2(cross(b2Vec,vVec))

! ######################################################################
! #      TO BE CHECKED       ###########################################
! ######################################################################
    ! Find the elements belonging to each section ---------------
    ! allocate and initialise
    allocate(box_secloads(0:n_sect+1))
    do is = 0 , n_sect+1 ! 0,n_sect+1 to add dummy ghost extrem sections
      box_secloads(is)%nelems = 0
      allocate(box_secloads(is)%elems(size(comps(id_comp)%el)))
      box_secloads(is)%elems = -333
      allocate(box_secloads(is)%fracs(size(comps(id_comp)%el)))
      box_secloads(is)%fracs = -333.3_wp
      allocate(box_secloads(is)%cen(3,size(comps(id_comp)%el)))
      box_secloads(is)%cen   = -333.3_wp
    end do


    do ie = 1 , size(comps(id_comp)%el) ! loop over elems

      ! Reset some useful variables -----------------------------
      nCoordVert = 0.0_wp
      secVert = 0
      interSectPoints = 0.0_wp
      interSectAreas  = 0.0_wp
      nInterSect = 0
      nNodeInt = 0
      nver = comps(id_comp)%el(ie)%n_ver

      ! Compute cen and area of the elements --------------------
      ! since no prepare or load geom has been called yet.
      ! TODO: move this part in the time loop, first timestep
      ! compute the centre and the area of the element
      comps(id_comp)%el(ie)%cen =  &
        sum(comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver - &
                                      comps(id_comp)%i_points(1)+1),2) / &
                                      real(comps(id_comp)%el(ie)%n_ver,wp)
      comps(id_comp)%el(ie)%area =  0.5_wp * norm2( &
        cross ( comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(3)) &
              - comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(1)) &
            ,   comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(2)) &
              - comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(nver))&
              ) )

      ! Compute the distance of the centre of the elem from -----
      ! lateral faces of the box
      do i1 = 1 , 4
        distance(i1) = sum ( &
                  (comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,i1) ) * &
                                                    normalLateralFaces(:,i1) )
      end do

      ! if the centre of the element belongs to the box ->
      !  compute the slice contributions
      if ( all( distance .gt. 0.0_wp ) ) then

        ! nCoord of the cen of the elem (useless??)
        nCoord = sum ( &
                (comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,1)) * nVec )

        ! nCoord of the vertices of the elem ( (r-refPoi(1))\cdot nVec )
        do iv = 1 , nver
          nCoordVert(iv) = &
              sum ( ( comps(id_comp)%loc_points(:, &
                                      comps(id_comp)%el(ie)%i_ver(iv) ) &
                     - refPoiLateralFaces(:,1) ) * nVec )
        end do

        ! elems with at least one point s.t. nCoor \in (nCoorMix,nCoorMax)
        if ( ( maxval(nCoordVert) .ge. 0.0_wp ) .and. &
            ( minval(nCoordVert) .le. nCoordMax-nCoordMin ) ) then

            ! find the sections where the ver of the elem belong -------------
            ! !!! if the one vert is outside the box, secVert keep 0 value !!!
            do iv = 1 , nver
              if ( nCoordVert(iv) - nCoordSec(n_sect+1) .gt. 0.0_wp  ) then
                  secVert(iv) = n_sect + 1
              elseif ( nCoordVert(iv) - nCoordSec(1) .le. 0.0_wp  ) then
                  secVert(iv) = 0
              else
                do is = 1 , n_sect !TODO: improve this loop; exit when found
                  if ( ( nCoordVert(iv) - nCoordSec(is)   .gt. 0.0_wp ) .and. &
                        ( nCoordVert(iv) - nCoordSec(is+1) .le. 0.0_wp ) ) then
                    secVert(iv) = is
                  end if
                end do
              end if
            end do

            ! If the elements belong to: -------------------------------------
            ! (a) .gt. 2 sections ----> error.
            ! (b) .eq. 1 section  ----> easy
            ! (c) .eq. 2 section  ----> find intersections and partial contributions
            !                          to the sections

            ! (a) check if the elements belong to less than three sections --> otherwise ERROR
            if ( maxval(secVert(1:nver) ) &
                - minval(secVert(1:nver) ) .gt. 1 ) then
                write(msg,'(A)') 'The analyzed component '&
                &//comps(id_comp)%comp_name//' has an element which spans more &
                &than two sections, this cannot be handled. This generally means &
                &that an element is significantly wider than the section spacing,&
                & try reducing the number of sections'
              call error(this_mod_name,this_sub_name,msg)
            end if
            ! (b)
            if ( all(  secVert .eq. secVert(1) ) ) then
              box_secloads(secVert(1))%nelems = box_secloads(secVert(1))%nelems + 1
              box_secloads(secVert(1))%elems(box_secloads(secVert(1))%nelems) = ie
              box_secloads(secVert(1))%fracs(box_secloads(secVert(1))%nelems) = 1.0_wp
              box_secloads(secVert(1))%cen(:,box_secloads(secVert(1))%nelems) = &
                                                              comps(id_comp)%el(ie)%cen
            else !(c)

              !TODO: treat elements partially belonging to the box.  if ( ... )

              ! Find intersections with the plane delimiting the sections ----
              do iv = 1 , nver

                index2 = mod(iv,comps(id_comp)%el(ie)%n_ver)+1 ! following node

                ! if two consecutive nodes belong to different sections
                !  ---> find intersection
                if ( secVert( index2 ) .ne. secVert(iv) ) then
                  nInterSect = nInterSect + 1
                  node1 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(iv) )
                  node2 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(index2) )
                  interSectPoints(:,nInterSect) = node1 + (node2-node1) * &
                        (nCoordSec( max( secVert(iv), secVert(index2) ) ) - nCoordVert(iv) ) / &
                        (nCoordVert(index2) - nCoordVert(iv) )  ! ^---- max(.,.): id. of the section

                  if ( secVert(iv) .eq. secVert(1) ) then ! sec.1
                    nNodeInt(1,nInterSect) = iv     ; nNodeInt(2,nInterSect) = index2
                  else ! sec.2
                    nNodeInt(1,nInterSect) = index2 ; nNodeInt(2,nInterSect) = iv
                  end if
                end if
              end do

              ! Count the nodes belonging to the splitted elem ---------------
              iSec1 = secVert(1) ; sec1_nVer = 1 ; sec2_nVer = 0
              do iv = 2 , comps(id_comp)%el(ie)%n_ver
                if ( secVert(iv) .ne. iSec1 ) then
                  iSec2 = secVert(iv) ; sec2_nVer = sec2_nVer + 1
                else
                  sec1_nVer = sec1_nVer + 1
                end if
              end do

              ! Compute the area of the subelem with the minimum n.of points -
              !  (since it could be only TRI or QUAD) and compute the area of
              !  the othe elem as a difference.
              if ( sec1_nVer .le. sec2_nVer ) then
                if ( sec1_nVer .gt. 2 ) then
                    call error(this_sub_name, this_mod_name, 'Generic sectional &
                    & loads error')
                end if

                interSectAreas(1) = 0.5_wp * norm2( &
                    cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver( nNodeInt(1,1) ) ) , &
                            interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver( nNodeInt(1,2) ) ) ) ) ! nNodes(2,iSec1) ) )
                interSectAreas(2) = comps(id_comp)%el(ie)%area - interSectAreas(1)

                if ( sec1_nVer .eq. 1 ) then
                  interSectCen(:,1) = ( interSectPoints(:,1) + interSectPoints(:,2) + &
                        comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,1)) ) ) / 3.0_wp
                elseif ( sec1_nVer .eq. 2 ) then
                  interSectCen(:,1) = ( interSectPoints(:,1) + interSectPoints(:,2) + &
                        comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,1)) ) + &
                        comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,2)) ) ) / 4.0_wp
                else
                end if
                interSectCen(:,2) = ( comps(id_comp)%el(ie)%area * comps(id_comp)%el(ie)%cen - &
                       interSectAreas(1) * interSectCen(:,1) ) / interSectAreas(2)

              else
                if ( sec2_nVer .gt. 2 ) then
                    call error(this_sub_name, this_mod_name, 'Generic sectional &
                    & loads error')
                end if

              interSectAreas(2) = 0.5_wp * norm2( &
                      cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver( nNodeInt(2,1) ) ) , &
                          interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver( nNodeInt(2,2) ) ) ) ) ! nNodes(2,iSec1) ) )
              interSectAreas(1) = comps(id_comp)%el(ie)%area - interSectAreas(2)

              if ( sec2_nVer .eq. 1 ) then
                interSectCen(:,2) = ( interSectPoints(:,1) + interSectPoints(:,2) + &
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,1)) ) ) / 3.0_wp
              elseif ( sec2_nVer .eq. 2 ) then
                interSectCen(:,2) = ( interSectPoints(:,1) + interSectPoints(:,2) + &
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,1)) ) + &
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,2)) ) ) / 4.0_wp
              else
              end if
              interSectCen(:,1) = ( comps(id_comp)%el(ie)%area * comps(id_comp)%el(ie)%cen - &
                      interSectAreas(2) * interSectCen(:,2) ) / interSectAreas(1)

            end if

            ! Update structures
            box_secloads(isec1)%nelems = box_secloads(isec1)%nelems + 1
            box_secloads(isec1)%elems(box_secloads(isec1)%nelems) = ie
            box_secloads(isec1)%fracs(box_secloads(isec1)%nelems) = interSectAreas(1) / comps(id_comp)%el(ie)%area
            box_secloads(isec1)%cen(:,box_secloads(isec1)%nelems) = interSectCen(:,1)
            box_secloads(isec2)%nelems = box_secloads(isec2)%nelems + 1
            box_secloads(isec2)%elems(box_secloads(isec2)%nelems) = ie
            box_secloads(isec2)%fracs(box_secloads(isec2)%nelems) = interSectAreas(2) / comps(id_comp)%el(ie)%area
            box_secloads(isec2)%cen(:,box_secloads(isec2)%nelems) = interSectCen(:,2)


          end if

        end if ! elems with at least one point s.t. nCoor in (nCoorMix,nCoorMax)

      end if ! centre of the elems between the lateral faces of the box

    end do ! loop over elems

    ! Allocations and time loop ++++++++++++++++++++++++++++++++++++++++

    ! Find the id of the reference where the loads must be projected ---
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
    call open_hdf5_file(trim(filename),floc)
    call load_refs(floc,refs_R,refs_off) ! ,refs_G,refs_f,refs_tag)
    call close_hdf5_file(floc)

    ref_id = comps(id_comp)%ref_id
    write(msg,'(A,I0,A)') '   Employing the local reference frame for moment &
      &and force projection.'//nl//'   Reference tag: '&
      //trim(comps(id_comp)%ref_tag)//'  (ID: ',ref_id,')'
    call printout(trim(msg))

    ! allocate tmp array to store the results --------------
    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( sec_loads(n_time,n_sect,n_loads) ) ; sec_loads = -333.0_wp ! initialisation for DEBUG
    allocate( time(n_time) ) ; time = -333.0_wp
    allocate( ref_mat(n_time,9) , off_mat(n_time,3) )
    ires = 0
    
    do it = an_start, an_end, an_step

      ires = ires + 1

      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)

      ! Load the references and move the points ---
      call load_refs(floc,refs_R,refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, points_virtual, refs_R, refs_off)
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)

      call close_hdf5_file(floc)

      do is = 1 , n_sect

        ! from local coordinates to coordinates in the base ref.frame
        r_axis_bas(:,is) = matmul( refs_R(:,:,ref_id) , r_axis(:,is) ) + &
                          refs_off(:,ref_id)

        ! force
        F_bas = 0.0_wp ; M_bas = 0.0_wp
        do ie = 1 , box_secloads(is)%nelems

          !from local coordinates to coordinates in the base ref.frame
          box_secloads_cen = matmul( refs_R(:,:,ref_id) , &
                                    box_secloads(is)%cen(:,ie) )

          F_bas1 = comps(id_comp)%el( box_secloads(is)%elems(ie) )%dforce * &
                                      box_secloads(is)%fracs(ie)
          F_bas  = F_bas + F_bas1
          !TODO: add moment (it requires axis computation)
          M_bas  = M_bas + cross( box_secloads_cen -          &
                                  r_axis_bas(:,is) , F_bas1 ) &
                        + comps(id_comp)%el(                 &
                            box_secloads(is)%elems(ie) )%dmom ! updated 2018-07-12
        end do

        ! From global to local coordinates of forces and moments
        sec_loads(ires,is,1:3) = matmul( &
            transpose( refs_R(:,:, ref_id) ) , F_bas ) / y_span(is)
        ! TODO: moment
        ! moment ( only the component around the <ax_coord> )
        M_bas = matmul( &
            transpose( refs_R(:,:, ref_id) ) , M_bas )
        sec_loads(ires,is,4) = M_bas(2) / y_span(is)
      end do

      ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
      off_mat(ires,:) = refs_off(:,ref_id)
      time(ires) = t

    end do

    if(average) then
      allocate(sec_loads_ave(1,size(sec_loads,2),size(sec_loads,3)))
      sec_loads_ave(1,:,:) = sum(sec_loads, 1)/real(size(sec_loads,1),wp)
      select case(trim(out_frmt))
        case('dat')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)
          call dat_out_sectional ( filename, components_names(1), y_cen, &
                                    y_span, chord, time(1:1), &
                            sec_loads_ave, ref_mat , off_mat, average )
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_ave.plt'
          call tec_out_sectional (filename, time(1:1), sec_loads_ave, y_cen, &
                                                                      y_span , chord)
      end select
      deallocate(sec_loads_ave)
    else
      select case(trim(out_frmt))
        case('dat')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)
          call dat_out_sectional ( filename, components_names(1), y_cen, &
                                y_span, chord, time, sec_loads, &
                                ref_mat, off_mat, average)
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
          call tec_out_sectional ( filename, time, sec_loads, y_cen, y_span, chord)
      end select
    endif

    ! destroy box_secloads structure
    do is = lbound(box_secloads,1) , ubound(box_secloads,1)
      deallocate(box_secloads(is)%elems)
      deallocate(box_secloads(is)%fracs)
    end do

    deallocate(r_axis,r_axis_bas)
    deallocate(nCoordSec,nCoordCen)
    deallocate(box_coord_tmp,box_coord)
    deallocate(y_cen,y_span)

  case default
    call error(this_sub_name, this_mod_name, 'No valid kind of generation &
    &method for the geometry was found in geometry file. This is highly &
    &unlikely to happen, there might be troubles with the geometry file &
    &employed')

  end select

  call destroy_elements(comps)
  deallocate(comps,components_names)

  write(msg,'(A,I0,A)') nl//'++++++++++ Sectional loads done'//nl
  call printout(trim(msg))

end subroutine post_sectional

! ----------------------------------------------------------------------




! ----------------------------------------------------------------------

end module mod_post_sectional
