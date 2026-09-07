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
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================

!> Module containing the subroutines to perform chordwise loads
!! analysis during postprocessing
module mod_post_chordwise

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
    dat_out_chordwise

  use mod_tecplot_out, only: &
    tec_out_chordwise

  use mod_math, only: &
    cross, linear_interp, vec2mat

  implicit none

  public :: post_chordwise

  private

  character(len=*), parameter :: this_mod_name = 'mod_post_chordwise'
  character(len=max_char_len) :: msg

  contains


!TODO (to be tested):
! - parametric components aligned with y-axis
! - tol_y_cen is hard-coded. write tol_y_cen as an input
! - preCICE coupled components 
subroutine post_chordwise(sbprms, basename, data_basename, an_name, &
                          ia, out_frmt, components_names, all_comp, &
                          an_start, an_end, an_step, average )
  type(t_parse), pointer                                  :: sbprms
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
  integer(h5loc)                                          :: floc, gloc, cloc, ploc 
  real(wp), allocatable                                   :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                   :: vort(:), cp(:)
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

  !> chordwise load: paramteric components
  real(wp)                                                :: axis_dir(3) , axis_nod(3)
  integer                                                 :: is                                 
  integer, parameter                                      :: n_loads = 3   ! 
  real(wp), allocatable                                   :: y_cen(:) , y_span(:)
  real(wp), parameter                                     :: tol_y_cen = 1.0e-3_wp
  
  integer                                                 :: ie , ic , it , i1, ista, idir, ii
  integer                                                 :: ipan_minus, ipan_plus 
  integer                                                 :: n_station 
  real(wp), allocatable                                   :: span_station(:)
  real(wp), allocatable                                   :: y_cen_tras(:) 
  real(wp), allocatable                                   :: chord(:) 
  integer                                                 :: chord_start, chord_end
  integer,  allocatable                                   :: id_minus(:), id_plus(:) 
  real(wp)                                                :: u_inf(3), P_inf, rho 
  real(wp)                                                :: R_cen_minus(3,3), R_cen_plus(3,3)
  !> field to interpolate  
  real(wp)                                                :: force_minus(3) = 0.0_wp
  real(wp)                                                :: force_plus(3)  = 0.0_wp 
  real(wp)                                                :: tang_minus(3)  = 0.0_wp   
  real(wp)                                                :: tang_plus(3)   = 0.0_wp
  real(wp)                                                :: nor_minus(3)   = 0.0_wp   
  real(wp)                                                :: nor_plus(3)    = 0.0_wp
  real(wp)                                                :: cen_minus(3)   = 0.0_wp   
  real(wp)                                                :: cen_plus(3)    = 0.0_wp
  real(wp)                                                :: pres_minus     = 0.0_wp     
  real(wp)                                                :: pres_plus      = 0.0_wp
  real(wp)                                                :: cp_minus       = 0.0_wp
  real(wp)                                                :: cp_plus        = 0.0_wp
  real(wp)                                                :: chord_minus    = 0.0_wp
  real(wp)                                                :: chord_plus     = 0.0_wp
  !> interpolated field
  real(wp), allocatable                                   :: force_int(:,:,:,:)   
  real(wp), allocatable                                   :: tang_int(:,:,:,:)    
  real(wp), allocatable                                   :: nor_int(:,:,:,:)    
  real(wp), allocatable                                   :: cen_int(:,:,:,:) 
  real(wp), allocatable                                   :: pres_int(:,:,:)      
  real(wp), allocatable                                   :: cp_int(:,:,:)      
  real(wp), allocatable                                   :: chord_int(:)
  !> average field 
  real(wp), allocatable                                   :: force_ave(:,:,:,:)   
  real(wp), allocatable                                   :: tang_ave(:,:,:,:)    
  real(wp), allocatable                                   :: nor_ave(:,:,:,:)    
  real(wp), allocatable                                   :: cen_ave(:,:,:,:) 
  real(wp), allocatable                                   :: pres_ave(:,:,:)      
  real(wp), allocatable                                   :: cp_ave(:,:,:)      
  logical :: bracket_found   ! Flag to indicate if chord bracketing is found

  character(len=max_char_len)                             :: filename
  character(len=max_char_len)                             :: comp_input

  character(len=*), parameter :: this_sub_name = 'post_chordwise' 

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' chordwise loads'//nl
  call printout(trim(msg))


  ! Some warnings and errors -------------------------------
  ! WARNING: sectional loads are computed for one components only
  !  at time. If no component is specified --> return
  n_comp = countoption(sbprms,'component')
  if ( n_comp .le. 0 ) then
    call warning(this_mod_name, this_sub_name, 'No component specified for &
                              &chordwise_loads analysis. Skipped analysis.')
    return
  else if ( n_comp .ge. 2 ) then
    call warning(this_mod_name, this_sub_name, &
        'More than one component specified &
        &for ''chordwise_loads'' analysis: just the first one is considered. &
        &Please run another ''chordwise_analysis'' if you need it on more than &
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

  n_station = getint(sbprms,'n_station')
  span_station = getrealarray(sbprms,'span_station', n_station)
  comp_input = trim(comps(1)%comp_input ) 
  select case( trim(comp_input) )

  case ( 'parametric', 'pointwise')
    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
    ax_coor = 2   ! 2. the parametric elements is defined
                  !    along the y-axis 

    nelem_span = comps(id_comp)%parametric_nelems_span
    nelem_chor = comps(id_comp)%parametric_nelems_chor
    n_sect = nelem_span

    do i1 = 1, size(comps(id_comp)%loc_points, 2) 
      comps(id_comp)%loc_points(:, i1) = &
      matmul(comps(id_comp)%coupling_node_rot,comps(id_comp)%loc_points(:,i1))
    end do

    do i1 = 1, size(comps(id_comp)%el) 
      comps(id_comp)%el(i1)%cen = &
      matmul(comps(id_comp)%coupling_node_rot,comps(id_comp)%el(i1)%cen)
    end do

    
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

    allocate(id_minus(n_station))
    allocate(id_plus(n_station))
    allocate(y_cen_tras(n_sect)) 
    allocate(chord_int(n_station)) 

    !> Translate centers on reference station
    do ista = 1, n_station
      y_cen_tras = y_cen - span_station(ista)
      !> Get index of the stripes across the span stripe
      bracket_found = .false.
      do is = 1, n_sect - 1
        if (y_cen_tras(is) <= (1e-6_wp) .and. y_cen_tras(is + 1) >= (1e-6_wp)) then
          id_minus(ista) = is
          id_plus(ista) = is + 1
          bracket_found = .true.
          exit ! exit loop once bracket is found
        elseif (y_cen_tras(is) <= (1e-6_wp) .and. y_cen_tras(is + 1) >= (1e-6_wp)) then
          id_minus(ista) = is - 1
          id_plus(ista) = is
          bracket_found = .true.
          exit ! exit loop once bracket is found
        endif
      end do

      if (.not. bracket_found) then
        ! Chord bracketing not found for current station
        call error(this_sub_name,this_mod_name,'Error: Chord bracketing not found for station')
      elseif (id_minus(ista) < 1 .or. id_plus(ista) > n_sect) then
        ! Chord segments out of range
        call error(this_sub_name,this_mod_name,'Error: Chord segments out of range for station')
      endif

      chord_minus = chord(id_minus(ista))
      chord_plus = chord(id_plus(ista))

      !> Interpolate chord at current station
      call linear_interp((/chord_minus, chord_plus/), &
                      (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),&
                      span_station(ista), chord_int(ista))
    end do


    ! Find the coordinate of the reference points on the axis --------------
    !  ( with coord. y_cen )
    if ( abs(axis_dir(2)) .lt. 1e-6_wp ) then
      call error('dust_post','','Wrong definition of the axis in&
            & sectional_loads analysis: abs(axis_dir(2)) .lt. 1e-6.&
            & STOP')
    end if

    ! Find the ref_id or the reference frame where loads are projected -----
    ref_id = comps(id_comp)%ref_id
    write(msg,'(A,I0,A)') '   Employing the local reference frame for moment &
      &and force projection.'//nl//'   Reference tag: '&
      //trim(comps(id_comp)%ref_tag)//'  (ID: ',ref_id,')'
    call printout(trim(msg))

    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( time(n_time) ) ; time = -333.0_wp
    ires = 0

    allocate(force_int(n_time,n_station,nelem_chor,3)); force_int = 0.0_wp 
    allocate(tang_int(n_time,n_station,nelem_chor,3));  tang_int = 0.0_wp 
    allocate(nor_int(n_time,n_station,nelem_chor,3));   nor_int = 0.0_wp 
    allocate(cen_int(n_time,n_station,nelem_chor,3));   cen_int = 0.0_wp 
    allocate(pres_int(n_time,n_station,nelem_chor));    pres_int = 0.0_wp 
    allocate(cp_int(n_time,n_station,nelem_chor));      cp_int = 0.0_wp 
    
    do it = an_start, an_end, an_step
      ires = ires + 1

      !> Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)
      !> Load free-stream parameters
      call open_hdf5_group(floc,'Parameters',ploc)
      call read_hdf5(u_inf,'u_inf',ploc)
      call read_hdf5(P_inf,'P_inf',ploc)
      call read_hdf5(rho,'rho_inf',ploc)
      call close_hdf5_group(ploc) 

      ! Load the references and move the points ---
      call load_refs(floc, refs_R, refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                  filen = trim(filename) )
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)
      
      call close_hdf5_file(floc)

      
      ! compute chordwise loads. Loop over the panels of each section
      do ista = 1, n_station 
        ipan_minus = 0
        ipan_plus = 0
        ipan_minus = (id_minus(ista)- 1)*nelem_chor
        ipan_plus = (id_plus(ista)- 1)*nelem_chor 

        do ic = 1, nelem_chor ! loop over chord
          ipan_minus = ipan_minus + 1  
          !> 4d matrix: time, station, chord, [x, y, z] 
          force_minus = comps(id_comp)%el(ipan_minus)%dforce
          tang_minus = comps(id_comp)%el(ipan_minus)%tang(:,2)/&
                        norm2(comps(id_comp)%el(ipan_minus)%tang(:,2))
          nor_minus = comps(id_comp)%el(ipan_minus)%nor/&
                        norm2(comps(id_comp)%el(ipan_minus)%nor)
          cen_minus = comps(id_comp)%el(ipan_minus)%cen 
          !> 3d matrix: time, station, chord(pres) 
          pres_minus = comps(id_comp)%el(ipan_minus)%pres 
          cp_minus = (pres_minus - P_inf)/&
                                    (0.5_wp*rho*norm2(u_inf)**2) 

          ipan_plus = ipan_plus + 1 
          !> 4d matrix: time, station, chord, [x, y, z] 
          force_plus = comps(id_comp)%el(ipan_plus)%dforce
          tang_plus = comps(id_comp)%el(ipan_plus)%tang(:,2)/&
                      norm2(comps(id_comp)%el(ipan_plus)%tang(:,2))
          nor_plus = comps(id_comp)%el(ipan_plus)%nor/&
                      norm2(comps(id_comp)%el(ipan_plus)%nor)
          cen_plus = comps(id_comp)%el(ipan_plus)%cen 
          !> 3d matrix: time, station, chord(pres) 
          pres_plus = comps(id_comp)%el(ipan_plus)%pres 
          cp_plus  = (pres_plus - P_inf)/& 
                                    (0.5_wp*rho*norm2(u_inf)**2) 
#if USE_PRECICE
          call vec2mat(comps(id_comp)%el(ipan_plus)%ori, R_cen_minus)
          R_cen_minus = matmul(transpose(comps(id_comp)%coupling_node_rot), R_cen_minus) 
          call vec2mat(comps(id_comp)%el(ipan_plus)%ori, R_cen_plus)
          R_cen_plus = matmul(transpose(comps(id_comp)%coupling_node_rot), R_cen_plus) 
#else     
          R_cen_minus = refs_R(:,:, ref_id)
          R_cen_plus = refs_R(:,:, ref_id)
#endif
          ! From global to local coordinates of forces and moments
          force_minus = matmul(transpose(R_cen_minus), force_minus) & 
                                          / y_span(id_minus(ista)) 
          tang_minus = matmul(transpose(R_cen_minus), tang_minus)                                                   
          nor_minus = matmul(transpose(R_cen_minus), nor_minus)                                                   
          cen_minus = matmul(transpose(R_cen_minus), cen_minus) 

          force_plus = matmul(transpose(R_cen_plus), force_plus) & 
                                          / y_span(id_plus(ista)) 
          tang_plus = matmul(transpose( R_cen_plus), tang_plus)                                                   
          nor_plus = matmul(transpose( R_cen_plus), nor_plus)                                                   
          cen_plus = matmul(transpose(R_cen_plus), cen_plus) 
          !> interpolate  
          do idir = 1,3
            call linear_interp((/force_minus(idir), force_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), force_int(ires,ista,ic,idir)) 
            call linear_interp((/tang_minus(idir), tang_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), tang_int(ires,ista,ic,idir)) 
            call linear_interp((/nor_minus(idir), nor_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), nor_int(ires,ista,ic,idir)) 
            call linear_interp((/cen_minus(idir), cen_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), cen_int(ires,ista,ic,idir)) 
          enddo 
          
          call linear_interp((/pres_minus, pres_plus/), & 
                            (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                            span_station(ista), pres_int(ires,ista,ic)) 
          call linear_interp((/cp_minus, cp_plus/), & 
                            (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                            span_station(ista), cp_int(ires,ista,ic)) 
          
          
        end do ! loop over chord
        
      end do ! loop over sections

      time(ires) = t

    end do ! loop over time steps

    if(average) then 
      allocate(force_ave(1,n_station,nelem_chor,3)); force_ave = 0.0_wp 
      allocate(tang_ave(1,n_station,nelem_chor,3));  tang_ave = 0.0_wp 
      allocate(nor_ave(1,n_station,nelem_chor,3));   nor_ave = 0.0_wp 
      allocate(cen_ave(1,n_station,nelem_chor,3));   cen_ave = 0.0_wp 
      allocate(pres_ave(1,n_station,nelem_chor));    pres_ave = 0.0_wp 
      allocate(cp_ave(1,n_station,nelem_chor));      cp_ave = 0.0_wp 

      force_ave(1,:,:,:) = sum(force_int,1)/real(size(force_int,1),wp)
      tang_ave(1,:,:,:)  = sum(tang_int,1)/real(size(tang_int,1),wp)
      nor_ave(1,:,:,:)   = sum(nor_int,1)/real(size(nor_int,1),wp)
      cen_ave(1,:,:,:)   = sum(cen_int,1)/real(size(cen_int,1),wp)
      pres_ave(1,:,:)    = sum(pres_int,1)/real(size(pres_int,1),wp)
      cp_ave(1,:,:)      = sum(cp_int,1)/real(size(cp_int,1),wp)
      select case(trim(out_frmt))
        case('dat')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)
          call dat_out_chordwise (filename, components_names(1), time(1:1), &
                        force_ave, tang_ave, nor_ave, cen_ave, pres_ave, & 
                        cp_ave, average, n_station, span_station, chord_int)
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_ave.plt' 
          call tec_out_chordwise(filename, time(1:1), &
                              force_ave, tang_ave, nor_ave, cen_ave, pres_ave, & 
                              cp_ave, n_station, span_station, chord_int)
      end select 
    else 
      select case(trim(out_frmt))
        case('dat')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)
          call dat_out_chordwise (filename, components_names(1), time, &
                        force_int, tang_int, nor_int, cen_int, pres_int, & 
                        cp_int, average, n_station, span_station, chord_int)
        case('tecplot')
          write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt' 
          call tec_out_chordwise(filename, time, &
                              force_int, tang_int, nor_int, cen_int, pres_int, & 
                              cp_int, n_station, span_station, chord_int)
      end select 
    endif 

  case('cgns') 
    call error(this_sub_name, this_mod_name, 'case cgns not implemented so far')
  case default
    call error(this_sub_name, this_mod_name, 'type not known')   
  end select 
  
    
  
end subroutine 
end module 
