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


!> Module to treat the most simple input-output from ascii formatted data
!! files
module mod_parametric_io

use mod_param, only: &
  wp, max_char_len, nl, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption

use mod_spline, only: & 
  hermite_spline_profile, catmull_rom_chain

use mod_math, only: & 
  linear_interp, unique, unique_row, linspace, sort_vector_real, & 
  geoseries, geoseries_both, geoseries_hinge

use mod_hinges, only: &
  t_hinge, t_hinge_input 

use, intrinsic :: ieee_arithmetic
!----------------------------------------------------------------------

implicit none

public :: read_mesh_parametric, read_actuatordisk_parametric , &
          define_section , define_division, geoseries, geoseries_both
private

character(len=*), parameter :: this_mod_name = 'mod_parametric_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

subroutine read_mesh_parametric(mesh_file, ee, rr, y_fountain, nelem_chord, ref_chord_fraction, ElType, &
                    npoints_chord_tot, nelem_span_tot, hinges, n_hinges, mesh_mirror, mesh_symmetry, &
                    nelem_span_list, airfoil_list_actual, i_airfoil_e, normalised_coord_e, & 
                    aero_table_out, thickness)

  character(len=*), intent(in)              :: mesh_file
  integer  , allocatable, intent(out)       :: ee(:,:)
  real(wp) , allocatable, intent(out)       :: rr(:,:)
  integer  , intent(in)                     :: nelem_chord
  integer  , intent(out), optional          :: npoints_chord_tot, nelem_span_tot 
  integer  , intent(in)                     :: n_hinges
  logical  , intent(in)                     :: mesh_mirror, mesh_symmetry

  real(wp) , intent(out)                    :: y_fountain      
  type(t_hinge_input), allocatable, intent(inout) :: hinges(:)
  type(t_parse)                             :: pmesh_prs
  integer                                   :: ee_size , rr_size
  logical                                   :: twist_linear_interp
  logical, intent(out), optional            :: aero_table_out
  logical                                   :: aero_table
  real(wp), allocatable , intent(out), optional ::  thickness(:,:)
  real(wp), allocatable                     :: thickness_section1(:), thickness_section2(:) 
  real(wp)                                  :: thickness_section
  integer                                   :: nelem_chord_tot
  integer                                   :: npoint_chord_tot, npoint_span_tot
  integer                                   :: nRegions, nSections
  integer                                   :: iRegion, iSection 
  integer                                   :: iChord, iSpan 
  integer                                   :: iElement, iPoint
  real(wp), intent(in)                      :: ref_chord_fraction
  real(wp), allocatable                     :: ref_point(:)
  !> Data read from file
  !> Sections 
  real(wp), allocatable                     :: chord_list(:) 
  real(wp), allocatable                     :: twist_list(:)
  integer                                   :: iAirfoil
  character(len=max_char_len), allocatable  :: airfoil_list(:)
  character(len=max_char_len), allocatable  :: airfoil_table_list(:)
  character(len=max_char_len), allocatable , intent(out), optional :: airfoil_list_actual(:)
  integer,  allocatable, intent(out), optional :: i_airfoil_e(:,:)
  real(wp), allocatable, intent(out), optional :: normalised_coord_e(:,:)
  !> geometric series parameters
  real(wp)                                  :: r, r_le, r_te, r_le_fix, r_te_fix, r_le_moving, r_te_moving
  real(wp)                                  :: x_refinement
  !> hinge mesh 
  real(wp), allocatable                     :: rrv_le(:,:), rrv_te(:,:), ac_line(:,:) 
  integer                                   :: ih, ia
  real(wp), allocatable                     :: csi_hinge_not_unique(:), csi_hinge(:)
  real(wp), allocatable                     :: delta_x(:), delta_x_no_off(:), csi_adim(:) 
  integer, allocatable                      :: point_region(:)
  real(wp)                                  :: merge_tol = 0.0_wp
  logical                                   :: adaptive_mesh 
  !> Regions  
  integer , allocatable, intent(out)        :: nelem_span_list(:)
  real(wp)                                  :: rr_sta
  real(wp), allocatable                     :: rr_adim(:)
  real(wp), allocatable                     :: span_list(:) 
  real(wp), allocatable                     :: sweep_list(:) 
  real(wp), allocatable                     :: dihed_list(:)
  real(wp), allocatable                     :: division(:), divisionIB(:), divisionOB(:) 
  !> geometric series parameters (region)
  real(wp), allocatable                     :: r_in_list(:), r_ob_list(:), r_span_list(:)
  real(wp), allocatable                     :: y_ref_list(:)

  logical                                   :: adjust_xz = .false.
  character(len=max_char_len), allocatable  :: type_span_list(:)
  integer                                   :: n_type_span
  character, intent(in)                     :: ElType
  !> Sections 1. 2.
  real(wp), allocatable :: xySection1(:,:), xySection2(:,:), xyAirfoil1(:,:), xyAirfoil2(:,:)
  real(wp), allocatable :: rrSection1(:,:), rrSection2(:,:)
  real(wp)                                  :: dx_ref , dy_ref , dz_ref
  integer                                   :: ista, iend, i, j

  !> Linear interpolation of the twist angle
  real(wp), allocatable                     :: rr_tw(:,:) , rr_tw_1(:,:) , rr_tw_2(:,:)
  real(wp)                                  :: dx_ref_1, dy_ref_1, dz_ref_1
  real(wp)                                  :: interp_weight
  real(wp), allocatable                     :: chord_fraction(:)
  real(wp), allocatable                     :: span_fraction(:)
  character(len=max_char_len)               :: type_chord
  integer                                   :: i1
  integer                                   :: iSec
  real(wp)                                  :: dy_actual_airfoils
  real(wp)                                  :: dy_sections, csi
  real(wp)                                  :: twist_rad
  integer                                   :: i_aero1, i_aero2, nAirfoils
  character(len=*), parameter               :: this_sub_name = 'read_mesh_parametric'
  aero_table = .false.

  !Prepare all the parameters to be read in the file
  ! Global parameters
  
  call pmesh_prs%CreateLogicalOption('airfoil_table_correction', &
                'include presence of aerodynamic .c81 for corrections', &
                'F', &
                multiple=.false.);
  call pmesh_prs%CreateStringOption('type_chord',&
                'type of chord-wise division: uniform, cosine, cosineLE, cosineTE',&
                'uniform', &
                multiple=.false.);
  call pmesh_prs%CreateRealArrayOption('starting_point',&
                'Starting point (inboard TE), (x, y, z)', &
                '(/0.0, 0.0, 0.0/)',&
                multiple=.false.);
  call pmesh_prs%CreateLogicalOption('twist_linear_interpolation',&
                'Linear interpolation of the twist angle, for the whole component',&
                'F',&
                multiple=.false.); 

  
  ! TODO: do we really need to allow the user to define the parameter above?
  ! the reference chord fraction is a number in the range [0,1] that defines
  ! - where the reference point is located in the root chord
  ! - the line related to the sweep angle
  ! - the point around which the sections are rotated to impose twist


  !> Section parameters
  call pmesh_prs%CreateRealOption(     'chord', 'section chord', &
                multiple=.true.)
  call pmesh_prs%CreateRealOption(     'twist', 'section twist angle', &
                multiple=.true.)
  call pmesh_prs%CreateStringOption( 'airfoil', 'section airfoil', &
                multiple=.true.)                
  call pmesh_prs%CreateStringOption( 'airfoil_table', 'airfoil table path', &
                multiple=.true.)
  !> geoseries parameters 
  call pmesh_prs%CreateRealOption( 'r', 'growth ratio of the elements at edge', '0.125', &
                multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_le', 'growth ratio of the elements at leading edge', &
                '0.143', multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_te', 'growth ratio of the elements at trailing edge', &
                '0.066', multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_le_fix', 'growth ratio of the elements at leading edge &
                fixed part', '0.125', multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_te_fix', 'growth ratio of the elements at trailing edge &
                fixed part', '0.143', multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_le_moving', 'growth ratio of the elements at leading edge &
                moving part', '0.143', multiple=.false.)
  call pmesh_prs%CreateRealOption( 'r_te_moving', 'growth ratio of the elements at trailing edge & 
                moving part', '0.1', multiple=.false.) 
  call pmesh_prs%CreateRealOption( 'x_refinement', 'chordwise station to which the refinement start', &
                '0.5', multiple=.false.)

  !> Region parameters
  call pmesh_prs%CreateRealOption(    'span', 'region span',&
                multiple=.true.)
  call pmesh_prs%CreateRealOption(   'sweep', 'region sweep angle [degrees]',&
                multiple=.true.)
  call pmesh_prs%CreateRealOption(   'dihed', 'region dihedral angle [degrees]',&
                multiple=.true.)
  call pmesh_prs%CreateIntOption( 'nelem_span', 'number of span-wise elements in the region',&
                multiple=.true.)
  call pmesh_prs%CreateStringOption('type_span', 'type of span-wise division: &
                &uniform, cosine, cosineIB, cosineOB, equalarea, geoseries, geoseries_ob, geoseries_ib',&
                multiple=.true.)

  !> geoseries parameters (region)
  call pmesh_prs%CreateRealOption('r_ib', 'growth ratio of the elements inboard', &
                multiple=.true.)
  call pmesh_prs%CreateRealOption('r_ob', 'growth ratio of the elements at outboard', &
                multiple=.true.)
  call pmesh_prs%CreateRealOption('y_refinement', 'spanwise station to which the refinement start', &
                multiple=.true.) 
  
  !> y_fountain 
  call pmesh_prs%CreateRealOption('y_fountain', 'Upper y-coordinates beyond which the wake particles are not released', & 
                                  '0.0',  multiple=.false.)
  
  !> Read the parameters
  call pmesh_prs%read_options(mesh_file,printout_val=.true.)

  nSections  = countoption(pmesh_prs, 'chord')
  nRegions   = countoption(pmesh_prs, 'span')
  aero_table = getlogical(pmesh_prs, 'airfoil_table_correction')

  !> geometric series parameters 
  r = getreal(pmesh_prs, 'r')
  r_le = getreal(pmesh_prs, 'r_le')
  r_te = getreal(pmesh_prs, 'r_te')
  r_le_fix = getreal(pmesh_prs, 'r_le_fix')
  r_te_fix = getreal(pmesh_prs, 'r_te_fix')
  r_le_moving = getreal(pmesh_prs, 'r_le_moving')
  r_te_moving = getreal(pmesh_prs, 'r_te_moving')
  x_refinement = getreal(pmesh_prs, 'x_refinement')

  !> y_fountain 
  y_fountain = getreal(pmesh_prs, 'y_fountain') 
  
  ! Check that nSections = nRegion + 1
  if ( nSections .ne. nRegions + 1 ) then
    call error(this_sub_name, this_mod_name, 'Unconsistent input: &
          &nSections .ne. nRegions. Stop.')
  end if

  ref_point          = getrealarray(pmesh_prs,'starting_point',3)

  twist_linear_interp= getlogical(pmesh_prs,'twist_linear_interpolation')

  !Check the number of inputs
  if(countoption(pmesh_prs,'nelem_span') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &number of "nelem_span" different from number of regions.')
  if(countoption(pmesh_prs,'dihed') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &number of "dihed" different from number of regions.')
  if(countoption(pmesh_prs,'sweep') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &number of "sweep" different from number of regions.')

  if(countoption(pmesh_prs,'twist') .ne. nSections ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &number of "twist" different from number of sections.')
  if(countoption(pmesh_prs,'airfoil') .ne. nSections ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &number of "airfoil" different from number of sections.')

  ! Get total number of elements and initialize arrays
  allocate(nelem_span_list(nRegions));   nelem_span_list = 0
  allocate(      span_list(nRegions));         span_list = 0.0_wp
  allocate(     sweep_list(nRegions));        sweep_list = 0.0_wp
  allocate(     dihed_list(nRegions));        dihed_list = 0.0_wp
  allocate(    r_span_list(nRegions));       r_span_list = 0.0_wp
  allocate(     y_ref_list(nRegions));        y_ref_list = 0.0_wp
  allocate(      r_in_list(nRegions));         r_in_list = 0.0_wp
  allocate(      r_ob_list(nRegions));         r_ob_list = 0.0_wp 

  nelem_span_tot = 0
  do iRegion = 1,nRegions
    nelem_span_list(iRegion) = getint(pmesh_prs,'nelem_span')
    span_list(iRegion)       = getreal(pmesh_prs,    'span' )
    sweep_list(iRegion)      = getreal(pmesh_prs,    'sweep')
    dihed_list(iRegion)      = getreal(pmesh_prs,    'dihed')
    nelem_span_tot = nelem_span_tot + nelem_span_list(iRegion)
  enddo

  ! type_span
  allocate( type_span_list(nRegions));
  n_type_span = countoption(pmesh_prs,'type_span')
  if ( n_type_span .eq. 0 ) then ! default is 'uniform'
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = 'uniform'
    end do
  else if ( n_type_span .eq. nRegions ) then
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = getstr(pmesh_prs,'type_span')
      !> check geometry series type 
      if (trim(type_span_list(iRegion)) .eq. 'geoseries') then 
        y_ref_list(iRegion)      = getreal(pmesh_prs,   'y_refinement')
        r_in_list(iRegion)       = getreal(pmesh_prs,   'r_ib')
        r_ob_list(iRegion)       = getreal(pmesh_prs,   'r_ob') 
      elseif (trim(type_span_list(iRegion)) .eq. 'geoseries_ib') then 
        r_in_list(iRegion)       = getreal(pmesh_prs,   'r_ib')
      elseif (trim(type_span_list(iRegion)) .eq. 'geoseries_ob') then
        r_ob_list(iRegion)       = getreal(pmesh_prs,   'r_ob')
      endif
    end do
  else
    write(*,*) ' Mesh file   : ' , trim(mesh_file)
    write(*,*) ' Number of span element distribution type definitions &
                &(type_span): ' , n_type_span
    write(*,*) ' Number of regions : ' , nRegions
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
            & the number of type of span element distribution defined is &
            &different from the number of span regions')
  end if

  allocate(chord_list  (nSections))  ; chord_list = 0.0_wp
  allocate(twist_list  (nSections))  ; twist_list = 0.0_wp
  allocate(airfoil_list(nSections))
  allocate(airfoil_table_list(nSections))

  allocate(thickness_section1(nRegions)); thickness_section1 = 0.0_wp
  allocate(thickness_section2(nRegions)); thickness_section2 = 0.0_wp

  do iSection= 1, nSections
    chord_list(iSection)   = getreal(pmesh_prs,'chord')
    twist_list(iSection)   = getreal(pmesh_prs,'twist')
    airfoil_list(iSection) = getstr(pmesh_prs,'airfoil')

    if (ElType .eq.'v' .and. aero_table) then
      airfoil_table_list(iSection) = getstr(pmesh_prs,'airfoil_table')
    endif   

  enddo
  type_chord = getstr(pmesh_prs,'type_chord')

  if (ElType.eq.'p') then
    nelem_chord_tot = 2*nelem_chord
  else
    nelem_chord_tot = nelem_chord
  endif

  npoint_chord_tot = nelem_chord_tot + 1
  npoint_span_tot  = nelem_span_tot  + 1

  ee_size = nelem_chord_tot * nelem_span_tot
  rr_size = npoint_chord_tot * npoint_span_tot

  ! Directly build connectivity matrix
  allocate(ee(4,ee_size)) ; ee = 0
  do iChord = 1,nelem_chord_tot
    do iSpan = 1,nelem_span_tot
      iElement =  nelem_chord_tot*(iSpan-1) + iChord
      iPoint   = npoint_chord_tot*(iSpan-1) + iChord
      ee(1,iElement) = iPoint + npoint_chord_tot         ! iPoint
      ee(2,iElement) = iPoint                            ! iPoint + 1
      ee(3,iElement) = iPoint + 1                        ! iPoint + npoint_chord_tot + 1
      ee(4,iElement) = iPoint + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
    enddo
  enddo

  allocate(rr(3,rr_size)) ; rr = 0.0_wp
  ! get chordwise division
  allocate(chord_fraction(nelem_chord+1))

  adaptive_mesh = .false. 
  do ih = 1, n_hinges
    if (hinges(ih)%adaptive_mesh) then
      adaptive_mesh = .true.
    endif
  enddo

  !> build adaptive mesh for hinged case   
  if ((n_hinges .ge. 1) .and. adaptive_mesh) then 
    
    !> allocate surface points for the component 
    allocate(rrv_le(2,nSections));  rrv_le = 0.0_wp 
    allocate(rrv_te(2,nSections));  rrv_te = 0.0_wp
    !> allocate aerodynamic center line for the component 
    allocate(ac_line(2,nSections));  ac_line = 0.0_wp 
    !> initialize first section 
    if (mesh_mirror) then 
      rrv_le (1,nRegions + 1) = -chord_list(nRegions + 1)*ref_chord_fraction
      rrv_te (1,nRegions + 1) = chord_list(nRegions + 1)*(1 - ref_chord_fraction)
    else
      rrv_le (1,1) = -chord_list(1)*ref_chord_fraction
      rrv_te (1,1) = chord_list(1)*(1 - ref_chord_fraction)
    endif     
    do iRegion = 1, nRegions
      !> aerodynamic center line      
      if (mesh_mirror) then 
        ac_line(1, iRegion + 1) = span_list(iRegion)*sin(sweep_list(iRegion)*pi/180.0_wp) + ac_line(1, iRegion)
        ac_line(2, iRegion + 1) = -(span_list(iRegion) + ac_line(2, iRegion))
      else 
        ac_line(1, iRegion + 1) = span_list(iRegion)*sin(sweep_list(iRegion)*pi/180.0_wp) + ac_line(1, iRegion)
        ac_line(2, iRegion + 1) = +(span_list(iRegion) + ac_line(2, iRegion))
      endif 
    enddo

    if (mesh_mirror) then  
      ac_line(:, nRegions+1:1:-1) = ac_line
      chord_list(nRegions+1:1:-1) = chord_list

      do iRegion = 1, nRegions + 1
      !> leading edge and trailing edge points in wind axis  
        rrv_le(1, iRegion) = ac_line(1, iRegion) - chord_list(iRegion)*ref_chord_fraction 
        rrv_te(1, iRegion) = ac_line(1, iRegion) + chord_list(iRegion)*(1 - ref_chord_fraction) 
        rrv_le(2, iRegion) = ac_line(2, iRegion)
        rrv_te(2, iRegion) = ac_line(2, iRegion)
      end do 
      chord_list(nRegions+1:1:-1) = chord_list
    else 
    
      do iRegion = 1, nRegions
        !> leading edge and trailing edge points in wind axis  
        rrv_le(1, iRegion + 1) = ac_line(1, iRegion + 1) - chord_list(iRegion + 1)*ref_chord_fraction 
        rrv_te(1, iRegion + 1) = ac_line(1, iRegion + 1) + chord_list(iRegion + 1)*(1 - ref_chord_fraction) 
        rrv_le(2, iRegion + 1) = ac_line(2, iRegion + 1)
        rrv_te(2, iRegion + 1) = ac_line(2, iRegion + 1)
      end do 
    
    endif 
    
    !> add ref_point (starting_point) to rrv and ac_line 
    rrv_le(1,:)   = ref_point(1) + rrv_le(1,:) 
    rrv_te(1,:)   = ref_point(1) + rrv_te(1,:) 
    ac_line(1,:)  = ref_point(1) + ac_line(1,:) 

    if (mesh_mirror) then
      rrv_le(2,:)   = -ref_point(2) + rrv_le(2,:) 
      rrv_te(2,:)   = -ref_point(2) + rrv_te(2,:) 
      ac_line(2,:)  = -ref_point(2) + ac_line(2,:)       
    else
      rrv_le(2,:)   = ref_point(2) + rrv_le(2,:) 
      rrv_te(2,:)   = ref_point(2) + rrv_te(2,:) 
      ac_line(2,:)  = ref_point(2) + ac_line(2,:)    
    endif

    do iRegion = 1, nRegions
      do ih = 1, n_hinges
        if (mesh_mirror .and. hinges(ih)%adaptive_mesh ) then 
          if ((hinges(ih)%node2(2) .le. ac_line(2, iRegion + 1)) .and. & 
              (hinges(ih)%node2(2) .ge. ac_line(2, iRegion))) then 

              !> interpolate at hinge station to get the leading edge point
              call linear_interp((/ rrv_le(1, iRegion), rrv_le(1, iRegion + 1)/), &
                                  (/rrv_le(2, iRegion), rrv_le(2, iRegion + 1)/), & 
                                  hinges(ih)%node2(2), hinges(ih)%le2(1))
              !> interpolate at hinge station to get the trailing edge point
              call linear_interp((/ rrv_te(1, iRegion), rrv_te(1, iRegion + 1)/), &
                                  (/rrv_te(2, iRegion), rrv_te(2, iRegion + 1)/), & 
                                  hinges(ih)%node2(2), hinges(ih)%te2(1))
              hinges(ih)%le2(2) = hinges(ih)%node2(2) 
              hinges(ih)%te2(2) = hinges(ih)%node2(2)
              !> hinge chord  
              hinges(ih)%chord2 = abs(hinges(ih)%le2(1) - hinges(ih)%te2(1))
              !> hinge node adimensional location along chord
              hinges(ih)%csi2 = (hinges(ih)%node2(1) - hinges(ih)%le2(1)) / hinges(ih)%chord2 
          endif    

          if ((hinges(ih)%node1(2) .le. ac_line(2, iRegion + 1)) .and. & 
              (hinges(ih)%node1(2) .ge. ac_line(2, iRegion))) then 
              !> interpolate at hinge station to get the leading edge point
              call linear_interp((/rrv_le(1, iRegion), rrv_le(1, iRegion + 1)/), &
                                (/ rrv_le(2, iRegion), rrv_le(2, iRegion + 1)/), & 
                                hinges(ih)%node1(2), hinges(ih)%le1(1))
              !> interpolate at hinge station to get the trailing edge point
              call linear_interp((/rrv_te(1, iRegion), rrv_te(1, iRegion + 1)/), &
                                (/ rrv_te(2, iRegion), rrv_te(2, iRegion + 1)/), & 
                                hinges(ih)%node1(2), hinges(ih)%te1(1))
              hinges(ih)%le1(2) = hinges(ih)%node1(2) 
              hinges(ih)%te1(2) = hinges(ih)%node1(2)
              !> hinge chord
              hinges(ih)%chord1 = abs(hinges(ih)%le1(1) - hinges(ih)%te1(1))
              !> hinge node adimensional location along chord
              hinges(ih)%csi1 = (hinges(ih)%node1(1) - hinges(ih)%le1(1)) / hinges(ih)%chord1 
          endif 
        elseif ( hinges(ih)%adaptive_mesh) then
          if ((hinges(ih)%node1(2) .ge. ac_line(2, iRegion)) .and. & 
              (hinges(ih)%node1(2) .le. ac_line(2, iRegion + 1))) then 
              
              !> interpolate at hinge station to get the leading edge point
              call linear_interp((/ rrv_le(1, iRegion), rrv_le(1, iRegion + 1)/), &
                                  (/rrv_le(2, iRegion), rrv_le(2, iRegion + 1)/), & 
                                  hinges(ih)%node1(2), hinges(ih)%le1(1))

              !> interpolate at hinge station to get the trailing edge point
              call linear_interp((/ rrv_te(1, iRegion), rrv_te(1, iRegion + 1)/), &
                                  (/rrv_te(2, iRegion), rrv_te(2, iRegion + 1)/), & 
                                  hinges(ih)%node1(2), hinges(ih)%te1(1))
              
              hinges(ih)%le1(2) = hinges(ih)%node1(2) 
              hinges(ih)%te1(2) = hinges(ih)%node1(2)
              !> hinge chord  
              hinges(ih)%chord1 = abs(hinges(ih)%le1(1) - hinges(ih)%te1(1))
              !> hinge node adimensional location along chord
              hinges(ih)%csi1 = (hinges(ih)%node1(1) - hinges(ih)%le1(1)) / hinges(ih)%chord1 
          endif    

          if ((hinges(ih)%node2(2) .ge. ac_line(2, iRegion)) .and. & 
                  (hinges(ih)%node2(2) .le. ac_line(2, iRegion + 1))) then 
              !> interpolate at hinge station to get the leading edge point
              call linear_interp((/ rrv_le(1, iRegion), rrv_le(1, iRegion + 1)/), &
                                (/rrv_le(2, iRegion), rrv_le(2, iRegion + 1)/), & 
                                hinges(ih)%node2(2), hinges(ih)%le2(1))
              !> interpolate at hinge station to get the trailing edge point
              call linear_interp((/ rrv_te(1, iRegion), rrv_te(1, iRegion + 1)/), &
                                (/rrv_te(2, iRegion), rrv_te(2, iRegion + 1)/), & 
                                hinges(ih)%node2(2), hinges(ih)%te2(1))
              hinges(ih)%le2(2) = hinges(ih)%node2(2) 
              hinges(ih)%te2(2) = hinges(ih)%node2(2)
              !> hinge chord
              hinges(ih)%chord2 = abs(hinges(ih)%le2(1) - hinges(ih)%te2(1))
              !> hinge node adimensional location along chord
              hinges(ih)%csi2 = (hinges(ih)%node2(1) - hinges(ih)%le2(1)) / hinges(ih)%chord2 
          endif          
        endif 
      enddo
    enddo
    
    !> cast all the adimensional location into a single vector  
    allocate(csi_hinge_not_unique(2*size(hinges)+1)); csi_hinge_not_unique = 0.0_wp 
    do ih = 1, size(hinges)
      csi_hinge_not_unique(ih) = hinges(ih)%csi1
      csi_hinge_not_unique(ih+2) = hinges(ih)%csi2      
      merge_tol = merge_tol + hinges(ih)%merge_tol 
    enddo
    !> take average and avoid numerical issues
    merge_tol = merge_tol/real(size(hinges),wp) + 1e-16_wp

    !> delete doubled nodes or merge them in single one if the distance is lower the 1% the chord lenght 
    call unique(csi_hinge_not_unique, csi_hinge, merge_tol)  
    
    allocate(point_region(size(csi_hinge) + 1))
    allocate(delta_x(size(csi_hinge) + 1)); delta_x = 0.0_wp
    allocate(delta_x_no_off(size(csi_hinge) + 1)); delta_x_no_off = 0.0_wp ! maybe useless 
    
    do ia = 1, size(csi_hinge) + 1
      !> first leading edge region 
      if (ia .eq. 1) then 
        delta_x(ia) = csi_hinge(ia)
        point_region(ia) = floor(real(nelem_chord,wp)*delta_x(ia))
        allocate(csi_adim(point_region(ia) + 1)); csi_adim = 0.0_wp
        call define_division(type_chord, point_region(ia), csi_adim, & 
                            r, r_le, r_te, r_le_fix, r_le_moving, r_te_fix, r_te_moving, x_refinement)

        chord_fraction(1:point_region(ia) + 1) = delta_x(ia)*csi_adim 
        deallocate(csi_adim)

      !> trailing edge region 
      elseif (ia .eq. (size(csi_hinge) + 1)) then 
        delta_x_no_off(ia) = 1.0_wp - csi_hinge(ia-1) 
        delta_x(ia) = delta_x(ia - 1) + delta_x_no_off(ia) 
        point_region(ia) = nelem_chord - sum(point_region(1:ia - 1)) 
        
        allocate(csi_adim(point_region(ia) + 1)); csi_adim = 0.0_wp
        call define_division(type_chord, point_region(ia), csi_adim , &
                            r, r_le, r_te, r_le_fix, r_le_moving, r_te_fix, r_te_moving, x_refinement)

        chord_fraction((sum(point_region(1:ia - 1)) + 1) : (sum(point_region) + 1) ) = & 
            (delta_x(ia) - delta_x(ia - 1))*csi_adim + delta_x(ia - 1)  
        deallocate(csi_adim) 
        
      !> regions between two hinge points 
      else 
        delta_x_no_off(ia) = csi_hinge(ia) - csi_hinge(ia-1) 
        delta_x(ia) = delta_x(ia - 1) + delta_x_no_off(ia) 
        point_region(ia) = floor(real(nelem_chord,wp)*delta_x_no_off(ia)) 
        
        allocate(csi_adim(point_region(ia) + 1)); csi_adim = 0.0_wp
        call define_division(type_chord, point_region(ia), csi_adim, & 
                            r, r_le, r_te, r_le_fix, r_le_moving, r_te_fix, r_te_moving, x_refinement) 

        chord_fraction((sum(point_region(1:ia - 1)) + 1) : (sum(point_region(1:ia)) + 1) ) = & 
            (delta_x(ia) - delta_x(ia - 1))*csi_adim + delta_x(ia - 1)   
      
        deallocate(csi_adim)
      endif 
    enddo 
    deallocate(point_region, ac_line, rrv_le, rrv_te)
  else
    call define_division(type_chord, nelem_chord, chord_fraction, & 
                        r, r_le, r_te, r_le_fix, r_le_moving, r_te_fix, r_te_moving, x_refinement)
  endif 

  if (allocated(rrSection1))  deallocate(rrSection1)
  if (allocated(rrSection2)) deallocate(rrSection2) 
  ! Initialize the span division to the maximum dimension
  allocate(span_fraction(maxval(nelem_span_list))) ; span_fraction = 0.0_wp
  allocate(rrSection1(3,npoint_chord_tot)) ; rrSection1 = 0.0_wp
  allocate(rrSection2(3,npoint_chord_tot)) ; rrSection2 = 0.0_wp
  ! Initialise dr_ref for the definition of the actual reference_point
  !  of each bay (by updating)
  dx_ref = ref_point(1)       ! 0.0_wp
  dy_ref = ref_point(2)       ! 0.0_wp
  dz_ref = ref_point(3)       ! 0.0_wp

  ! === new-2019-02-06 ===
  ! check airfoil_list input ----
  if ( trim(airfoil_list(1)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The first "airfoil"&
          & cannot be set as "interp".')
  end if
  if ( trim(airfoil_list(nSections)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The last "airfoil"&
          & cannot be set as "interp".')
  end if
  
  ista = 1
  iend = npoint_chord_tot
  ! Loop over regions
  do iRegion = 1, nRegions
    if ( iRegion .gt. 1 ) then  ! first section = last section of the previous region
      rrSection1 = rrSection2
      thickness_section1(iRegion) = thickness_section2(iRegion-1)      
    else   ! first section                      ! build points
      call define_section(chord_list(iRegion), trim(adjustl(airfoil_list(iRegion))), &
                          twist_list(iRegion), ElType, nelem_chord,              &
                          type_chord , chord_fraction, ref_chord_fraction,       &
                          ref_point, xySection1, thickness_section)
        rrSection1(1,:) = xySection1(1,:) + ref_point(1)
        rrSection1(2,:) = 0.0_wp          + ref_point(2)     ! <--- read from region structure
        rrSection1(3,:) = xySection1(2,:) + ref_point(3)
      ! Update rr
      rr(:,ista:iend) = rrSection1
      thickness_section1(iRegion) = thickness_section
    end if
    
    ! It is possible to define the airfoils on some of the sections only.
    !  When the shape of the airfoil is not defined on a section, it is interpolated
    if ( trim(adjustl(airfoil_list(iRegion + 1))) .ne. 'interp' ) then  ! read the field 'airfoil'

      call define_section( chord_list(iRegion+1), trim(adjustl(airfoil_list(iRegion+1))), &
                            twist_list(iRegion+1), ElType, nelem_chord,                    &
                            type_chord , chord_fraction, ref_chord_fraction,               &
                            ref_point, xySection2, thickness_section)
      thickness_section2(iRegion)         = thickness_section ! second section 
    else ! interpolation
      do iSec = iRegion + 1 , nRegions+1
        if ( airfoil_list(iSec) .ne. 'interp' ) then
          dy_actual_airfoils = sum( abs(span_list(iRegion:iSec-1)) )
        exit
        end if
      end do
      dy_sections = abs(span_list(iRegion))
      csi = dy_sections / dy_actual_airfoils ! adimensional "coord" for interpolation

      call define_section( 1.0_wp , trim(adjustl(airfoil_list(iSec))), &
                            0.0_wp , ElType, nelem_chord,               &
                            type_chord , chord_fraction, 0.0_wp,        &
                            ref_point, xyAirfoil2)

      ! Compute the coordinates xySection2(), after removing the offset
      if ( allocated(xySection2) ) deallocate(xySection2)
      
      allocate(xySection2(size(xyAirfoil2,1),size(xyAirfoil2,2)))

      if ( allocated(xyAirfoil1) ) deallocate(xyAirfoil1)

      allocate(xyAirfoil1(size(xyAirfoil2,1),size(xyAirfoil2,2)))

      ! coords of prev section, back in the local normalized frame 0->1
      xyAirfoil1(1,:) = ( rrSection1(1,:) - dx_ref ) / chord_list(iRegion)
      xyAirfoil1(2,:) = ( rrSection1(3,:) - dz_ref ) / chord_list(iRegion)

      twist_rad = twist_list(iRegion) * 4.0_wp * atan(1.0_wp) / 180.0_wp
      xyAirfoil1 = matmul( &
                  reshape( (/ cos(twist_rad), sin(twist_rad) , &
                            -sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , &
                                                              xyAirfoil1 )
      
      xyAirfoil1(1,:) = xyAirfoil1(1,:)+ ref_chord_fraction
      
      ! linear interpolation, casting back to ref_chord frame
      xySection2 = ( (1-csi) * xyAirfoil1 + &
                        csi  * xyAirfoil2 ) * chord_list(iRegion+1)
      xySection2(1,:) = xySection2(1,:) - ref_chord_fraction*chord_list(iRegion+1)

      twist_rad = twist_list(iRegion+1) * 4.0_wp * atan(1.0_wp) / 180.0_wp
      xySection2 = matmul( &
                  reshape( (/ cos(twist_rad),-sin(twist_rad) , &
                              sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , &
                                                          xySection2 )

    end if

    if ( abs( sweep_list(iRegion) ) .gt. 60.0_wp ) then
      call warning(this_sub_name, this_mod_name, 'Requested a sweep angle of &
        &more than 60 degrees. Are you sure of this input?')
    end if
    if ( abs( dihed_list(iRegion) ) .gt. 60.0_wp ) then
      call warning(this_sub_name, this_mod_name, 'Requested a dihedral angle of &
        &more than 60 degrees. Are you sure of this input?')
    end if

    !> save before update
    dx_ref_1 = dx_ref;  dy_ref_1 = dy_ref;  dz_ref_1 = dz_ref

    dx_ref = span_list(iRegion) * tan( sweep_list(iRegion)* pi / 180.0_wp ) + dx_ref
    dy_ref = span_list(iRegion)                                             + dy_ref
    dz_ref = span_list(iRegion) * tan( dihed_list(iRegion)* pi / 180.0_wp ) + dz_ref

    rrSection2(1,:) = xySection2(1,:) + dx_ref
    rrSection2(2,:) = 0.0_wp          + dy_ref  ! <--- read from region structure
    rrSection2(3,:) = xySection2(2,:) + dz_ref

    ! Interpolation of the nodes of the region i (between sections i and i+1)
    do i1 = 1 , nelem_span_list(iRegion)

      ista = iend + 1
      iend = iend + npoint_chord_tot

      if ( .not. twist_linear_interp ) then

        if ( trim(type_span_list(iRegion)) .eq. 'uniform' ) then
          ! uniform spacing in span
          rr(:,ista:iend) = rrSection1 + real(i1,wp) / &
                            real(nelem_span_list(iRegion),wp) * &
                            ( rrSection2 - rrSection1 )
          adjust_xz = .false.

        elseif ( trim(type_span_list(iRegion)) .eq. 'cosine' ) then
          ! cosine  spacing in span
          rr(:,ista:iend) = 0.5_wp * ( rrSection1 + rrSection2 ) - &
                            0.5_wp * ( rrSection2 - rrSection1 ) * &
                      cos( real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) )
          adjust_xz = .false.

        elseif ( trim(type_span_list(iRegion)) .eq. 'cosineOB' ) then
          ! cosine  spacing in span: outboard refinement
          rr(:,ista:iend) = rrSection1 + &
                          ( rrSection2 - rrSection1 ) * &
              sin( 0.5_wp*real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) )
          adjust_xz = .false. 

        elseif ( trim(type_span_list(iRegion)) .eq. 'cosineIB' ) then
          ! cosine  spacing in span: inboard refinement
          rr(:,ista:iend) = rrSection2 - &
                          ( rrSection2 - rrSection1 ) * &
              cos( 0.5_wp*real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) )
          adjust_xz = .false.

        elseif ( trim(type_span_list(iRegion)) .eq. 'equalarea' ) then 
          rr(2, ista:iend) = sqrt(real(i1,wp)/real(nelem_span_list(iRegion),wp))* & 
                            (rrSection2(2,:) - rrSection1(2,:)) + rrSection1(2,:) 
          adjust_xz = .true.

        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries' ) then
          do j = 1, npoint_chord_tot
            call geoseries_both(rrSection1(2,j), rrSection2(2,j), nelem_span_list(iRegion), & 
                                y_ref_list(iRegion), r_in_list(iRegion), r_ob_list(iRegion), division)  
            rr(2, ista + j - 1) = division(i1 + 1)
          enddo 
          adjust_xz = .true.

        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries_ob' ) then
          do j = 1, npoint_chord_tot
            call geoseries(rrSection1(2,j), rrSection2(2,j), nelem_span_list(iRegion), & 
                          r_ob_list(iRegion), divisionIB, divisionOB)  
            rr(2, ista + j - 1) = divisionOB(i1 + 1)
          enddo 
          adjust_xz = .true.

        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries_ib' ) then
          do j = 1, npoint_chord_tot
            call geoseries(rrSection1(2,j), rrSection2(2,j), nelem_span_list(iRegion), &
                          r_in_list(iRegion), divisionIB, divisionOB)  
            rr(2, ista + j - 1) = divisionIB(i1 + 1)
          enddo 
          adjust_xz = .true.
        
        else 
          write(*,*) ' Mesh file   : ' , trim(mesh_file)
          write(*,*) ' type_span   : ' , trim(type_span_list(iRegion))
          call error(this_sub_name, this_mod_name, 'Incorrect input: &
                & type_span must be equal to uniform, cosine, cosineIB, cosineOB, equalarea, &
                & geoseries, geoseries_ib, geoseries_ob.')
        end if
        !> adjust the x and z coordinate 
        if (adjust_xz) then
          rr(1,ista:iend) = rrSection1(1,:) + (rr(2, ista:iend) - rrSection1(2,:))*&
                            (rrSection2(1,:) - rrSection1(1,:))/(rrSection2(2,:) - rrSection1(2,:))
          rr(3,ista:iend) = rrSection1(3,:) + (rr(2, ista:iend) - rrSection1(2,:))*&
                            (rrSection2(3,:) - rrSection1(3,:))/(rrSection2(2,:) - rrSection1(2,:))
        endif  
      else !-> linear interpolation of the twist angle

        allocate(rr_tw(  2,npoint_chord_tot))
        allocate(rr_tw_1(2,npoint_chord_tot))
        allocate(rr_tw_2(2,npoint_chord_tot))

        ! === Transform sections back to local reference frames ===
        ! Section 1
        !> remove offset in x,z
        rr_tw_1(1,:) = rrSection1(1,:) - dx_ref_1
        rr_tw_1(2,:) = rrSection1(3,:) - dz_ref_1
        !> rotate section back ( coord. in the local ref. frame )
        twist_rad = twist_list(iRegion) * pi/180.0_wp
        rr_tw_1 = matmul( &
                  reshape( (/ cos(twist_rad), sin(twist_rad) , &
                              -sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , &
                                                                rr_tw_1 )
        ! Section 2
        !> remove offset in x,z
        rr_tw_2(1,:) = rrSection2(1,:) - dx_ref
        rr_tw_2(2,:) = rrSection2(3,:) - dz_ref
        !> rotate section back ( coord. in the local ref. frame )
        twist_rad = twist_list(iRegion+1) * pi/180.0_wp
        rr_tw_2 = matmul( &
                  reshape( (/ cos(twist_rad), sin(twist_rad), &
                            -sin(twist_rad), cos(twist_rad)/), (/2,2/)), rr_tw_2 )
        !> Interpolation weight 
        if ( trim(type_span_list(iRegion)) .eq. 'uniform' ) then
          interp_weight = real(i1,wp) / real(nelem_span_list(iRegion),wp)
        else if ( trim(type_span_list(iRegion)) .eq. 'cosine' ) then
          interp_weight = 0.5_wp * ( &
                        1.0_wp - cos( real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) ) )
        else if ( trim(type_span_list(iRegion)) .eq. 'cosineOB' ) then
          interp_weight = sin( 0.5_wp*real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) )
        else if ( trim(type_span_list(iRegion)) .eq. 'cosineIB' ) then
          interp_weight = 1.0_wp - cos( 0.5_wp*real(i1,wp)*pi/ real(nelem_span_list(iRegion),wp) )
        else if ( trim(type_span_list(iRegion)) .eq. 'equalarea' ) then
          interp_weight = sqrt(real(i1,wp)/real(nelem_span_list(iRegion),wp))
        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries' ) then
          call geoseries_both(0.0_wp, 1.0_wp, nelem_span_list(iRegion), &
                              y_ref_list(iRegion), r_in_list(iRegion), r_ob_list(iRegion), division)
          interp_weight = division(i1 + 1) 
        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries_ob' ) then
          call geoseries(0.0_wp, 1.0_wp, nelem_span_list(iRegion), & 
            r_ob_list(iRegion), divisionIB, divisionOB)
          interp_weight = divisionOB(i1 + 1) 
        elseif ( trim(type_span_list(iRegion)) .eq. 'geoseries_ib' ) then
          call geoseries(0.0_wp, 1.0_wp, nelem_span_list(iRegion), r_in_list(iRegion), divisionIB, divisionOB)
          interp_weight = divisionIB(i1 + 1) 
        else
          write(*,*) ' Mesh file   : ' , trim(mesh_file)
          write(*,*) ' type_span   : ' , trim(type_span_list(iRegion))
          call error(this_sub_name, this_mod_name, 'Incorrect input: &
                & type_span must be equal to uniform, cosine, cosineIB, cosineOB, equalarea, &
                & geoseries, geoseries_ib, geoseries_ob.')
        end if

        ! === x,z coordinates ===
        rr_tw = rr_tw_1 + ( rr_tw_2 - rr_tw_1 ) * interp_weight
        !> rotation ( linear interpolation of the twist angle )
        twist_rad = pi/180.0_wp * ( twist_list(iRegion) + &
            ( twist_list(iRegion+1)-twist_list(iRegion) ) * interp_weight )
        rr_tw = matmul( &
                reshape( (/ cos(twist_rad),-sin(twist_rad) , &
                          sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , &
                                                                  rr_tw )
        rr(1,ista:iend) = rr_tw(1,:) + dx_ref_1 + &
                            ( dx_ref - dx_ref_1 ) * interp_weight
        rr(3,ista:iend) = rr_tw(2,:) + dz_ref_1 + &
                            ( dz_ref - dz_ref_1 ) * interp_weight
        ! === y coordinate ===
        rr(2,ista:iend) = rrSection1(2,:) + &
                        ( rrSection2(2,:) - rrSection1(2,:) ) * interp_weight

        deallocate(rr_tw, rr_tw_1, rr_tw_2)

      end if

    end do ! iRegion 

  enddo ! iRegion 
  
  ! lots of deallocation missing causing memory leakage 
  if ( allocated(xySection1) ) deallocate(xySection1)
  if ( allocated(xySection2) ) deallocate(xySection2)
  

  !> Interpolation of airfoil table for corrected vl only
  if (ElType .eq.'v' .and. aero_table) then
    allocate(i_airfoil_e(2,nelem_span_tot))
    i_airfoil_e = 0
    
    allocate(normalised_coord_e(2,nelem_span_tot))
    normalised_coord_e = 0.0_wp
    allocate(thickness(2,nelem_span_tot)); 
    thickness = 0.0_wp
      !> Check airfoil_list()
    if ( trim(airfoil_table_list(1)) .eq. 'interp' ) then
      call error(this_sub_name, this_mod_name, 'The first "airfoil"&
            & cannot be set as "interp".')
    end if
    if ( trim(airfoil_table_list(nSections)) .eq. 'interp' ) then
      call error(this_sub_name, this_mod_name, 'The last "airfoil"&
            & cannot be set as "interp".')
    end if

    nAirfoils = 0
    do i = 1 , nSections
      if ( trim(airfoil_table_list(i)) .ne. 'interp' ) nAirfoils = nAirfoils + 1
    end do

    allocate(airfoil_list_actual(nAirfoils))

    iAirfoil = 0
    
    do i = 1 , nSections
      if ( trim(airfoil_table_list(i)) .ne. 'interp' ) then
        iAirfoil = iAirfoil + 1
        airfoil_list_actual(iAirfoil) = trim(airfoil_table_list(i))
        call check_file_exists(airfoil_list_actual(iAirfoil), this_sub_name, &
              this_mod_name)
      end if
    end do

    i_aero2 = 1 
    ista = 1 
    iAirfoil = 0
    
    do iAirfoil = 1, nAirfoils - 1

      i_aero1 = i_aero2
      do i = i_aero1 + 1, nSections
        if (airfoil_table_list(i) .ne. 'interp') then
          i_aero2 = i 
          exit
        end if
      end do
      
      iend = sum(nelem_span_list(1:i_aero2-1))
      
      do iSpan = ista, iend
        i_airfoil_e(1, iSpan) = iAirfoil
        i_airfoil_e(2, iSpan) = iAirfoil + 1
        
        normalised_coord_e(2,iSpan) = (rr(2,(iSpan+1)*npoint_chord_tot) - rr(2,ista*2)) / &
                                      (rr(2,(iend +1)*npoint_chord_tot) - rr(2,ista*2))
        !> curvature interpolation 
        !> inboard
        if (iSpan .eq. ista) then 
          thickness(1,iSpan) = thickness_section1(iAirfoil) 
        else          
          thickness(1,iSpan) = thickness(2,iSpan - 1)
        endif 
        !> outboard
        thickness(2,iSpan) = thickness_section1(iAirfoil)*(1.0_wp - normalised_coord_e(2,iSpan)) + &
                              thickness_section2(iAirfoil)*normalised_coord_e(2,iSpan)

        if ( iSpan .ne. ista ) then 
          normalised_coord_e(1,iSpan) = normalised_coord_e(2,iSpan-1)          
        end if

      end do

      ista = iend + 1  

    end do
  endif

  ! optional output ----
  npoints_chord_tot = npoint_chord_tot
  if (present(aero_table_out)) then
    aero_table_out = aero_table
    !thickness_out = thickness
  end if
  ! optional output ----

end subroutine read_mesh_parametric
!-------------------------------------------------------------------------------

subroutine define_section(chord, airfoil, twist, ElType, nelem_chord, &
                          type_chord , chord_fraction, reference_chord_fraction,&
                          reference_point, point_list, thickness_out)

  real(wp), allocatable , intent(out)               :: point_list(:,:)
  real(wp), intent(in)                              :: reference_point(:) 
  real(wp), intent(inout)                           :: chord_fraction(:)
  character(len=*) , intent(in)                     :: type_chord
  real(wp), intent(in)                              :: reference_chord_fraction, twist, chord
  integer, intent(in)                               :: nelem_chord
  character, intent(in)                             :: ElType
  character(len=*) , intent(in)                     :: airfoil
  real(wp),  intent(out), optional                  :: thickness_out
  real(wp)                                          :: thickness
  real(wp)                                          :: twist_rad
  real(wp), allocatable                             :: points_mean_line(:,:)
  integer                                           :: i, ascii
  character(len=*), parameter                       :: this_sub_name='define_section'

  twist_rad = twist * 4.0_wp * atan(1.0_wp) / 180.0_wp

  ! Airfoil geometry 
  ! - read coordinate from file: if <airfoil> is 'xxxxxx.dat'
  ! - build as a member of a airfoil family:
  !   - NACA 4-digit: NACAmpss
  !   - NACA 5-digit: NACAxxxxx
  !   - ...
  if ( airfoil(len_trim(airfoil)-3 : len_trim(airfoil)) .eq. '.dat' ) then

    call check_file_exists(airfoil, this_sub_name, this_mod_name)
    call read_airfoil ( airfoil , ElType , nelem_chord , chord_fraction, point_list, thickness)

  else if ( airfoil(1:4) .eq. 'NACA' ) then
  
    ! Check on NACA syntax, only digits and spaces are allowed
    do i= 5, len(airfoil) 
      ascii = iachar(airfoil(i:i))
      if ( (ascii .gt. 57 .or. ascii .lt. 48) .and. ascii .ne. 32 ) then
        call error(this_sub_name, this_mod_name, ' only 4-digit and some 5-digit &
      &NACA airfoils implemented. Check your input syntax or provide the coordinates&
      & of the airfoil as a .dat file ')
      end if
    end do
  
    if ( len_trim(airfoil) .eq. 8 ) then      ! NACA 4-digit -------
      call naca4digits(airfoil(5:8), nelem_chord, chord_fraction, &
                        points_mean_line , point_list, thickness)

    elseif ( len_trim(airfoil) .eq. 9 ) then
      call naca5digits(airfoil(5:9), nelem_chord, chord_fraction, &
                        points_mean_line , point_list, thickness)     
    else
      call error(this_sub_name, this_mod_name, ' only 4-digit and some 5-digit &
      &NACA airfoils implemented. Check your input syntax or provide the coordinates&
      & of the airfoil as a .dat file ')    
    end if
    if (ElType .eq. 'v') then 
      point_list = points_mean_line
    end if
  else
    call error(this_sub_name, this_mod_name, ' only 4-digit and some 5-digit &
      &NACA airfoils implemented. Check your input syntax or provide the coordinates&
      & of the airfoil as a .dat file ')
  end if

  ! Geometric transformation +++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 1. translation : reference_chord_fraction (defined for a unit-chord airfoil)
  point_list(1,:) = point_list(1,:) - reference_chord_fraction
  ! 2. scaling     : chord
  point_list = point_list * chord
  ! 3. rotation    : twist
  point_list        = matmul( reshape( (/ cos(twist_rad),-sin(twist_rad) , &
                                          sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , &
                                                                    point_list )
  thickness = thickness * chord 
  !> optional output
  if (present(thickness_out)) then
    thickness_out = thickness
  end if
end subroutine define_section

!-------------------------------------------------------------------------------

subroutine naca4digits(airfoil_name, nelem_chord,&
                      chord_fraction, points_mean_line , points, thick)
  character(len=*), intent(in) :: airfoil_name
  integer , intent(in)  :: nelem_chord
  real(wp), intent(in)  :: chord_fraction(:)
  real(wp), allocatable , intent(out) :: points_mean_line(:,:), points(:,:)
  real(wp), intent(out) :: thick 
  real(wp), allocatable :: points_upper(:,:), points_lower(:,:)

  integer :: mm, pp, ss, iPoint
  real(wp) :: m,p,s, xa, xac, ml, theta, thickness
  character :: str1
  character(len=2) :: str2
  integer :: ierr

  str1 = airfoil_name(1:1)
  read(str1,*,iostat=ierr) mm
  str1 = airfoil_name(2:2)
  read(str1,*,iostat=ierr) pp
  str2 = airfoil_name(3:4)
  read(str2,*,iostat=ierr) ss

  m = real(mm,wp)/100.0_wp
  p = real(pp,wp)/10.0_wp
  s = real(ss,wp)/100.0_wp

  allocate(points_upper(2,nelem_chord+1)) ; points_upper = 0.0_wp
  allocate(points_lower(2,nelem_chord+1)) ; points_lower = 0.0_wp

  if ( allocated(points_mean_line) )  deallocate(points_mean_line)
  allocate(points_mean_line(2,nelem_chord+1)) ; points_mean_line = 0.0_wp
  
  do iPoint = 1,nelem_chord+1

    xa = chord_fraction(iPoint)
    ! Define mean line ml and local slope theta
    ml = 0.0_wp
    theta = 0.0_wp
    if (p>0) then
      if (xa <= p) then
        ml = m/p**2 * (2.0_wp*p*xa - xa**2)
        theta = 2.0_wp*m/p**2 * (p - xa)
      else
        ml = m/(1.0_wp-p)**2 * (1.0_wp-2.0_wp*p + 2.0_wp*p*xa - xa**2)
        theta = 2.0_wp*m/(1.0_wp-p)**2 * (p - xa)
      endif
    endif

    ! Thickness
    thickness = 5.0_wp*s*( 0.2969_wp*sqrt(max(xa, 1e-16_wp))  &
                          - 0.1260_wp* xa       &
                          - 0.3516_wp*(xa**2)   &
                          + 0.2843_wp*(xa**3)   &
                          - 0.1036_wp*(xa**4))     ! closed TE

    points_mean_line(1,iPoint) = xa
    points_mean_line(2,iPoint) = ml

    points_upper(1,iPoint) = xa - thickness*sin(theta)
    points_upper(2,iPoint) = ml + thickness*cos(theta)

    points_lower(1,iPoint) = xa + thickness*sin(theta)
    points_lower(2,iPoint) = ml - thickness*cos(theta)

  enddo

  thick = s

  if ( allocated(points) ) deallocate(points)
  allocate(points(2,2*nelem_chord+1))
  points(:,            1:  nelem_chord  ) = points_lower(:,nelem_chord+1:2:-1)
  points(:,nelem_chord+1:2*nelem_chord+1) = points_upper

endsubroutine naca4digits

!----------------------------------------------------------------------

subroutine naca5digits(airfoil_name, nelem_chord,&
                        chord_fraction, points_mean_line , points, thick)
  character(len=*), intent(in)          :: airfoil_name
  integer , intent(in)                  :: nelem_chord
  real(wp), intent(in)                  :: chord_fraction(:)
  real(wp), allocatable , intent(out)   :: points_mean_line(:,:), points(:,:)
  real(wp),               intent(out)   :: thick
  real(wp), allocatable                 :: points_upper(:,:), points_lower(:,:)
  integer :: ss, iPoint, L, Q, P
  real(wp) :: s, xa, xac, ml, theta, thickness, mult, r, k1
  character :: str1
  character(len=2) :: str2
  integer :: ierr
  character(len=*), parameter :: this_sub_name = 'naca5digits'


  str1 = airfoil_name(1:1)
  read(str1,*,iostat=ierr) L
  str1 = airfoil_name(2:2)
  read(str1,*,iostat=ierr) P
  str1 = airfoil_name(3:3)
  read(str1,*,iostat=ierr) Q
  str2 = airfoil_name(4:5)
  read(str2,*,iostat=ierr) ss

  s = real(ss,wp)/100.0_wp

  !A limited number of airfoils have been implemented, perform some checks...
  if(Q .ne. 0) call error(this_sub_name, this_mod_name, &
                          'Reverse 5 digits naca profiles not yet implemented')
  select case(P)
    case(1)
      r = 0.0580_wp
      k1 = 361.400_wp
    case(2)
      r = 0.1260_wp
      k1 = 51.640_wp
    case(3)
      r = 0.2025_wp
      k1 = 15.957_wp
    case(4)
      r = 0.2900_wp
      k1 = 6.643_wp
    case(5)
      r = 0.0580_wp
      k1 = 3.230_wp
    case default
      call error(this_sub_name, this_mod_name, &
                          '5 digit naca not valid')
  end select

  !camber is created for 0.3 of design Cl (L=2). All others are generated by
  !linearly scaling the camber
  mult = (real(L,wp)*3.0_wp/20.0_wp) / 0.3_wp

  allocate(points_upper(2,nelem_chord+1)) ; points_upper = 0.0_wp
  allocate(points_lower(2,nelem_chord+1)) ; points_lower = 0.0_wp

  if ( allocated(points_mean_line) )  deallocate(points_mean_line)
  allocate(points_mean_line(2,nelem_chord+1)) ; points_mean_line = 0.0_wp

  do iPoint = 1,nelem_chord+1

    xa = chord_fraction(iPoint)

    ! Define mean line ml and local slope theta
    ml = 0.0_wp
    theta = 0.0_wp
    if (xa <= r) then
      ml = mult*k1/6.0_wp*(xa**3 -3.0_wp*r*xa**2+r**2*(3.0_wp-r)*xa)
      theta = mult*k1/6.0_wp*(3.0_wp*xa**2 -6.0_wp*r*xa+r**2*(3.0_wp-r))
    else
      ml = mult*k1*r**3/6.0_wp*(1-xa)
      theta = -mult*k1*r**3/6.0_wp
    endif

    ! Thickness
    thickness = 5.0_wp*s*(0.2969_wp*sqrt(max(xa,1e-16_wp)) - 0.1260_wp*xa - &
                0.3516_wp*(xa**2) + 0.2843_wp*(xa**3) - 0.1015_wp*(xa**4))

    points_mean_line(1,iPoint) = xa
    points_mean_line(2,iPoint) = ml

    points_upper(1,iPoint) = xa - thickness*sin(theta)
    points_upper(2,iPoint) = ml + thickness*cos(theta)

    points_lower(1,iPoint) = xa + thickness*sin(theta)
    points_lower(2,iPoint) = ml - thickness*cos(theta)

  enddo

  thick = s 

  if ( allocated(points) ) deallocate(points)
  allocate(points(2,2*nelem_chord+1))
  points(:,            1:  nelem_chord  ) = points_lower(:,nelem_chord+1:2:-1)
  points(:,nelem_chord+1:2*nelem_chord+1) = points_upper


endsubroutine naca5digits

!-------------------------------------------------------------------------------

subroutine read_airfoil (filen, ElType , nelems_chord , csi_half, rr, thickness )

  character(len=*), intent(in)                :: filen
  character(len=*), intent(in)                :: ElType
  integer         , intent(in)                :: nelems_chord
  real(wp)        , intent(in)                :: csi_half(:)
  real(wp), allocatable , intent(out)         :: rr(:,:)
  real(wp)        , intent(out)               :: thickness

  real(wp) , allocatable      :: rr_dat(:,:) 
  real(wp) , allocatable      :: rr_up(:,:), rr_low(:,:)
  real(wp) , allocatable      :: rr_mean(:,:)
  integer                     :: np_dat, min_X, idx_m, n
  real(wp)                    :: m, p
  real(Wp)                    :: mean_id_min, tang_le_low, tang_le_up
  real(wp) , allocatable      :: point_up_updated(:,:), point_low_updated(:,:), point_le_updated(:,:)
  real(wp) , allocatable      :: xq_mean(:), yq_mean(:) 
  real(wp) , allocatable      :: rr2mesh_up(:,:), rr2mesh_low(:,:), rr2mesh_mean(:,:), rr2mesh(:,:)
  real(wp) , allocatable      :: alpha_up(:), alpha_low(:), alpha_mean(:)
  real(wp) , allocatable      :: dl_up(:), dl_low(:), dl_mean(:) 
  real(wp) , allocatable      :: rr_geo_up(:,:), rr_geo_low(:,:), rr_geo_mean(:,:)
  real(wp) , allocatable      :: rr_sort(:,:), x_sort(:), x_unique(:)
  real(wp) , allocatable      :: csi_up(:), csi_low(:), csi_mean(:) 
  real(wp) , allocatable      :: rr_unique_up(:,:), rr_unique_low(:,:), rr_mean_dat(:,:)
  real(wp) , allocatable      :: xq_up(:), yq_up(:), xq_le(:), yq_le(:), yq_low(:), xq_low(:)
  real(wp)                    :: tang_te_up, tang_te_low, tang_mean_le, tang_mean_te, tang_mean_le_new
  real(wp)                    :: signed_area
  integer  , parameter        :: n_query = 20000 
  integer  , parameter        :: maxiter = 100
  integer  , parameter        :: num_points = 10000 !> density of points for interpolation 
  real(wp) , parameter        :: max_tang_te = 2.5_wp !> max allowed slope at TE (~68 deg)
  integer                     :: fid, ierr
  integer                     :: i1 , id_le, i, j, ind_x 
  integer  , allocatable      :: ind_sort(:)
  character(len=*), parameter :: this_sub_name = 'read_airfoil' 
  
  !> Read coordinates from dat file 
  call new_file_unit(fid, ierr)
  open(unit=fid,file=trim(adjustl(filen)) )
  read(fid,*) np_dat

  allocate(rr_dat(2,np_dat)); rr_dat = 0.0_wp 
  
  do i1 = 1 , np_dat
    read(fid,*) rr_dat(:,i1)
  end do
  close(fid)  
  
  !> retro-compatibility with old vl format 
  if ((ElType .eq. 'v') .and. (rr_dat(1, 1) .eq. 0.0_wp) .and. (rr_dat(1, np_dat) .eq. 1.0_wp)) then 
    allocate(rr_mean(2, np_dat));                rr_mean = 0.0_wp 
    rr_mean(1,:) = rr_dat(1, :)
    rr_mean(2,:) = rr_dat(2, :) 
    m = maxval(rr_mean(2,:))
    idx_m = maxloc(rr_mean(2,:),1)  ! y_mean line max thickness
    p = rr_mean(1, idx_m)           ! x_mean line max thickness 
    tang_mean_le = 2.0_wp*m/p ! need to get the center of the circle -> linearized approach
    tang_mean_te = (rr_mean(2, size(rr_mean,2) - 1) - rr_mean(2, size(rr_mean,2)))/ &
                   (rr_mean(1, size(rr_mean,2) - 1) - rr_mean(1, size(rr_mean,2)))
  
    allocate(xq_mean(n_query)); xq_mean = 0.0_wp
    allocate(yq_mean(n_query)); yq_mean = 0.0_wp 
    call linspace(0.0_wp, 1.0_wp, xq_mean)  
    !> start iterative loop to get a more precise value of tang_mean_le 
    do i = 1, maxiter
      call hermite_spline_profile(rr_mean(1,:), rr_mean(2,:), xq_mean, &
                                  tang_mean_le, tang_mean_te, yq_mean)
      m = maxval(yq_mean)
      idx_m = maxloc(yq_mean,1)   ! y_mean line max thickness
      p = xq_mean(idx_m)           ! x_mean line max thickness
      tang_mean_le_new = 2.0_wp*m/p
      if (abs(tang_mean_le_new - tang_mean_le) .lt. 1e-6_wp) exit
      tang_mean_le = tang_mean_le_new 
    enddo 
    !> mean line 
    allocate(rr2mesh_mean(2,size(xq_mean))); rr2mesh_mean = 0.0_wp
    rr2mesh_mean(1,:) = xq_mean
    rr2mesh_mean(2,:) = yq_mean
    !> get length and derivative (alpha) along the mean line
    allocate(alpha_mean(size(rr2mesh_mean(1,:)))); alpha_mean = 0.0_wp
    alpha_mean(1) = tang_mean_le
    allocate(dl_mean(size(rr2mesh_mean(1,:))));     dl_mean = 0.0_wp
    dl_mean(1) = 0.0_wp
    do i = 2, size(rr2mesh_mean(1,:))
      dl_mean(i) = dl_mean(i - 1) + sqrt((rr2mesh_mean(1,i) - rr2mesh_mean(1,i-1))**2.0_wp + &
                                         (rr2mesh_mean(2,i) - rr2mesh_mean(2,i-1))**2.0_wp)
      alpha_mean(i) = atan2(rr2mesh_mean(2,i) - rr2mesh_mean(2,i-1), &
                            rr2mesh_mean(1,i) - rr2mesh_mean(1,i-1))
    end do
    !> interpolation of the mean line
    allocate(csi_mean(size(csi_half))); csi_mean = 0.0_wp
    csi_mean = csi_half*dl_mean(size(dl_mean))
    allocate(rr_geo_mean(2, size(csi_mean))); rr_geo_mean = 0.0_wp
    do i = 1, size(csi_mean)
      do j = 1, size(dl_mean) - 1
        if ((csi_mean(i) .ge. dl_mean(j)) .and. (csi_mean(i) .le. dl_mean(j + 1))) then
          rr_geo_mean(1,i) = rr2mesh_mean(1, j) + cos(alpha_mean(j))*(csi_mean(i) - dl_mean(j))
          rr_geo_mean(2,i) = rr2mesh_mean(2, j) + sin(alpha_mean(j))*(csi_mean(i) - dl_mean(j))
        endif
      enddo
    enddo
    thickness = 0.12_wp 
  else 

    !> leading edge index
    id_le = minloc(rr_dat(1,:), 1)

  !---------------------------------------------------------------
  ! Determine contour orientation using the signed polygon area.
  !  signed_area > 0 : counter-clockwise  (Selig)
  !  signed_area < 0 : clockwise          (Old DUST)
  !---------------------------------------------------------------
    n = size(rr_dat, 2)
    signed_area = 0.0_wp

    do i = 1, n-1
       signed_area = signed_area + &
          rr_dat(1,i) * rr_dat(2,i+1) - &
          rr_dat(1,i+1) * rr_dat(2,i)
    end do

    signed_area = signed_area + &
       rr_dat(1,n) * rr_dat(2,1) - &
       rr_dat(1,1) * rr_dat(2,n)

    signed_area = 0.5_wp * signed_area

    if (signed_area > 0.0_wp) then
       ! Selig format: TE -> upper -> LE -> lower -> TE
       allocate(rr_up(2, id_le));                    rr_up = 0.0_wp
       allocate(rr_low(2, size(rr_dat,2)-id_le+1));  rr_low = 0.0_wp

       rr_up  = rr_dat(:,1:id_le)
       rr_up  = rr_up(:, size(rr_up,2):1:-1)       ! LE -> TE

       rr_low = rr_dat(:,id_le:size(rr_dat,2))
    else
       ! Original DUST format: TE -> lower -> LE -> upper -> TE
       allocate(rr_up(2, size(rr_dat,2)-id_le+1));  rr_up = 0.0_wp
       allocate(rr_low(2, id_le));                  rr_low = 0.0_wp

       rr_up  = rr_dat(:,id_le:size(rr_dat,2))

       rr_low = rr_dat(:,1:id_le)
       rr_low = rr_low(:, size(rr_low,2):1:-1)     ! LE -> TE
    endif

    !> x sort elements 
    allocate(rr_sort(2, size(rr_dat,2))); rr_sort = 0.0_wp 

    call sort_vector_real(rr_dat(1,:), size(rr_dat,2), x_sort, ind_sort) 

    rr_sort(1,:) = x_sort
    rr_sort(2,:) = rr_dat(2,ind_sort) 

    !> remove duplicate elements 
    call unique(rr_sort(1,:), x_unique, 1e-6_wp)

    !> linear interpolation of upper and lower surfaces
    allocate(rr_unique_up(2, size(x_unique))); rr_unique_up = 0.0_wp
    allocate(rr_unique_low(2, size(x_unique))); rr_unique_low = 0.0_wp

    rr_unique_up(1,:) = x_unique
    rr_unique_low(1,:) = x_unique
    
    do i = 1, size(x_unique)
      call linear_interp(rr_up(2,:), rr_up(1,:),   rr_unique_up(1,i), rr_unique_up(2,i)) 
      call linear_interp(rr_low(2,:), rr_low(1,:), rr_unique_low(1,i), rr_unique_low(2,i)) 
    end do

    !> build middle line 
    allocate(rr_mean_dat(2, size(rr_unique_up,2))); rr_mean_dat = 0.0_wp
    rr_mean_dat(1,:) = rr_unique_low(1,:)
    rr_mean_dat(2,:) = (rr_unique_up(2,:) + rr_unique_low(2,:))/2.0_wp

    allocate(rr_mean(2, 50)); rr_mean = 0.0_wp 
    call linspace(0.0_wp, 1.0_wp, rr_mean(1,:)) 
    do i = 1, size(rr_mean, 2)
      call linear_interp(rr_mean_dat(2,:), rr_mean_dat(1,:), rr_mean(1,i), rr_mean(2,i))
    enddo  

    !> --- FIX TE SLOPES ---
    n = size(rr_unique_up, 2)
    tang_te_up = (rr_unique_up(2, n) - rr_unique_up(2, n-1)) / &
                 max(rr_unique_up(1, n) - rr_unique_up(1, n-1), 1e-6_wp)

    n = size(rr_unique_low, 2)
    tang_te_low = (rr_unique_low(2, n) - rr_unique_low(2, n-1)) / &
                  max(rr_unique_low(1, n) - rr_unique_low(1, n-1), 1e-6_wp)

    tang_te_up  = max(min(tang_te_up, max_tang_te), -max_tang_te)
    tang_te_low = max(min(tang_te_low, max_tang_te), -max_tang_te)

    !> --- FIX LE SLOPES ---
    tang_le_up = (rr_unique_up(2, 2) - rr_unique_up(2, 1)) / &
                 max(rr_unique_up(1, 2) - rr_unique_up(1, 1), 1e-6_wp)
    tang_le_low = (rr_unique_low(2, 2) - rr_unique_low(2, 1)) / &
                  max(rr_unique_low(1, 2) - rr_unique_low(1, 1), 1e-6_wp)

    tang_le_up  = max(min(tang_le_up, 10.0_wp), -10.0_wp)
    tang_le_low = max(min(tang_le_low, 10.0_wp), -10.0_wp)

    !> tangent at mean line
    tang_mean_le = (rr_mean(2,2)- rr_mean(2,1)) / (rr_mean(1,2)-rr_mean(1,1))
    tang_mean_te = (rr_mean(2, size(rr_mean,2))- rr_mean(2, size(rr_mean,2)-1)) / &
                    (rr_mean(1, size(rr_mean,2))-rr_mean(1, size(rr_mean,2)-1)) 
    
    allocate(point_up_updated(2, size(rr_unique_up,2) - 2)); point_up_updated = 0.0_wp
    allocate(point_low_updated(2, size(rr_unique_low,2) -2)); point_low_updated = 0.0_wp
    allocate(point_le_updated(2, 5)); point_le_updated = 0.0_wp
    point_up_updated = rr_unique_up(:, 3:size(rr_unique_up,2)) 
    point_low_updated = rr_unique_low(:, 3:size(rr_unique_low,2))
    
    !> switch xy coordinates for LE
    point_le_updated(1,3:5) = rr_unique_up(2,1:3) 
    point_le_updated(2,3:5) = rr_unique_up(1,1:3)
    point_le_updated(1,1)   = rr_unique_low(2,3)
    point_le_updated(2,1)   = rr_unique_low(1,3)
    point_le_updated(1,2)   = rr_unique_low(2,2)
    point_le_updated(2,2)   = rr_unique_low(1,2)
    
    !> upper part
    allocate(xq_up(num_points)); xq_up = 0.0_wp 
    call linspace(rr_unique_up(1,3), rr_unique_up(1,size(rr_unique_up,2)), xq_up) 
    allocate(yq_up(num_points)); yq_up = 0.0_wp
    call hermite_spline_profile(point_up_updated(1,:), point_up_updated(2,:), & 
                                xq_up, tang_le_up, tang_te_up, yq_up)
    
    !> lower part
    allocate(xq_low(num_points)); xq_low = 0.0_wp
    call linspace(rr_unique_low(1,3), rr_unique_low(1,size(rr_unique_low,2)), xq_low)
    allocate(yq_low(num_points)); yq_low = 0.0_wp 
    call hermite_spline_profile(point_low_updated(1,:), point_low_updated(2,:), & 
                                xq_low, tang_le_low, tang_te_low, yq_low)
    
    !> flip lower part
    yq_low = yq_low(size(yq_low):1:-1)
    xq_low = xq_low(size(xq_low):1:-1)

    !> leading edge part
    allocate(xq_le(num_points)); xq_le = 0.0_wp
    call linspace(point_le_updated(1,1), point_le_updated(1,5), xq_le)
    allocate(yq_le(num_points)); yq_le = 0.0_wp
    call hermite_spline_profile(point_le_updated(1,:), point_le_updated(2,:), & 
                                xq_le, 1.0_wp/tang_le_low, 1.0_wp/tang_le_up, yq_le)
    
    !> mean line
    allocate(xq_mean(num_points)); xq_mean = 0.0_wp
    call linspace(0.0_wp, 1.0_wp, xq_mean)
    allocate(yq_mean(num_points)); yq_mean = 0.0_wp
    call hermite_spline_profile(rr_mean(1,:), rr_mean(2,:), &
                                xq_mean, tang_mean_le, tang_mean_te, yq_mean)
    
    !> concatenate upper and lower part 
    allocate(rr2mesh(2, size(xq_up) + size(xq_low) + size(xq_le) - 2)); rr2mesh = 0.0_wp
    rr2mesh(1,:) = (/xq_low, yq_le(2:size(yq_le)-1), xq_up/)
    rr2mesh(2,:) = (/yq_low, xq_le(2:size(xq_le)-1), yq_up/) 
    
    min_x = minloc(rr2mesh(1,:),1)

    !> upper part
    allocate(rr2mesh_up(2, size(rr2mesh,2) - min_x + 1)); rr2mesh_up = 0.0_wp
    rr2mesh_up(1,:) = rr2mesh(1,min_x:size(rr2mesh,2))
    rr2mesh_up(2,:) = rr2mesh(2,min_x:size(rr2mesh,2)) 

    !> get length and derivative (alpha) along upper line
    allocate(alpha_up(size(rr2mesh_up, 2))); alpha_up = 0.0_wp
    allocate(dl_up(size(rr2mesh_up, 2)));    dl_up = 0.0_wp
    alpha_up(1) = atan2(rr2mesh_up(2,2) - rr2mesh_up(2,1), &
                        rr2mesh_up(1,2) - rr2mesh_up(1,1))
    
    do i = 2, size(rr2mesh_up, 2)
      dl_up(i) = dl_up(i - 1) + sqrt((rr2mesh_up(1,i) - rr2mesh_up(1,i-1))**2.0_wp + &
                                     (rr2mesh_up(2,i) - rr2mesh_up(2,i-1))**2.0_wp)
      alpha_up(i) = atan2(rr2mesh_up(2,i) - rr2mesh_up(2,i-1), &
                          rr2mesh_up(1,i) - rr2mesh_up(1,i-1))  
    end do
    
    !> lower part 
    allocate(rr2mesh_low(2, min_x)); rr2mesh_low = 0.0_wp
    
    rr2mesh_low(1,:) = rr2mesh(1,1:min_x)
    rr2mesh_low(2,:) = rr2mesh(2,1:min_x)
    
    !> flip lower part [0 to 1]
    rr2mesh_low(1,:) = rr2mesh_low(1,size(rr2mesh_low,2):1:-1)
    rr2mesh_low(2,:) = rr2mesh_low(2,size(rr2mesh_low,2):1:-1)
    
    allocate(alpha_low(size(rr2mesh_low,2))); alpha_low = 0.0_wp
    allocate(dl_low(size(rr2mesh_low,2)));    dl_low = 0.0_wp
    alpha_low(1) = atan2(rr2mesh_low(2,2) - rr2mesh_low(2, 1), &
                         rr2mesh_low(1,2) - rr2mesh_low(1, 1))
    
    do i = 2, size(rr2mesh_low,2)
      dl_low(i) = dl_low(i - 1) + sqrt((rr2mesh_low(1,i) - rr2mesh_low(1,i-1))**2.0_wp + &
                                       (rr2mesh_low(2,i) - rr2mesh_low(2,i-1))**2.0_wp)
      alpha_low(i) = atan2(rr2mesh_low(2,i) - rr2mesh_low(2,i-1), &
                           rr2mesh_low(1,i) - rr2mesh_low(1,i-1))
    end do
    
    !> mean line 
    allocate(rr2mesh_mean(2, size(xq_mean))); rr2mesh_mean = 0.0_wp
    rr2mesh_mean(1,:) = xq_mean
    rr2mesh_mean(2,:) = yq_mean
    
    allocate(alpha_mean(size(rr2mesh_mean,2))); alpha_mean = 0.0_wp
    allocate(dl_mean(size(rr2mesh_mean,2)));    dl_mean = 0.0_wp
    alpha_mean(1) = atan(tang_mean_le) 
    do i = 2, size(rr2mesh_mean(1,:))
      dl_mean(i) = dl_mean(i - 1) + sqrt((rr2mesh_mean(1,i) - rr2mesh_mean(1,i-1))**2.0_wp + &
                                         (rr2mesh_mean(2,i) - rr2mesh_mean(2,i-1))**2.0_wp)
      alpha_mean(i) = atan2(rr2mesh_mean(2,i) - rr2mesh_mean(2,i-1), &
                            rr2mesh_mean(1,i) - rr2mesh_mean(1,i-1))
    end do

    !> interpolation upper part 
    allocate(csi_up(size(csi_half))); csi_up = 0.0_wp
    csi_up = csi_half*dl_up(size(dl_up)) 
    allocate(rr_geo_up(2, size(csi_up))); rr_geo_up = 0.0_wp 
    
    do i = 1, size(csi_up) 
      do j = 1, size(dl_up) - 1 
        if ((csi_up(i) .ge. dl_up(j)) .and. (csi_up(i) .le. dl_up(j + 1))) then 
          rr_geo_up(1,i) = rr2mesh_up(1, j) + cos(alpha_up(j))*(csi_up(i) - dl_up(j)) 
          rr_geo_up(2,i) = rr2mesh_up(2, j) + sin(alpha_up(j))*(csi_up(i) - dl_up(j)) 
        endif   
      enddo
    enddo
    
    !> interpolation lower part 
    allocate(csi_low(size(csi_half))); csi_low = 0.0_wp
    csi_low = csi_half*dl_low(size(dl_low))
    allocate(rr_geo_low(2, size(csi_low))); rr_geo_low = 0.0_wp
    do i = 1, size(csi_low)
      do j = 1, size(dl_low) - 1
        if ((csi_low(i) .ge. dl_low(j)) .and. (csi_low(i) .le. dl_low(j + 1))) then
          rr_geo_low(1,i) = rr2mesh_low(1, j) + cos(alpha_low(j))*(csi_low(i) - dl_low(j))
          rr_geo_low(2,i) = rr2mesh_low(2, j) + sin(alpha_low(j))*(csi_low(i) - dl_low(j))
        endif
      enddo
    enddo

    !> flip lower part [1 to 0]
    rr_geo_low(1,:) = rr_geo_low(1,size(rr_geo_low,2):1:-1) 
    rr_geo_low(2,:) = rr_geo_low(2,size(rr_geo_low,2):1:-1)

    !> interpolation mean line
    allocate(csi_mean(size(csi_half))); csi_mean = 0.0_wp
    csi_mean = csi_half*dl_mean(size(dl_mean))
    allocate(rr_geo_mean(2, size(csi_mean))); rr_geo_mean = 0.0_wp
    
    do i = 1, size(csi_mean)
      do j = 1, size(dl_mean) - 1
        if ((csi_mean(i) .ge. dl_mean(j)) .and. (csi_mean(i) .le. dl_mean(j + 1))) then
          rr_geo_mean(1,i) = rr2mesh_mean(1, j) + cos(alpha_mean(j))*(csi_mean(i) - dl_mean(j))
          rr_geo_mean(2,i) = rr2mesh_mean(2, j) + sin(alpha_mean(j))*(csi_mean(i) - dl_mean(j))
        endif
      enddo
    enddo

    if (ElType .eq. 'p') then 
      allocate(rr(2,(size(rr_geo_up, 2) + size(rr_geo_low, 2) - 1))); rr = 0.0_wp 
      rr(1,:) = (/rr_geo_low(1,1:size(rr_geo_low,2)-1), rr_geo_up(1,1:size(rr_geo_up,2))/) 
      rr(2,:) = (/rr_geo_low(2,1:size(rr_geo_low,2)-1), rr_geo_up(2,1:size(rr_geo_up,2))/)
    elseif (ElType .eq. 'v') then 
      allocate(rr(2,size(rr_geo_mean, 2))); rr = 0.0_wp 
      rr(1,:) = rr_geo_mean(1,:)
      rr(2,:) = rr_geo_mean(2,:) 
    endif 
    thickness = maxval(rr(2,:)) - minval(rr(2,:))
  end if

  !> cleanup
  if (allocated(rr_dat))            deallocate(rr_dat) 
  if (allocated(rr_up))             deallocate(rr_up)
  if (allocated(rr_low))            deallocate(rr_low)
  if (allocated(rr_sort))           deallocate(rr_sort)
  if (allocated(rr_unique_up))      deallocate(rr_unique_up)
  if (allocated(rr_unique_low))     deallocate(rr_unique_low)
  if (allocated(rr_mean_dat))       deallocate(rr_mean_dat)
  if (allocated(rr_mean))           deallocate(rr_mean)
  if (allocated(point_up_updated))  deallocate(point_up_updated)
  if (allocated(point_low_updated)) deallocate(point_low_updated)
  if (allocated(point_le_updated))  deallocate(point_le_updated)
  if (allocated(xq_up))             deallocate(xq_up)
  if (allocated(yq_up))             deallocate(yq_up)
  if (allocated(xq_low))            deallocate(xq_low)
  if (allocated(yq_low))            deallocate(yq_low)
  if (allocated(xq_le))             deallocate(xq_le)
  if (allocated(yq_le))             deallocate(yq_le)
  if (allocated(xq_mean))           deallocate(xq_mean)
  if (allocated(yq_mean))           deallocate(yq_mean)
  if (allocated(rr2mesh))           deallocate(rr2mesh)
  if (allocated(rr2mesh_up))        deallocate(rr2mesh_up)
  if (allocated(rr2mesh_low))       deallocate(rr2mesh_low)
  if (allocated(rr2mesh_mean))      deallocate(rr2mesh_mean)
  if (allocated(alpha_up))          deallocate(alpha_up)
  if (allocated(alpha_low))         deallocate(alpha_low)
  if (allocated(alpha_mean))        deallocate(alpha_mean)
  if (allocated(dl_up))             deallocate(dl_up)
  if (allocated(dl_low))            deallocate(dl_low)
  if (allocated(dl_mean))           deallocate(dl_mean)
  if (allocated(csi_up))            deallocate(csi_up)
  if (allocated(csi_low))           deallocate(csi_low)
  if (allocated(csi_mean))          deallocate(csi_mean)
  if (allocated(rr_geo_up))         deallocate(rr_geo_up)
  if (allocated(rr_geo_low))        deallocate(rr_geo_low)
  if (allocated(rr_geo_mean))       deallocate(rr_geo_mean)
  
end subroutine read_airfoil

!-------------------------------------------------------------------------------
subroutine define_division(type_mesh, nelem, division, &
                          r, r_le, r_te, r_le_aft, r_te_aft, r_le_fore, r_te_fore, xh)

  integer, intent(in)                :: nelem
  real(wp), intent(in)               :: r, r_le, r_te, r_le_aft
  real(Wp), intent(in)               :: r_te_aft, r_le_fore, r_te_fore, xh
  character(len=*), intent(in)       :: type_mesh
  real(wp), allocatable, intent(out) :: division(:)

  real(wp),        allocatable       :: division_(:)
  real(wp)                           :: step
  integer                            :: iPoint
  character(len=*), parameter        :: this_sub_name = 'define_division'
  
  if(.not. allocated(division)) allocate(division(nelem+1));  division = 0.0_wp
  step = 1.0_wp/(real(nelem,wp) + 1e-12_wp)

  select case (trim(type_mesh))
  case ("uniform")
    do iPoint = 1,nelem+1
      division(iPoint) = (real(iPoint-1,wp))*step
    enddo
  case ("cosine")
    do iPoint = 1,nelem+1
      division(iPoint) = (1.0_wp - cos(pi*(real(iPoint-1,wp))*step))/2.0_wp
    enddo
  case ("cosine_le", "cosine_ib")
    do iPoint = 1,nelem+1
      division(iPoint) = 1.0_wp - cos(pi/2.0_wp*(real(iPoint-1,wp))*step)
    enddo
  case ("cosine_te", "cosine_ob")
    do iPoint = 1,nelem+1
      division(iPoint) = sin(pi/2.0_wp*((real(iPoint-1,wp))*step))
    enddo
  case ("geoseries_le", "geoseries_ib")
    call geoseries(0.0_wp, 1.0_wp, nelem, r, division, division_) 
  case ("geoseries_te", "geoseries_ob")
    call geoseries(0.0_wp, 1.0_wp, nelem, r, division_, division) 
  case ("geoseries")
    call geoseries_both(0.0_wp, 1.0_wp, nelem, xh, r_le, r_te, division)
  case ("geoseries_hi", "geoseries_as")
    call geoseries_hinge(0.0_wp, 1.0_wp, nelem, xh, r_le_aft, r_te_aft, r_le_fore, r_te_fore, division)
  case default
    call error(this_sub_name, this_mod_name, 'Incorrect input: &
        & type_chord must be equal to uniform, cosine, cosine_le, cosine_te, geoseries_le, geoseries_te, &
        & geoseries, geoseries_hi')
  end select

end subroutine define_division


!-------------------------------------------------------------------------------
subroutine define_thickness(rr_geo, thickness)
  real(wp), intent(out) :: thickness
  real(wp), intent(inout) :: rr_geo(:,:)

  real(wp), allocatable :: y_up_interp(:), y_low_interp(:), rr_geo_tmp(:,:)
  real(wp), allocatable :: x_coor(:), x_geo_shift(:) 
  real(wp), allocatable :: mean_line(:), curv_per(:)
  real(wp), allocatable :: curv(:), ds(:), m_up(:), m_low(:)
  real(wp), allocatable :: x_mid_point(:), y_mid_point(:)
  real(wp), allocatable :: vec_intersect_up(:,:), vec_intersect_low(:,:), thick_vec(:)
  real(wp), allocatable :: rr_up(:,:), rr_low(:,:)
  real(wp)              :: m_1, m_2, m_3
  real(wp)              :: mat_up(2,2), mat_low(2,2)
  real(wp)              :: vec_up(2), vec_low(2)
  integer               :: i, id_change, len_up, len_low

  allocate(rr_geo_tmp(size(rr_geo,1),2))
  rr_geo_tmp = rr_geo 
  call unique(rr_geo_tmp(1,:), x_coor, 1e-10_wp)
  allocate(y_up_interp(size(x_coor))); y_up_interp = 0.0_wp
  allocate(y_low_interp(size(x_coor))); y_low_interp = 0.0_wp

  !> split upper and lower part
  allocate(x_geo_shift(size(rr_geo,2)-1)); x_geo_shift = 0.0_wp 
  do i = 1,size(rr_geo,2) - 1
    x_geo_shift(i) = rr_geo(1,i+1) - rr_geo(1,i)
  end do 

  id_change = 0
  !> get id where the sign change 
  do i = 1, size(x_geo_shift) - 1
    if (x_geo_shift(i) .lt. 1e-10_wp .and. x_geo_shift(i+1) .gt. 1e-10_wp) then 
      id_change = i 
    endif
  enddo 

  len_low = id_change + 1 
  len_up = size(rr_geo,2)-id_change + 1

  allocate(rr_low(2,len_low)); rr_low = 0.0_wp
  allocate(rr_up(2,len_up)); rr_up = 0.0_wp
  
  rr_low = rr_geo(:,1:len_low)
  rr_up = rr_geo(:,len_low:size(rr_geo,2))
  
  do i = 2, size(x_coor) 
    call linear_interp(rr_up(2,:), rr_up(1,:), x_coor(i), y_up_interp(i))
    call linear_interp(rr_low(2,len_low:1:-1), rr_low(1,len_low:1:-1), x_coor(i), y_low_interp(i))
  enddo 
  
  !> evaluate mean line
  allocate(mean_line(size(y_up_interp))); mean_line = 0.0_wp
  mean_line = (y_up_interp + y_low_interp)/2
  !> length of each mean line segment
  allocate(ds(size(mean_line)-1));    ds = 0.0_wp 
  !> curvature of mean line, up and lower side
  allocate(curv(size(mean_line)-1));  curv = 0.0_wp 
  allocate(m_up(size(mean_line)-1));  m_up = 0.0_wp 
  allocate(m_low(size(mean_line)-1)); m_low = 0.0_wp 
  !> mid point 
  allocate(x_mid_point(size(mean_line)-1)); x_mid_point = 0.0_wp 
  allocate(y_mid_point(size(mean_line)-1)); y_mid_point = 0.0_wp 

  do i = 1, (size(mean_line)-1)
    ds(i) = x_coor(i+1)-x_coor(i)
    !> calculate the derivative 
    curv(i) = (mean_line(i+1) - mean_line(i))/ds(i)
    m_up(i) = (y_up_interp(i+1) - y_up_interp(i))/ds(i)
    m_low(i) = (y_low_interp(i+1) - y_low_interp(i))/ds(i)
    !> calculate mid point 
    y_mid_point(i) = (mean_line(i+1)+mean_line(i))/2
    x_mid_point(i) = (x_coor(i+1) + x_coor(i))/2
  end do

  !> take the perpendicular
  allocate(curv_per(size(mean_line)-1)); curv_per = 0.0_wp
  curv_per = -1/curv

  !> calculate the intersection between perpendicular lines and airfoil 
  allocate(vec_intersect_up(2,(size(mean_line)-1)));  vec_intersect_up = 0.0_wp
  allocate(vec_intersect_low(2,(size(mean_line)-1))); vec_intersect_low = 0.0_wp
  allocate(thick_vec(size(mean_line)-1)); thick_vec = 0.0_wp

  do i = 1, size(mean_line)-1
    m_1 = curv_per(i) 
    m_2 = m_up(i) 
    m_3 = m_low(i) 

    !> check infinity
    if (abs(m_1) .gt. 1e+16_wp) then
      if (abs(m_2) .gt. 1e+16_wp) then
          vec_intersect_up(1,i) = x_mid_point(i)
          vec_intersect_up(2,i) = y_mid_point(i)
      elseif (abs(m_3) .gt. 1e+16_wp) then
          vec_intersect_low(1,i) = x_mid_point(i)
          vec_intersect_low(2,i) = y_mid_point(i)
      else           
          vec_intersect_up(1,i) = x_mid_point(i)
          vec_intersect_up(2,i) = m_2*(x_mid_point(i) - x_coor(i)) + y_up_interp(i)
          vec_intersect_low(1,i) = x_mid_point(i)
          vec_intersect_low(2,i) = m_3*(x_mid_point(i) - x_coor(i)) + y_low_interp(i)
      end if
    else  
      mat_up(1,:) = 1/(-m_1 + m_2)*(/1.0_wp, -1.0_wp/)
      mat_up(2,:) = 1/(-m_1 + m_2)*(/m_2, -m_1/)

      vec_up(1) = y_mid_point(i) - m_1*x_mid_point(i)
      vec_up(2) = y_up_interp(i) - m_2*x_coor(i)

      vec_intersect_up(:,i) = matmul(mat_up, vec_up)

      mat_low(1,:) = 1/(-m_1 + m_3)*(/1.0_wp, -1.0_wp/)
      mat_low(2,:) = 1/(-m_1 + m_3)*(/m_3, -m_1/)

      vec_low(1) = y_mid_point(i) - m_1*x_mid_point(i)
      vec_low(2) = y_low_interp(i) - m_3*x_coor(i)

      vec_intersect_low(:,i) = matmul(mat_low, vec_low)
    end if
    !> calculate length of each segment 
    thick_vec(i) = sqrt((vec_intersect_low(1,i)-vec_intersect_up(1,i))**2 + &
                    (vec_intersect_low(2,i)-vec_intersect_up(2,i))**2)    
  end do

  !> finally the max value is the airfoil thickness 
  thickness = maxval(thick_vec) 
  

end subroutine define_thickness

!-------------------------------------------------------------------------------

subroutine read_actuatordisk_parametric(mesh_file,ee,rr)
  character(len=*), intent(in) :: mesh_file
  integer  , allocatable, intent(out) :: ee(:,:)
  real(wp) , allocatable, intent(out) :: rr(:,:)

  real(wp), allocatable :: x(:), y(:)
  type(t_parse) :: pmesh_prs
  character :: ElType
  integer :: nstep, ax, ip
  real(wp) :: r, theta
  integer :: ind1, ind2
  character(len=*), parameter ::  this_sub_name = 'read_actuatordisk_parametric' 

  call pmesh_prs%CreateStringOption('el_type', &
                'element type (temporary) p panel v vortex ring &
                & l lifting line a actuator disk')
  call pmesh_prs%CreateRealOption('radius', 'Radius of the actuator disk')
  call pmesh_prs%CreateIntOption('nstep','Number of subdivisions')
  call pmesh_prs%CreateIntOption('Axis','Which axis to align the disk')

  !read the parameters
  call pmesh_prs%read_options(trim(mesh_file),printout_val=.true.)

  ElType = getstr(pmesh_prs,'el_type')

  if(trim(ElType) .ne. 'a') call error(this_sub_name, this_mod_name, &
    'This should have not happened, a team of professionals is under way to &
    &remove the evidence')

  r = getreal(pmesh_prs,'radius')
  nstep = getint(pmesh_prs,'nstep')
  ax = getint(pmesh_prs,'Axis')


  allocate(x(nstep), y(nstep))

  do ip = 1, nstep
    theta = real(ip-1,wp) * 2.0_wp*pi/real(nstep,wp)
    x(ip) = cos(theta)*r
    y(ip) = sin(theta)*r
  enddo

  ind1 = 1+mod(ax,3)
  ind2 = 1+mod(ind1,3)
  allocate(rr   (3,nstep)) ; rr = 0.0_wp
  rr(ind1,:) = x
  rr(ind2,:) = y

  allocate(ee(nstep,1)) ; ee = 0
  do ip = 1,nstep
    ee(ip,1) = ip
  enddo

  deallocate(x,y)

end subroutine read_actuatordisk_parametric

!-------------------------------------------------------------------------------

end module mod_parametric_io
