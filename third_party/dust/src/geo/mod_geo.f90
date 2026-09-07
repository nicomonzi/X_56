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
!!          Alberto Savino
!!=========================================================================

!> Module to treat the geometry of the solid bodies

module mod_geometry

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

use mod_sim_param, only: &
  sim_param

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption &
  , finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc, &
  check_file_exists

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_cgns_io, only: &
  read_mesh_cgns

use mod_parametric_io, only: &
  read_mesh_parametric

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_stripe, only: &
  t_stripe 

use mod_liftlin, only: &
  t_liftlin, t_liftlin_p

use mod_actuatordisk, only: &
  t_actdisk

use mod_virtual, only: &
  c_elem_virtual, t_elem_virtual_p

use mod_c81, only: &
  t_aero_tab , read_c81_table , interp_aero_coeff

use mod_math, only: &
  cross

use mod_reference, only: &
  t_ref, build_references, update_all_references, destroy_references  

use mod_hinges, only: &
  t_hinge, initialize_hinge_config

use mod_hdf5_io, only: &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  write_hdf5_attr, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5

#if USE_PRECICE
use mod_precice_rbf, only: &
  t_precice_rbf
#endif

!----------------------------------------------------------------------

implicit none

public :: t_geo, t_geo_component, t_tedge, &
          create_geometry, update_geometry, destroy_geometry, &
          calc_geo_data_pan, calc_geo_data_ll, calc_geo_data_ad,&
          calc_node_vel , calc_geo_vel, calc_geo_vel_virtual, destroy_elements

private

!----------------------------------------------------------------------

!> Geometry component
!!
!! A geometry component is a single geometrical component of the model,
!! such as a wing, an aileron, a rotor blade etc.
!!
!! It is used to load the mesh and link it to the proper refernce frame.
!!
!! It also contains the vector of all the aerodynamic elements generated
!! by the grid
!!
!! It also contains some temporary arrays employed during geometry generation
type :: t_geo_component

  !> Element type
  character(len=max_char_len) :: comp_el_type

  !> Name of the component
  character(len=max_char_len) :: comp_name

  !> Type of input: basic, cgns, parametric
  character(len=max_char_len) :: comp_input

  !> Coupling w/ external structural code: .true., .false.
  logical :: coupling

  !> Coupling type: 'll', 'rigid', 'rbf'
  character(len=max_char_len) :: coupling_type

  !> Coupling node: local coordinates (in the component local ref.frame)
  ! of the coupling node (why initilized?)
  real(wp) :: coupling_node(3) = 0.0_wp

  !> Coupling node orientation: orientation of the coupling node, w.r.t.
  ! the local ref.rame of the geometrical component (why initialized?)
  real(wp) :: coupling_node_rot(3,3) = 0.0_wp

#if USE_PRECICE
  type(t_precice_rbf) :: rbf
#endif

  !> Id of the component (warning: not always defined)
  integer :: comp_id

  !> Number of elements of the component
  integer :: nelems

  !> Number of virtual elements of the component
  integer :: nelems_virtual

  !> Number of elements in chord
  integer :: n_c
  
  !> Number of virtual elements in chord
  integer :: n_c_virtual

  !> Number of elements in span
  integer :: n_s


  !> Elements of the group
  !! The main memory allocation happens here, they will be pointed from
  !! somewhere else
  class(c_pot_elem), allocatable     :: el(:)

  class(c_elem_virtual), allocatable :: el_virtual(:)

  !> Global indexes of the points in the component
  integer, allocatable :: i_points(:)

  !> Global indexes of the virtual points in the component
  integer, allocatable :: i_points_virtual(:)

#if USE_PRECICE
  !> Global PreCICE indices of the points for coupling
  integer, allocatable :: i_points_precice(:)



#endif

  !> Reference frame tag
  character(len=max_char_len) :: ref_tag

  ! Reference frame id
  integer :: ref_id

  !> Points in local reference frame
  real(wp), allocatable :: loc_points(:,:), loc_points_virtual(:,:)  
  real(wp), allocatable :: loc_cen(:,:), loc_cen_virtual(:,:)
  real(wp), allocatable :: loc_ctr_pt(:,:) 
  
  !> Number of surface panels in the component
  integer :: nSurfPan
  !> Number of vortex rings in the component
  integer :: nVortRin
  !> Number of lifting line elements in the component
  integer :: nLiftLin
  !> Number of lifting line ele
  integer :: nActDisk

  !> Is the component moving?
  logical :: moving

  !> Arrays for lifting line components
  character(len=max_char_len), allocatable :: airfoil_list(:)
  !> Number of ll elements in each region. Size(nelem_span_list) = nRegions
  integer, allocatable :: nelem_span_list(:)
  !> Normalised coordinate of sections in each wing region.
  !! Size(...) = nRegions+1
  real(wp),allocatable :: normalised_coord_e(:,:)
  !> twist angle for flat elements
  real(wp),allocatable :: theta_e(:)
  !> Id of the airfoil elements (index in airfoil_list char array)
  integer ,allocatable :: i_airfoil_e(:,:)
  character(len=5) :: aero_correction
  !> airfoil thickness for dynamic stall 
  real(wp), allocatable :: thickness(:,:)
  type(t_stripe), allocatable :: stripe(:)

#if USE_PRECICE
  !> Chord vector reference, required to assign the motion ot the
  ! TE of a deformable LL element, knowing the LE motion
  real(wp), allocatable :: c_ref_p(:,:)
  !> "Chord" vector reference, referred to the center of the elements
  real(wp), allocatable :: c_ref_c(:,:)
  !> x_ac offset for LL
  real(wp), allocatable :: xac(:)
#endif

  !> Hinges
  integer :: n_hinges
  type(t_hinge), allocatable :: hinge(:)


  !> Dimensions of parametric elements only
  integer :: parametric_nelems_span , parametric_nelems_chor, parametric_nelems_chor_virtual

end type  t_geo_component

!-----------------------------------


!> Geometry type
!!
!! The type contains the whole geometry, with all the components and reference
!! frames and all the points of the mesh
type :: t_geo

  !> Total number of implicit elements
  integer :: nelem_impl

  !> Number of explicit elements
  integer :: nelem_expl

  !> Number of virtual elements
  integer :: nelem_virtual

  !> Number of surface panel elements
  integer :: nSurfPan
  !> Number of vortex ring elements
  integer :: nVortLatt
  !> Number of lifting line elements
  integer :: nLiftLin
  !> Number of Actuator disk elements
  integer :: nActDisk


  ! Only for implicit elements (for slicing)
  !> Global id of surface panel elements
  integer, allocatable :: idSurfPan(:)
  !> Global id of vortex ring elements
  integer, allocatable :: idVortLatt(:)
  !> Surfpan id of global surface panel elements (G2L: global to local)
  integer, allocatable :: idSurfPanG2L(:)

  !> Number of statical implicit elements
  integer :: nstatic_impl
  !> Number of moving implicit elements
  integer :: nmoving_impl

  !> Number of statical surfpan elements (for slicing)
  integer :: nstatic_SurfPan

  !> Number of static or moving lifting lines
  integer :: nstatic_expl, nmoving_expl

  !> All the components of the geometry
  type(t_geo_component), allocatable :: components(:)

  !> Points (element vertexes)
  real(wp), allocatable :: points(:,:)

  !> Virtual points (element vertexes)
  real(wp), allocatable :: points_virtual(:,:)

#if USE_PRECICE
  !> Velocity of the points (element vertexes), to be used by
  ! coupled components
  real(wp), allocatable :: points_vel(:,:), points_vel_virtual(:,:)
#endif

  !> All the reference frames of the geometry
  type(t_ref), allocatable :: refs(:)

end type t_geo

!-----------------------------------

!> Trailing edge type
!!
!! The type contains the whole trailing edge info

type t_tedge

  !> Global id of the elements at the TE
  type(t_pot_elem_p), allocatable :: e(:,:)

  !> Global id of the elements at the TE belonging to the fountain 
  type(t_pot_elem_p), allocatable :: e_fountain(:,:) 

  !> Global id of the nodes on the TE
  integer , allocatable           :: i(:,:)

  !> TE id of the nodes of the TE elements
  integer , allocatable           :: ii(:,:)

  !> TE id of neighboring TE elements
  integer , allocatable           :: neigh(:,:)

  !> Relative orientation of the neighboring TE elements
  integer , allocatable           :: o(:,:)

  !> Unit vector at TE nodes
  real(wp), allocatable           :: t(:,:)
  real(wp), allocatable           :: t_hinged(:,:) ! considering hinge deflection
  integer, allocatable            :: is_hinged(:)

  !> Reference frame of the TE nodes
  integer , allocatable           :: ref(:)

  !> Component of the TE
  integer , allocatable           :: icomp(:)

  !> Individual scaling for each component
  real(wp), allocatable           :: scaling(:)

  !> Number of TE elements only for surface panels
  integer                         :: nte_surfpan 
  integer, allocatable            :: id_lifting_line(:) 
  
end type t_tedge


!----------------------------------------------------------------------

interface destroy_elements

  module procedure destroy_elements_geo
  module procedure destroy_elements_comps

end interface


!----------------------------------------------------------------------

character(len=*), parameter :: this_mod_name = 'mod_geo'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Create all the goemetry
!!
!! The geometry creation proceeds in this way:
!! -# The references frames are read from the file and built
!! -# The components are built by loading them from the geometry file
!! -# The number of elements of each kind are counted and the geometrical
!!    quantities of each element are calculated
!! -# The pointer arrays are sorted with the static elements before and the
!!    moving elements after, and in the total elements surface panels
!!    and vortex rings before, then lifting lines and finally actuator disks
subroutine create_geometry(geo_file_name, ref_file_name, in_file_name,  geo, &
                      te, elems_impl, elems_expl, elems_ad, elems_ll, elems_virtual, elems_corr, &
                      elems_non_corr, elems_tot, airfoil_data, target_file, run_id)
  character(len=*),    intent(in)                    :: geo_file_name
  character(len=*),    intent(inout)                 :: ref_file_name
  character(len=*),    intent(in)                    :: in_file_name
  type(t_geo),         intent(out), target           :: geo
  type(t_impl_elem_p), allocatable, intent(out)      :: elems_impl(:)
  type(t_expl_elem_p), allocatable, intent(out)      :: elems_expl(:)
  type(t_expl_elem_p), allocatable, intent(out)      :: elems_ad(:)
  type(t_liftlin_p),   allocatable, intent(out)      :: elems_ll(:)
  type(t_elem_virtual_p), allocatable, intent(out)   :: elems_virtual(:)
  type(t_pot_elem_p),  allocatable, intent(out)      :: elems_tot(:)
  type(t_pot_elem_p),  allocatable, intent(out)      :: elems_corr(:)
  type(t_pot_elem_p),  allocatable, intent(out)      :: elems_non_corr(:)
  
  type(t_tedge), intent(out)                         :: te
  type(t_aero_tab) , allocatable, intent(out)        :: airfoil_data(:)
  character(len=*), intent(in)                       :: target_file
  integer, intent(in)                                :: run_id(10)
  real(wp)                                           :: tstart
  real(wp), allocatable                              :: cen(:,:), cen_virtual(:,:)
  real(wp), allocatable                              :: ctr_pt(:,:)

  integer :: i, j, is, im,  i_comp, i_ll, i_ad
  integer :: i_non_corr, i_corr, i_tot, i_expl, i_impl, i_virtual
  type(t_impl_elem_p), allocatable :: temp_static_i(:), temp_moving_i(:)
  type(t_expl_elem_p), allocatable :: temp_static_e(:), temp_moving_e(:)
  integer(h5loc)                                     :: floc
  integer                                            :: iSurfPan, iVortLatt
  character(len=max_char_len)                        :: msg
  real(wp)                                           :: ac_ll(3)
  character(len=*), parameter                        :: this_sub_name='create_geometry'

  tstart = sim_param%t0

  !> Build the reference frames
  ! read which is the reference frame file, if default the main input file is
  ! employed
  if (trim(ref_file_name) .eq. 'no_set') then
    ref_file_name = trim(in_file_name)
  endif
  call check_file_exists(trim(ref_file_name), this_sub_name, this_mod_name)
  call build_references(geo%refs, trim(ref_file_name))

! -----------------------------------------------------------------------------------------
  !> Create the output geometry file (if it is not restarting with the same name)
  if(trim(geo_file_name).ne.trim(target_file)) then
    call new_hdf5_file(trim(target_file), floc)
    call write_hdf5_attr(run_id, 'run_id', floc)
    call close_hdf5_file(floc)
  endif

  !> Load the components from the file created by the preprocessor
  call check_preproc(geo_file_name)
  call load_components(geo, trim(geo_file_name), trim(target_file), te)

  call import_aero_tab(geo, airfoil_data)

  ! Initialisation
  geo%nelem_impl    = 0
  geo%nelem_expl    = 0
  geo%nelem_virtual = 0
  geo%nstatic_impl  = 0
  geo%nmoving_impl  = 0
  geo%nstatic_expl  = 0
  geo%nmoving_expl  = 0
  geo%nstatic_SurfPan = 0
  geo%nSurfPan   = 0
  geo%nVortLatt  = 0
  geo%nLiftLin   = 0
  geo%nActDisk   = 0

  ! count the elements
  do i_comp = 1,size(geo%components)

    select case(trim(geo%components(i_comp)%comp_el_type))
      case('p')
        if(geo%components(i_comp)%moving) then
          geo%nmoving_impl = geo%nmoving_impl + geo%components(i_comp)%nelems
        else
          geo%nstatic_impl = geo%nstatic_impl + geo%components(i_comp)%nelems
          geo%nstatic_SurfPan = geo%nstatic_SurfPan + geo%components(i_comp)%nelems
        endif
        geo%nSurfPan      = geo%nSurfPan + geo%components(i_comp)%nelems
        geo%nelem_impl    = geo%nelem_impl + geo%components(i_comp)%nelems
        geo%nelem_virtual = geo%nelem_virtual + geo%components(i_comp)%nelems_virtual

      case('v')
        if(geo%components(i_comp)%moving) then
          geo%nmoving_impl = geo%nmoving_impl + geo%components(i_comp)%nelems
        else
          geo%nstatic_impl = geo%nstatic_impl + geo%components(i_comp)%nelems
        endif
        geo%nVortLatt = geo%nVortLatt + geo%components(i_comp)%nelems
        geo%nelem_impl = geo%nelem_impl + geo%components(i_comp)%nelems
        geo%nelem_virtual = geo%nelem_virtual + geo%components(i_comp)%nelems_virtual

      case('l')
        if(geo%components(i_comp)%moving) then
          geo%nmoving_expl = geo%nmoving_expl + geo%components(i_comp)%nelems
        else
          geo%nstatic_expl = geo%nstatic_expl + geo%components(i_comp)%nelems
        endif
        geo%nLiftLin = geo%nLiftLin + geo%components(i_comp)%nelems
        geo%nelem_expl = geo%nelem_expl + geo%components(i_comp)%nelems
        geo%nelem_virtual = geo%nelem_virtual + geo%components(i_comp)%nelems_virtual

      case('a')
        if(geo%components(i_comp)%moving) then
          geo%nmoving_expl = geo%nmoving_expl + geo%components(i_comp)%nelems
        else
          geo%nstatic_expl = geo%nstatic_expl + geo%components(i_comp)%nelems
        endif
        geo%nActDisk = geo%nActDisk + geo%components(i_comp)%nelems
        geo%nelem_expl = geo%nelem_expl + geo%components(i_comp)%nelems

      case default
        call error (this_sub_name, this_mod_name, 'Unknow type of element: '&
                  //trim(geo%components(i_comp)%comp_el_type))
    end select

  enddo

  !> Calculate the geometric quantities
  !  already update the geometry for the first time to get the right
  !  starting geometrical condition
  call prepare_geometry(geo)
  call update_geometry(geo, te, tstart, update_static=.true., time_cycle = .false. )

  if(sim_param%debug_level .ge. 3) then
    call printout(nl//' Geometry details:' )
    write(msg,'(A,I9)') '  number of total elements:    ', &
                                                geo%nelem_impl + geo%nelem_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of implicit elements: ' ,geo%nelem_impl
    call printout(msg)
    write(msg,'(A,I9)') '  number of explicit elements: ' ,geo%nelem_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of static elements:   ' , &
                                            geo%nstatic_impl + geo%nstatic_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of moving elements:   ' , &
                                            geo%nmoving_impl + geo%nmoving_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of surface panels:    ' ,geo%nSurfpan
    call printout(msg)
    write(msg,'(A,I9)') '  number of vortex lattices:   ' ,geo%nVortLatt
    call printout(msg)
    write(msg,'(A,I9)') '  number of lifting lines:     ' ,geo%nLiftLin
    call printout(msg)
    write(msg,'(A,I9)') '  number of actuator disks:    ' ,geo%nActDisk
    call printout(msg)

  endif

  !Create the vector of pointers to all the elements
  allocate(elems_impl(geo%nelem_impl), elems_expl(geo%nelem_expl),  &
            elems_tot(geo%nelem_impl+geo%nelem_expl))
  allocate(elems_ad(geo%nactdisk), elems_ll(geo%nLiftLin))
  allocate(elems_virtual(geo%nelem_virtual))

  i_impl=0; i_expl=0; i_ad=0; i_ll=0; i_non_corr = 0; i_corr = 0;  i_tot=0; i_virtual=0;

  !> First count the corrected elements 
  do i_comp = 1,size(geo%components)
    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'v' .and. & 
      trim(geo%components(i_comp)%aero_correction) .eq. 'true') then
      do j = 1,size(geo%components(i_comp)%el)
        i_corr = i_corr + 1
      end do 
    end if 
  end do 
  allocate(elems_corr(i_corr)) 
  allocate(elems_non_corr(geo%nelem_impl+geo%nelem_expl - i_corr))
  i_corr = 0 

  do i_comp = 1,size(geo%components)
    
    if (trim(geo%components(i_comp)%comp_el_type) .ne. 'a') then ! virtual elements
      do j = 1,size(geo%components(i_comp)%el_virtual)
        i_virtual = i_virtual + 1
        elems_virtual(i_virtual)%p => geo%components(i_comp)%el_virtual(j) 
      enddo
    end if

    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p') then ! panels, implicit and not corrected
      do j = 1,size(geo%components(i_comp)%el)
        i_tot = i_tot+1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        i_impl = i_impl+1
        select type(el => geo%components(i_comp)%el(j)); class is(c_impl_elem)
          elems_impl(i_impl)%p => el
        end select
        i_non_corr = i_non_corr + 1 
        elems_non_corr(i_non_corr)%p => geo%components(i_comp)%el(j)
      enddo
    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then ! vortex lattice, implicit, can be corrected or not

      do j = 1,size(geo%components(i_comp)%el)
        i_tot = i_tot+1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        i_impl = i_impl+1
        select type(el => geo%components(i_comp)%el(j)); class is(c_impl_elem)
          elems_impl(i_impl)%p => el
        end select
        if (trim(geo%components(i_comp)%aero_correction) .eq. 'false') then ! not corrected
          i_non_corr = i_non_corr + 1 
          elems_non_corr(i_non_corr)%p => geo%components(i_comp)%el(j)
        elseif (trim(geo%components(i_comp)%aero_correction) .eq. 'true') then ! corrected
          i_corr = i_corr + 1 
          elems_corr(i_corr)%p => geo%components(i_comp)%el(j) 
        end if
      end do
    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then ! lifting line, explicit, not corrected
      do j = 1,size(geo%components(i_comp)%el)
        i_tot = i_tot + 1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        i_expl = i_expl + 1
        i_ll = i_ll + 1
        select type(el => geo%components(i_comp)%el(j)); class is(t_liftlin)
          elems_ll(i_ll)%p => el
          elems_expl(i_expl)%p => el
        end select
        i_non_corr = i_non_corr + 1
        elems_non_corr(i_non_corr)%p => geo%components(i_comp)%el(j)
      enddo
    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'a') then ! actuator disk, explicit, not corrected
      do j = 1,size(geo%components(i_comp)%el)
        i_tot = i_tot + 1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        i_expl = i_expl + 1
        i_ad = i_ad + 1
        select type(el => geo%components(i_comp)%el(j)); class is(c_expl_elem)
          elems_ad(i_ad)%p => el
          elems_expl(i_expl)%p => el
        end select
        i_non_corr = i_non_corr + 1
        elems_non_corr(i_non_corr)%p => geo%components(i_comp)%el(j)
      enddo
    end if
  enddo

  !> Sort implicit elements: first static, then moving 
  !  fill in the two temporaries
  allocate(temp_static_i(geo%nstatic_impl), temp_moving_i(geo%nmoving_impl))
  is = 0; im = 0;
  do i = 1,geo%nelem_impl
    if(elems_impl(i)%p%moving) then
      im = im+1
      temp_moving_i(im) = elems_impl(i)
    else
      is = is+1
      temp_static_i(is) = elems_impl(i)
    endif
  enddo

  !Now might be more bombproof to deallocate and allocate, but for the moment..
  elems_impl(1:geo%nstatic_impl) = temp_static_i
  elems_impl(geo%nstatic_impl+1:geo%nelem_impl) = temp_moving_i

  deallocate(temp_static_i, temp_moving_i)

  !Update the indexing since we re-ordered the vector
  do i = 1,geo%nelem_impl
    elems_impl(i)%p%id = i
  end do

  ! Save global indices of implicit elements for slicing
  allocate(geo%idSurfPan(geo%nSurfPan),geo%idVortLatt(geo%nVortLatt))
  allocate(geo%idSurfPanG2L(geo%nelem_impl)) ; geo%idSurfPanG2L = 0
  iSurfPan = 0 ; iVortLatt = 0
  do i = 1,geo%nelem_impl
    select type( el => elems_impl(i)%p ); class is(t_surfpan)
      iSurfPan = iSurfPan + 1
      geo%idSurfPan(iSurfPan) = el%id
      geo%idSurfPanG2L(i    ) = iSurfPan
    end select
    select type( el => elems_impl(i)%p ); class is(t_vortlatt)
      iVortLatt= iVortLatt + 1
      geo%idVortLatt(iVortLatt) = el%id ! get from here the vortex lattice id
    end select
  end do

  !> Now re-order the explicit elements
  allocate(temp_static_e(geo%nstatic_expl), temp_moving_e(geo%nmoving_expl))
  is = 0; im = 0;
  do i = 1,geo%nelem_expl
    if(elems_expl(i)%p%moving) then
      im = im+1
      temp_moving_e(im) = elems_expl(i)
    else
      is = is+1
      temp_static_e(is) = elems_expl(i)
    endif
  enddo
  if(geo%nelem_expl .gt. 0) then
    elems_expl(1:geo%nstatic_expl) = temp_static_e
    elems_expl(geo%nstatic_expl+1:geo%nelem_expl) = temp_moving_e
  endif
  deallocate(temp_static_e, temp_moving_e)

  !> Update the indexing since we re-ordered the vector
  do i = 1,geo%nelem_expl
    elems_expl(i)%p%id = i+geo%nelem_impl
  end do

  call create_strip_connectivity(geo)

  !> Initialisation of the field surf_vel to zero (for surfpan only)
  do i=1, geo%nelem_impl
    select type(el=>elems_impl(i)%p) ; class is(t_surfpan)
      el%surf_vel = 0.0_wp
#if USE_PRECICE
    if (geo%components(el%comp_id)%coupling) then
      call el%create_local_velocity_stencil( &    
              geo%components(el%comp_id)%coupling_node_rot)
      !> chtls stencil
      call el%create_chtls_stencil( &             
              geo%components(el%comp_id)%coupling_node_rot)
    else 
      call el%create_local_velocity_stencil( &    
                geo%refs(geo%components(el%comp_id)%ref_id)%R_g)
      !> chtls stencil
      call el%create_chtls_stencil( &             
              geo%refs(geo%components(el%comp_id)%ref_id)%R_g)
    endif 
#else 
    call el%create_local_velocity_stencil( &    
            geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
    !> chtls stencil
    call el%create_chtls_stencil( &             
            geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
#endif 
    end select
  end do

  !> Update connectivity matrix for cen and ac_stripe for  
#if USE_PRECICE

  do i_comp = 1, size(geo%components)
    if ((geo%components(i_comp)%coupling) .and. (trim(geo%components(i_comp)%coupling_type) .eq. 'rbf')) then 
      !> build connectivity for the panel center
      !> CAUTION: ll%cen is at 50% chord, but forces must be applied at 25% (cfr. mod_post_integral l. 297)
      allocate(cen(3, size(geo%components(i_comp)%el))); cen = 0.0_wp 

      if (geo%components(i_comp)%comp_el_type .eq. 'l') then 
        do i = 1, size(geo%components(i_comp)%el)
          !> el(i)%cen could be in a dust-defined reference frame, convert it back to base reference
          ! to interface with coupling
          ac_ll = sum ( geo%components(i_comp)%el(i)%ver(:,1:2),2 ) / 2.0_wp
          cen(:,i) = matmul(transpose(geo%refs(geo%components(i_comp)%ref_id)%R_g),&
            ac_ll -geo%refs(geo%components(i_comp)%ref_id)%of_g) 
        end do 
      else
        do i = 1, size(geo%components(i_comp)%el)
          !> el(i)%cen could be in a dust-defined reference frame, convert it back to base reference
          ! to interface with coupling
          cen(:,i) = matmul(transpose(geo%refs(geo%components(i_comp)%ref_id)%R_g),&
            geo%components(i_comp)%el(i)%cen -geo%refs(geo%components(i_comp)%ref_id)%of_g) 
        end do 
      endif
      
      call geo%components(i_comp)%rbf%build_connectivity(cen, geo%components(i_comp)%coupling_node_rot)
      !> transfer index and weight matrix 
      geo%components(i_comp)%rbf%cen%ind = geo%components(i_comp)%rbf%point%ind
      geo%components(i_comp)%rbf%cen%wei = geo%components(i_comp)%rbf%point%wei
      
      !> store the read points into the local cen  
      allocate(geo%components(i_comp)%loc_cen(3,size(cen,2)))
      geo%components(i_comp)%loc_cen = cen 

      !> cleanup 
      deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)
      deallocate(cen)

      !> build connectivity for the virtual panel center 
      allocate(cen_virtual(3, size(geo%components(i_comp)%el_virtual))); cen_virtual = 0.0_wp 

      do i = 1, size(geo%components(i_comp)%el_virtual)
        !> el(i)%cen could be in a dust-defined reference frame, convert it back to base reference
        ! to interface with coupling
        cen_virtual(:,i) = matmul(transpose(geo%refs(geo%components(i_comp)%ref_id)%R_g),&
          geo%components(i_comp)%el_virtual(i)%cen -geo%refs(geo%components(i_comp)%ref_id)%of_g) 
      end do 
      
      call geo%components(i_comp)%rbf%build_connectivity(cen_virtual, geo%components(i_comp)%coupling_node_rot)
      !> transfer index and weight matrix 
      geo%components(i_comp)%rbf%cen_virtual%ind = geo%components(i_comp)%rbf%point%ind
      geo%components(i_comp)%rbf%cen_virtual%wei = geo%components(i_comp)%rbf%point%wei
      
      !> store the read points into the local cen  
      allocate(geo%components(i_comp)%loc_cen_virtual(3,size(cen_virtual,2)))
      geo%components(i_comp)%loc_cen_virtual = cen_virtual 

      !> cleanup 
      deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)
      deallocate(cen_virtual)
      
      !> build connectivity for the stripe (vl_corrected) for velocity evaluation 
      if (trim(geo%components(i_comp)%comp_el_type) .eq. 'v' .and. &
          trim(geo%components(i_comp)%aero_correction) .eq. 'true') then 

        allocate(ctr_pt(3, size(geo%components(i_comp)%stripe))); ctr_pt = 0.0_wp 
          
        do i = 1, size(geo%components(i_comp)%stripe) 
          
          ctr_pt(:,i) = geo%components(i_comp)%stripe(i)%ctr_pt
          !> store the read points into the local control point          
          geo%components(i_comp)%stripe(i)%loc_ctr_pt = ctr_pt(:,i)

        end do       

        call geo%components(i_comp)%rbf%build_connectivity(ctr_pt, geo%components(i_comp)%coupling_node_rot)
        !> transfer index and weight matrix 
        geo%components(i_comp)%rbf%ctr_pt%ind = geo%components(i_comp)%rbf%point%ind
        geo%components(i_comp)%rbf%ctr_pt%wei = geo%components(i_comp)%rbf%point%wei 
        
        !> cleanup 
        deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)
        deallocate(ctr_pt) 
      
      !> build connectivity for lifting line control points for velocity evaluation 
      elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then 

          allocate(ctr_pt(3, size(geo%components(i_comp)%el))); ctr_pt = 0.0_wp 

          do i = 1, size(geo%components(i_comp)%el)
            select type(el => geo%components(i_comp)%el(i))
              type is(t_liftlin) 
              ctr_pt(:,i) = el%cen + el%tang_cen * el%chord / 2.0_wp ! control point of lifting line 
              !> store the read points into the local control point  
              el%loc_ctr_pt = ctr_pt(:,i)
            end select    
          end do       

          call geo%components(i_comp)%rbf%build_connectivity(ctr_pt, geo%components(i_comp)%coupling_node_rot)
          !> transfer index and weight matrix 
          geo%components(i_comp)%rbf%ctr_pt%ind = geo%components(i_comp)%rbf%point%ind
          geo%components(i_comp)%rbf%ctr_pt%wei = geo%components(i_comp)%rbf%point%wei
    
          !> cleanup 
          deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)
          deallocate(ctr_pt)
      endif 

    endif 
  end do 
#endif 

end subroutine create_geometry

!----------------------------------------------------------------------

!> Load the components mesh
! 
!  Here all the components from the geometry file are loaded.
!  The components linked to multiple reference frames are multiplied and
!  attached to their relative reference frames.
!  The points and the elements are then genereated.
!  Finally the trailing edge is created
subroutine load_components(geo, in_file, out_file, te)
  type(t_geo), intent(inout), target    :: geo
  character(len=*), intent(in)          :: in_file
  character(len=*), intent(in)          :: out_file
  type(t_tedge), intent(out)            :: te
  integer                               :: i1, i2, i3
  integer, allocatable                  :: ee(:,:), ee_virtual(:,:)
  real(wp), allocatable                 :: rr(:,:), rr_virtual(:,:)
  character(len=max_char_len)           :: comp_el_type, comp_name, comp_input
  integer                               :: points_offset, points_offset_virtual, n_vert 
  integer                               :: elems_offset, elems_offset_virtual
  real(wp), allocatable                 :: points_tmp(:,:), points_virtual_tmp(:,:)
  character(len=max_char_len)           :: ref_tag, ref_tag_m
  integer                               :: ref_id, iref
  character(len=max_char_len)           :: msg, cname, cname_write
  integer(h5loc)                        :: floc, gloc, cloc , geo_loc , te_loc, cloc2
  integer(h5loc)                        :: floc_out, gloc_out
  integer                               :: n_comp, i_comp 
  integer                               :: n_comp_input, i_comp_input, n_comp_write
  integer                               :: n_mult, i_mult
  logical                               :: mult
  !> Connectivity and te structures
  integer , allocatable                 :: neigh(:,:)
  !> Aerotable correction for vl 
  character(len=5)                      :: aero_table
  real(wp), allocatable                 :: thickness(:,:)
  !> Lifting Line elements
  real(wp), allocatable                 :: normalised_coord_e(:,:), theta_e(:)
  integer                 , allocatable :: i_airfoil_e(:,:)
  character(max_char_len) , allocatable :: airfoil_list(:)
  integer                 , allocatable :: nelem_span_list(:)
  !> Coupling
  character(len=max_char_len)           :: comp_coupling_str
  character(len=max_char_len)           :: comp_coupling_type
  real(wp) :: comp_coupling_node(3)       = 0.0_wp    ! <- Initialization(?)
  real(wp) :: comp_coupling_node_rot(3,3) = 0.0_wp    ! <- Initialization(?)
  logical :: comp_coupling

#if USE_PRECICE
  real(wp), allocatable                 :: c_ref_p(:,:)
  real(wp), allocatable                 :: c_ref_c(:,:)
  real(wp), allocatable                 :: comp_coupling_nodes(:,:)
  integer                               :: points_offset_precice, np_precice, np_precice_tot
  integer                               :: comp_ind
  integer, allocatable                  :: hinge_ind(:)
  integer, allocatable                  :: ind_coupling(:)
  integer                               :: n_nodes_coupling_hinges
  integer                               :: n_nodes_coupling_hinge_1
#endif

  real(wp)                              :: coupling_node_rot(3,3) = 0.0_wp
  !> Hinges
  integer                               :: n_hinges, ih
  character(len=2)                      :: hinge_id_str
  integer(h5loc)                        :: hloc, hiloc, hloc2
  !> Parametric elements
  integer                               :: par_nelems_span , par_nelems_chor, par_nelems_chor_virtual
  !> Trailing edge 
  integer , allocatable                 :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
  integer , allocatable                 :: neigh_te(:,:) , o_te(:,:), e_te_fountain(:,:)
  real(wp), allocatable                 :: t_te(:,:)
  real(wp)                              :: scale_te
  integer                               :: ne_te , nn_te, n_e_te_panel 
  !> Tmp arrays
  type(t_pot_elem_p) , allocatable      :: e_te_tmp(:,:), e_te_fountain_tmp(:,:) 
  integer, allocatable                  :: i_te_tmp(:,:), ii_te_tmp(:,:)
  integer , allocatable                 :: neigh_te_tmp(:,:), o_te_tmp(:,:)
  real(wp), allocatable                 :: t_te_tmp(:,:), t_hinged_te_tmp(:,:)
  integer , allocatable                 :: t_is_hinged_tmp(:)
  integer , allocatable                 :: ref_te_tmp(:)
  integer , allocatable                 :: icomp_te_tmp(:)
  real(wp), allocatable                 :: scale_te_tmp(:)
  !> # n. elements and nodes at TE ( of the prev. comps)
  integer                               :: ne_te_prev , nn_te_prev
  real(wp)                              :: trac, rad
  logical                               :: rewrite_geo
  logical                               :: is_hinge
  integer, allocatable                  :: ind_h(:)
  character(len=*), parameter           :: this_sub_name = 'load_components'

  !>Re-write the geometry only if it is not restarting from a previous run
  ! with the same name
  rewrite_geo = .true.
  if(trim(in_file).eq.trim(out_file)) rewrite_geo = .false.

  call printout('Reading geometry from file "'//trim(in_file)//'"')
  call open_hdf5_file(trim(in_file),floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp,'NComponents',gloc)

  !> First it is necessary to count how many actual components will
  !  be present, including the multiples. Due to pointers pointing to
  !  allocatables the components array cannot be expanded later
  n_comp_input = n_comp

  do i_comp_input = 1, n_comp_input
    write(cname,'(A,I3.3)') 'Comp',i_comp_input
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(ref_tag,'RefTag',cloc)
    !Look for the reference frame of the component
    ref_id = -1
    do iref = 0,ubound(geo%refs,1)
      if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag)) then
        !set id
        ref_id = iref
      endif
    enddo
    !if not found the reference
    if (ref_id .lt. 0) then
      write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                    ' a reference with tag ',trim(ref_tag),' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif
    mult = .false.
    n_mult = 1
    mult = geo%refs(ref_id)%multiple
    if(mult) then
      n_mult = geo%refs(ref_id)%n_mult
      n_comp = n_comp+n_mult-1 !(one was already counted)
    endif
    call close_hdf5_group(cloc)
  enddo

  !> Allocate the components of the right full size
  allocate(geo%components(n_comp))

  !> Open and prepare the output file (if not restarting with the same name)
  if(rewrite_geo) then
    call open_hdf5_file(trim(out_file),floc_out)
    call new_hdf5_group(floc_out,'Components',gloc_out)
    call write_hdf5(n_comp,'NComponents',gloc_out)
  endif

  elems_offset = 0
  elems_offset_virtual = 0
  n_comp_write = n_comp

#if USE_PRECICE
  points_offset_precice = 0
#endif

  i_comp = 1
  do i_comp_input = 1,n_comp_input

    write(cname,'(A,I3.3)') 'Comp',i_comp_input
    call open_hdf5_group(gloc,trim(cname),cloc)


    ! ====== REFERENCES ======

    call read_hdf5(ref_tag,'RefTag',cloc)

    !> Look for the reference frame of the component
    ref_id = -1
    do iref = 0,ubound(geo%refs,1)
      if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag)) then
        ! set id
        ref_id = iref
      endif
    enddo
    !> If not found the reference
    if (ref_id .lt. 0) then
      write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                    ' a reference with tag ',trim(ref_tag),' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif

    !> Check if it is a multiple reference
    mult = .false.
    n_mult = 1
    mult = geo%refs(ref_id)%multiple

    ! ====== READING =====
    call read_hdf5(comp_el_type,'ElType',cloc)
    call read_hdf5(comp_name,'CompName',cloc)
    call read_hdf5(comp_input,'CompInput',cloc)

    !> Coupling
    call read_hdf5(comp_coupling_str     ,'Coupled',cloc)
    call read_hdf5(comp_coupling_type    ,'CouplingType',cloc)
    call read_hdf5(comp_coupling_node    ,'CouplingNode',cloc)
    call read_hdf5(comp_coupling_node_rot,'CouplingNodeOrientation',cloc)
    comp_coupling = .false.
#if USE_PRECICE
    if ( trim(comp_coupling_str) .eq. 'true' ) comp_coupling = .true.
#else
    if ( trim(comp_coupling_str) .eq. 'true' ) then
      call error (this_sub_name, this_mod_name, &
          ' Coupled = T for component'//trim(comp_name)// &
          ', but no coupled simulation is expected, since dust &
          &has been compiled without #USE_PRECICE option. Stop'//nl)
    endif
#endif

    !> Some errors (or to do implementation?) for multiple components
    !> Multiple component and Coupling
    if ( mult .and. comp_coupling ) then
      call error (this_sub_name, this_mod_name, &
          ' Coupled = T for "Muliple" component'//trim(comp_name)// &
          ', but this is not allowed (at least so far). Stop'//nl)
    end if

    !> Hinges
    call open_hdf5_group(cloc, 'Hinges', hloc)
    call read_hdf5(n_hinges, 'n_hinges', hloc)
    !> Some errors (or to do implementation?) for multiple components
    !> Multiple component and Hinges
    if ( mult .and. ( n_hinges .gt. 0 ) ) then
      call error (this_sub_name, this_mod_name, &
          ' n_hinges .gt. 0 for "Multiple" component'//trim(comp_name)// &
          ', but this is not allowed (so far, at least). Stop'//nl)
    end if

    if(mult) then
      n_mult = geo%refs(ref_id)%n_mult
    endif

    !> Cycle on all the possible multiple references, if it is not multiple
    !  it will be passed just once
    !  TODO: consider wrapping in another subroutine
    do i_mult = 1,n_mult
      ref_tag_m = ref_tag
      !> set the component id
      geo%components(i_comp)%comp_id = i_comp
      if(mult)  then
        !> look again for the multiple reference
        write(ref_tag_m,'(A,I2.2)') trim(ref_tag)//'__',i_mult
        ref_id = -1
        do iref = 0,ubound(geo%refs,1)
          if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag_m)) then
            !> set id
            ref_id = iref
          endif
        enddo
        !> If not found the reference
        if (ref_id .lt. 0) then
          write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                      ' a reference with tag ',trim(ref_tag_m),' was not found'
          call error(this_sub_name, this_mod_name, msg)
        endif
      endif

      geo%components(i_comp)%ref_id  = ref_id
      geo%components(i_comp)%ref_tag = trim(ref_tag_m)
      geo%components(i_comp)%moving  = geo%refs(ref_id)%moving

      geo%components(i_comp)%coupling          =      comp_coupling
      geo%components(i_comp)%coupling_type     = trim(comp_coupling_type)
      geo%components(i_comp)%coupling_node     =      comp_coupling_node
      geo%components(i_comp)%coupling_node_rot =      comp_coupling_node_rot

#if USE_PRECICE
      coupling_node_rot = geo%components(i_comp)%coupling_node_rot
#else
      coupling_node_rot(1,:) = (/ 1._wp, 0._wp, 0._wp/)
      coupling_node_rot(2,:) = (/ 0._wp, 1._wp, 0._wp/)
      coupling_node_rot(3,:) = (/ 0._wp, 0._wp, 1._wp/)
#endif

      !> Overwrite moving for coupled components
      if ( comp_coupling )  geo%components(i_comp)%moving = .true.

      ! ====== READING =====
      geo%components(i_comp)%comp_el_type = trim(comp_el_type)
      geo%components(i_comp)%comp_input   = trim(comp_input  )

      if(mult) then
        write(geo%components(i_comp)%comp_name,'(A,I2.2)') trim(comp_name)&
                                                            &//'__',i_mult
      else
        geo%components(i_comp)%comp_name = trim(comp_name)
      endif

      ! Geometry --------------------------
      call open_hdf5_group(cloc,'Geometry',geo_loc)
      call read_hdf5_al(ee   ,'ee'   ,geo_loc)
      call read_hdf5_al(rr   ,'rr'   ,geo_loc)
      call read_hdf5_al(ee_virtual   ,'ee_virtual'   ,geo_loc)
      call read_hdf5_al(rr_virtual   ,'rr_virtual'   ,geo_loc)
      call read_hdf5_al(neigh,'neigh',geo_loc)

      !> element-specific reads
      if ( comp_el_type(1:1) .eq. 'l' ) then

        call read_hdf5_al(airfoil_list      ,'airfoil_list'      ,geo_loc)
        call read_hdf5_al(nelem_span_list   ,'nelem_span_list'   ,geo_loc)
        call read_hdf5_al(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc)
        call read_hdf5_al(normalised_coord_e,'normalised_coord_e',geo_loc)
        call read_hdf5_al(theta_e,           'theta_e'           ,geo_loc)

        allocate(geo%components(i_comp)%airfoil_list(size(airfoil_list)))
        geo%components(i_comp)%airfoil_list = airfoil_list

        allocate(geo%components(i_comp)%nelem_span_list(size(nelem_span_list)))
        geo%components(i_comp)%nelem_span_list = nelem_span_list

        allocate(geo%components(i_comp)%i_airfoil_e( &
              size(i_airfoil_e,1),size(i_airfoil_e,2)) )
        geo%components(i_comp)%i_airfoil_e = i_airfoil_e

        allocate(geo%components(i_comp)%normalised_coord_e( &
              size(normalised_coord_e,1),size(normalised_coord_e,2)))
        geo%components(i_comp)%normalised_coord_e = normalised_coord_e

        allocate(geo%components(i_comp)%theta_e(size(theta_e)))
        geo%components(i_comp)%theta_e = theta_e

      elseif (comp_el_type(1:1) .eq. 'v') then
        
        !> Additional fields for vl corrections 
        call read_hdf5(aero_table,  'aero_table',  geo_loc)

        if (trim(aero_table) .eq. 'true') then 

          call read_hdf5_al(airfoil_list      ,'airfoil_list'      ,geo_loc)
          call read_hdf5_al(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc)
          call read_hdf5_al(normalised_coord_e,'normalised_coord_e',geo_loc)
          call read_hdf5_al(thickness          ,'thickness'       ,geo_loc)
          
          allocate(geo%components(i_comp)%airfoil_list(size(airfoil_list)))
          geo%components(i_comp)%airfoil_list = airfoil_list
          
          allocate(geo%components(i_comp)%i_airfoil_e( &
                size(i_airfoil_e,1),size(i_airfoil_e,2)) )
          geo%components(i_comp)%i_airfoil_e = i_airfoil_e
          
          allocate(geo%components(i_comp)%normalised_coord_e( &
                size(normalised_coord_e,1),size(normalised_coord_e,2)))
          geo%components(i_comp)%normalised_coord_e = normalised_coord_e
          
          allocate(geo%components(i_comp)%thickness( &
                    size(thickness,1),size(thickness,2)))
          geo%components(i_comp)%thickness = thickness
          
          geo%components(i_comp)%aero_correction = trim(aero_table)
          sim_param%vl_correction = .true. 
        else
          geo%components(i_comp)%aero_correction = 'false'
        endif

      else if (comp_el_type(1:1) .eq. 'a') then
        call read_hdf5(trac,'Traction',cloc)
        call read_hdf5(rad,'Radius',cloc)
      end if

      !> Hinges 
      geo%components(i_comp)%n_hinges = n_hinges
      allocate( geo%components(i_comp)%hinge(n_hinges) )
      do ih = 1, n_hinges

        !> Open hinge group
        write(hinge_id_str,'(I2.2)') ih
        call open_hdf5_group(hloc, 'Hinge_'//hinge_id_str, hiloc)

        !> Read input and fill component%hinge fields
        call read_hdf5(geo%components(i_comp)%hinge(ih)%nodes_input,   'Nodes_Input', hiloc)
        call read_hdf5(geo%components(i_comp)%hinge(ih)%offset,        'Offset', hiloc)
        call read_hdf5(geo%components(i_comp)%hinge(ih)%span_blending, 'Spanwise_Blending', hiloc)
        call read_hdf5(geo%components(i_comp)%hinge(ih)%ref_dir,       'Ref_Dir', hiloc)
        call read_hdf5(geo%components(i_comp)%hinge(ih)%tag,           'Tag', hiloc)
        call read_hdf5_al(geo%components(i_comp)%hinge(ih)%ref%rr,     'rr', hiloc)

        geo%components(i_comp)%hinge(ih)%n_nodes = size( geo%components(i_comp)%hinge(ih)%ref%rr, 2)

        call read_hdf5( geo%components(i_comp)%hinge(ih)%input_type, 'Hinge_Rotation_Input'    , hiloc)

        !> if component is not already moving, and hinge is not constant, make component moving
        !> if component is moving, do nothing (keep it moving)
        !> if component is not moving, and hinge is constant, do nothing (keep it not moving)
        if ( .not. geo%components(i_comp)%moving .and.  &
           & .not. trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'function:const' ) then
          geo%components(i_comp)%moving = .true.
        endif

        !> Actual input only for input_type = function:..., otherwise dummy inputs
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_ampl , 'Hinge_Rotation_Amplitude', hiloc)
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_omega, 'Hinge_Rotation_Omega', hiloc)
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_phase, 'Hinge_Rotation_Phase', hiloc)
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_ampl_init , 'Hinge_Rotation_Amplitude_initial', hiloc)
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_cosine_cycl, 'Hinge_Rotation_cosine_cycles', hiloc)   
        call read_hdf5( geo%components(i_comp)%hinge(ih)%f_initial_time, 'Hinge_Rotation_Initial_Time', hiloc)   

        if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'coupling' ) then
          call read_hdf5_al( geo%components(i_comp)%hinge(ih)%i_coupling_nodes, &
                            'Coupling_Nodes', hiloc )
        end if

        allocate( geo%components(i_comp)%hinge(ih)%theta( &
                  geo%components(i_comp)%hinge(ih)%n_nodes ) )
        allocate( geo%components(i_comp)%hinge(ih)%theta_old( &
                  geo%components(i_comp)%hinge(ih)%n_nodes ) )

        !> Initialize theta: set the values of theta and theta_old fields
        call geo%components(i_comp)%hinge(ih)%init_theta( t=0.0_wp )

        call close_hdf5_group( hiloc )

        !> Initialize reference and actual configuration
        call initialize_hinge_config( geo%components(i_comp)%hinge(ih)%ref , &
                                      geo%components(i_comp)%hinge(ih) )
        call initialize_hinge_config( geo%components(i_comp)%hinge(ih)%act , &
                                      geo%components(i_comp)%hinge(ih) )

        !> Allocate and initialize hinge node coords in the actual configuration
        allocate( geo%components(i_comp)%hinge(ih)%act%rr( &
                3,geo%components(i_comp)%hinge(ih)%n_nodes ) )

        !> *** to do *** Only for non-coupled hinges?
        geo%components(i_comp)%hinge(ih)%act%rr = &
                                      geo%components(i_comp)%hinge(ih)%ref%rr
        !> Build hinge connectivity and weights, for grid nodes
        call geo%components(i_comp)%hinge(ih)%build_connectivity( rr, coupling_node_rot )
        !> and for cell centers
        call geo%components(i_comp)%hinge(ih)%build_connectivity_cen( rr, ee, coupling_node_rot)

      end do

#if USE_PRECICE
      !> PreCICE coupling 
      if ( trim(comp_coupling_str) .eq. 'true' ) then

        if (trim(comp_coupling_type) .eq. 'rigid' ) then
          call read_hdf5_al(c_ref_p, 'c_ref_p', geo_loc)
          call read_hdf5_al(c_ref_c, 'c_ref_c', geo_loc)
          allocate(geo%components(i_comp)%c_ref_p( size(c_ref_p,1) , &
                                                    size(c_ref_p,2) ) )
          allocate(geo%components(i_comp)%c_ref_c( size(c_ref_c,1) , &
                                                    size(c_ref_c,2) ) )
          geo%components(i_comp)%c_ref_p = c_ref_p
          geo%components(i_comp)%c_ref_c = c_ref_c

          !> For rigid coupling, the motion of all the nodes of the components
          !  is defined through the motion of the coupling_node
          np_precice = 1

          !> Allocate and fill i_points_precice array containing the
          !  connectivity between dust with PreCICE nodes
          allocate(geo%components(i_comp)%i_points_precice( np_precice ))
          geo%components(i_comp)%i_points_precice = &
                          (/((i3),i3=points_offset_precice+1, &
                                      points_offset_precice+np_precice)/)
          points_offset_precice = points_offset_precice + np_precice

        elseif ( trim(comp_coupling_type) .eq. 'rbf' ) then

          !> Load all the coupling nodes, both structural and hinges
          call read_hdf5_al( comp_coupling_nodes, 'CouplingNodes', geo_loc )

          !> Find n. nodes of coupled hinges
          n_nodes_coupling_hinges = 0
          do ih = 1, n_hinges
            if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'coupling' ) then
              !> N. nodes of the actual hinge and update overall count
              n_nodes_coupling_hinge_1 = size(geo%components(i_comp)%hinge(ih)%i_coupling_nodes,1)
              n_nodes_coupling_hinges = n_nodes_coupling_hinges + n_nodes_coupling_hinge_1

              allocate( geo%components(i_comp)%hinge(ih)%i_points_precice( n_nodes_coupling_hinge_1) )
            end if
          end do

          !> N. of coupling nodes, overall and structural only
          np_precice_tot = size(comp_coupling_nodes,2)
          np_precice     = np_precice_tot - n_nodes_coupling_hinges

          !> Allocate comp%i_points_precice()
          allocate(geo%components(i_comp)%i_points_precice(np_precice))

          !> Evaluate precice/dust connectivity for structural and hinge nodes
          allocate( ind_coupling( np_precice ) )
          allocate(hinge_ind(n_hinges)) ; hinge_ind = 0 ; comp_ind = 0

          !> Global numbering in i_points_precice: add points_offset_precice
          if ( n_hinges .gt. 0 ) then
            do i3 = 1, np_precice_tot
              is_hinge = .false.
              do ih = 1, n_hinges
                if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'coupling' ) then
                  if ( any(i3 .eq. geo%components(i_comp)%hinge(ih)%i_coupling_nodes) ) then
                    is_hinge = .true.
                    exit
                  end if
                end if
              end do
              if (is_hinge) then
                !> Hinge node, ih-th hinge
                hinge_ind(ih) = hinge_ind(ih) + 1
                geo%components(i_comp)%hinge(ih)%i_points_precice( hinge_ind(ih) ) = &
                                                  i3 + points_offset_precice
              else                                                    
                !> Structural node
                comp_ind = comp_ind + 1
                ind_coupling( comp_ind ) = i3
                geo%components(i_comp)%i_points_precice( comp_ind ) = i3 + points_offset_precice
              end if
            end do
          else
            ind_coupling = (/ ( i3, i3=1, size(ind_coupling,1) ) /)
            geo%components(i_comp)%i_points_precice = ind_coupling + points_offset_precice
          end if

          !> Define %rbf%nodes to compute struct-aero connectivity, with build_connectivity()
          ! %rbf%nodes: coordinates of the coupling nodes in the local reference frame
          ! %rbf%rrb  : coordinates of the coupling nodes in the global reference frame,
          !             - here only allocated and itialized to a meaningless value,
          !             - to be updated at each timestep, with data read from precice,
          !               in t_precice%update_elems

          allocate( geo%components(i_comp)%rbf%rrb  (3, np_precice) , &
                    geo%components(i_comp)%rbf%rrb_rot  (3, np_precice))
          geo%components(i_comp)%rbf%nodes = comp_coupling_nodes(:,ind_coupling)
          geo%components(i_comp)%rbf%rrb   = -333.3_wp

          allocate(ind_h(n_nodes_coupling_hinges))
          i3 = 1  
          do ih = 1, n_hinges
            if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'coupling' ) then
              !> N. nodes of the actual hinge and update overall count
              n_nodes_coupling_hinge_1 = size(geo%components(i_comp)%hinge(ih)%i_coupling_nodes,1)

              ind_h(i3:i3+n_nodes_coupling_hinge_1-1) = geo%components(i_comp)%hinge(ih)%i_points_precice - points_offset_precice
              i3 = i3 +n_nodes_coupling_hinge_1
            end if
          end do  

          do ih = 1, n_hinges
            if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. 'coupling' ) then

              !> Hinge node coordinates in the coupling reference space
              allocate( geo%components(i_comp)%hinge(ih)%nodes(3,hinge_ind(ih)) )
              geo%components(i_comp)%hinge(ih)%nodes = &
                  comp_coupling_nodes(:, geo%components(i_comp)%hinge(ih)%i_points_precice - &
                                        points_offset_precice )

              !> Connectivity between hinge nodes and other structural nodes
              call geo%components(i_comp)%hinge(ih)%build_connectivity_hin( &
                  comp_coupling_nodes, ind_h) 
            end if
          end do
          
          !> build connectivity for the rr nodes 
          call geo%components(i_comp)%rbf%build_connectivity(rr, coupling_node_rot)        

          !> transfer index and weight matrix 
          geo%components(i_comp)%rbf%nod%ind = geo%components(i_comp)%rbf%point%ind
          geo%components(i_comp)%rbf%nod%wei = geo%components(i_comp)%rbf%point%wei                 
          deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)

          !> build connectivity for the rr_virtual nodes 
          call geo%components(i_comp)%rbf%build_connectivity(rr_virtual, coupling_node_rot)        

          !> transfer index and weight matrix 
          geo%components(i_comp)%rbf%nod_virtual%ind = geo%components(i_comp)%rbf%point%ind
          geo%components(i_comp)%rbf%nod_virtual%wei = geo%components(i_comp)%rbf%point%wei                 
          deallocate(geo%components(i_comp)%rbf%point%ind, geo%components(i_comp)%rbf%point%wei)

          !> Update offset of precice/dust coupling nodes
          points_offset_precice = points_offset_precice + np_precice_tot

          allocate(geo%components(i_comp)%c_ref_c(0,0))
          allocate(geo%components(i_comp)%c_ref_p(0,0))

          deallocate(ind_coupling, hinge_ind, ind_h)  

        endif
      else
        np_precice = 0
        allocate(geo%components(i_comp)%i_points_precice( np_precice ))
      end if
#endif


      !> For PARAMETRIC elements only:
      !  parametric_nelems_span , parametric_nelems_chor
      if ( trim(comp_input) .eq. 'parametric' .or. &
          trim(comp_input) .eq. 'pointwise') then
        call read_hdf5(par_nelems_span,'parametric_nelems_span',geo_loc)
        call read_hdf5(par_nelems_chor,'parametric_nelems_chor',geo_loc)
        call read_hdf5(par_nelems_chor_virtual,'parametric_nelems_chor_virtual',geo_loc)
        geo%components(i_comp)%n_s = par_nelems_span
        geo%components(i_comp)%n_c = par_nelems_chor         
        geo%components(i_comp)%n_c_virtual = par_nelems_chor_virtual   
      end if

      call close_hdf5_group(geo_loc)


      ! Trailing Edge (not all elements build the trailing edge)
      if( comp_el_type(1:1) .eq. 'p' .or. &
          comp_el_type(1:1) .eq. 'v' .or. &
          comp_el_type(1:1) .eq. 'l') then
        call open_hdf5_group(      cloc,  'Trailing_Edge', te_loc)
        call read_hdf5_al(         e_te,           'e_te', te_loc)
        call read_hdf5_al(e_te_fountain,  'e_te_fountain', te_loc)
        call read_hdf5(    n_e_te_panel,   'n_e_te_panel', te_loc)
        call read_hdf5_al(         i_te,           'i_te', te_loc)
        call read_hdf5_al(        ii_te,          'ii_te', te_loc)
        call read_hdf5_al(     neigh_te,       'neigh_te', te_loc)
        call read_hdf5_al(         o_te,           'o_te', te_loc)
        call read_hdf5_al(         t_te,           't_te', te_loc)
        !> debug 
        !do i1 = 1, size(e_te,2)
        !  write(*,*) 'e_te         ', e_te(1,i1) 
        !  write(*,*) 'e_te_fountain', e_te_fountain(1,i1)
        !enddo 
        if(check_dset_hdf5('scale_te',te_loc)) then
          call read_hdf5(scale_te, 'scale_te',te_loc)
        else
          scale_te = 1.0_wp
        endif
        call close_hdf5_group(te_loc)
      else
        allocate(e_te(0,0), i_te(2,0), ii_te(2,0), neigh_te(2,0), o_te(2,0),&
                t_te(3,0), e_te_fountain(2,0))
      endif
      ! Re-write all that was read in the new output file (cannot just copy
      ! the read file since the components order will be changed when emplying
      ! the multiple components). If it is restarting with the same name
      ! avoid write

      write(cname_write,'(A,I3.3)') 'Comp',geo%components(i_comp)%comp_id
      if(rewrite_geo) then
        call new_hdf5_group(gloc_out,trim(cname_write),cloc2)
        call write_hdf5(ref_id,'RefId',cloc2)
        call write_hdf5(trim(ref_tag_m),'RefTag',cloc2)
        call write_hdf5(trim(comp_el_type),'ElType',cloc2)
        call write_hdf5(trim(geo%components(i_comp)%comp_name), &
                                                            'CompName',cloc2)
        call write_hdf5(trim(geo%components(i_comp)%comp_input),&
                                                            'CompInput',cloc2)
        if ( geo%components(i_comp)%coupling ) then
          call write_hdf5('true','Coupled',cloc2)
          call write_hdf5(trim(geo%components(i_comp)%coupling_type), &
                          'CouplingType',cloc2)
          call write_hdf5(geo%components(i_comp)%coupling_node, &
                          'CouplingNode',cloc2)
          call write_hdf5(geo%components(i_comp)%coupling_node_rot, &
                          'CouplingNodeOrientation',cloc2)

        else
          call write_hdf5('false','Coupled',cloc2)
          call write_hdf5('none', 'CouplingType',cloc2)
          call write_hdf5((/0.0_wp,0.0_wp,0.0_wp/), 'CouplingNode',cloc2)
          call write_hdf5(reshape( (/ 1._wp, 0._wp, 0._wp, &
                                      0._wp, 1._wp, 0._wp, &
                                      0._wp, 0._wp, 1._wp/),(/3,3/)), &
                          'CouplingNodeOrientation',cloc2)
        end if
        call new_hdf5_group(cloc2,'Geometry',geo_loc)

        call write_hdf5(ee   ,'ee'   ,geo_loc)
        call write_hdf5(rr   ,'rr'   ,geo_loc)
        call write_hdf5(ee_virtual   ,'ee_virtual'   ,geo_loc)
        call write_hdf5(rr_virtual   ,'rr_virtual'   ,geo_loc)
        call write_hdf5(neigh,'neigh',geo_loc)

#if USE_PRECICE
        call write_hdf5(geo%components(i_comp)%i_points_precice, &
                        'i_points_PreCICE', geo_loc)
#endif
        if ( comp_el_type(1:1) .eq. 'l' ) then
          call write_hdf5(airfoil_list       ,'airfoil_list'        ,geo_loc)
          call write_hdf5(nelem_span_list    ,'nelem_span_list'     ,geo_loc)
          call write_hdf5(i_airfoil_e        ,'i_airfoil_e'         ,geo_loc)
          call write_hdf5(normalised_coord_e ,'normalised_coord_e'  ,geo_loc)
          call write_hdf5(theta_e            ,'theta_e'             ,geo_loc)

        elseif (comp_el_type(1:1) .eq. 'v' ) then
          call write_hdf5(trim(aero_table)     ,'aero_table'        ,geo_loc)          
          if (trim(aero_table) .eq. 'true') then
            call write_hdf5(airfoil_list       ,'airfoil_list'      ,geo_loc)
            call write_hdf5(i_airfoil_e        ,'i_airfoil_e'       ,geo_loc)
            call write_hdf5(normalised_coord_e ,'normalised_coord_e',geo_loc)
            call write_hdf5(thickness          ,'thickness'         ,geo_loc)
          endif
        else if (comp_el_type(1:1) .eq. 'a') then
          call write_hdf5(trac, 'Traction',cloc2)
          call write_hdf5(rad,  'Radius'  ,cloc2)
        endif

        if (trim(comp_input) .eq. 'parametric' .or. & 
            trim(comp_input) .eq. 'pointwise') then
          !> Write HDF5 fields
          call write_hdf5(par_nelems_span        ,'parametric_nelems_span'        ,geo_loc)
          call write_hdf5(par_nelems_chor        ,'parametric_nelems_chor'        ,geo_loc)
          call write_hdf5(par_nelems_chor_virtual,'parametric_nelems_chor_virtual',geo_loc)
        end if

        call close_hdf5_group(geo_loc)

        !> Hinges
        call new_hdf5_group(cloc2, 'Hinges' , hloc2)
        call write_hdf5(n_hinges, 'n_hinges', hloc2)

        do ih = 1, n_hinges
          !> Open hinge group
          write(hinge_id_str,'(I2.2)') ih
          call new_hdf5_group(hloc2, 'Hinge_'//hinge_id_str, hiloc)

          !> Read input and fill component%hinge fields
          call write_hdf5( geo%components(i_comp)%hinge(ih)%nodes_input, &
                                                            'Nodes_Input', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%offset, &
                                                            'Offset', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%span_blending, &
                                                  'Spanwise_Blending', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%ref_dir, &
                                                            'Ref_Dir', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%tag, 'Tag', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%ref%rr, 'rr', hiloc)

          call write_hdf5( geo%components(i_comp)%hinge(ih)%input_type, &
                                                'Hinge_Rotation_Input'    , hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_ampl , &
                                                'Hinge_Rotation_Amplitude', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_omega, &
                                                'Hinge_Rotation_Omega', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_phase, &
                                                'Hinge_Rotation_Phase', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_ampl_init, &
                                                'Hinge_Rotation_Amplitude_initial', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_cosine_cycl, &
                                                'Hinge_Rotation_cosine_cycles', hiloc)
          call write_hdf5( geo%components(i_comp)%hinge(ih)%f_initial_time, &
                                                'Hinge_Rotation_Initial_Time', hiloc)
          call close_hdf5_group(hiloc)

        end do

        call close_hdf5_group(hloc2)

        if( comp_el_type(1:1) .eq. 'p' .or. &
            comp_el_type(1:1) .eq. 'v' .or. &
            comp_el_type(1:1) .eq. 'l') then
          call new_hdf5_group(cloc2,    'Trailing_Edge', te_loc)
          call write_hdf5(    e_te,              'e_te', te_loc)
          call write_hdf5(e_te_fountain,'e_te_fountain', te_loc)
          call write_hdf5(n_e_te_panel,  'n_e_te_panel', te_loc)
          call write_hdf5(    i_te,              'i_te', te_loc)
          call write_hdf5(   ii_te,             'ii_te', te_loc)
          call write_hdf5(neigh_te,          'neigh_te', te_loc)
          call write_hdf5(    o_te,              'o_te', te_loc)
          call write_hdf5(    t_te,              't_te', te_loc)
          call write_hdf5(scale_te,          'scale_te', te_loc)
          call close_hdf5_group(te_loc)
        endif

        call close_hdf5_group(cloc2)
      endif

      ! ======= CREATING ELEMENTS ======

      !> Treat the points 
      if(allocated(geo%points)) then
        points_offset = size(geo%points,2)
      else
        points_offset = 0
      endif

      !> Treat the virtual points 
      if(allocated(geo%points_virtual)) then
        points_offset_virtual = size(geo%points_virtual,2)
      else
        points_offset_virtual = 0
      endif


      !> store the read points into the local points
      allocate(geo%components(i_comp)%loc_points(3,size(rr,2)))
      allocate(geo%components(i_comp)%loc_points_virtual(3,size(rr_virtual,2)))
      geo%components(i_comp)%loc_points = rr
      geo%components(i_comp)%loc_points_virtual = rr_virtual


      !> Now for the moments the points are stored here without moving them,
      !  will be moved later, consider not storing them here at all
      allocate(points_tmp(3,size(rr,2)+points_offset))
      if (points_offset .gt. 0) points_tmp(:,1:points_offset) = geo%points
      points_tmp(:,points_offset+1:points_offset+size(rr,2)) = rr
      call move_alloc(points_tmp, geo%points)
      allocate(geo%components(i_comp)%i_points(size(rr,2)))
      geo%components(i_comp)%i_points = &
                        (/((i3),i3=points_offset+1,points_offset+size(rr,2))/)

      allocate(points_virtual_tmp(3,size(rr_virtual,2) + points_offset_virtual))
      if (points_offset_virtual .gt. 0) points_virtual_tmp(:,1:points_offset_virtual) = geo%points_virtual
      points_virtual_tmp(:,points_offset_virtual+1:points_offset_virtual+size(rr_virtual,2)) = rr_virtual
      call move_alloc(points_virtual_tmp, geo%points_virtual)
      allocate(geo%components(i_comp)%i_points_virtual(size(rr_virtual,2)))
      geo%components(i_comp)%i_points_virtual = &
                        (/((i3),i3=points_offset_virtual+1,points_offset_virtual+size(rr_virtual,2))/)                      

      !> Treat the elements 
      !allocate the elements of the component of the right kind
      geo%components(i_comp)%nelems = size(ee,2)
      geo%components(i_comp)%nelems_virtual = size(ee_virtual,2)
      allocate(geo%components(i_comp)%el_virtual(size(ee_virtual,2)))
      
      select case(trim(geo%components(i_comp)%comp_el_type))
        case('p')
          allocate(t_surfpan::geo%components(i_comp)%el(size(ee,2)))
        case('v')
          allocate(t_vortlatt::geo%components(i_comp)%el(size(ee,2)))
        case('l')
          allocate(t_liftlin::geo%components(i_comp)%el(size(ee,2)))
        case('a')
          allocate(t_actdisk::geo%components(i_comp)%el(size(ee,2)))
        case default
          call error(this_sub_name, this_mod_name, &
            'Unknown type of element: '//geo%components(i_comp)%comp_el_type)
      end select
      
      !> fill (some) of the real elements fields
      do i2=1,size(ee,2)
        !> vertices
        n_vert = count(ee(:,i2).ne.0)
        allocate(geo%components(i_comp)%el(i2)%i_ver(n_vert))
        allocate(geo%components(i_comp)%el(i2)%neigh(n_vert))
        geo%components(i_comp)%el(i2)%n_ver = n_vert
        geo%components(i_comp)%el(i2)%i_ver(1:n_vert) = &
                                              ee(1:n_vert,i2) + points_offset
        do i3=1,n_vert
          if ( neigh(i3,i2) .ne. 0 ) then
            geo%components(i_comp)%el(i2)%neigh(i3)%p => &
                                    geo%components(i_comp)%el(neigh(i3,i2))
          else
            !> Do nothing, keep the neighbour pointer not associated
            geo%components(i_comp)%el(i2)%neigh(i3)%p => null()
          end if
        end do

        !> Motion
        geo%components(i_comp)%el(i2)%moving = geo%components(i_comp)%moving

      enddo
      n_vert = 0

      !> fill (some) of the virtual elements fields
      do i2=1,size(ee_virtual,2)

        !> vertices
        n_vert = count(ee_virtual(:,i2).ne.0)
        allocate(geo%components(i_comp)%el_virtual(i2)%i_ver(n_vert))
        geo%components(i_comp)%el_virtual(i2)%n_ver = n_vert
        geo%components(i_comp)%el_virtual(i2)%i_ver(1:n_vert) = &
                                              ee_virtual(1:n_vert,i2) + points_offset_virtual
        !> Motion
        geo%components(i_comp)%el_virtual(i2)%moving = geo%components(i_comp)%moving

      enddo

      geo%components(i_comp)%nSurfPan = 0; geo%components(i_comp)%nVortRin = 0;
      if(comp_el_type(1:1) .eq. 'p') &
                                  geo%components(i_comp)%nSurfPan = size(ee,2)
      if(comp_el_type(1:1) .eq. 'v') &
                                  geo%components(i_comp)%nVortRin = size(ee,2)

      !> If it is an actuator disk read the traction
      if(geo%components(i_comp)%comp_el_type(1:1) .eq. 'a') then
        select type (el=>geo%components(i_comp)%el)
        type is(t_actdisk)
          do i2 = 1,size(el)
            el(i2)%traction = trac
            el(i2)%radius = rad
          enddo
        end select
      endif

      !> Trailing Edge 
      ne_te = 0; nn_te=0
      if(allocated(e_te)) ne_te = size(e_te,2)
      if(allocated(i_te)) nn_te = size(i_te,2)
      if (sim_param%debug_level .ge. 3) then
        write(msg,'(A,I0,A,I0)') ' Trailing edge: elements: ', &
                                                        ne_te, ' nodes ', nn_te
        call printout(msg)
      endif
      if (.not.allocated(te%e)) then ! it should be enough
        allocate(te%e         (2,ne_te))
        allocate(te%e_fountain(2,ne_te))
        !do i1 = 1, ne_te
        !  write(*,*) 'e_te(1,i1) = ', e_te(1,i1) 
        !  write(*,*) 'e_te_fountain(1,i1) = ', e_te_fountain(1,i1) 
        !enddo 
        !> number of panels in the trailing edge
        te%nte_surfpan = n_e_te_panel   
        do i1 = 1, ne_te
          te%e(1,i1)%p => null()
          te%e(2,i1)%p => null()
          te%e_fountain(1,i1)%p => null()
          te%e_fountain(2,i1)%p => null() 
          te%e(1,i1)%p  => geo%components(i_comp)%el(e_te(1,i1))
          if (e_te_fountain(1,i1) .gt. 0) then  
            te%e_fountain(1,i1)%p  => geo%components(i_comp)%el(e_te_fountain(1,i1))
          endif 
          if(e_te(2,i1) .gt. 0) &
            te%e(2,i1)%p  => geo%components(i_comp)%el(e_te(2,i1))            
          if(e_te_fountain(2, i1) .gt. 0) &
            te%e_fountain(2,i1)%p  => geo%components(i_comp)%el(e_te_fountain(2,i1)) 
        enddo

        allocate(te%i    (2,nn_te) ) ; te%i     =     i_te
        allocate(te%ii   (2,ne_te) ) ; te%ii    =    ii_te
        allocate(te%neigh(2,ne_te) ) ; te%neigh = neigh_te
        allocate(te%o    (2,ne_te) ) ; te%o     =     o_te
        allocate(te%t    (3,nn_te) ) ; te%t     =     t_te
        allocate(te%t_hinged    (3,nn_te) ) ; te%t_hinged     =     te%t
        allocate(te%is_hinged    (nn_te) ) ; te%is_hinged     =     -1
        allocate(te%ref  (  nn_te) ) ; te%ref   = geo%components(i_comp)%ref_id
        allocate(te%icomp(  nn_te) ) ; te%icomp = i_comp
        allocate(te%scaling(  nn_te) ) 
        
        te%scaling = scale_te*sim_param%first_panel_scaling

      elseif (ne_te .gt. 0) then
        nn_te_prev = size(te%i,2)
        ne_te_prev = size(te%e,2)

        allocate(e_te_tmp(2,size(te%e,2)+ne_te))
        allocate(e_te_fountain_tmp(2,size(te%e,2)+ne_te)) 
        e_te_tmp(:,             1:size(te%e,2)    ) = te%e
        e_te_fountain_tmp(:,    1:size(te%e,2)    ) = te%e_fountain 
        do i1 = 1,ne_te 
          e_te_tmp(1,size(te%e,2)+i1)%p => null()
          e_te_tmp(2,size(te%e,2)+i1)%p => null()
          e_te_fountain_tmp(1,size(te%e,2)+i1)%p => null() 
          e_te_fountain_tmp(2,size(te%e,2)+i1)%p => null() 
          e_te_tmp(1,size(te%e,2)+i1)%p  => &
                                      geo%components(i_comp)%el(e_te(1,i1))
          if (e_te_fountain(1,i1) .gt. 0) then
            e_te_fountain_tmp(1,size(te%e,2)+i1)%p  => &
                                      geo%components(i_comp)%el(e_te_fountain(1,i1))
          endif
          if(e_te(2,i1) .gt. 0) then
            e_te_tmp(2,size(te%e,2)+i1)%p  => &
                                      geo%components(i_comp)%el(e_te(2,i1))
          endif 
          if(e_te_fountain(2,i1) .gt. 0) then
            e_te_fountain_tmp(2,size(te%e,2)+i1)%p  => &
                                      geo%components(i_comp)%el(e_te_fountain(2,i1)) 
          endif
        enddo
        call move_alloc(e_te_tmp, te%e)
        call move_alloc(e_te_fountain_tmp, te%e_fountain) 

        allocate(i_te_tmp(2,size(te%i,2)+nn_te))
        i_te_tmp(:,             1:size(te%i,2)    ) = te%i
        i_te_tmp(:,size(te%i,2)+1:size(i_te_tmp,2)) = i_te + points_offset
        call move_alloc(i_te_tmp,te%i)
      
        allocate(ii_te_tmp(2,size(te%ii,2)+ne_te))
        ii_te_tmp(:,              1:size(te%ii,2)    ) = te%ii
        ii_te_tmp(:,size(te%ii,2)+1:size(ii_te_tmp,2)) = ii_te + nn_te_prev
        call move_alloc(ii_te_tmp,te%ii)

        allocate(neigh_te_tmp(2,size(te%neigh,2)+ne_te))
        neigh_te_tmp(:,                 1:size(te%neigh,2)    ) = te%neigh
        where ( neigh_te .ne. 0 ) neigh_te = neigh_te + ne_te_prev
        neigh_te_tmp(:,size(te%neigh,2)+1:size(neigh_te_tmp,2)) = neigh_te
        call move_alloc(neigh_te_tmp,te%neigh)

        allocate( o_te_tmp(2,size(te%o ,2)+ne_te))
        o_te_tmp(:,              1:size(te%o ,2)    ) = te%o
        o_te_tmp(:,size(te%o ,2)+1:size( o_te_tmp,2)) =  o_te
        call move_alloc(o_te_tmp,te%o )

        allocate( t_te_tmp(3,size(te%t ,2)+nn_te))
        t_te_tmp(:,              1:size(te%t ,2)    ) = te%t
        t_te_tmp(:,size(te%t ,2)+1:size( t_te_tmp,2)) =  t_te
        call move_alloc(t_te_tmp,te%t )

        allocate( t_hinged_te_tmp(3,size(te%t_hinged ,2)+nn_te))
        t_hinged_te_tmp(:,              1:size(te%t_hinged ,2)    ) = te%t_hinged
        t_hinged_te_tmp(:,size(te%t_hinged ,2)+1:size( t_hinged_te_tmp,2)) = t_te
        call move_alloc(t_hinged_te_tmp,te%t_hinged )

        allocate( t_is_hinged_tmp(size(te%is_hinged,1)+nn_te))
        t_is_hinged_tmp(1:size(te%is_hinged ,1)) = te%is_hinged
        t_is_hinged_tmp(size(te%is_hinged ,1)+1:size( t_is_hinged_tmp,1)) =  -1
        call move_alloc(t_is_hinged_tmp,te%is_hinged )

        allocate(ref_te_tmp(size(te%ref   )+nn_te))
        ref_te_tmp(                 1:size(te%ref   )   ) =  te%ref
        ref_te_tmp(  size(te%ref  )+1:size(ref_te_tmp  )) = &
                                                  geo%components(i_comp)%ref_id
        call move_alloc(ref_te_tmp,te%ref)

        allocate(icomp_te_tmp(size(te%icomp)+nn_te))
        icomp_te_tmp(               1:size(te%icomp    )) = te%icomp
        icomp_te_tmp(size(te%icomp)+1:size(icomp_te_tmp)) = i_comp
        call move_alloc(icomp_te_tmp,te%icomp)

        allocate( scale_te_tmp(size(te%scaling)+nn_te))
        scale_te_tmp(1:size(te%scaling )    ) = te%scaling
        scale_te_tmp(size(te%scaling )+1:size( scale_te_tmp)) = &
                                        scale_te*sim_param%first_panel_scaling
        call move_alloc(scale_te_tmp,te%scaling )
      end if

      ! 2018-07-5. Deassociate neighboring elements through the TE
      if(trim(geo%components(i_comp)%comp_el_type) .eq. 'p') then
      do i1 = 1,ne_te
        do i2 = 1,geo%components(i_comp)%el(e_te(1,i1))%n_ver
          if ( neigh(i2,e_te(1,i1)) .eq. e_te(2,i1) ) then
            geo%components(i_comp)%el(e_te(1,i1))%neigh(i2)%p => null()
          end if
        end do
        do i2 = 1,geo%components(i_comp)%el(e_te(2,i1))%n_ver
          if ( neigh(i2,e_te(2,i1)) .eq. e_te(1,i1) ) then
            geo%components(i_comp)%el(e_te(2,i1))%neigh(i2)%p => null()
          end if
        end do
      end do
      endif

      deallocate(e_te, i_te, ii_te, neigh_te, o_te, t_te, e_te_fountain)

      !> Update elems_offset for the next component
      elems_offset = elems_offset + size(ee,2)
      elems_offset_virtual = elems_offset_virtual + size(ee_virtual,2)
      !> Cleanup
      deallocate(ee, ee_virtual, rr, rr_virtual, neigh)
      i_comp = i_comp + 1
    enddo !i_mult

    call close_hdf5_group(hloc)
    call close_hdf5_group(cloc)

  enddo !i_comp

  !> Update the total number of components
  !  call write_hdf5(i_comp-1,'NComponents',gloc)
  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)
  if(rewrite_geo) then
    call close_hdf5_group(gloc_out)
    call close_hdf5_file(floc_out)
  endif

  !> Define comp_id for all the elems
  do i_comp = 1 , size(geo%components)
    do i1 = 1 , size(geo%components(i_comp)%el)
      geo%components(i_comp)%el(i1)%comp_id = i_comp
      geo%components(i_comp)%el_virtual(i1)%comp_id = i_comp
    end do
  end do

#if USE_PRECICE
  allocate(geo%points_vel(3, size(geo%points,2))); geo%points_vel = 0.0_wp
  allocate(geo%points_vel_virtual(3, size(geo%points_virtual,2))); geo%points_vel_virtual = 0.0_wp
#endif


end subroutine load_components

!----------------------------------------------------------------------

subroutine import_aero_tab(geo,coeff)
  type(t_geo), intent(inout), target             :: geo
  type(t_aero_tab) , allocatable , intent(inout) :: coeff(:)
  integer                                        :: n_tmp , n_tmp2
  character(len=max_char_len) , allocatable      :: list_tmp(:)
  character(len=max_char_len) , allocatable      :: list_tmp_tmp(:)
  integer , allocatable                          :: i_airfoil_e_tmp (:,:)
  integer                                        :: i_c, n_c, i_a, n_a, i_l

  n_tmp = 30
  allocate(list_tmp(n_tmp))  

  ! Count # of different airfoil
  n_c = size(geo%components)
  n_a = 0
  do i_c = 1, n_c

    if ( trim(geo%components(i_c)%comp_el_type) .eq. 'l' .or. &
        (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
        trim(geo%components(i_c)%aero_correction) .eq. 'true')) then 

      allocate( i_airfoil_e_tmp( &
          size(geo%components(i_c)%i_airfoil_e,1), &
          size(geo%components(i_c)%i_airfoil_e,2) ) )
      i_airfoil_e_tmp = 0

      do i_a = 1 , size(geo%components(i_c)%airfoil_list)

        if (all( geo%components(i_c)%airfoil_list(i_a) .ne. &
                                                      list_tmp(1:n_a))) then
          ! new airfoil ----
          n_a = n_a + 1

          where (geo%components(i_c)%i_airfoil_e .eq. i_a) i_airfoil_e_tmp = n_a

          ! if n_a > n_tmp --> movalloc
          if ( n_a .gt. n_tmp ) then
            n_tmp2 = n_tmp + n_tmp
            allocate(list_tmp_tmp(n_tmp2))
            list_tmp_tmp(1:n_tmp) = list_tmp
            deallocate(list_tmp)
            call move_alloc(list_tmp_tmp,list_tmp)
            n_tmp = n_tmp2
          end if

          list_tmp(n_a) = geo%components(i_c)%airfoil_list(i_a)

        else
          !> Airfoil already used: find the element and replace the global index
          do i_l = 1 , n_a
            if ( geo%components(i_c)%airfoil_list(i_a) .eq. list_tmp(i_l) ) exit
          end do

          where (geo%components(i_c)%i_airfoil_e .eq. i_a) i_airfoil_e_tmp = i_l

        end if

      end do

    geo%components(i_c)%i_airfoil_e = i_airfoil_e_tmp

    deallocate(i_airfoil_e_tmp)

    end if

  end do

  !> Read tables and fill coeff structure
  allocate(coeff(n_a))
  do i_a = 1 , n_a
    call read_c81_table( list_tmp(i_a) , coeff(i_a) )
  end do


end subroutine import_aero_tab

!----------------------------------------------------------------------

!> Prepare the geometry allocating all the relevant fields for each
!  kind of element
subroutine prepare_geometry(geo)
  type(t_geo), intent(inout), target :: geo
  integer                            :: i_comp, ie
  integer                            :: nsides
  class(c_pot_elem), pointer         :: elem
  class(c_elem_virtual), pointer     :: elem_virtual
  character(len=*), parameter        :: this_sub_name = 'prepare_geometry'

  do i_comp = 1,size(geo%components)

    do ie = 1,size(geo%components(i_comp)%el)
      
      elem => geo%components(i_comp)%el(ie)

      nsides = size(elem%i_ver)

      !> Fields common to each element
      allocate(elem%ver(3,nsides))
      allocate(elem%edge_vec(3,nsides))
      allocate(elem%edge_len(nsides))
      allocate(elem%edge_uni(3,nsides))

      !> Type-specific fields
      select type(elem)
        !> Surface panel
        class is(t_surfpan)
          allocate(elem%verp(3,nsides))
          allocate(elem%cosTi(nsides))
          allocate(elem%sinTi(nsides))

        class is(t_liftlin)

          allocate(elem%tang_cen(3))
          allocate(elem%bnorm_cen(3))  

         elem%csi_cen = 0.5_wp * &
                        sum(geo%components(i_comp)%normalised_coord_e(:,ie))
          elem%i_airfoil =  geo%components(i_comp)%i_airfoil_e(:,ie)
          elem%twist     =  geo%components(i_comp)%theta_e(ie)
        class is(t_vortlatt)

          
        class is(t_actdisk)

        class default
          call error(this_sub_name, this_mod_name, 'Unknown element type')
      end select

      !> Position and velocity increment, due to hinge motion
      elem%dcen_h     = 0.0_wp
      elem%dcen_h_old = 0.0_wp
      elem%dvel_h     = 0.0_wp
    enddo

    do ie = 1,size(geo%components(i_comp)%el_virtual)
      
      elem_virtual => geo%components(i_comp)%el_virtual(ie)

      nsides = size(elem_virtual%i_ver)

      !> Fields common to each element
      allocate(elem_virtual%ver(3,nsides))
      allocate(elem_virtual%edge_vec(3,nsides))
      allocate(elem_virtual%edge_len(nsides))
      allocate(elem_virtual%edge_uni(3,nsides))

      !> Position and velocity increment, due to hinge motion
      elem_virtual%dcen_h     = 0.0_wp
      elem_virtual%dcen_h_old = 0.0_wp
      elem_virtual%dvel_h     = 0.0_wp
    enddo
    
  enddo

end subroutine prepare_geometry

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! panel element (vortex ring or surface panel)
subroutine calc_geo_data_pan(elem, vert)
  class(c_pot_elem), intent(inout) :: elem
  real(wp), intent(in)             :: vert(:,:)
  integer                          :: nsides, is
  real(wp)                         :: nor(3), tanl(3)

  nsides = size(vert,2)

  elem%n_ver = nsides

  !> Vertices
  elem%ver = vert

  !> Center
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)
  
  !> Unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross(vert(:,3) - vert(:,1) , &
                vert(:,4) - vert(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross(vert(:,3) - vert(:,2) , &
                vert(:,1) - vert(:,2)     )
  end if

  elem%area = 0.5_wp * norm2(nor)
  elem%nor = nor / norm2(nor)

  !> Local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( vert(:,nsides) + vert(:,1) ) - elem%cen

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )

  !> Vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_tri(is)) - vert(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_qua(is)) - vert(:,is)
    end do
  end if

  !> Edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  !> Unit vector
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  !> Surface panels own fields
  select type(elem)
  type is(t_surfpan)
    do is = 1 , nsides
      !> cosTi , sinTi
      elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) )
      elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) )
      !> Projection of the vertices on the mean plane
      elem%verp(:,is) = vert(:,is) - elem%nor * &
                      sum( (vert(:,is) - elem%cen ) * elem%nor )
    end do
  end select

  !> Initialise %dforce
  elem%dforce = 0.0_wp

  !> Initialise %dmom
  elem%dmom   = 0.0_wp

end subroutine calc_geo_data_pan


!> Calculate the geometrical quantities of a lifting line element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! lifting line element
subroutine calc_geo_data_ll(elem,vert)
  class(c_pot_elem), intent(inout) :: elem
  real(wp), intent(in) :: vert(:,:)
  integer :: is, nsides
  real(wp):: nor(3), tanl(3)

  nsides = size(vert,2)

  elem%n_ver = nsides

  !> Vertices
  elem%ver = vert

  !> Center, for the lifting line is the mid-point
  elem%cen =  sum ( vert(:,1:2),2 ) / 2.0_wp
  !> Unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross(vert(:,3) - vert(:,1) , &
                vert(:,4) - vert(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross(vert(:,3) - vert(:,2) , &
                vert(:,1) - vert(:,2)     )
  end if

  elem%area = 0.5_wp * norm2(nor)
  elem%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( vert(:,nsides) + vert(:,1) ) - elem%cen

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )


  !> Vector connecting two consecutive vertices:
  !> edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_tri(is)) - vert(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_qua(is)) - vert(:,is)
    end do
  end if

  !> Edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  !> Unit vector
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  ! ll-specific fields
  select type(elem)
    type is(t_liftlin)
    elem%tang_cen = elem%edge_uni(:,2) - elem%edge_uni(:,4)
    elem%tang_cen = elem%tang_cen / norm2(elem%tang_cen)

    elem%bnorm_cen = cross(elem%tang_cen, elem%nor)
    elem%bnorm_cen = elem%bnorm_cen / norm2(elem%bnorm_cen)
    elem%chord = sum(elem%edge_len((/2,4/)))*0.5_wp
  end select

  !> Initialise %dforce
  elem%dforce = 0.0_wp

  !> Initialise %dmom
  elem%dmom   = 0.0_wp

end subroutine calc_geo_data_ll


!> Calculate the geometrical quantities of an actuator disk
!!
!! The subroutine calculates all the relevant geometrical quantities of an
!! actuator disk
subroutine calc_geo_data_ad(elem,vert)
  class(c_pot_elem), intent(inout)  :: elem
  real(wp), intent(in)              :: vert(:,:)
  integer                           :: nsides, is
  real(wp)                          :: nor(3), tanl(3)
  integer                           :: nxt

  nsides = size(vert,2)

  elem%n_ver = nsides

  !> Vertices
  elem%ver = vert

  !> Center, for the lifting line is the mid-point
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)

  elem%area = 0.0_wp; elem%nor = 0.0_wp
  do is = 1, nsides
    nxt = 1+mod(is,nsides)
    nor = cross(vert(:,is) - elem%cen,&
                vert(:,nxt) - elem%cen )
    elem%area = elem%area + 0.5_wp * norm2(nor)
    elem%nor = elem%nor + nor/norm2(nor)
  enddo
    elem%nor = elem%nor/real(nsides,wp)

  !> Local tangent unit vector: aligned with first node, normal to n
  tanl = (vert(:,1)-elem%cen)-sum((vert(:,1)-elem%cen)*elem%nor)*elem%nor

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )

  !> Vector connecting two consecutive vertices:
  do is = 1 , nsides
    nxt = 1+mod(is,nsides)
    elem%edge_vec(:,is) = vert(:,nxt) - vert(:,is)
  end do

  !> Edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  !> Unit vector
  do is = 1 , nsides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  !> Initialise %dforce
  elem%dforce = 0.0_wp

  !> Initialise %dmom
  elem%dmom   = 0.0_wp

end subroutine calc_geo_data_ad


!> Calculate the local velocity on the panels to then enforce the
!! boundary condition
subroutine calc_geo_vel(elem, G, f)
  class(c_pot_elem), intent(inout)  :: elem
  real(wp), intent(in)              :: f(3), G(3,3)

  elem%ub = f + matmul(G,elem%cen)

end subroutine calc_geo_vel

subroutine calc_geo_vel_virtual(elem, G, f)
  class(c_elem_virtual), intent(inout)  :: elem
  real(wp), intent(in)              :: f(3), G(3,3)

  elem%ub = f + matmul(G, elem%cen)

end subroutine calc_geo_vel_virtual

!> Calculate velocity of a point whose coordinate is rr
!! boundary condition
!!
subroutine calc_node_vel( r, G, f, v)
  real(wp), intent(in)  :: r(3)
  real(wp), intent(in)  :: f(3), G(3,3)
  real(wp), intent(out) :: v(3)

  v = f + matmul(G,r)

end subroutine calc_node_vel

!----------------------------------------------------------------------

!> Updatad node-element connectivity ( ee array )
! It is assumed that nodes and elements are sorted as follows
!
!  1-----5-----9-----13----17----21   <--- LE
!  |  1  |  4  |  7  |  10 |  13 |
!  2-----6-----10----14----18----22
!  |  2  |  5  |  8  |  11 |  14 |
!  3-----7-----11----15----19----23
!  |  3  |  6  |  9  |  12 |  15 |
!  4-----8-----12----16----20----24   <--- TE
! and the ee array is built as  5, 1, 2, 6
!                               6, 2, 3, 7
!                               7, 3, 4, 8
!                               9, 5, 6,10
!                                , ...
!                              23,19,20,24
subroutine create_strip_connectivity(geo)
  type(t_geo),  target, intent(inout)  :: geo
  real(wp)                             :: h2
  integer                              :: io_te , io_tip
  integer                              :: n_el , ie_ind
  integer                              :: i_comp , i_el
  integer                              :: i_c , i_s , n_s , n_c , i_c2, n_pan, i
  real(wp)                             :: nor(3), tanl(3)
  character(len=*), parameter          :: this_sub_name = 'create_strip_connectivity'

  do i_comp = 1 , size(geo%components)
    associate(comp => geo%components(i_comp))
      
    !> Connectivity built for lifting surface only
    if ( ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' .or. &
            geo%components(i_comp)%comp_el_type(1:1) .eq. 'l' )  ) then
      
      
      n_el = size(comp%el)
            
      !> Some checks added: the first element must be at the LE,
      !  the last at the TE
      if ( associated(comp%el(  1 )%neigh(1)%p)) then
        call error(this_sub_name, this_mod_name, &
                  'First element must be at the LE')
      end if
      if ( associated(comp%el(n_el)%neigh(3)%p) ) then
        call error(this_sub_name, this_mod_name, &
                  'Last element must be at the TE')
      end if
    
      !> Initialisation
      io_tip = 1 ; io_te = 0 ; n_s = 0
      do i_el = 1 , n_el
        ie_ind = comp%el(i_el)%id
        if ( .not. associated(comp%el(i_el)%neigh(1)%p) ) then
          comp%el(i_el)%stripe_1%p => null()
          n_s = n_s + 1
        else
          comp%el(i_el)%stripe_1%p => comp%el(i_el)%neigh(1)%p
        end if
      
        comp%el(i_el)%dy = sum( cross( comp%el(i_el)%edge_uni(:,2) ,    &
                                      comp%el(i_el)%edge_vec(:,3)  ) * &
                                comp%el(i_el)%nor )
      
        h2 = sum( cross( comp%el(i_el)%edge_uni(:,2) ,                        &
                        comp%el(i_el)%ver(:,1) - comp%el(i_el)%ver(:,3) ) * &
                  comp%el(i_el)%nor )
      
        if ( abs( h2 - comp%el(i_el)%dy ) .gt. 1.0e-4_wp .and. &
            sim_param%debug_level .ge. 7) then
          call warning(this_sub_name, this_mod_name, 'non parallel stripe edges (rough check)')
        end if
      
      end do
    
      !> Allocate and fill comp%strip_elem array
      if ( mod(n_el, n_s) .ne. 0 ) then
        call error(this_sub_name, this_mod_name, ' The number of elements of&
            & a parametric element is not an exact multiple of the number of&
            & spanwise stripes. There is something wrong in the geometry input&
            & file')
      end if
      
      n_c = n_el / n_s  ! integer division to find number of chord panels
      
      comp%n_s = n_s
      comp%n_c = n_c  

      allocate(comp%stripe(n_s))

      do i_s = 1, n_s

        allocate(comp%stripe(i_s)%panels(n_c))
        comp%stripe(i_s)%area = 0.0_wp
        allocate(comp%stripe(i_s)%ver(3,4))
        comp%stripe(i_s)%ver = 0.0_wp 
        allocate(comp%stripe(i_s)%edge_vec(3,4))
        comp%stripe(i_s)%edge_vec = 0.0_wp
        allocate(comp%stripe(i_s)%edge_len(4))
        comp%stripe(i_s)%edge_len = 0.0_wp
        allocate(comp%stripe(i_s)%edge_uni(3,4))
        comp%stripe(i_s)%edge_uni = 0.0_wp
        
        
        do i_c = 1, n_c
          comp%stripe(i_s)%panels(i_c)%p => comp%el(i_c+(i_s-1)*n_c)
          comp%stripe(i_s)%area = comp%stripe(i_s)%area + comp%el(i_c+(i_s-1)*n_c)%area
        end do   
          
          if (trim(comp%comp_el_type) .eq. 'v' .and. &
              trim(comp%aero_correction) .eq. 'true') then

            comp%stripe(i_s)%n_ver = 4 ! hardcoded, but stripe should be always a quadrilateral element
            comp%stripe(i_s)%csi_cen = 0.5_wp * sum(comp%normalised_coord_e(:,i_s)) 
            comp%stripe(i_s)%i_airfoil =  comp%i_airfoil_e(:,i_s)
            !> stripe verteces 
            comp%stripe(i_s)%ver(:,1) = comp%el(1+(i_s-1)*n_c)%ver(:,1) 
            comp%stripe(i_s)%ver(:,2) = comp%el(1+(i_s-1)*n_c)%ver(:,2)
            comp%stripe(i_s)%ver(:,3) = comp%el(n_c+(i_s-1)*n_c)%ver(:,3)
            comp%stripe(i_s)%ver(:,4) = comp%el(n_c+(i_s-1)*n_c)%ver(:,4) 

            ! > get the control point as the point at 0.5 c of the chord (!)
            ! It is updated at every timestep to account for curvature
            comp%stripe(i_s)%cen = sum(comp%stripe(i_s)%ver(:,1:2),2)/2.0_wp 
            
            nor = cross(comp%stripe(i_s)%ver(:,3) - comp%stripe(i_s)%ver(:,1),&
                        comp%stripe(i_s)%ver(:,4) - comp%stripe(i_s)%ver(:,2))
            
            comp%stripe(i_s)%nor = nor/norm2(nor) 
            
            tanl = 0.5_wp * ( comp%stripe(i_s)%ver(:,4) + comp%stripe(i_s)%ver(:,1) ) - &
                            comp%stripe(i_s)%cen
            
            comp%stripe(i_s)%tang(:,1) = tanl / norm2(tanl)
            comp%stripe(i_s)%tang(:,2) = cross(comp%stripe(i_s)%nor, comp%stripe(i_s)%tang(:,1))
            
            !> Vector connecting two consecutive vertices:
            do i = 1, 4
              comp%stripe(i_s)%edge_vec(:,i) = comp%stripe(i_s)%ver(:,next_qua(i)) - comp%stripe(i_s)%ver(:,i)
            end do

            !> Edge: edge_len(:)
            do i = 1, 4
              comp%stripe(i_s)%edge_len(i) = norm2(comp%stripe(i_s)%edge_vec(:,i))
            end do
          
            !> Unit vector
            do i = 1,4
              comp%stripe(i_s)%edge_uni(:,i) = comp%stripe(i_s)%edge_vec(:,i) / comp%stripe(i_s)%edge_len(i)
            end do

            !> Area 
            n_pan = size(comp%stripe(i_s)%panels)
            comp%stripe(i_s)%area = 0.0_wp 
            do i = 1, n_pan
              comp%stripe(i_s)%area = comp%stripe(i_s)%area + comp%stripe(i_s)%panels(i)%p%area 
            end do

            comp%stripe(i_s)%tang_cen = comp%stripe(i_s)%edge_uni(:,2) - comp%stripe(i_s)%edge_uni(:,4)
            comp%stripe(i_s)%tang_cen = comp%stripe(i_s)%tang_cen / norm2(comp%stripe(i_s)%tang_cen)
  
            comp%stripe(i_s)%bnorm_cen = cross(comp%stripe(i_s)%tang_cen, comp%stripe(i_s)%nor)
            comp%stripe(i_s)%bnorm_cen = comp%stripe(i_s)%bnorm_cen / norm2(comp%stripe(i_s)%bnorm_cen)
            comp%stripe(i_s)%chord = sum(comp%stripe(i_s)%edge_len((/2,4/)))*0.5_wp
            
            comp%stripe(i_s)%thickness = 0.12*comp%stripe(i_s)%chord !> workaround to make it compile 
            
            comp%stripe(i_s)%ctr_pt = comp%stripe(i_s)%cen +  & 
                                      comp%stripe(i_s)%tang_cen * comp%stripe(i_s)%chord / 2.0_wp
            
            !> Update velocity 
            call calc_node_vel(comp%stripe(i_s)%ctr_pt, geo%refs(comp%ref_id)%G_g, &
                                geo%refs(comp%ref_id)%f_g, comp%stripe(i_s)%ub)
            
          endif    
          
      end do 

      do i_s = 1 , n_s
        do i_c = 1 , n_c
          allocate( comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c) )
          do i_c2 = 1 , i_c
            comp%el(i_c+(i_s - 1)*n_c)%stripe_elem(i_c2)%p => &
                                          comp%el(i_c2+(i_s-1)*n_c)
            comp%el(i_c+(i_s - 1)*n_c)%n_c = n_c
            comp%el(i_c+(i_s - 1)*n_c)%n_s = n_s
          end do
        end do
        
      end do
    
    end if    
    end associate
  end do
  
  

end subroutine create_strip_connectivity


function move_points(pp, R, of)  result(rot_pp)
  real(wp), intent(in)    :: pp(:,:)
  real(wp), intent(in)    :: R(:,:)
  real(wp), intent(in)    :: of(:)
  real(wp)                :: rot_pp(size(pp,1),size(pp,2))

  rot_pp = matmul(R,pp)
  rot_pp(1,:) = rot_pp(1,:) + of(1)
  rot_pp(2,:) = rot_pp(2,:) + of(2)
  rot_pp(3,:) = rot_pp(3,:) + of(3)

end function move_points


!> Update all the geometry
!!
!! The geometry is updated at the beginning of the simulation, or after
!! a movement of the geometry. At the beginning the geometrical data of
!! all the elements is calculated, after only the moving elements are
!! updated.
!!
!! Also the velocity of the centerpoint of the elements is calculated
subroutine update_geometry(geo, te, t, update_static, time_cycle)
  type(t_geo), intent(inout) :: geo
  type(t_tedge), intent(inout) ::te
  real(wp), intent(in) :: t
  logical, intent(in) :: update_static
  logical, intent(in) :: time_cycle

  real(wp), allocatable :: rr_hinge_contig(:,:)
  real(wp), allocatable :: rr_hinge_contig_virtual(:,:)
  integer :: i_comp, ie, ih, i_el

  !> Update all the references
  call update_all_references(geo%refs,t)

  do i_comp = 1,size(geo%components)
    associate(comp => geo%components(i_comp))
      
      !> Update only rigid components, or update at first timestep
      ! *** to do *** quite a dirty implementation
      if ( ( .not. comp%coupling ) .or. update_static ) then
        
        if (comp%moving .or. update_static) then
          
          !> store %nor at previous time step, for moving, used few lines
          ! below to evaluate unit normal time derivative dn_dt
          if ( .not. update_static ) then
            do ie = 1 , size(comp%el)
              comp%el(ie)%nor_old = comp%el(ie)%nor
            end do
            do ie = 1 , size(comp%el_virtual)
              comp%el_virtual(ie)%nor_old = comp%el_virtual(ie)%nor
            end do
          end if

          !> Move points of the component
          geo%points(:,comp%i_points) = move_points(comp%loc_points, &
                                        geo%refs(comp%ref_id)%R_g, &
                                        geo%refs(comp%ref_id)%of_g)
          geo%points_virtual(:,comp%i_points_virtual) = move_points(comp%loc_points_virtual, &
                                        geo%refs(comp%ref_id)%R_g, &
                                        geo%refs(comp%ref_id)%of_g)                                       

          !> Update geometrical data of the elements, after geometry motion
          do ie = 1,size(comp%el)
            call comp%el(ie)%calc_geo_data(geo%points(:,comp%el(ie)%i_ver))
          enddo
          do ie = 1,size(comp%el_virtual)
            call comp%el_virtual(ie)%calc_geo_data_virtual(geo%points_virtual(:,comp%el_virtual(ie)%i_ver))
          enddo

          !> Unit normal vector time derivative, dn_dt
          ! set %nor_old = %nor for static elements and first time step,
          ! this if statement completes the if statement found few lines
          ! above, avoiding initialization issues (dirty implementation?)
          if ( update_static ) then
            do ie = 1,size(comp%el)
              comp%el(ie)%nor_old = comp%el(ie)%nor
            end do
            do ie = 1,size(comp%el_virtual)
              comp%el_virtual(ie)%nor_old = comp%el_virtual(ie)%nor
            end do            
          end if

          do ie = 1 , size(comp%el)
            comp%el(ie)%dn_dt = (comp%el(ie)%nor - comp%el(ie)%nor_old)/sim_param % dt
          end do
          do ie = 1 , size(comp%el)
            comp%el_virtual(ie)%dn_dt = (comp%el_virtual(ie)%nor - comp%el_virtual(ie)%nor_old)/sim_param % dt
          end do          

          !> Calculate the velocity of the centers to assing b.c.
          do ie = 1,size(comp%el)
            call calc_geo_vel(comp%el(ie), geo%refs(comp%ref_id)%G_g, &
                                          geo%refs(comp%ref_id)%f_g)
          enddo
          do ie = 1,size(comp%el_virtual)
            call calc_geo_vel_virtual(comp%el_virtual(ie), geo%refs(comp%ref_id)%G_g, &
                                          geo%refs(comp%ref_id)%f_g)
          enddo
          !> Update stripe velocity 
          if (time_cycle .and. & 
              trim(comp%comp_el_type) .eq. 'v' .and. &
              trim(comp%aero_correction) .eq. 'true') then
              
            do ie = 1, size(comp%stripe)              
              call calc_node_vel(comp%stripe(ie)%ctr_pt, geo%refs(comp%ref_id)%G_g, &
                                geo%refs(comp%ref_id)%f_g, comp%stripe(ie)%ub)
            end do  
          endif 
        end if

      elseif ( comp%coupling .or. update_static ) then

        do ie = 1,size(comp%el)

          comp%el(ie)%dn_dt = (comp%el(ie)%nor - comp%el(ie)%nor_old)/sim_param%dt
        end do
        do ie = 1,size(comp%el_virtual)

          comp%el_virtual(ie)%dn_dt = (comp%el_virtual(ie)%nor - comp%el_virtual(ie)%nor_old)/sim_param%dt
        end do

      end if  ! if ( .not. comp%coupling )

    !> Hinges
    !> Get trailing edge direction to be rotated for hinge deflection
      
    do ih = 1, comp%n_hinges
      !> Update:
      ! - hinge node coordinates, act%rr
      ! - hinge node orientation (unit vectors h,v,n)
      ! - hinge rotation angle, theta
      ! for non-coupled components only (so far)
      if ( (.not. comp%coupling) ) then
        
        !> hinge nodes, points and orientation
        call comp%hinge(ih)%update_hinge_nodes( geo%refs(comp%ref_id)%R_g, &
                                                geo%refs(comp%ref_id)%of_g )
        !> hinge rotation, theta
        call comp%hinge(ih)%update_theta( t )
        
        !> Allocating contiguous array to pass to %hinge_deflection procedure
        allocate(rr_hinge_contig(3,size(comp%i_points)))
        rr_hinge_contig = geo%points(:, comp%i_points)
        
        if (comp%moving .or. update_static) then 
          call comp%hinge(ih)%hinge_deflection(comp%i_points, rr_hinge_contig,  t, te%i, te%t_hinged )
        endif
        
        geo%points(:, comp%i_points) = rr_hinge_contig
        
        deallocate(rr_hinge_contig)

      else

        !> Updated in mod_precice/update_elems

      end if


    end do

    do ih = 1, comp%n_hinges
      !> Update:
      ! - hinge node coordinates, act%rr
      ! - hinge node orientation (unit vectors h,v,n)
      ! - hinge rotation angle, theta
      ! for non-coupled components only (so far)
      if ( (.not. comp%coupling) ) then
        
        !> hinge nodes, points and orientation
        call comp%hinge(ih)%update_hinge_nodes( geo%refs(comp%ref_id)%R_g, &
                                                geo%refs(comp%ref_id)%of_g )
        !> hinge rotation, theta
        call comp%hinge(ih)%update_theta( t )
        
        !> Allocating contiguous array to pass to %hinge_deflection procedure
        allocate(rr_hinge_contig_virtual(3,size(comp%i_points_virtual)))
        rr_hinge_contig_virtual = geo%points_virtual(:, comp%i_points_virtual)
        
        if (comp%moving .or. update_static) then 
          call comp%hinge(ih)%hinge_deflection(comp%i_points_virtual, rr_hinge_contig_virtual,  t, te%i, te%t_hinged )
        endif
        
        geo%points_virtual(:, comp%i_points_virtual) = rr_hinge_contig_virtual
        
        deallocate(rr_hinge_contig_virtual)

      else

        !> Updated in mod_precice/update_elems

      end if


    end do

    ! *** to do ***
    ! Update surface velocity, considering both hinges with prescribed motion
    ! and those coupled with the structural components
    ! *** to do ***
    if ( comp%n_hinges .gt. 0 ) then
      !> Then update geometrical data, position and velocity of the elem centers
      do ie = 1,size(comp%el)

        !> Save dcen_h at previous timestep
        comp%el(ie)%dcen_h_old = comp%el(ie)%dcen_h

        !> Save the position of the centre, before hinge rotation
        comp%el(ie)%dcen_h = comp%el(ie)%cen

        !> Update new position of the centers, taking into account hinge motions
        call comp%el(ie)%calc_geo_data(geo%points(:,comp%el(ie)%i_ver))

        !> Evaluate dcen_h, delta position due to hinge motion
        comp%el(ie)%dcen_h = comp%el(ie)%cen - comp%el(ie)%dcen_h

        !> Evaluate dvel_h, delta velocity due to hinge motion, with first
        ! order finite difference
        comp%el(ie)%dvel_h = ( comp%el(ie)%dcen_h - comp%el(ie)%dcen_h_old ) / &
                                                                  sim_param % dt
        !> Update on-body velocity
        comp%el(ie)%ub = comp%el(ie)%ub + comp%el(ie)%dvel_h

        !> Update time derivative of the unit normal vector
        comp%el(ie)%dn_dt = comp%el(ie)%dn_dt + &
                            ( comp%el(ie)%nor - comp%el(ie)%nor_old ) / &
                                                                sim_param % dt
      end do
      do ie = 1,size(comp%el_virtual)

        !> Save dcen_h at previous timestep
        comp%el_virtual(ie)%dcen_h_old = comp%el_virtual(ie)%dcen_h

        !> Save the position of the centre, before hinge rotation
        comp%el_virtual(ie)%dcen_h = comp%el_virtual(ie)%cen

        !> Update new position of the centers, taking into account hinge motions
        call comp%el_virtual(ie)%calc_geo_data_virtual(geo%points_virtual(:,comp%el_virtual(ie)%i_ver))

        !> Evaluate dcen_h, delta position due to hinge motion
        comp%el_virtual(ie)%dcen_h = comp%el_virtual(ie)%cen - comp%el_virtual(ie)%dcen_h

        !> Evaluate dvel_h, delta velocity due to hinge motion, with first
        ! order finite difference
        comp%el_virtual(ie)%dvel_h = ( comp%el_virtual(ie)%dcen_h - comp%el_virtual(ie)%dcen_h_old ) / &
                                                                  sim_param % dt
        !> Update on-body velocity
        comp%el_virtual(ie)%ub = comp%el_virtual(ie)%ub + comp%el_virtual(ie)%dvel_h

        !> Update time derivative of the unit normal vector
        comp%el_virtual(ie)%dn_dt = comp%el_virtual(ie)%dn_dt + &
                            ( comp%el_virtual(ie)%nor - comp%el_virtual(ie)%nor_old ) / &
                                                                sim_param % dt
      end do

    end if

    end associate
  enddo


end subroutine update_geometry

!----------------------------------------------------------------------

!> Destroy the pointers array of the components
!!
!! This subroutine is needed just for compatibility with gcc 4.8
!! Newer versions of the compiler and ifort deallocate also the
!! vector of pointers alongside the geometry in destroy_geometry,
!! however gcc 4.8 gives a SIGSEGV and needs and explicit deallocation
subroutine destroy_elements_geo(geo)
  type(t_geo), intent(inout) :: geo
  integer                    :: ic

  do ic=1,size(geo%components)
    if(allocated(geo%components(ic)%el)) deallocate(geo%components(ic)%el)
  enddo

end subroutine destroy_elements_geo

!----------------------------------------------------------------------

!> Destroy the pointers array of the components
!!
!! This subroutine is needed just for compatibility with gcc 4.8
!! Newer versions of the compiler and ifort deallocate also the
!! vector of pointers alongside the geometry in destroy_geometry,
!! however gcc 4.8 gives a SIGSEGV and needs and explicit deallocation
subroutine destroy_elements_comps(components)
  type(t_geo_component), intent(inout) :: components(:)
  integer                              :: ic

  do ic=1,size(components)
    if(allocated(components(ic)%el)) deallocate(components(ic)%el)
  enddo

end subroutine destroy_elements_comps

!> Destroy all the contents of a geometry type
!
subroutine destroy_geometry(geo, elems)
  type(t_geo), intent(out)                     :: geo
  type(t_pot_elem_p), allocatable, intent(out) :: elems(:)
  integer                                      :: i
  !call destroy_references(geo%refs)

  !> Dummy to avoid warnings
  i = size(elems)

end subroutine destroy_geometry

end module mod_geometry
