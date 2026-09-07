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


!> This is a module dedicated to the handling of the geometry during the
!! postprocessing
!!
!! It uses the same data structures of the postprocessing as well as
!! performing similar tasks, but much less cumbersome

module mod_geo_postpro

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

use mod_sim_param, only: &
  t_sim_param

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption &
  , finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

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

use mod_liftlin, only: &
  t_liftlin

use mod_virtual, only: &
  c_elem_virtual, t_elem_virtual_p

use mod_actuatordisk, only: &
  t_actdisk

use mod_c81, only: &
  t_aero_tab , read_c81_table , interp_aero_coeff

use mod_math, only: &
  cross

use mod_reference, only: &
  t_ref, build_references, update_all_references, destroy_references

use mod_hdf5_io, only: &
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

use mod_hinges, only: &
  initialize_hinge_config

use mod_geometry, only: &
  t_geo, t_geo_component , calc_geo_data_pan , calc_geo_vel, calc_geo_vel_virtual

use mod_stringtools, only: &
  LowCase, IsInList, strip_mult_appendix

!----------------------------------------------------------------------

implicit none

public :: load_components_postpro, update_points_postpro , &
          prepare_geometry_postpro, expand_actdisk_postpro

private

!----------------------------------------------------------------------

character(len=*), parameter :: this_mod_name = 'mod_geo_postpro'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
                                    components_names, all_comp)
  type(t_geo_component), allocatable, intent(inout) :: comps(:)
  real(wp), allocatable, intent(out)                :: points(:,:), points_virtual(:,:)
  integer, intent(out)                              :: nelem, nelem_virtual
  integer(h5loc), intent(in)                        :: floc
  character(len=*), allocatable, intent(inout)      :: components_names(:)
  logical, intent(in)                               :: all_comp

  type(t_geo_component), allocatable                :: comp_temp(:)
  integer                                           :: i1 , i2, i3
  integer, allocatable                              :: ee(:,:), ee_virtual(:,:)
  real(wp), allocatable                             :: rr(:,:), rr_virtual(:,:)
  real(wp), allocatable                             :: ori(:,:)
  character(len=max_char_len)                       :: comp_el_type, comp_name
  character(len=max_char_len)                       :: comp_input, comp_name_stripped
  integer                                           :: points_offset, points_offset_virtual, n_vert
  real(wp), allocatable                             :: points_tmp(:,:)
  character(len=max_char_len)                       :: ref_tag
  integer                                           :: ref_id
  character(len=max_char_len)                       :: cname 
  integer(h5loc)                                    :: gloc, cloc, geo_loc
  integer                                           :: n_comp, i_comp, n_comp_tot 
  integer                                           :: i_comp_tot, i_comp_tmp
  integer                                           :: parametric_nelems_span
  integer                                           :: parametric_nelems_chor, parametric_nelems_chor_virtual
  real(wp)                                          :: coupling_node_rot(3,3) = 0.0_wp
  
  !> Hinges
  integer(h5loc)                                    :: hiloc, hloc
  integer                                           :: ih, n_hinges
  real(wp)                                          :: rotation_amplitude
  character(len=2)                                  :: hinge_id_str

  !> vl corrected 
  character(len=5)                                  :: aero_table

  !> Some structure to handle multiple components
  character(len=max_char_len), allocatable          :: components(:) , components_tmp(:)
  character(len=max_char_len)                       :: component_stripped

  character(len=max_char_len)                       :: comp_coupling_str

  character(len=*), parameter                       :: this_sub_name = 'load_components_postpro'

  ! Read all the components
  call open_hdf5_group(floc,'Components', gloc)
  call read_hdf5(n_comp_tot,'NComponents',gloc)

  allocate(components(n_comp_tot))
  do i_comp_tot = 1 , n_comp_tot
    write(cname,'(A,I3.3)') 'Comp',i_comp_tot
    call open_hdf5_group(gloc, trim(cname),cloc)
    call read_hdf5(components(i_comp_tot),'CompName',cloc)
    call read_hdf5(comp_input,'CompInput',cloc)
    call close_hdf5_group(cloc)
  end do
  

! RE-BUILD components_names() to host multiple components ++++++++++++++++++
! components_names: input for post-processing
! comp_name: all the components of the model (blades: Hub01__01, ..., Hub01__Nb)
! comp_name_stripped: Hub01__01, ..., Hub01__Nb  ---> Hub01
! multiple components: if components_names(i1) is <Hub>, then add all the blades to the output

! TODO: check user inputs ( to avoid rotorll and rotorll__01 to be considered twice )
! TODO: check multiple components ( double if statementes ... )
  allocate(components_tmp(n_comp_tot))
  i_comp_tmp = 0
  if ( allocated(components_names) ) then
    do i1 = 1 , size(components_names)

      do i2 = 1 , n_comp_tot

        call strip_mult_appendix(components(i2), component_stripped, '__')


        ! CASE #1. Ex.: components_names(i1) .eq. rotorll
        if ( trim(components_names(i1)) .eq. trim(component_stripped) ) then
          i_comp_tmp = i_comp_tmp + 1
          components_tmp(i_comp_tmp) = trim(components(i2))

        else

! **** Encapsulated if stat **** for debug
          if ( trim(components_names(i1)) .eq. trim(components(i2)) ) then
            i_comp_tmp = i_comp_tmp + 1
            components_tmp(i_comp_tmp) = trim(components(i2))
          end if

        end if

      end do

    end do

    deallocate(components_names)
    allocate(components_names(i_comp_tmp))
    components_names = components_tmp(1:i_comp_tmp)

  ! If components_names was not allocated, all the components are required
  else ! ( .not. allocated(components_names) ) then

    allocate(components_names(n_comp_tot))
    do i_comp_tot = 1 , n_comp_tot
      write(cname,'(A,I3.3)') 'Comp',i_comp_tot
      call open_hdf5_group(gloc,trim(cname),cloc)
      call read_hdf5(components_names(i_comp_tot),'CompName',cloc)
      call close_hdf5_group(cloc)
    end do

  end if

  !TODO: check if this is necessary but not enough
  ! to avoid Components that do not exist as an input in dust_post.in
  if ( .not. allocated(components_names) ) then
    call error(this_sub_name, this_mod_name, &
                  'components_names .not. allocated. Something strange &
                  &and unexpected happened. It could be a bug')
  end if
  if ( size(components_names) .le. 0 ) then
    call error(this_sub_name, this_mod_name, &
                  'No component found corresponding to the requested one(s) &
                  &for this analysis')
  end if

!  allocate(comps(n_comp))
  nelem = 0; nelem_virtual = 0
  i_comp = 0; n_comp = 0
  do i_comp_tot = 1,n_comp_tot

    write(cname,'(A,I3.3)') 'Comp',i_comp_tot
    call open_hdf5_group(gloc,trim(cname),cloc)

    call read_hdf5(comp_name,'CompName',cloc)
    call read_hdf5(comp_input,'CompInput',cloc)

    !Strip the appendix of multiple components, to load all the multiple
    !components at once
    call strip_mult_appendix(comp_name, comp_name_stripped, '__')



    if(IsInList(comp_name, components_names) .or. all_comp) then

      i_comp = i_comp+1; n_comp = n_comp+1
      allocate(comp_temp(n_comp))
      if(allocated(comps)) comp_temp(1:n_comp-1) = comps
      call move_alloc(comp_temp, comps)

      comps(i_comp)%comp_id = i_comp_tot
      call read_hdf5(ref_id,'RefId',cloc)
      call read_hdf5(ref_tag,'RefTag',cloc)

      comps(i_comp)%ref_id  = ref_id
      comps(i_comp)%ref_tag = trim(ref_tag)

      ! ====== READING =====
      call read_hdf5(comp_el_type,'ElType',cloc)
      comps(i_comp)%comp_el_type = trim(comp_el_type)

      comps(i_comp)%comp_name = trim(comp_name)
      comps(i_comp)%comp_input = trim(comp_input)

      call read_hdf5(comp_coupling_str,'Coupled',cloc)
      if ( trim(comp_coupling_str) .eq. 'true' ) then
        comps(i_comp)%coupling = .true.
        call read_hdf5(coupling_node_rot,'CouplingNodeOrientation',cloc)
        comps(i_comp)%coupling_node_rot = coupling_node_rot
      else
        comps(i_comp)%coupling = .false.
        comps(i_comp)%coupling_node_rot(1,:) = (/ 1._wp, 0._wp, 0._wp/)
        comps(i_comp)%coupling_node_rot(2,:) = (/ 0._wp, 1._wp, 0._wp/)
        comps(i_comp)%coupling_node_rot(3,:) = (/ 0._wp, 0._wp, 1._wp/)
      end if

      ! Geometry --------------------------
      call open_hdf5_group(cloc,'Geometry',geo_loc)
      call read_hdf5_al(ee   ,'ee'   ,geo_loc)
      call read_hdf5_al(rr   ,'rr'   ,geo_loc)
      call read_hdf5_al(ee_virtual   ,'ee_virtual'   ,geo_loc)
      call read_hdf5_al(rr_virtual   ,'rr_virtual'   ,geo_loc)
      if ( trim(comps(i_comp)%comp_input) .eq. 'parametric' .or. & 
          trim(comps(i_comp)%comp_input) .eq. 'pointwise') then
        call read_hdf5( parametric_nelems_span ,'parametric_nelems_span',geo_loc)
        call read_hdf5( parametric_nelems_chor ,'parametric_nelems_chor',geo_loc)
        call read_hdf5( parametric_nelems_chor_virtual ,'parametric_nelems_chor_virtual',geo_loc)
        comps(i_comp)%parametric_nelems_span = parametric_nelems_span
        comps(i_comp)%parametric_nelems_chor = parametric_nelems_chor
        comps(i_comp)%parametric_nelems_chor_virtual = parametric_nelems_chor_virtual
      end if

      if ( trim(comps(i_comp)%comp_el_type) .eq. 'v' ) then
        call read_hdf5(aero_table,  'aero_table', geo_loc)
        comps(i_comp)%aero_correction = trim(aero_table)
      end if 
  
      call close_hdf5_group(geo_loc)

      ! ======= CREATING ELEMENTS ======

      ! --- treat the points ---
      if(allocated(points)) then
        points_offset = size(points,2)
      else
        points_offset = 0
      endif
      if(allocated(points_virtual)) then
        points_offset_virtual = size(points_virtual,2)
      else
        points_offset_virtual = 0
      endif

      !store the read points into the local points
      allocate(comps(i_comp)%loc_points(3,size(rr,2)))
      comps(i_comp)%loc_points = rr
      allocate(comps(i_comp)%loc_points_virtual(3,size(rr_virtual,2)))
      comps(i_comp)%loc_points_virtual = rr_virtual
      !Now for the moments the points are stored here without moving them,
      !will be moved later, consider not storing them here at all
      allocate(points_tmp(3,size(rr,2)+points_offset))
      if (points_offset .gt. 0) points_tmp(:,1:points_offset) = points
      points_tmp(:,points_offset+1:points_offset+size(rr,2)) = rr
      call move_alloc(points_tmp, points)
      allocate(comps(i_comp)%i_points(size(rr,2)))
      comps(i_comp)%i_points = (/((i3),i3=points_offset+1,points_offset+size(rr,2))/)

      allocate(points_tmp(3,size(rr_virtual,2)+points_offset_virtual))
      if (points_offset_virtual .gt. 0) points_tmp(:,1:points_offset_virtual) = points_virtual
      points_tmp(:,points_offset_virtual+1:points_offset_virtual+size(rr_virtual,2)) = rr_virtual
      call move_alloc(points_tmp, points_virtual)
      allocate(comps(i_comp)%i_points_virtual(size(rr_virtual,2)))
      comps(i_comp)%i_points_virtual = (/((i3),i3=points_offset_virtual+1,points_offset_virtual+size(rr_virtual,2))/)

      ! --- treat the elements ---

      !allocate the elements of the component of the right kind
      comps(i_comp)%nelems = size(ee,2)
      comps(i_comp)%nelems_virtual = size(ee_virtual,2)
      nelem = nelem + comps(i_comp)%nelems
      nelem_virtual = nelem_virtual + comps(i_comp)%nelems_virtual
      select case(trim(comps(i_comp)%comp_el_type))
        case('p')
          allocate(t_surfpan::comps(i_comp)%el(size(ee,2)))
        case('v')
          allocate(t_vortlatt::comps(i_comp)%el(size(ee,2)))
        case('l')
          allocate(t_liftlin::comps(i_comp)%el(size(ee,2)))
        case('a')
          allocate(t_actdisk::comps(i_comp)%el(size(ee,2)))
        case default
          call error(this_sub_name, this_mod_name, &
                  'Unknown type of element: '//comps(i_comp)%comp_el_type)
      end select
      allocate(comps(i_comp)%el_virtual(size(ee_virtual,2)))

      !> fill (some) of the elements fields
      do i2=1,size(ee,2)

        !> vertices
        n_vert = count(ee(:,i2).ne.0)
        allocate(comps(i_comp)%el(i2)%i_ver(n_vert))
        comps(i_comp)%el(i2)%n_ver = n_vert
        comps(i_comp)%el(i2)%i_ver(1:n_vert) = &
                                              ee(1:n_vert,i2) + points_offset
      enddo

      n_vert = 0

      !> fill (some) of the virtual elements fields
      do i2=1,size(ee_virtual,2)

        !> vertices
        n_vert = count(ee_virtual(:,i2).ne.0)
        allocate(comps(i_comp)%el_virtual(i2)%i_ver(n_vert))
        comps(i_comp)%el_virtual(i2)%n_ver = n_vert
        comps(i_comp)%el_virtual(i2)%i_ver(1:n_vert) = &
                                              ee_virtual(1:n_vert,i2) + points_offset_virtual
      enddo


      ! Update elems_offset for the next component
      !elems_offset = elems_offset + size(ee,2)

      !> Hinges ---
      call open_hdf5_group(cloc,'Hinges', hloc)
      call read_hdf5(n_hinges, 'n_hinges', hloc)
      comps(i_comp)%n_hinges = n_hinges
      allocate( comps(i_comp)%hinge(n_hinges) )
      do ih = 1, n_hinges

        !> Open hinge group
        write(hinge_id_str,'(I2.2)') ih
        call open_hdf5_group(hloc, 'Hinge_'//hinge_id_str, hiloc)

        !> read input and fill component%hinge fields
        call read_hdf5( comps(i_comp)%hinge(ih)%nodes_input, &
                                              'Nodes_Input', hiloc)
        call read_hdf5( comps(i_comp)%hinge(ih)%offset, &
                                              'Offset', hiloc)
        call read_hdf5( comps(i_comp)%hinge(ih)%span_blending, &
                                              'Spanwise_Blending', hiloc)
        call read_hdf5( comps(i_comp)%hinge(ih)%ref_dir, &
                                              'Ref_Dir', hiloc)
        call read_hdf5( comps(i_comp)%hinge(ih)%tag, &
                                              'Tag', hiloc)
        call read_hdf5_al( comps(i_comp)%hinge(ih)%ref%rr, 'rr', hiloc)
        comps(i_comp)%hinge(ih)%n_nodes = size( &
          comps(i_comp)%hinge(ih)%ref%rr, 2)

        call read_hdf5( comps(i_comp)%hinge(ih)%input_type, &
                                            'Hinge_Rotation_Input'    , hiloc)
        call read_hdf5( rotation_amplitude, 'Hinge_Rotation_Amplitude', hiloc)
        allocate( comps(i_comp)%hinge(ih)%theta( &
                  comps(i_comp)%hinge(ih)%n_nodes ) )
        allocate( comps(i_comp)%hinge(ih)%theta_old( &
                  comps(i_comp)%hinge(ih)%n_nodes ) )
        comps(i_comp)%hinge(ih)%theta     = rotation_amplitude
        comps(i_comp)%hinge(ih)%theta_old = 0.0_wp

        call close_hdf5_group( hiloc )

        !> Initialize reference and actual configuration
        call initialize_hinge_config( comps(i_comp)%hinge(ih)%ref , &
                                      comps(i_comp)%hinge(ih) )
        call initialize_hinge_config( comps(i_comp)%hinge(ih)%act , &
                                      comps(i_comp)%hinge(ih) )
        ! Allocate and initialize hinge node coords in the actual configuration
        allocate( comps(i_comp)%hinge(ih)%act%rr( 3,comps(i_comp)%hinge(ih)%n_nodes ) )
        
        comps(i_comp)%hinge(ih)%act%rr = comps(i_comp)%hinge(ih)%ref%rr
#if USE_PRECICE
        coupling_node_rot = comps(i_comp)%coupling_node_rot
#else
        coupling_node_rot(1,:) = (/ 1._wp, 0._wp, 0._wp/)
        coupling_node_rot(2,:) = (/ 0._wp, 1._wp, 0._wp/)
        coupling_node_rot(3,:) = (/ 0._wp, 0._wp, 1._wp/)
#endif
        call comps(i_comp)%hinge(ih)%build_connectivity( rr, coupling_node_rot)
        call comps(i_comp)%hinge(ih)%build_connectivity_cen( rr, ee, coupling_node_rot)

      end do

      call close_hdf5_group( hloc )


      !cleanup
      deallocate(ee, rr, ee_virtual, rr_virtual)

      if (allocated(ori)) then
        deallocate(ori) 
      endif 

    endif !load the element because in list

    call close_hdf5_group(cloc)

  enddo !i_comp
  call close_hdf5_group(gloc)

end subroutine load_components_postpro

!----------------------------------------------------------------------

!> Prepare the geometry allocating all the relevant fields for each
!! kind of element
subroutine prepare_geometry_postpro(comps)
  type(t_geo_component), intent(inout), target :: comps(:)

  integer                         :: i_comp, ie
  integer                         :: nsides
  class(c_pot_elem), pointer      :: elem
  class(c_elem_virtual), pointer  :: elem_virtual
  character(len=*), parameter     :: this_sub_name = 'prepare_geometry_postpro'


  do i_comp = 1,size(comps)
    do ie = 1,size(comps(i_comp)%el)
      elem => comps(i_comp)%el(ie)

      nsides = size(elem%i_ver)

      !Fields common to each element
      allocate(elem%ver(3,nsides))

      select type(elem)
        class is(t_surfpan)
          allocate(elem%verp(3,nsides))
          allocate(elem%edge_vec(3,nsides))
          allocate(elem%edge_len(nsides))
          allocate(elem%edge_uni(3,nsides))
          allocate(elem%cosTi(nsides))
          allocate(elem%sinTi(nsides))

        class is(t_vortlatt)
          allocate(elem%edge_vec(3,nsides))
          allocate(elem%edge_len(nsides))
          allocate(elem%edge_uni(3,nsides))

        class is(t_liftlin)
          allocate(elem%edge_vec(3,nsides))
          allocate(elem%edge_len(nsides))
          allocate(elem%edge_uni(3,nsides))

        class is(t_actdisk)
          allocate(elem%edge_vec(3,nsides))
          allocate(elem%edge_len(nsides))
          allocate(elem%edge_uni(3,nsides))

        class default
          call error(this_sub_name, this_mod_name, 'Unknown element type')
      end select

    enddo
    do ie = 1,size(comps(i_comp)%el_virtual)
      elem_virtual => comps(i_comp)%el_virtual(ie)
      nsides = size(elem_virtual%i_ver)
      allocate(elem_virtual%edge_vec(3,nsides))
      allocate(elem_virtual%edge_len(nsides))
      allocate(elem_virtual%edge_uni(3,nsides))
    enddo
  enddo

end subroutine prepare_geometry_postpro

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! panel element (vortex ring or surface panel)
subroutine calc_geo_data_postpro(elem,vert)
  class(c_pot_elem), intent(inout) :: elem
  real(wp), intent(in) :: vert(:,:)

  integer :: nsides, is
  real(wp):: nor(3), tanl(3)
  integer :: nxt

  nsides = size(vert,2)

  elem%n_ver = nsides

  ! vertices
  elem%ver = vert

  ! center
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)

  select type(elem)

    class default

      ! unit normal and area
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

    class is(t_actdisk)

      elem%area = 0.0_wp; elem%nor = 0.0_wp
      do is = 1, nsides
        nxt = 1+mod(is,nsides)
        nor = cross(vert(:,is) - elem%cen,&
                    vert(:,nxt) - elem%cen )
        elem%area = elem%area + 0.5_wp * norm2(nor)
        elem%nor = elem%nor + nor/norm2(nor)
      enddo
        elem%nor = elem%nor/real(nsides,wp)

      ! local tangent unit vector: aligned with first node, normal to n
      tanl = (vert(:,1)-elem%cen)-sum((vert(:,1)-elem%cen)*elem%nor)*elem%nor

      elem%tang(:,1) = tanl / norm2(tanl)
      elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )

  end select

  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  do is = 1 , nsides
    nxt = 1+mod(is,nsides)
    elem%edge_vec(:,is) = vert(:,nxt) - vert(:,is)
  end do

  ! edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  ! variables for surfpan elements only
  select type(elem)
    class is(t_surfpan)
      do is = 1 , nsides
        elem%verp(:,is) = vert(:,is) - elem%nor * &
                          sum( (vert(:,is) - elem%cen ) * elem%nor )
        elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) )
        elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) )
      end do

    class default

  end select


end subroutine calc_geo_data_postpro

!----------------------------------------------------------------------

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

!----------------------------------------------------------------------
subroutine update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                  refs_G , refs_f, filen)
  type(t_geo_component), intent(inout)  :: comps(:)
  real(wp), intent(inout)               :: points(:,:), points_virtual(:,:)
  real(wp), intent(in)                  :: refs_R(:,:,0:)
  real(wp), intent(in)                  :: refs_off(:,0:)
  real(wp), optional , intent(in)       :: refs_G(:,:,0:)
  real(wp), optional , intent(in)       :: refs_f(:,0:)

#if USE_PRECICE
  real(wp), allocatable                  :: rr(:,:), rr_virtual(:,:)
  real(wp), allocatable                  :: ori(:,:)
#endif
  character(len=*), optional, intent(in) :: filen
  character(max_char_len)                :: cname
  integer(h5loc)                         :: floc, gloc, cloc, rloc

  real(wp)                               :: time_todo = 0.0_wp  ! *** to do *** pass time as an input
  integer(h5loc)                         :: hiloc, hloc
  real(wp), allocatable                  :: rr_hinge_contig(:,:), rr_hinge_contig_virtual(:,:)
  integer                                :: i_comp, ie, ih
  character(len=2)                       :: hinge_id_str


  do i_comp = 1,size(comps)
    associate(comp => comps(i_comp))
#if USE_PRECICE
  if ( .not. comp%coupling ) then
#endif
    !> Move points of a rigid component, not coupled with an external software
    points(:,comp%i_points) = move_points(comp%loc_points, &
                              refs_R(:,:,comp%ref_id), &
                              refs_off(:,comp%ref_id))
    points_virtual(:,comp%i_points_virtual) = move_points(comp%loc_points_virtual, &
                              refs_R(:,:,comp%ref_id), &
                              refs_off(:,comp%ref_id))
#if USE_PRECICE
  else
    !> Read points of a coupled component, from result files

    !> Open result hdf5 file and read points coordinates
    call open_hdf5_file( trim(filen), floc )
    call open_hdf5_group( floc, 'Components', gloc )
    write(cname,'(A,I3.3)') 'Comp', comp%comp_id
    call open_hdf5_group( gloc, trim(cname), cloc )
    call open_hdf5_group( cloc, 'Geometry', rloc )
    call read_hdf5_al(rr,  'rr', rloc)
    call read_hdf5_al(rr_virtual,  'rr_virtual', rloc)
    call read_hdf5_al(ori, 'ori', rloc)
    
    points(:,comp%i_points) = rr
    points_virtual(:,comp%i_points_virtual) = rr_virtual
    do ie = 1, size(comp%el)
      comp%el(ie)%ori = ori(ie,:)
    enddo 
    !> Quite dirty: open and close for each component
    call close_hdf5_group(rloc)
    call close_hdf5_group(cloc)
    call close_hdf5_group(gloc)
    call close_hdf5_file(floc)
  endif
#endif


  !> Hinges 
  !> Re-open and close result hdf5 file and groups
  call open_hdf5_file ( trim(filen), floc )
  call open_hdf5_group( floc, 'Components', gloc )
  write(cname,'(A,I3.3)') 'Comp', comp%comp_id !i_comp
  call open_hdf5_group( gloc, trim(cname), cloc )
  call open_hdf5_group( cloc, 'Hinges', hloc )
  do ih = 1, comp%n_hinges

    !> Read hinge node position, orientation and theta
    write(hinge_id_str,'(I2.2)') ih
    call open_hdf5_group( hloc, 'Hinge_'//hinge_id_str, hiloc)
    call read_hdf5_al( comp%hinge(ih)%act%rr, 'act_rr', hiloc)
    call read_hdf5_al( comp%hinge(ih)%act%v , 'act_v ', hiloc)
    call read_hdf5_al( comp%hinge(ih)%act%h , 'act_h ', hiloc)
    call read_hdf5_al( comp%hinge(ih)%act%n , 'act_n ', hiloc)
    call read_hdf5_al( comp%hinge(ih)%theta    , 'theta'    , hiloc)
    call read_hdf5_al( comp%hinge(ih)%theta_old, 'theta_old', hiloc)

    !> Allocating contiguous array
    allocate(rr_hinge_contig(3,size(comp%i_points)))
    rr_hinge_contig = points(:, comp%i_points)
    call comp%hinge(ih)%hinge_deflection( comp%i_points, &
                        rr_hinge_contig, time_todo, postpro=.true. )
    points(:, comp%i_points) = rr_hinge_contig

    deallocate(rr_hinge_contig)

    !> Allocating virtual contiguous array
    allocate(rr_hinge_contig_virtual(3,size(comp%i_points_virtual)))
    rr_hinge_contig_virtual = points_virtual(:, comp%i_points_virtual)
    call comp%hinge(ih)%hinge_deflection( comp%i_points_virtual, &
                        rr_hinge_contig_virtual, time_todo, postpro=.true. )
    points_virtual(:, comp%i_points_virtual) = rr_hinge_contig_virtual

    deallocate(rr_hinge_contig_virtual)


    call close_hdf5_group(hiloc)

  end do

  call close_hdf5_group(hloc)
  call close_hdf5_group(cloc)
  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)
    !> Hinges 
    do ie = 1,size(comp%el)
      call calc_geo_data_postpro(comp%el(ie),points(:,comp%el(ie)%i_ver))

      if ( present(refs_G) .and. present(refs_f) ) then
        !> Calculate the velocity of the centers to impose the boundary condition
        call calc_geo_vel(comp%el(ie), refs_G(:,:,comp%ref_id) , &
                                        refs_f(:,comp%ref_id) )
      end if
    enddo

    do ie = 1,size(comp%el_virtual)
      call comp%el_virtual(ie)%calc_geo_data_virtual(points_virtual(:,comp%el_virtual(ie)%i_ver))

      if ( present(refs_G) .and. present(refs_f) ) then
        !> Calculate virtual panel center velocity for postprocessing output.
        call calc_geo_vel_virtual(comp%el_virtual(ie), refs_G(:,:,comp%ref_id), &
                                  refs_f(:,comp%ref_id))
      end if
    enddo

    end associate
  enddo

end subroutine update_points_postpro

!----------------------------------------------------------------------

subroutine expand_actdisk_postpro(comps, points, points_virtual, points_exp, & 
                                  points_virtual_exp, elems, elems_virtual)

  type(t_geo_component), intent(in)  :: comps(:)
  real(wp), intent(in)               :: points(:,:), points_virtual(:,:)
  real(wp), allocatable, intent(out) :: points_exp(:,:), points_virtual_exp(:,:)
  integer, allocatable, intent(out)  :: elems(:,:), elems_virtual(:,:)

  real(wp), allocatable              :: pt_tmp(:,:), pt_tmp_virtual(:,:)
  integer, allocatable               :: ee_tmp(:,:), ee_tmp_virtual(:,:)
  integer                            :: i_comp, ie, extra_offset, iv, ipt, ipt1
  integer                            :: start_pts, start_cen

  extra_offset = 0
  allocate(points_exp(3,0), elems(4,0))
  allocate(points_virtual_exp(3,0), elems_virtual(4,0))
  do i_comp = 1,size(comps)
    associate(cmp=>comps(i_comp))
    select type(el => cmp%el)
      type is(t_actdisk)
        !> make space also for the centers
        allocate(pt_tmp(3,size(points_exp,2) + &
                  size(cmp%loc_points,2) + cmp%nelems))
        pt_tmp(:,1:size(points_exp,2)) = points_exp
        start_pts = size(points_exp,2)
        start_cen = size(points_exp,2) + size(cmp%i_points)
        pt_tmp(:,start_pts+1 : start_cen) = points(:,cmp%i_points)
        do ie = 1,cmp%nelems
          pt_tmp(:,start_cen+ie) = el(ie)%cen
        enddo
        call move_alloc(pt_tmp, points_exp)
        extra_offset = extra_offset + cmp%nelems

        allocate(ee_tmp(4,size(elems,2)+size(cmp%i_points)))
        ee_tmp(:,1:size(elems,2)) = elems
        ee_tmp(:,size(elems,2)+1:size(ee_tmp,2)) = 0

        ipt = 1
        do ie = 1,cmp%nelems
          ipt1 = ipt
          do iv = 1,el(ie)%n_ver-1
            ee_tmp(1,size(elems,2)+ipt) = start_pts+ipt
            ee_tmp(2,size(elems,2)+ipt) = start_pts+ipt+1
            ee_tmp(3,size(elems,2)+ipt) = start_cen+ie
            ipt = ipt+1
          enddo
          !> last element
            ee_tmp(1,size(elems,2)+ipt) = start_pts+ipt
            ee_tmp(2,size(elems,2)+ipt) = start_pts+ipt1
            ee_tmp(3,size(elems,2)+ipt) = start_cen+ie
            ipt = ipt+1
        enddo
        call move_alloc(ee_tmp, elems)

      class default
        allocate(pt_tmp(3,size(points_exp,2)+size(cmp%loc_points,2)))
        pt_tmp(:,1:size(points_exp,2)) = points_exp
        pt_tmp(:,size(points_exp,2)+1:size(pt_tmp,2)) = points(:,cmp%i_points)
        call move_alloc(pt_tmp, points_exp)

        allocate(ee_tmp(4,size(elems,2)+cmp%nelems))
        ee_tmp(:,1:size(elems,2)) = elems
        ee_tmp(:,size(elems,2)+1:size(ee_tmp,2)) = 0
        do ie = 1,cmp%nelems
          ee_tmp(1:el(ie)%n_ver,size(elems,2)+ie) =  &
                                                  el(ie)%i_ver+extra_offset
        enddo
        call move_alloc(ee_tmp, elems)

        allocate(pt_tmp_virtual(3,size(points_virtual_exp,2)+size(cmp%loc_points_virtual,2)))
        pt_tmp_virtual(:,1:size(points_virtual_exp,2)) = points_virtual_exp
        pt_tmp_virtual(:,size(points_virtual_exp,2)+1:size(pt_tmp_virtual,2)) = points_virtual(:,cmp%i_points_virtual)
        call move_alloc(pt_tmp_virtual, points_virtual_exp)

        allocate(ee_tmp_virtual(4,size(elems_virtual,2)+cmp%nelems_virtual))
        ee_tmp_virtual(:,1:size(elems_virtual,2)) = elems_virtual
        ee_tmp_virtual(:,size(elems_virtual,2)+1:size(ee_tmp_virtual,2)) = 0
        do ie = 1,cmp%nelems_virtual
          ee_tmp_virtual(1:cmp%el_virtual(ie)%n_ver,size(elems_virtual,2)+ie) =  &
                                                  cmp%el_virtual(ie)%i_ver + extra_offset
        enddo
        call move_alloc(ee_tmp_virtual, elems_virtual)

      end select
    end associate
  enddo
end subroutine
!----------------------------------------------------------------------

end module mod_geo_postpro
