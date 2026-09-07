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
!!          Andrea Colli
!!          Alessandro Cocco
!!          Alberto Savino
!!          Federico Gentile
!!          Matteo Dall'Ora
!!=========================================================================

!> Module to generate the geometry from different kinds of inputs, from mesh
!! files or parametric input

module mod_build_geo

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, getsuboption, &
  countoption, finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_file_exists

use mod_basic_io, only: &
  read_real_array_from_file, read_mesh_basic, write_basic

use mod_cgns_io, only: &
  read_mesh_cgns

use mod_parametric_io, only: &
  read_mesh_parametric, read_actuatordisk_parametric

use mod_pointwise_io, only: &
  read_mesh_pointwise , read_mesh_pointwise_ll

use mod_stringtools, only : &
  isempty

use mod_ll_io, only: &
  read_mesh_ll, read_xac_offset

use mod_hinges, only: &
  t_hinge, t_hinge_input, build_hinges, hinge_input_parser

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

use mod_math, only: &
  cross, sort_vector_real

!----------------------------------------------------------------------

implicit none

public :: build_geometry

private

!----------------------------------------------------------------------

character(len=max_char_len) :: msg
character(len=*), parameter :: this_mod_name = 'mod_build_geo'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Builds the whole geometry from the inputs provided
!!
!! Namely calls build_component() for each component
subroutine build_geometry(geo_files, ref_tags, comp_names, output_file, &
                          global_tol_sew , global_inner_prod_te )
  character(len=*), intent(in) :: geo_files(:)
  character(len=*), intent(in) :: ref_tags(:)
  character(len=*), intent(in) :: comp_names(:)
  character(len=*), intent(in) :: output_file
  real(wp),intent(in)          :: global_tol_sew , global_inner_prod_te
  integer                      :: n_geo, i_geo
  integer(h5loc)               :: file_loc, group_loc

  n_geo = size(geo_files)

  call new_hdf5_file(output_file, file_loc)
  call new_hdf5_group(file_loc, 'Components', group_loc)
  call write_hdf5(n_geo,'NComponents',group_loc)
  do i_geo = 1,n_geo
    call build_component(group_loc, trim(geo_files(i_geo)), &
                        trim(ref_tags(i_geo)), trim(comp_names(i_geo)), i_geo , &
                        global_tol_sew , global_inner_prod_te )
  enddo

  call close_hdf5_group(group_loc)
  call close_hdf5_file(file_loc)

end subroutine build_geometry

!-----------------------------------------------------------------------

!> Builds a single geometrical component
!!
!! The subroutine reads the input in the input file for the specific
!! component, then builds the component by either loading the mesh
!! from the mesh file or generating a parametric component
subroutine build_component(gloc, geo_file, ref_tag, comp_tag, comp_id, &
                            global_tol_sew , global_inner_prod_te)
  
  integer(h5loc), intent(in)   :: gloc
  character(len=*), intent(in) :: geo_file
  character(len=*), intent(in) :: ref_tag
  character(len=*), intent(in) :: comp_tag
  integer, intent(in)          :: comp_id
  real(wp),intent(in)          :: global_tol_sew , global_inner_prod_te
  type(t_parse)                :: geo_prs
  character(len=max_char_len)  :: mesh_file
  integer, allocatable         :: ee(:,:), ee_virtual(:,:)
  real(wp), allocatable        :: rr(:,:), rr_virtual(:,:)
  real(wp)                     :: y_fountain
  real(wp), allocatable        :: rr_sym(:,:), rr_sym_virtual(:,:)
  real(wp), allocatable        :: theta_e(:), theta_p(:), chord_p(:)
  character(len=max_char_len)  :: comp_el_type
  character                    :: ElType
  character(len=max_char_len)  :: mesh_file_type
  !> Hinge 
  integer                      :: n_hinges
  type(t_parse), pointer       :: hinge_prs
  type(t_parse), pointer       :: fun_prs, file_prs, coupling_prs
  type(t_hinge_input), allocatable :: hinges(:)
  character(len=max_char_len)  :: hinge_str
  character(len=2)             :: hinge_id
  integer(h5loc)               :: hinge_loc, hinge_i_loc
  
  !> Coupling 
  logical :: coupled_comp
  character(len=max_char_len) :: coupled_str  ! there is no hdf_write for logical
  character(len=max_char_len) :: coupling_type
  real(wp)                    :: coupling_node(3)
  real(wp)                    :: coupling_node_rot(3,3)
#if USE_PRECICE
  character(len=max_char_len) :: coupling_node_file
  real(wp), allocatable       :: coupling_nodes(:,:)
  integer                     :: n_coupling_nodes, io
#endif

  !> Symmetry and mirror 
  logical                     :: mesh_symmetry, mesh_mirror
  real(wp)                    :: symmetry_point(3), symmetry_normal(3)
  real(wp)                    :: mirror_point(3),   mirror_normal(3)
  character(len=max_char_len) :: comp_name
  integer(h5loc)              :: comp_loc, geo_loc, te_loc

  !> Span 
  character(len=max_char_len), allocatable :: airfoil_list(:)
  integer , allocatable                    :: nelem_span_list(:)
  integer , allocatable                    :: i_airfoil_e(:,:)
  real(wp), allocatable                    :: normalised_coord_e(:,:)
  real(wp)                                 :: trac, radius, length
  logical                                  :: aero_table
  real(wp), allocatable                    :: thickness(:,:)
  !> Section names for CGNS
  integer                                  :: nSections, iSection
  character(len=max_char_len), allocatable :: sectionNamesCGNS(:)
  
  !> Offset and scaling factor for CGNS
  real(wp)         :: offset(3)
  real(wp)         :: scaling
  integer          :: npoints_chord_tot, npoints_chord_tot_virtual, nelems_span, nelems_span_tot

  !> Chord elements for parametric meshing
  integer          :: nelem_chord, nelem_chord_virtual
  real(wp)         :: ref_chord_fraction

  !> Connectivity and te structures
  integer , allocatable                    :: neigh(:,:)

  !> Parameters for gap sewing and te identification
  real(wp)         :: tol_sewing , inner_product_threshold
  integer          :: i_count
  
  !> Projection of the unit tangent exiting from te
  logical                     :: te_proj_logical
  character(len=max_char_len) :: te_proj_dir
  real(wp) , allocatable      :: te_proj_vec(:)
  logical                     :: suppress_te
  real(wp)                    :: scale_te

  !> Trailing edge 
  integer , allocatable       :: e_te(:,:), i_te(:,:), ii_te(:,:)
  integer , allocatable       :: neigh_te(:,:), o_te(:,:)
  real(wp), allocatable       :: rr_te(:,:), t_te(:,:)
  integer                     :: n_e_te_panel
  integer, allocatable        :: e_te_fountain(:, :) 
  integer :: i

#if USE_PRECICE
  real(wp), allocatable :: c_ref_p(:,:)
  real(wp), allocatable :: c_ref_c(:,:)
  integer               :: n_non_zero, j
#endif
  !> ll offset
  character(len=max_char_len) :: xac_offset_filen
  logical                     :: xac_offset_io
  real(wp), allocatable       :: xac_p(:)
  
  integer, allocatable :: dummy_nelem_span_list(:)

  character(len=*), parameter :: this_sub_name = 'build_component'


  !===== Define the parameters ====

  !> Mesh file
  call geo_prs%CreateStringOption('mesh_file','Mesh file definition')
  call geo_prs%CreateStringOption('mesh_file_type','Mesh file type')
  
  !> Element types
  call geo_prs%CreateStringOption('el_type', &
              'element type (temporary) p:panel, v:vortex ring, l:lifting line')
  
  !> Coupled simulation component
  call geo_prs%CreateLogicalOption('coupled', &
              'Component of a coupled simulation, w/ a structural solver?', 'F')
  call geo_prs%CreateStringOption('coupling_type', &
              'Type of the coupling: &
              &ll   : ll/beam coupling, &
              &rigid: rigid component/node, &
              &rbf: rigid component/node, &
              &none', 'none')
  call geo_prs%CreateRealArrayOption('coupling_node', &
              'Node for rigid coupling in the reference configuration (x, y, z)', &
              '(/0.0, 0.0, 0.0/)')
  call geo_prs%CreateStringOption('coupling_node_file', &
              'File containing the nodes for FSI (fluid structure interaction)')
  call geo_prs%CreateRealArrayOption('coupling_node_orientation', &
              'Orientation of the node for rigid coupling. This array contains the &
              &local components (in the local reference frame of the geometrical &
              &component) of the unit vectors of the coupling node ref.frame: &
              &A = (/ I_nod_1^loc, I_nod_2^loc, I_nod_3^loc /)', &
              '(/ 1.,0.,0., 0.,1.,0., 0.,0.,1. /)')
  call geo_prs%CreateRealArrayOption('CouplingNodeOrientationMirror', &
              'Orientation of the node for rigid coupling. This array contains the &
              &local components (in the local reference frame of the geometrical &
              &component) of the unit vectors of the coupling node ref.frame: &
              &A = (/ I_nod_1^loc, I_nod_2^loc, I_nod_3^loc /)', &
              '(/ 1.,0.,0., 0.,1.,0., 0.,0.,1. /)')
  !> ll xac offset (meant for coupled simulations)
  call geo_prs%CreateLogicalOption('Offset_xac', 'Offset in the chordwise &
              &direction of the LL elem w.r.t. the structural nodes', 'F')
  call geo_prs%CreateStringOption('Offset_xac_file','Filename of the &
              &chordwise input')
  !hinges
  call hinge_input_parser( geo_prs, hinge_prs, fun_prs, file_prs, coupling_prs )

  !> Symmetry
  call geo_prs%CreateLogicalOption('mesh_symmetry',&
              'Reflect and double the geometry?', 'F')
  call geo_prs%CreateRealArrayOption('symmetry_point', &
              'Center point of symmetry plane, (x, y, z)', &
              '(/0.0, 0.0, 0.0/)')
  call geo_prs%CreateRealArrayOption('symmetry_normal', &
              'Normal of symmetry plane, (xn, yn, zn)', &
              '(/0.0, 1.0, 0.0/)')

  !> Mirroring
  call geo_prs%CreateLogicalOption('mesh_mirror',&
              'Reflect and discard the original geometry', 'F')
  call geo_prs%CreateRealArrayOption('mirror_point', &
              'Center point of mirror plane, (x, y, z)', &
              '(/0.0, 0.0, 0.0/)')
  call geo_prs%CreateRealArrayOption('mirror_normal', &
              'Normal of mirror plane, (xn, yn, zn)', &
              '(/1.0, 0.0, 0.0/)')
  
  !> Parameters for the actuator disks
  call geo_prs%CreateRealOption('traction', &
              'Traction of the rotor')
  call geo_prs%CreateRealOption('radius', &
              'Radius of the rotor')

  !> Parameters for gap sewing and edge identification
  call geo_prs%CreateRealOption('tol_se_wing', &
              'Dimension of the gaps that will be filled, to close TE', &
              '0.001')  ! very rough default value
  call geo_prs%CreateRealOption('inner_product_te', &
              'Threshold of the inner product of adiacent elements &
              &across the TE','-0.5')  ! very rough default value

  !> Parameters for te unit tangent vector generation
  call geo_prs%CreateLogicalOption('proj_te','Remove some component &
              &from te vectors.')
  call geo_prs%CreateStringOption('proj_te_dir','Parallel or normal to &
              &ProjTeVector direction.','parallel')
  call geo_prs%CreateRealArrayOption('proj_te_vector','Vector used for &
              &the te projection.')
  call geo_prs%CreateLogicalOption('suppress_te','Suppress the trailing edge &
              &from the component','F')
  call geo_prs%CreateRealOption('scale_te','Scale the t.e. individually in &
              &each component','1.0')

  !> Section name from CGNS file
  call geo_prs%CreateStringOption('section_name', &
              'Section name from CGNS file', multiple=.true.)

  !> Scaling factor to scale CGNS files
  call geo_prs%CreateRealArrayOption('offset',&
              'Offset the points coordinates (before scaling)','(/0.0, 0.0, 0.0/)')
  call geo_prs%CreateRealOption('scaling_factor', &
                                'Scaling of the points coordinates.', '1.0')

  !> Body of revolution
  call geo_prs%CreateRealOption('rev_length',&
              'Length of the body of revolution measured from nose to nose')
  call geo_prs%CreateRealOption('rev_radius', &
              'Radius of the body of revolution section.')
  call geo_prs%CreateRealOption('rev_nose_radius', &
              'Radius of the body of revolution nose.')
  call geo_prs%CreateIntOption('rev_nelem_long', &
              'Number of elements along the length of the body revolution.')
  call geo_prs%CreateIntOption('rev_nelem_rev', &
          'Number of elements around the body of revolution circumference.')

  !> Parametric mesh
  call geo_prs%CreateIntOption('nelem_chord',&
          'number of chord-wise elements', &
          multiple=.false.);
  call geo_prs%CreateIntOption('nelem_chord_virtual',&
          'number of virtual chord-wise elements', &
          multiple=.false.);
  call geo_prs%CreateRealOption('reference_chord_fraction',&
          'Reference chord fraction', &
          '0.25',&
          multiple=.false.);

  !=====
  !===== Read the parameters ====
  !=====
  call check_file_exists(geo_file,this_sub_name,this_mod_name)
  call geo_prs%read_options(geo_file,printout_val=.false.)

  mesh_file_type  = getstr(geo_prs,'mesh_file_type')
  
  !> Symmetry
  mesh_symmetry   = getlogical(geo_prs, 'mesh_symmetry')
  symmetry_point  = getrealarray(geo_prs, 'symmetry_point',3)
  symmetry_normal = getrealarray(geo_prs, 'symmetry_normal',3)
  
  !> Mirror 
  mesh_mirror   = getlogical(geo_prs, 'mesh_mirror')
  mirror_point  = getrealarray(geo_prs, 'mirror_point',3)
  mirror_normal = getrealarray(geo_prs, 'mirror_normal',3)
  
  ! Trailing edge options 
  suppress_te   = getlogical(geo_prs, 'suppress_te')
  scale_te      = getreal(geo_prs, 'scale_te')

  comp_el_type  = getstr(geo_prs,'el_type')
  ElType        = comp_el_type(1:1)

  if ((trim(ElType) .ne. 'p') .and. ( trim(ElType) .ne. 'v') .and. &
      (trim(ElType) .ne. 'l') .and. ( trim(ElType) .ne. 'a')) then

    call error (this_sub_name, this_mod_name, 'Wrong ElType ('// &
        trim(ElType)//') for component'//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//'. ElType must be:'// &
        ' p,v,l,a.')
  end if
  
  if ( trim(mesh_file_type) .eq. 'cgns' .and. trim(ElType) .ne. 'p' ) then
    call error(this_sub_name, this_mod_name, &
        'Element type for a cgns file can only be panel. ElType = ''p'' ' )
  end if

  !> Hinges ---------------------------------------------------------------
  n_hinges = getint(geo_prs, 'n_hinges')

  if ( n_hinges .gt. 0 ) then ! Read Hinges suboptions
    call build_hinges( geo_prs, n_hinges, hinges )
  end if

  !> Coupling -------------------------------------------------------------
  coupled_comp = getlogical(geo_prs, 'coupled')
  coupled_str = 'false'
#if USE_PRECICE
  if ( coupled_comp ) then

    coupled_str   = 'true'
    coupling_type = getstr(geo_prs, 'coupling_type')
    coupling_node = getrealarray(geo_prs, 'coupling_node', 3)
    
    if ( countoption(geo_prs, 'coupling_node_file') .ne. 0 ) &
        coupling_node_file = getstr(geo_prs, 'coupling_node_file')
    
    coupling_node_rot = reshape( &
                        getrealarray(geo_prs, 'CouplingNodeOrientation', 9), &
                        (/3,3/) )

    if (  trim(coupling_type) .eq. 'rbf' ) then
      !> Open coupling_nodes_file and read nodes for FSI
      n_coupling_nodes = 0 
      io = 0

      open(unit=21, file=trim(coupling_node_file))

        do while ( io .eq. 0 ) 
          read(21,*,iostat=io) ; n_coupling_nodes = n_coupling_nodes + 1
        end do
        
        n_coupling_nodes = n_coupling_nodes - 1

      close(21)

      allocate(coupling_nodes(3,n_coupling_nodes)); coupling_nodes = 0.0_wp

      open(unit=21, file=trim(coupling_node_file))
      write(*,*) ' n_coupling_nodes: ', n_coupling_nodes
      do i = 1, n_coupling_nodes
        read(21,*) coupling_nodes(:,i)
      end do
      close(21)
    end if

    if ((trim(coupling_type) .ne. 'rigid' ) .and. &
        (trim(coupling_type) .ne. 'rbf'   ) ) then
      
      call error (this_sub_name, this_mod_name, &
        ' Coupled = T for component'//trim(comp_tag)// &
        ', but CouplingType is equal to: '//trim(coupling_type)// &
        ', it is not assigned, or it is assigned &
        &as ''none''. Assign a valid CouplingType:'//nl// &
        '  rigid: rigid component/node'//nl// &
        '  rbf  : generic component/node'//nl//'Stop'//nl)
    end if
  end if
#else
  if ( coupled_comp ) then
      call error (this_sub_name, this_mod_name, &
          ' Coupled = T for component'//trim(comp_tag)// &
          ', but no coupled simulation is expected, since dust &
          &has been compiled without #USE_PRECICE option. Stop'//nl)
  end if
#endif


  !> Parameters for gap sewing and edge identification
  !> TolSewing
  i_count = countoption(geo_prs,'tol_se_wing')
  
  if ( i_count .eq. 1 ) then
    tol_sewing = getreal(geo_prs,'tol_se_wing')
  else if ( i_count .eq. 0 ) then
    !> Inherit tol_sewing from general parameter
    tol_sewing = global_tol_sew
  else
    tol_sewing = getreal(geo_prs,'tol_se_wing')
    call warning(this_sub_name, this_mod_name, 'More than one &
        &TolSewing parameter defined for component'//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//'.')
  end if

  !> InnerProductTe
  i_count = countoption(geo_prs,'inner_product_te')

  if ( i_count .eq. 1 ) then
    inner_product_threshold = getreal(geo_prs,'inner_product_te')
  else if ( i_count .eq. 0 ) then
    !> Inherit inner_product_threshold from general parameter
    inner_product_threshold = global_inner_prod_te
  else
    inner_product_threshold = getreal(geo_prs,'inner_product_te')
    call warning(this_sub_name, this_mod_name, 'More than one &
        &InnerProductTe parameter defined for component'//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//'. First value used.')
  end if

  !> Trailing edge projection
  i_count = countoption(geo_prs,'proj_te')
  te_proj_logical = .false.
  if ( i_count .eq. 1 ) then
    te_proj_logical = getlogical(geo_prs,'proj_te')
  end if

  if ( te_proj_logical ) then
    i_count = countoption(geo_prs,'proj_te_dir' )
    if ( i_count .eq. 1 ) then
      te_proj_dir  = getstr(geo_prs, 'proj_te_dir')
    elseif ( i_count .lt. 1 ) then
      te_proj_dir = 'parallel'
      call warning (this_sub_name, this_mod_name, 'ProjTe = T, but &
        & no ProjTeVDir defined for component'//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//" Default: 'parallel'.")
    else
      te_proj_dir = getstr(geo_prs, 'proj_te_dir')
      call warning(this_sub_name, this_mod_name, 'ProjTe = T, and &
        & more than one ProjTeDir defined for component'//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//'. First value used.')
    end if
    if ( ( trim(te_proj_dir) .ne. 'parallel' ) .and. &
        ( trim(te_proj_dir) .ne. 'normal'   )         ) then
      call error (this_sub_name, this_mod_name, "Wrong 'ProjTeDir' &
        & input for component"//trim(comp_tag)// &
        ', defined in file: '//trim(geo_file)//'.')
    end if
  end if

  if ( te_proj_logical ) then
    i_count = countoption(geo_prs,'proj_te_vector' )
    if ( i_count .eq. 1 ) then
      te_proj_vec  = getrealarray(geo_prs, 'proj_te_vector',3)
    elseif ( i_count .lt. 1 ) then
      call error (this_sub_name, this_mod_name, 'ProjTe = T, but &
          & no ProjTeVector defined for component'//trim(comp_tag)// &
          ', defined in file: '//trim(geo_file)//'.')
    else
      te_proj_vec  = getrealarray(geo_prs, 'proj_te_vector',3)
      call warning(this_sub_name, this_mod_name, 'ProjTe = T, but &
          & more than one ProjTeVector defined for component'//trim(comp_tag)// &
          ', defined in file: '//trim(geo_file)//'. First value used.')
    end if
  end if

  !=====
  !===== Create the compontent group
  !=====
  write(comp_name,'(A,I3.3)')'Comp',comp_id
  call new_hdf5_group(gloc, trim(comp_name), comp_loc)
  call write_hdf5(comp_tag,'CompName',comp_loc)
  call write_hdf5(mesh_file_type,'CompInput',comp_loc)
  call write_hdf5(ref_tag,'RefTag',comp_loc)
  call write_hdf5(trim(comp_el_type),'ElType',comp_loc)
  call write_hdf5(trim(coupled_str),'Coupled',comp_loc)
  call write_hdf5(trim(coupling_type)    ,'CouplingType' ,comp_loc)
  call write_hdf5(     coupling_node     ,'CouplingNode' ,comp_loc)  
  call write_hdf5(     coupling_node_rot ,'CouplingNodeOrientation',comp_loc)
  call new_hdf5_group(comp_loc, 'Geometry', geo_loc)

  !=====
  !===== Read the input files
  !=====
  select case (trim(mesh_file_type))

  case('basic')
    mesh_file = getstr(geo_prs,'mesh_file')
    call read_mesh_basic(trim(mesh_file),ee, rr, ee_virtual, rr_virtual)

    
  case('cgns')
    mesh_file = getstr(geo_prs,'mesh_file')

    ! Check the selection of mesh sections
    nSections = countoption(geo_prs, 'section_name')

    allocate(sectionNamesCGNS(nSections))

    do iSection = 1, nSections
      sectionNamesCGNS(iSection) = trim(getstr(geo_prs, 'SectionName'))
    enddo

    call read_mesh_cgns(trim(mesh_file), sectionNamesCGNS,  ee, rr, ee_virtual, rr_virtual)

  case('revolution')

    if ( countoption(geo_prs,'mesh_file') .lt. 1 ) then

      !> Size of the body of revolution
      trac        = getreal(geo_prs, 'rev_nose_radius');
      length      = getreal(geo_prs, 'rev_length') - 2.0_wp * trac;
      radius      = getreal(geo_prs, 'rev_radius');
      nelems_span = getint (geo_prs, 'rev_nelem_long');

      if ( trac <= 0.0_wp ) then
        call error(this_sub_name, this_mod_name,  &
            'Input Rev_Nose_Radius lower than zero.')
      endif

      if ( radius <= 0.0_wp ) then
        call error(this_sub_name, this_mod_name,  &
            'Input Rev_Radius lower than zero.')
      endif

      if ( length <= 0.0_wp ) then
        call error(this_sub_name, this_mod_name,  &
            'Input Rev_Length is lower than 2*Rev_Nose_Radius.')
      endif

      if ( nelems_span < 1 ) then
        call error(this_sub_name, this_mod_name,  &
            'Input Rev_Nelem_long is lower than zero.')
      endif

      !> Generate 2D mesh
      allocate ( rr_te(nelems_span+1,2)  )
      call cigar2D(length,radius,trac,nelems_span,rr_te(:,1),rr_te(:,2))

    else

      mesh_file = getstr(geo_prs,'mesh_file')
      call read_real_array_from_file ( 2, mesh_file, rr_te )
      nelems_span = size(rr_te,1)-1

    endif

    !> Discretization of the body of revolution
    nSections = getint (geo_prs, 'rev_nelem_rev');

    if ( nSections < 1 ) then
      call error(this_sub_name, this_mod_name,  &
          'Input Rev_Nelem_rev is lower than zero.')
    endif

    ! 3D section
    allocate (rr (3,nSections*(nelems_span-1)+2), rr_virtual (3,nSections*(nelems_span-1)+2), &
                ee (4,nSections*nelems_span), ee_virtual (4,nSections*nelems_span))
    rr = 0; ee = 0; rr_virtual = 0; ee_virtual = 0
    call meshbyrev (rr_te(:,1),rr_te(:,2), nSections, rr, ee, rr_virtual, ee_virtual)

    deallocate (rr_te)

  case('parametric')

    if ((ElType .eq. 'v' .or. ElType .eq. 'p')) then
      nelem_chord = getint(geo_prs,'nelem_chord')  
      if (countoption(geo_prs, 'nelem_chord_virtual') .eq. 1) then
        nelem_chord_virtual = getint(geo_prs,'nelem_chord_virtual')
      else
        nelem_chord_virtual = nelem_chord
      end if
    end if
    ref_chord_fraction = getreal(geo_prs,'reference_chord_fraction')
    mesh_file          = geo_file

    if ((ElType .eq. 'v')) then

      !> TODO : actually it is possible to define the parameters in the GeoFile
      !> directly, find a good way to do this
    
      !> Construction of the virtual panels for collision avoidance
      call read_mesh_parametric(trim(mesh_file), ee_virtual, rr_virtual, y_fountain,       & 
                                nelem_chord_virtual, ref_chord_fraction, 'p',              &
                                npoints_chord_tot_virtual, nelems_span, hinges, n_hinges,  &
                                mesh_mirror, mesh_symmetry, nelem_span_list, airfoil_list, & 
                                i_airfoil_e, normalised_coord_e, aero_table, thickness) 
      
      !> Construction of the actual VL elements                         
      call read_mesh_parametric(trim(mesh_file), ee, rr, y_fountain,                       & 
                                nelem_chord, ref_chord_fraction, ElType,                   &
                                npoints_chord_tot, nelems_span, hinges, n_hinges,          & 
                                mesh_mirror, mesh_symmetry, nelem_span_list, airfoil_list, & 
                                i_airfoil_e, normalised_coord_e, aero_table, thickness)  
      
      !> Write additional fields for vl correction
      if (aero_table) then 

        if ( mesh_symmetry ) then
          call symmetry_update_vl_lists( nelem_span_list , &
                i_airfoil_e , normalised_coord_e, thickness)
        end if
        if ( mesh_mirror ) then
          call mirror_update_vl_lists( nelem_span_list , &
                  i_airfoil_e , normalised_coord_e, thickness )
        end if
  
        call write_hdf5(airfoil_list,       'airfoil_list',       geo_loc)
        call write_hdf5(i_airfoil_e,        'i_airfoil_e',        geo_loc)
        call write_hdf5(normalised_coord_e, 'normalised_coord_e', geo_loc)
        call write_hdf5('true',             'aero_table',         geo_loc)        
        call write_hdf5(thickness,          'thickness',         geo_loc)
      else
        call write_hdf5('false',            'aero_table',         geo_loc)
      endif

      !> Nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot =   nelems_span
        
    elseif (ElType .eq. 'p') then

      !> Construction of the virtual panels for collision avoidance
      call read_mesh_parametric(trim(mesh_file), ee_virtual, rr_virtual, y_fountain, &
                                nelem_chord_virtual, ref_chord_fraction, 'p', &
                                npoints_chord_tot_virtual, nelems_span, hinges, n_hinges, &
                                mesh_mirror, mesh_symmetry, nelem_span_list)

      !> Construction of the actual panels
      call read_mesh_parametric(trim(mesh_file), ee, rr, y_fountain, & 
                                nelem_chord, ref_chord_fraction, ElType, &
                                npoints_chord_tot, nelems_span, hinges, n_hinges, & 
                                mesh_mirror, mesh_symmetry, nelem_span_list)
                                
      !> Nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot =   nelems_span

    elseif(ElType .eq. 'l') then ! LIFTING LINE element

        if (countoption(geo_prs, 'nelem_chord_virtual') .eq. 1) then
          nelem_chord_virtual = getint(geo_prs,'nelem_chord_virtual')
        else
          nelem_chord_virtual = 20 ! Default value
        end if

      !> Construction of the virtual panels for collision avoidance
      call read_mesh_parametric(trim(mesh_file), ee_virtual, rr_virtual, y_fountain, &
                                nelem_chord_virtual, 0.25_wp, 'p', &
                                npoints_chord_tot_virtual, nelems_span, hinges, n_hinges,  &
                                mesh_mirror, mesh_symmetry, nelem_span_list, airfoil_list, & 
                                i_airfoil_e, normalised_coord_e, aero_table, thickness)

      !> Construction of the actual LL elements  
      call read_mesh_ll(trim(mesh_file), ee, rr, y_fountain, &
                        airfoil_list   , nelem_span_list   , &
                        i_airfoil_e    , normalised_coord_e, &
                        npoints_chord_tot, nelems_span, &
                        chord_p, theta_p, theta_e)
      
      ! nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot = nelems_span

      !> x_AC Offset
      xac_offset_io = getlogical(geo_prs,'Offset_xac')
      if ( xac_offset_io ) then
        xac_offset_filen = getstr(geo_prs,'Offset_xac_file')
        call read_xac_offset( xac_offset_filen, rr, &
                              xac_p )
      else
        allocate(xac_p(size(normalised_coord_e,2)+1)); xac_p = 0.0_wp
      end if
      call write_hdf5(xac_p,'x_ac_p',geo_loc)

      ! correction of the following list, if symmetry is required ---------
      if ( mesh_symmetry ) then
        call symmetry_update_ll_lists( nelem_span_list , &
              theta_e, theta_p , chord_p , i_airfoil_e , normalised_coord_e )
      end if
      if ( mesh_mirror ) then
        call mirror_update_ll_lists( nelem_span_list , &
              theta_e, theta_p , chord_p , i_airfoil_e , normalised_coord_e )
      end if

      ! -------------------------------------------------------------------

      call write_hdf5(airfoil_list   ,'airfoil_list'   ,geo_loc)
      call write_hdf5(nelem_span_list,'nelem_span_list',geo_loc)
      call write_hdf5(theta_p,'theta_p',geo_loc)
      call write_hdf5(chord_p,'chord_p',geo_loc)
      call write_hdf5(theta_e,'theta_e',geo_loc)
      call write_hdf5(i_airfoil_e,'i_airfoil_e',geo_loc)
      call write_hdf5(normalised_coord_e,'normalised_coord_e',geo_loc)

    elseif(ElType .eq. 'a') then ! ACTUATOR DISK
      call read_actuatordisk_parametric(trim(mesh_file),ee,rr)
      trac = getreal(geo_prs,'traction')
      call write_hdf5(trac,'Traction',comp_loc)
      radius = getreal(geo_prs,'radius')
      call write_hdf5(radius,'Radius', comp_loc)
    end if
    
  case('pointwise')
    
    if ((ElType .eq. 'v' .or. ElType .eq. 'p')) then
      nelem_chord         = getint(geo_prs,'nelem_chord')  
      if (countoption(geo_prs, 'nelem_chord_virtual') .eq. 1) then
        nelem_chord_virtual = getint(geo_prs,'nelem_chord_virtual')
      else
        nelem_chord_virtual = nelem_chord
      end if
    end if
    
    ref_chord_fraction = getreal(geo_prs,'reference_chord_fraction') 
    mesh_file = geo_file

    if ( ( ElType .eq. 'v' )) then  

      !> Construction of the virtual panels for collision avoidance
      call read_mesh_pointwise(trim(mesh_file), ee_virtual, rr_virtual, y_fountain, & 
                              nelem_chord_virtual, ref_chord_fraction, 'p',         &
                              npoints_chord_tot_virtual, nelems_span,               &
                              airfoil_list, i_airfoil_e, normalised_coord_e,        &
                              aero_table, thickness)  

      !> Construction of the actual VL elements
      call read_mesh_pointwise(trim(mesh_file), ee, rr, y_fountain,          & 
                              nelem_chord, ref_chord_fraction, ElType,       &
                              npoints_chord_tot, nelems_span,                &
                              airfoil_list, i_airfoil_e, normalised_coord_e, &
                              aero_table, thickness)  

      !> Write additional fields for vl correction
      if (aero_table) then 
      
      !> updates for symmetry and mirror must be done also for pointwise
      !  todo clean this up
      allocate(dummy_nelem_span_list(size(thickness)))
      if ( mesh_symmetry ) then
          call symmetry_update_vl_lists( dummy_nelem_span_list , &
                i_airfoil_e , normalised_coord_e, thickness)
        end if
        if ( mesh_mirror ) then
          call mirror_update_vl_lists( dummy_nelem_span_list , &
                  i_airfoil_e , normalised_coord_e, thickness )
        end if
      deallocate(dummy_nelem_span_list)
        
        call write_hdf5(airfoil_list,       'airfoil_list',       geo_loc)
        call write_hdf5(i_airfoil_e,        'i_airfoil_e',        geo_loc)
        call write_hdf5(normalised_coord_e, 'normalised_coord_e', geo_loc)        
        call write_hdf5('true',             'aero_table',         geo_loc)
        call write_hdf5(thickness,          'thickness',          geo_loc)        
      else
        call write_hdf5('false',            'aero_table',         geo_loc)
      endif

      !> Nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot = nelems_span
    
    elseif (( ElType .eq. 'p' ) ) then

      !> Construction of the virtual panels for collision avoidance
      call read_mesh_pointwise(trim(mesh_file), ee_virtual, rr_virtual, y_fountain, & 
                              nelem_chord_virtual, ref_chord_fraction, 'p', &
                              npoints_chord_tot_virtual, nelems_span )

      !> Construction of the actual surface panels
      call read_mesh_pointwise(trim(mesh_file), ee, rr, y_fountain,    & 
                              nelem_chord, ref_chord_fraction, ElType, &
                              npoints_chord_tot, nelems_span)
      nelems_span_tot = nelems_span

    elseif ( ElType .eq. 'l' ) then 

      if (countoption(geo_prs, 'nelem_chord_virtual') .eq. 1) then
        nelem_chord_virtual = getint(geo_prs,'nelem_chord_virtual')
      else
        nelem_chord_virtual = 20 ! Default value
      end if

      !> Construction of the virtual panels for collision avoidance
      call read_mesh_pointwise(trim(mesh_file), ee_virtual, rr_virtual, y_fountain, & 
                              nelem_chord_virtual, 0.25_wp, 'p',                    &
                              npoints_chord_tot_virtual, nelems_span)

      !> Construction of the actual LL elements
      call read_mesh_pointwise_ll(trim(mesh_file), ee, rr, y_fountain, 'l',         &
                                  airfoil_list, nelem_span_list,                    &
                                  i_airfoil_e, normalised_coord_e,                  &
                                  npoints_chord_tot, nelems_span,                   &
                                  chord_p, theta_p, theta_e )
                                  
      !> Nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot =   nelems_span

      !> Correction of the following list, if symmetry is required 
      if ( mesh_symmetry ) then
        call symmetry_update_ll_lists( nelem_span_list , &
              theta_e, theta_p , chord_p , i_airfoil_e , normalised_coord_e )
      end if

      !> Correction of the following list, if mirror is required 
      if ( mesh_mirror ) then
        call mirror_update_ll_lists( nelem_span_list , &
              theta_e, theta_p , chord_p , i_airfoil_e , normalised_coord_e )
      end if

      call write_hdf5(airfoil_list,       'airfoil_list',       geo_loc)
      call write_hdf5(nelem_span_list,    'nelem_span_list',    geo_loc)
      call write_hdf5(theta_p,            'theta_p',            geo_loc)
      call write_hdf5(chord_p,            'chord_p',            geo_loc)
      call write_hdf5(theta_e,            'theta_e',            geo_loc)
      call write_hdf5(i_airfoil_e,        'i_airfoil_e',        geo_loc)
      call write_hdf5(normalised_coord_e, 'normalised_coord_e', geo_loc)
    else
      call error(this_sub_name, this_mod_name, 'Requested an element of &
      &type '//trim(ElType)//' for the pointwise component '//trim(comp_tag)&
      &//' while only p (panels), v (vortex lattice) or l (lifting line) &
      &are aailable')
    end if

  case default
    call error(this_sub_name, this_mod_name, 'Unknown mesh file type')
  end select

  !> Scaling and offset of the mesh 

  offset   = getrealarray(geo_prs, 'offset',3)
  scaling  = getreal(geo_prs, 'scaling_factor')

  !> Apply offset 
  if (any(offset .ne. 0.0_wp)) then
    do i = 1,size(rr,2)
      rr(:,i) = rr(:,i) + offset
    enddo
    do i = 1,size(rr_virtual,2)
      rr_virtual(:,i) = rr_virtual(:,i) + offset
    enddo
  endif

  !> Apply scaling 
  if (scaling .ne. 1.0_wp) then
    rr = rr * scaling
    rr_virtual = rr_virtual * scaling
  endif

  !> coupling with preCICE 
#if USE_PRECICE
  if ( coupled_comp ) then

    if ( mesh_symmetry .or. mesh_mirror ) then
      
      if ( mesh_mirror ) then
        select case (trim(mesh_file_type))
          case('cgns', 'basic', 'revolution' )  ! TODO: check basic
            call mirror_mesh(ee_virtual, rr_virtual, mirror_point, mirror_normal)
            call mirror_mesh(ee, rr, mirror_point, mirror_normal)
            
          case('parametric','pointwise')
            call mirror_mesh_structured(ee_virtual, rr_virtual,  &
                                        npoints_chord_tot_virtual , nelems_span , &
                                        mirror_point, mirror_normal)
            call mirror_mesh_structured(ee, rr,  &
                                        npoints_chord_tot , nelems_span , &
                                        mirror_point, mirror_normal)

          case default
            call error(this_sub_name, this_mod_name,&
                'Mirror routines implemented for MeshFileType = &
                "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
                'MeshFileType = '//trim(mesh_file_type)//'. Stop.')
        end select
      end if
      
      if ( mesh_symmetry ) then
        select case (trim(mesh_file_type))
          case('cgns', 'basic', 'revolution' )  ! TODO: check basic
            call symmetry_mesh(ee_virtual, rr_virtual, symmetry_point, symmetry_normal)
            call symmetry_mesh(ee, rr, symmetry_point, symmetry_normal)
          case('parametric','pointwise')
            call symmetry_mesh_structured(ee_virtual, rr_virtual,  &
                                          npoints_chord_tot_virtual , nelems_span , &
                                          symmetry_point, symmetry_normal, rr_sym_virtual)
            call symmetry_mesh_structured(ee, rr,  &
                                          npoints_chord_tot , nelems_span , &
                                          symmetry_point, symmetry_normal, rr_sym)
            nelems_span_tot = 2*nelems_span
            ! TODO: fix mesh symmetry-> the rr before symmetry take the normal coupling nod, 
            ! instead the rr_sym take the coupling nod symmetry  
          case default
            call error(this_sub_name, this_mod_name,&
                      'Symmetry routines implemented for MeshFileType = &
                      & "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
                      'MeshFileType = '//trim(mesh_file_type))
        end select
      end if 
    end if

    if ( trim(coupling_type) .eq. 'rigid' ) then
      !> Rigid coupling between a rigid component and a "structural" node,
      ! defined as an input, coupling_node. This node represents the
      ! reference configuration for data communication between the aerodynamic
      ! and the structural solvers

      allocate(c_ref_p(3, size(rr,2))); c_ref_p = 0.0_wp
      do i =1, size(c_ref_p,2)
        !> Offset
        c_ref_p(:,i) = rr(:,i) - coupling_node
        !> Orientation
        c_ref_p(:,i) = matmul( transpose(coupling_node_rot), &
                              c_ref_p(:,i) )
      end do

      allocate(c_ref_c(3, size(ee,2))); c_ref_c = 0.0_wp
      do i =1, size(c_ref_c,2)
        n_non_zero = 0
        do j = 1, 4
          if ( ee(j,i) .ne. 0 ) then
            n_non_zero = n_non_zero + 1
            c_ref_c(:,i) = c_ref_c(:,i) + rr(:,ee(j,i))
          end if
        end do
        !> Offset
        c_ref_c(:,i) = c_ref_c(:,i)/dble(n_non_zero) - coupling_node
        !> Orientation
        c_ref_c(:,i) = matmul( transpose(coupling_node_rot), &
                              c_ref_c(:,i) )
      end do

      !> Write to hdf5 geo file
      call write_hdf5(c_ref_p,'c_ref_p',geo_loc)
      call write_hdf5(c_ref_c,'c_ref_c',geo_loc)

    elseif ( trim(coupling_type) .eq. 'rbf' ) then

      call write_hdf5( coupling_nodes,'CouplingNodes',geo_loc)
        rr = matmul( transpose(coupling_node_rot), rr )
        rr_virtual = matmul( transpose(coupling_node_rot), rr_virtual )
    end if

  else 
    !> non coupled_comp but compiled with preCICE 
    if ( mesh_symmetry ) then
      select case (trim(mesh_file_type))
        case('cgns', 'basic', 'revolution' )  ! TODO: check basic
          call symmetry_mesh(ee_virtual, rr_virtual, symmetry_point, symmetry_normal)
          call symmetry_mesh(ee, rr, symmetry_point, symmetry_normal)
          
        case('parametric','pointwise')
          call symmetry_mesh_structured(ee_virtual, rr_virtual,  &
                                          npoints_chord_tot_virtual , nelems_span , &
                                          symmetry_point, symmetry_normal, rr_sym_virtual)                                                   
          call symmetry_mesh_structured(ee, rr,  &
                                          npoints_chord_tot , nelems_span , &
                                          symmetry_point, symmetry_normal, rr_sym)
            nelems_span_tot = 2*nelems_span
        
        case default
          call error(this_sub_name, this_mod_name,&
                'Symmetry routines implemented for MeshFileType = &
                & "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
                'MeshFileType = '//trim(mesh_file_type))
      end select
    end if

    if ( mesh_mirror ) then

      select case (trim(mesh_file_type))
        case('cgns', 'basic', 'revolution' )  ! TODO: check basic
          call mirror_mesh(ee_virtual, rr_virtual, mirror_point, mirror_normal)
          call mirror_mesh(ee, rr, mirror_point, mirror_normal)
        case('parametric','pointwise')
          call mirror_mesh_structured(ee_virtual, rr_virtual,  &
                                        npoints_chord_tot_virtual , nelems_span , &
                                        mirror_point, mirror_normal)          
          call mirror_mesh_structured(ee, rr,  &
                                        npoints_chord_tot , nelems_span , &
                                        mirror_point, mirror_normal)
        case default
          call error(this_sub_name, this_mod_name,&
              'Mirror routines implemented for MeshFileType = &
              & "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
              'MeshFileType = '//trim(mesh_file_type)//'. Stop.')
      end select
    end if 
  endif  
#else 
  !> apply symmetry and mirror to the mesh 
  if ( mesh_symmetry ) then
    select case (trim(mesh_file_type))
      case('cgns', 'basic', 'revolution' )  ! TODO: check basic
        call symmetry_mesh(ee_virtual, rr_virtual, symmetry_point, symmetry_normal)
        call symmetry_mesh(ee, rr, symmetry_point, symmetry_normal)

      case('parametric','pointwise')
        call symmetry_mesh_structured(ee_virtual, rr_virtual,  &
                                        npoints_chord_tot_virtual , nelems_span , &
                                        symmetry_point, symmetry_normal, rr_sym_virtual)        
        call symmetry_mesh_structured(ee, rr,  &
                                        npoints_chord_tot , nelems_span , &
                                        symmetry_point, symmetry_normal, rr_sym)
          nelems_span_tot = 2*nelems_span

      case default
        call error(this_sub_name, this_mod_name,&
              'Symmetry routines implemented for MeshFileType = &
              & "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
              'MeshFileType = '//trim(mesh_file_type))
    end select
  end if

  if ( mesh_mirror ) then

    select case (trim(mesh_file_type))
      case('cgns', 'basic', 'revolution' )  ! TODO: check basic
        call mirror_mesh(ee_virtual, rr_virtual, mirror_point, mirror_normal)
        call mirror_mesh(ee, rr, mirror_point, mirror_normal)
      case('parametric','pointwise')
        call mirror_mesh_structured(ee_virtual, rr_virtual,  &
                                      npoints_chord_tot_virtual , nelems_span , &
                                      mirror_point, mirror_normal)        
        call mirror_mesh_structured(ee, rr,  &
                                      npoints_chord_tot , nelems_span , &
                                      mirror_point, mirror_normal)
      case default
        call error(this_sub_name, this_mod_name,&
            'Mirror routines implemented for MeshFileType = &
            & "cgns", "pointwise", "parametric", "basic", "revolution".'//nl// &
            'MeshFileType = '//trim(mesh_file_type)//'. Stop.')
    end select
  end if
#endif

  !=====
  !===== Write the geometry and connectivity
  !=====
  call write_hdf5(rr,'rr',geo_loc)
  call write_hdf5(ee,'ee',geo_loc)
  call write_hdf5(rr_virtual,'rr_virtual',geo_loc)
  call write_hdf5(ee_virtual,'ee_virtual',geo_loc)
  if (( trim(mesh_file_type) .eq. 'parametric' ) .or. &
      ( trim(mesh_file_type) .eq. 'pointwise'  ) ) then
    !write HDF5 fields
    call write_hdf5(  nelems_span_tot  ,'parametric_nelems_span',geo_loc)
    call write_hdf5(npoints_chord_tot-1,'parametric_nelems_chor',geo_loc)
    call write_hdf5(npoints_chord_tot_virtual-1,'parametric_nelems_chor_virtual',geo_loc)
  end if


  !=====
  !===== Build mesh specific fields:
  !===== * connectivity
  !===== * trailing edge
  !===== * velocity stencil
  !===== * stripe connectivity
  !=====
  selectcase(trim(mesh_file_type))
    case( 'basic', 'revolution' )

      call build_connectivity_general( ee , neigh )

      if ( mesh_symmetry ) call update_connectivity_symmetry( ee , rr , neigh )

      if ( ElType .eq. 'v' ) then
        write(msg,'(A)') 'Component named '//trim(comp_tag)//' introduced as &
          &"basic" with vortex ring elements: be sure that the inputs comply &
          &with the conventions of parametric component structures'
        call warning(this_sub_name, this_mod_name, msg)
        !TODO: IMPORTANT: what should be the behaviour here?
        if(.not. suppress_te) call build_te_parametric( ee , rr , y_fountain, ElType ,  &
            npoints_chord_tot , nelems_span_tot , &
            e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te, n_e_te_panel, e_te_fountain ) !te as an output
      else
        if(.not. suppress_te) call build_te_general ( ee , rr , y_fountain, ElType , &
                  tol_sewing , inner_product_threshold , &
                  te_proj_logical , te_proj_dir , te_proj_vec , &
                  e_te, e_te_fountain, i_te, rr_te, ii_te, neigh_te, o_te, t_te )
                                                            !te as an output
      end if

    case( 'cgns' )

      call build_connectivity_general( ee , neigh )

      if ( mesh_symmetry ) call update_connectivity_symmetry( ee , rr , neigh )

      if(.not. suppress_te) call build_te_general ( ee , rr , y_fountain, ElType , &
                tol_sewing , inner_product_threshold , &
                te_proj_logical , te_proj_dir , te_proj_vec , &
                e_te, e_te_fountain, i_te, rr_te, ii_te, neigh_te, o_te, t_te )
                                                          !te as an output

    case( 'parametric' , 'pointwise' )
      if ( ElType .eq. 'l' .or. ElType .eq. 'v' .or. ElType .eq. 'p' ) then
        call build_connectivity_parametric( ee ,     &
                      npoints_chord_tot , nelems_span_tot , &
                      neigh )
        if(.not. suppress_te) call build_te_parametric( ee , rr , y_fountain, ElType,   &
          npoints_chord_tot , nelems_span_tot , &
          e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te, n_e_te_panel, e_te_fountain ) !te as an output 
          !> debug output
          !do i = 1, size(rr, 2)
          !  write(*,*) 'rr', rr(:,i)
          !enddo
          !do i = 1, size(ee, 2)
          !  write(*,*) 'ee', i, ee(:,i)
          !enddo 
          !do i = 1, size(e_te_fountain,2) 
          !  write(*,*) 'e_te_fountain', e_te_fountain(:,i)
          !enddo
      elseif(ElType .eq. 'a') then
        allocate(neigh(size(ee,1), size(ee,2)))
        neigh = 0 !All actuator disk are isolated
      endif

#if USE_PRECICE
      if ( coupled_comp ) then
        ! *** to do *** check if the rotation is required only for 'rigid' coupling
        if ( ( trim(coupling_type) .eq. 'rigid' ) ) then
          !> Rotate t_te, if needed
          do i = 1, size(t_te,2)
            t_te(:,i) = matmul( transpose(coupling_node_rot), &
                                        t_te(:,i) )
          end do
        end if
      end if
#endif

    case default
      call error(this_sub_name, this_mod_name,&
            'Unknown option for building connectivity.')
  end select

  !=====
  !===== Write neighbours
  !=====
  call write_hdf5(neigh,'neigh',geo_loc)
  call close_hdf5_group(geo_loc)

  !=====
  !===== Write trailing edge
  !=====
  if (ElType .ne. 'a') then
    if(suppress_te) then
      !> Supppress the trailing edge for the current component
      allocate(e_te(2,0), i_te(2,0), rr_te(3,0), ii_te(2,0), neigh_te(2,0))
      allocate(o_te(2,0), t_te(3,0), e_te_fountain(2,0))
    endif
    call new_hdf5_group( comp_loc, 'Trailing_Edge', te_loc)
    call write_hdf5(         e_te,          'e_te',te_loc)
    call write_hdf5(e_te_fountain, 'e_te_fountain',te_loc)
    call write_hdf5( n_e_te_panel,  'n_e_te_panel',te_loc)
    call write_hdf5(         i_te,          'i_te',te_loc)
    call write_hdf5(        rr_te,         'rr_te',te_loc)
    call write_hdf5(        ii_te,         'ii_te',te_loc)
    call write_hdf5(     neigh_te,      'neigh_te',te_loc)
    call write_hdf5(         o_te,          'o_te',te_loc)
    call write_hdf5(         t_te,          't_te',te_loc)
    call write_hdf5(     scale_te,      'scale_te',te_loc)
    call close_hdf5_group(te_loc)
  endif

  !> Hinges 
  call new_hdf5_group(comp_loc, 'Hinges', hinge_loc)
  call write_hdf5(n_hinges, 'n_hinges', hinge_loc)
  do i = 1 , n_hinges

    !> Rotation
    if ( coupled_comp) then
      if ( trim(coupling_type) .eq. 'rbf' ) then
        hinges(i)%rr      = matmul( transpose(coupling_node_rot), hinges(i)%rr      )
        hinges(i)%ref_dir = matmul( transpose(coupling_node_rot), hinges(i)%ref_dir )
      end if
    end if

    write(hinge_id,'(I2.2)') i
    hinge_str = 'Hinge_'//hinge_id
    call new_hdf5_group(hinge_loc, trim(hinge_str), hinge_i_loc)
    call write_hdf5(hinges(i)%tag               , 'Tag'                     , hinge_i_loc)
    call write_hdf5(hinges(i)%nodes_input       , 'Nodes_Input'             , hinge_i_loc)
    call write_hdf5(hinges(i)%rr                , 'rr'                      , hinge_i_loc)
    call write_hdf5(hinges(i)%ref_dir           , 'Ref_Dir'                 , hinge_i_loc)
    call write_hdf5(hinges(i)%offset            , 'Offset'                  , hinge_i_loc)
    call write_hdf5(hinges(i)%span_blending     , 'Spanwise_Blending'       , hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_input    , 'Hinge_Rotation_Input'    , hinge_i_loc)
    !> Hinge_Rotation_Input = function:...
    ! (for other input types, dummy values, not to be used. Not a clean implementation)
    call write_hdf5(hinges(i)%rotation_amplitude, 'Hinge_Rotation_Amplitude', hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_omega    , 'Hinge_Rotation_Omega'    , hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_phase    , 'Hinge_Rotation_Phase'    , hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_amplitude_init, 'Hinge_Rotation_Amplitude_initial', hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_cosine_cycl, 'Hinge_Rotation_cosine_cycles', hinge_i_loc)
    call write_hdf5(hinges(i)%rotation_initial_time, 'Hinge_Rotation_Initial_Time', hinge_i_loc)
    !> Hinge_Rotation_Input = coupling
    if ( trim(hinges(i)%rotation_input) .eq. 'coupling' ) then
        call write_hdf5(hinges(i)%coupling_nodes, 'Coupling_Nodes', hinge_i_loc )
    endif
    call close_hdf5_group(hinge_i_loc)
  end do
  call close_hdf5_group(hinge_loc)

  call close_hdf5_group(comp_loc)

  !> cleanup
  deallocate(ee,rr)
  deallocate(ee_virtual, rr_virtual)
  if (allocated(e_te_fountain)) deallocate(e_te_fountain)
  !> hinges cleanup
  do i = 1 , n_hinges
    deallocate( hinges(i)%rr )
  end do
  
  if ( allocated(hinges) ) deallocate(hinges)

  !> Finalize parser
  call finalizeparameters(geo_prs)

end subroutine build_component


!> Subroutine used to double the mesh by reflecting it along a simmetry
!! plane
!!
!! Given a plane defined by a center point and a normal vector, the mesh
!! is doubled: all the points are reflected and new mirrored elements
!! introduced. The elements and points arrays are doubled.
subroutine symmetry_mesh(ee, rr, cent, norm)
  integer, allocatable, intent(inout)   :: ee(:,:)
  real(wp), allocatable, intent(inout)  :: rr(:,:)
  real(wp), intent(in)                  :: cent(3), norm(3)
  real(wp)                              :: n(3), d, l
  integer, allocatable                  :: ee_temp(:,:)
  real(wp), allocatable                 :: rr_temp(:,:)
  integer                               :: ip, np, ne
  integer                               :: ie, iv, nv
  character(len=*), parameter           :: this_sub_name = 'symmetry_mesh'
  
  ne = size(ee,2)
  np = size(rr,2)

  !> Enlarge size
  allocate(ee_temp(size(ee,1),2*ne))
  allocate(rr_temp(size(rr,1),2*np))

  !> First part equal
  ee_temp(:,1:ne) = ee
  rr_temp(:,1:np) = rr

  !> Second part of the elements: index incremented, need to rearrange the
  !  connectivity to preserve the normal
  ee_temp(:,ne+1:2*ne) = 0
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ee_temp(1,ne+ie) = np+ee(1,ie)
    do iv = 2,nv
      ee_temp(iv,ne+ie) = np+ee(nv-iv+2,ie)
    enddo
  enddo

  !> Calculate normal unit vector and distance from origin
  n = norm/norm2(norm)
  l = sum(cent * n)

  !> Now reflect the points
  do ip=1,np
    d = sum( rr(:,ip) * n) - l
    rr_temp(:,np+ip) = rr_temp(:,ip) - 2*d*n
  enddo

  !> Move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_temp, ee)

end subroutine symmetry_mesh


!> Subroutine used to double the mesh by reflecting it along a simmetry
!! plane for a parametric component.
!!
!! Given a plane defined by a center point and a normal vector, the mesh
!! is doubled: all the points are reflected and new mirrored elements
!! introduced. The elements and points arrays are doubled.
subroutine symmetry_mesh_structured( ee, rr, &
              npoints_chord_tot , nelems_span , cent, norm, rr_sym)

  integer, allocatable, intent(inout)   :: ee(:,:)
  real(wp), allocatable, intent(inout)  :: rr(:,:)
  real(wp), allocatable, intent(out)    :: rr_sym(:,:)
  real(wp), intent(in)                  :: cent(3), norm(3)
  integer , intent(in)                  :: npoints_chord_tot , nelems_span ! <- unused input
  integer                               :: nelems_chord
  real(wp)                              :: n(3), d, l
  integer , allocatable                 :: ee_temp(:,:) , ee_sort(:,:)
  real(wp), allocatable                 :: rr_temp(:,:)
  integer                               :: ip, np, ne
  integer                               :: ie, nv
  !> some checks  and warnings 
  real(wp)                              :: mabs  , m
  real(wp)                              :: minmaxPn
  real(wp), parameter                   :: eps = 1e-2_wp ! TODO: move it as an input
  integer                               :: imabs, i1, nsew
  logical                               :: sew_first_sec = .false. , sew_last_sec = .false.
  character(len=*), parameter           :: this_sub_name = 'symmetry_mesh_structured'


  !> Some checks 
  mabs = 0.0_wp
  do i1  = 1 , size(rr,2)
    if (      abs( sum((rr(:,i1)-cent)*norm) ) .gt. mabs ) then
      mabs  = abs( sum((rr(:,i1)-cent)*norm) ) ; imabs = i1
    end if
  end do

  m = sum( (rr(:,imabs) - cent) * norm )
  if ( m .gt. 0.0_wp ) then
    minmaxPn = sum( (rr(:,1)-cent)*norm )
    do i1  = 2 , size(rr,2)
      if (     sum( (rr(:,i1)-cent)*norm ) .lt. minmaxPn ) then
        minmaxPn = sum( (rr(:,i1)-cent)*norm )
      end if
    end do
    if ( minmaxPn .gt. eps ) then
      call error(this_sub_name, this_mod_name, 'Discontinuous component.')
    else if ( minmaxPn .lt. -eps ) then
      call error(this_sub_name, this_mod_name, 'Body compenetration.')
    end if
  else if ( m .lt. 0.0_wp ) then
    minmaxPn = sum( (rr(:,1)-cent)*norm )
    do i1  = 2 , size(rr,2)
      if (     sum( (rr(:,i1)-cent)*norm ) .gt. minmaxPn ) then
        minmaxPn = sum( (rr(:,i1)-cent)*norm )
      end if
    end do
    if ( minmaxPn .lt. -eps ) then
      call error(this_sub_name, this_mod_name, 'Discontinuous component.')
    else if ( minmaxPn .gt. eps ) then
      call error(this_sub_name, this_mod_name, 'Body compenetration.')
    end if
  end if

  !> A whole section must be sewed -
  nsew = 0
  do i1 = 1 , size(rr,2)
    if ( abs(sum( (rr(:,i1)-cent)*norm )) .lt. eps ) then
      nsew = nsew + 1

      if ( i1 .eq. 1          ) sew_first_sec = .true.
      if ( i1 .eq. size(rr,2) ) sew_last_sec = .true.

    end if
  end do

  ! === some checks ===
  !> no section to be sewed. (first, last) = (f , f)
  if ( ( .not. sew_first_sec ) .and. ( .not. sew_last_sec ) ) then
    call error(this_sub_name, this_mod_name, 'No section to be sewed. Stop')
  end if
  !> So far, avoid (first, last) = (f,t)
  if ( ( .not. sew_first_sec ) .and. ( sew_last_sec ) ) then
    call error(this_sub_name, this_mod_name, 'So far, the last section of a &
      &"parametric" of "pointwise" component can be sewed only if the first &
      &section is sewed. Stop.')
  end if
  !> (first, last) = (t , f) or (first, last) = (f , t)
  if ( ( ( nsew .ne. npoints_chord_tot ) .and. ( .not. sew_last_sec ) ) .or. &
        ( ( nsew .ne. npoints_chord_tot ) .and. ( .not. sew_first_sec ) )  ) then
    write(msg,'(A,I0,A,I0,A)') 'Wrong sewing of parametric component during &
          &symmetry reflection. Only ',nsew,' out of ',npoints_chord_tot,&
          'effectively sewed'
    call error(this_sub_name, this_mod_name, msg)
  end if
  !> (first, last) = (t , t)
  if ( ( nsew .ne. 2*npoints_chord_tot ) .and.  &
        ( ( sew_last_sec ) .and. ( sew_first_sec ) ) )  then
      write(msg,'(A,I0,A,I0,A)') 'Wrong sewing of parametric component during &
          &symmetry reflection. Only ',nsew,' out of ',npoints_chord_tot,&
          'effectively sewed'
    call error(this_sub_name, this_mod_name, msg)
  end if


  ne = size(ee,2); np = size(rr,2)

  ! enlarge size
  allocate(ee_temp(size(ee,1),2*ne))
! if ( ( .not. sew_first_sec ) .or. ( .not. sew_first_sec ) ) then ! only one sec to sew
  if ( ( .not. sew_last_sec ) ) then ! only one sec to sew
    allocate(rr_temp(size(rr,1),2*np-npoints_chord_tot))
  else ! two sections to sew
    allocate(rr_temp(size(rr,1),2*np-2*npoints_chord_tot))
  endif

  !first part of the rr,ee arrays remains the same
  ee_temp(:,1:ne) = ee
  rr_temp(:,1:np) = rr

  !second part of the elements: index incremented, need to rearrange the
  !connectivity to preserve the normal direction (outward pointing)
  ee_temp(:,ne+1:2*ne) = 0
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ! Quad elements only for structured meshes
    ee_temp(1,ne+ie) = np-npoints_chord_tot+ee(2,ie)
    ee_temp(2,ne+ie) = np-npoints_chord_tot+ee(1,ie)
    ee_temp(3,ne+ie) = np-npoints_chord_tot+ee(4,ie)
    ee_temp(4,ne+ie) = np-npoints_chord_tot+ee(3,ie)
  enddo

  ee_temp(1,ne+(/(i1,i1=1,npoints_chord_tot-1)/)) = (/(i1,i1=1,npoints_chord_tot-1)/)
  ee_temp(4,ne+(/(i1,i1=1,npoints_chord_tot-1)/)) = (/(i1,i1=2,npoints_chord_tot  )/)

  if ( sew_last_sec ) then
    ee_temp(2,2*ne-npoints_chord_tot+1+(/(i1,i1=1,npoints_chord_tot-1)/)) = &
        maxval(ee)-npoints_chord_tot + (/(i1,i1=1,npoints_chord_tot-1)/)
    ee_temp(3,2*ne-npoints_chord_tot+1+(/(i1,i1=1,npoints_chord_tot-1)/)) = &
        maxval(ee)-npoints_chord_tot + (/(i1,i1=2,npoints_chord_tot  )/)
  endif


  !calculate normal unit vector and distance from origin
  n = norm/norm2(norm)
  l = sum(cent * n)

  !now reflect the points
  if ( ( .not. sew_last_sec ) ) then ! only one sec to sew
    do ip=1,np-npoints_chord_tot
      d = sum( rr(:,ip+npoints_chord_tot) * n) - l
      rr_temp(:,np+ip) = rr_temp(:,npoints_chord_tot+ip) - 2*d*n
    enddo
  else
    do ip=1,np-2*npoints_chord_tot
      d = sum( rr(:,ip+npoints_chord_tot) * n) - l
      rr_temp(:,np+ip) = rr_temp(:,npoints_chord_tot+ip) - 2*d*n
    enddo
  endif

  ! sort ee array
  nelems_chord = npoints_chord_tot-1
  allocate(ee_sort(size(ee_temp,1),size(ee_temp,2)) ) ; ee_sort = 0
  do i1 = 1 , nelems_span
    ee_sort(:,1+(i1-1)*nelems_chord:i1*nelems_chord) = &
    ee_temp(:,2*ne-i1*nelems_chord+1:2*ne-(i1-1)*nelems_chord)
  end do
  do i1 = 1 , nelems_span
    ee_sort(:,1+(i1-1)*nelems_chord+ne:i1*nelems_chord+ne) = &
    ee_temp(:,1+(i1-1)*nelems_chord:i1*nelems_chord)
  end do
  rr_sym = rr_temp(:, np+1: np + (np-2*npoints_chord_tot))

! !move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_sort, ee)

  if ( allocated(ee_temp) )  deallocate(ee_temp)

end subroutine symmetry_mesh_structured

!----------------------------------------------------------------------

!> Subroutine used to reflect the mesh along a simmetry
!! plane.
!!
!! Given a plane defined by a center point and a normal vector, the mesh
!! is reflected: all the points are reflected and new mirrored elements
!! introduced. However the original mesh is not preserved and just
!! overwritten.
subroutine mirror_mesh(ee, rr, cent, norm)
  integer, allocatable, intent(inout)   :: ee(:,:)
  real(wp), allocatable, intent(inout)  :: rr(:,:)
  real(wp), intent(in)                  :: cent(3), norm(3)
  real(wp)                              :: n(3), d, l
  integer, allocatable                  :: ee_temp(:,:)
  real(wp), allocatable                 :: rr_temp(:,:)
  integer                               :: ip, np, ne
  integer                               :: ie, iv, nv
  character(len=*), parameter           :: this_sub_name = 'mirror_mesh'

  ne = size(ee,2) 
  np = size(rr,2)

  !> Enlarge size
  allocate(ee_temp(size(ee,1),ne)) ; ee_temp = 0
  allocate(rr_temp(size(rr,1),np)) ; rr_temp = rr

  !> Second part of the elements: index incremented, need to rearrange the
  !  connectivity to preserve the normal
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ee_temp(1,ie) = ee(1,ie)
    do iv = 2,nv
      ee_temp(iv,ie) = ee(nv-iv+2,ie)
    enddo
  enddo

  !> Calculate normal unit vector and distance from origin
  n = norm/norm2(norm)
  l = sum(cent*n)

  !> Now reflect the points
  do ip=1,np
    d = sum( rr(:,ip) * n) - l
    rr_temp(:,ip) = rr_temp(:,ip) - 2*d*n
  enddo

  !> Move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_temp, ee)

end subroutine mirror_mesh

!----------------------------------------------------------------------

!> Subroutine used to reflect the mesh along a simmetry
!! plane for parametric components.
!!
!! Given a plane defined by a center point and a normal vector, the mesh
!! is reflected: all the points are reflected and new mirrored elements
!! introduced. However the original mesh is not preserved and just
!! overwritten.
subroutine mirror_mesh_structured( ee, rr, &
              npoints_chord_tot , nelems_span , cent, norm )
  integer, allocatable, intent(inout)   :: ee(:,:)
  real(wp), allocatable, intent(inout)  :: rr(:,:)
  real(wp), intent(in)                  :: cent(3), norm(3)
  integer , intent(in)                  :: npoints_chord_tot
  integer , intent(in)                  :: nelems_span ! <- unused input
  integer                               :: nelems_chord
  real(wp)                              :: n(3), d, l
  integer , allocatable                 :: ee_temp(:,:) , ee_sort(:,:)
  real(wp), allocatable                 :: rr_temp(:,:)
  integer                               :: ip, np, ne
  integer                               :: ie, nv
  !> Some checks and warnings 
  real(wp), parameter                   :: eps = 1e-6_wp ! TODO: move it as an input
  integer                               :: i1
  character(len=*), parameter :: this_sub_name = 'mirror_mesh_structured'

  !> No check on sewing: the component is only mirrored,
  !  and the user must be carefully define it

  ne = size(ee, 2)
  np = size(rr, 2)

  !> Enlarge size
  allocate(ee_temp(size(ee,1),ne))
  allocate(rr_temp(size(rr,1),np))

  !> First part equal
  ee_temp(:,1:ne) = 0
  rr_temp(:,1:np) = rr

  !> Second part of the elements: index incremented, need to rearrange the
  !  connectivity to preserve the normal
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ! Quad elements only for structured meshes
    ee_temp(1,ie) = ee(2,ie)
    ee_temp(2,ie) = ee(1,ie)
    ee_temp(3,ie) = ee(4,ie)
    ee_temp(4,ie) = ee(3,ie)
  enddo

  !> Calculate normal unit vector and distance from origin
  n = norm/norm2(norm)
  l = sum(cent*n)

  !> Now reflect the points
  do ip=1,np
    d = sum( rr(:,ip) * n) - l
    rr_temp(:,ip) = rr_temp(:,ip) - 2*d*n
  enddo

  !> Sort ee array
  nelems_chord = npoints_chord_tot-1
  allocate(ee_sort(size(ee_temp,1),size(ee_temp,2)) ) ; ee_sort = 0
  do i1 = 1 , nelems_span
    ee_sort(:,1+(i1-1)*nelems_chord:i1*nelems_chord) = &
      ee_temp(:,ne-i1*nelems_chord+1:ne-(i1-1)*nelems_chord)
  end do

  !> Move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_sort, ee)

  if ( allocated(ee_temp) )  deallocate(ee_temp)

end subroutine mirror_mesh_structured


!> Build the connectivity of the elemens, general version
!!
!! Builds the neighbours connectivity of general elements, i.e. not
!! parametrically generated
!! TODO : connectivity is lost at the symmetry plane ---> fix it
subroutine build_connectivity_general ( ee , neigh )
  integer, allocatable, intent(in)  :: ee(:,:)
  integer, allocatable, intent(out) :: neigh(:,:)
  integer                           :: nelems , nverts
  integer                           :: nSides1 , nSides2
  integer                           :: vert1 , vert2
  integer                           :: ie1 , ie2 , iv1 , iv2
  character(len=*), parameter       :: this_sub_name = 'build_connectivity_general'

  nelems = size( ee , 2 )
  nverts = 4

  if ( size(ee,1) .ne. nverts ) then
    call error(this_sub_name, this_mod_name, 'Wrong size of the mesh&
          & connectivity. All the connectivity should be provided with &
          & 4 numbers with leading zeros for elements with fewer sides &
          &(i.e. triangles)')
  end if

  allocate( neigh(size(ee,1),size(ee,2)) ) ; neigh = 0

  do ie1 = 1 , nelems
    nSides1 = nverts - count( ee(:,ie1) .eq. 0 )

    do iv1 = 1 , nSides1

      ! if ( neigh(iv1,ie1) .ne. 0 ) cycle   ! <---- is it useful ?
      vert1 = ee(iv1,ie1)
      vert2 = ee(mod(iv1,nSides1)+1,ie1)

      do ie2 = ie1+1  , nelems

        if (any(ee(:,ie2) .eq. vert1 ) .and. &
            any(ee(:,ie2) .eq. vert2 ) ) then

          nSides2 = nverts - count( ee(:,ie2) .eq. 0 )

          do iv2 = 1 , nSides2

            if ( ee(iv2,ie2) .eq. vert2 ) then

              if ( ee(mod(iv2,nSides2)+1,ie2) .eq. vert1 ) then
                neigh(iv1,ie1) = ie2
                neigh(iv2,ie2) = ie1
              else
                write(*,*) ' first element : ' , ie1
                write(*,*) ' first element vertices : ' , ee(:,ie1)
                write(*,*) ' second element : ' , ie2
                write(*,*) ' second element vertices : ' , ee(:,ie2)
                call error(this_sub_name, this_mod_name, &
                  'Neighbouring elements with opposed normal orientation,&
                  &this is not allowed.')
              end if
            end if
          end do
        end if
      end do
    end do
  end do

end subroutine build_connectivity_general

!----------------------------------------------------------------------
!> Find if neighboring elements at the symmetry plane of the component:
! one elem belongs to the first half of the elem list, one belongs to the
! second half.
subroutine update_connectivity_symmetry( ee , rr , neigh )
  integer , allocatable, intent(in)    :: ee(:,:)
  real(wp), allocatable, intent(in)    :: rr(:,:)
  integer , allocatable, intent(inout) :: neigh(:,:)
  integer                              :: nelem, nelem_2  
  integer                              :: nsides1, nsides2
  integer                              :: i, j, k, l
  real(wp)                             :: tol = 1e-6_wp ! hardcoded

  nelem = size(ee,2)
  nelem_2 = nelem / 2

  do i = 1 , nelem_2

    ! QUAD or TRIA
    nsides1 = 4
    if ( ee(4,i) .eq. 0 ) nsides1 = 3

    do k = 1 , nsides1

      if ( neigh(k,i) .eq. 0 ) then ! missing neighbor

        do j = nelem_2+1 , nelem

          ! QUAD or TRIA
          nsides2 = 4
          if ( ee(4,j) .eq. 0 ) nsides2 = 3

            do l = 1 , nsides2

              if ( neigh(l,j) .eq. 0 ) then
                if ( sqrt( &
                    norm2( rr(:,ee(k,i)) - rr(:,ee(mod(l,nsides2)+1,j)) )**2.0_wp + &
                    norm2( rr(:,ee(l,j)) - rr(:,ee(mod(k,nsides1)+1,i)) )**2.0_wp ) .lt. tol ) then
                  !> Update neigh array
                  neigh(k,i) = j
                  neigh(l,j) = i
                end if
              end if
            end do
        end do
      end if
    end do
  end do

end subroutine update_connectivity_symmetry

!----------------------------------------------------------------------

!> Build the connectivity of the elements, parametric version
!!
!! Build the neighbours connectivity of the elements, for parametrically
!! generated quadrilateral elements
!!
!! WARNING: no fairing at the wing tip is allowed (up to now)
subroutine build_connectivity_parametric (ee, &
                  npoints_chord_tot, nelems_span, neigh)

  integer, allocatable , intent(in) :: ee(:,:)
  integer, intent(in)               :: npoints_chord_tot , nelems_span
  integer, allocatable , intent(out):: neigh(:,:)
  integer                           :: nelems_chord ! , nelems_span
  integer                           :: i1 , i2 , iel
  character(len=*), parameter       :: this_sub_name = 'build_component_parametric'

  if ( allocated(neigh) )   deallocate(neigh)
  allocate(neigh(4,size(ee,2))) ; neigh = 0

  nelems_chord = npoints_chord_tot - 1
  
  !> Inner strips 
  do i1 = 1 , nelems_span
    do i2 = 1 , nelems_chord
      iel = i2+(i1-1) * nelems_chord
      neigh(1,iel) = iel - 1
      neigh(2,iel) = iel - nelems_chord
      neigh(3,iel) = iel + 1
      neigh(4,iel) = iel + nelems_chord
    end do
  end do

  !> Correct tip elements
  neigh(2,(/ (    i1                          , i1 = 1,nelems_chord) /)) = 0
  neigh(4,(/ (i1+(nelems_span-1)*nelems_chord , i1 = 1,nelems_chord) /)) = 0
  
  !> Correct le-te elements
  neigh(1,(/ ( 1+(i1-1)*nelems_chord , i1 = 1,nelems_span) /)) = 0
  neigh(3,(/ (   (i1  )*nelems_chord , i1 = 1,nelems_span) /)) = 0

end subroutine build_connectivity_parametric

!----------------------------------------------------------------------

!> Detects the trailing edge and generates the relevant structures,
!! general version
!!
!! calls different subroutines to merge the neighbouring nodes (for
!! open trailing edges), build an updated connectivity, then
!! generate all the trailing edge structures
subroutine build_te_general ( ee , rr , y_fountain, ElType , &
                  tol_sewing , inner_prod_thresh ,  &
                  te_proj_logical , te_proj_dir , te_proj_vec , &
                  e_te, e_te_fountain, i_te, rr_te, ii_te, neigh_te, o_te, t_te)
                                                      
  integer   , intent(in)          :: ee(:,:)
  real(wp)  , intent(in)          :: rr(:,:)
  real(wp)  , intent(in)          :: y_fountain
  character , intent(in)          :: ElType
  real(wp)  , intent(in)          :: tol_sewing
  real(wp)  , intent(in)          :: inner_prod_thresh
  logical   , intent(in)          :: te_proj_logical
  character(len=*), intent(in)    :: te_proj_dir
  real(wp)  , intent(in)          :: te_proj_vec(:)
  !> TE structures
  integer , allocatable           :: e_te(:,:), i_te(:,:) , ii_te(:,:)
  integer , allocatable           :: neigh_te(:,:), o_te(:,:)
  integer , allocatable           :: e_te_fountain(:,:) 
  real(wp), allocatable           :: rr_te(:,:), t_te(:,:)
  real(wp), allocatable           :: rr_m(:,:)
  integer , allocatable           :: ee_m(:,:), i_m(:,:)
  !> 'closed-te' connectivity
  integer, allocatable            :: neigh_m(:,:)
  character(len=*), parameter     :: this_sub_name = 'build_te_general'


  call merge_nodes_general ( rr , ee , tol_sewing , rr_m , ee_m , i_m  )
  call build_connectivity_general ( ee_m , neigh_m )

  call find_te_general ( rr , ee , neigh_m , inner_prod_thresh , &
                      te_proj_logical , te_proj_dir , te_proj_vec , &
                      e_te, e_te_fountain, i_te, rr_te, ii_te, neigh_te, o_te, t_te )

  contains



subroutine merge_nodes_general (ri, ei, tol, rr, ee, im)
  real(wp), intent(in)               :: ri(:,:)
  integer , intent(in)               :: ei(:,:)
  real(wp), intent(in)               :: tol
  real(wp), allocatable, intent(out) :: rr(:,:)
  integer , allocatable, intent(out) :: ee(:,:) , im(:,:)
  integer, allocatable               :: im_tmp(:,:)
  integer                            :: nSides1, nSides2
  integer                            :: n_nodes, n_elems
  integer                            :: n_merge
  integer                            :: i1, i2, i_e1, i_e2, i_v1, i_v2
  integer                            :: ind1, inext2, iprev2, i_n1
  character(len=*), parameter :: this_sub_name = 'merge_nodes_general'

  n_nodes = size(ri,2)
  n_elems = size(ei,2)

  allocate(rr(size(ri,1),size(ri,2))) ; rr = ri
  allocate(ee(size(ei,1),size(ei,2))) ; ee = ei

  allocate(im_tmp(2,n_nodes)) ; im_tmp = 0

  n_merge = 0
  do i1 = 1 , n_nodes
    do i2 = i1 + 1 , n_nodes
      if ( norm2(ri(:,i2)-ri(:,i1)) .lt. tol ) then ! check if two nodes are very close
        do i_e2 = 1 , size(ei,2)
          ! 0 elem. assumed to be possible only in ei(4,:)
          nSides2 = count( ei(:,i_e2) .ne. 0 )
          do i_v2 = 1 , nSides2
            if ( ei(i_v2,i_e2) .eq. i2 ) then ! find elems2 where i2 belongs to
              ! avoid merging nodes belonging to the same element
              if ( all(ei(:,i_e2) .ne. i1 ) ) then
                inext2 = mod(i_v2,nSides2) + 1
                iprev2 = mod(i_v2+nSides2-2,nSides2) + 1
                ! find the elem i_e2 where i2 belongs to, and
                ! find if there is another couple of nodes to be sewed (in another iter.)
                do i_e1 = 1 , size(ei,2) ! find elems1 where i1 belongs to
                  nSides1 = count( ei(:,i_e1) .ne. 0 )
                  do i_v1 = 1 , nSides1
                    if ( ei(i_v1,i_e1) .eq. i1 ) then ! find elems1 where i1 belongs to
                      do i_n1 = 1 , nSides1-1
                        ind1 = mod( i_v1 + i_n1 - 1 , nSides1 ) + 1
                        if ( ( ( norm2( ri(:,ei(ind1  ,i_e1))   - &
                                        ri(:,ei(inext2,i_e2)) ) .lt. tol ) .or. &
                                ( norm2( ri(:,ei(ind1  ,i_e1))   - &
                                        ri(:,ei(iprev2,i_e2)) ) .lt. tol ) ) .and. &
                              ( all( ei(:,i_e1) .ne. ei(i_v2,i_e2) ) ) ) then
                              n_merge = n_merge + 1
                              rr(:,i1) = 0.5_wp * (ri(:,i1)+ri(:,i2))
                              rr(:,i2) = rr(:,i1)
                              ee(i_v2,i_e2) = i1
                              im_tmp(:,n_merge) = (/ i1 , i2 /)
                        end if
                      end do
                    end if
                  end do
                end do
              end if
            end if
          end do
        end do
      end if
    end do
  end do

! ! Misleading message: the number n_merge could be not so easy to be interpreted ---
! write(msg,'(A,I0,A)') 'During merging of nodes for trailing edge &
!       &detection ',n_merge,' nodes were found close enough to be &
!       &merged'
! call printout(msg)
! ! Misleading message: the number n_merge could be not so easy to be interpreted ---

  allocate(im(2,n_merge))  
  im = im_tmp(:,1:n_merge)

  deallocate(im_tmp)

end subroutine merge_nodes_general


subroutine find_te_general ( rr , ee , neigh_m , inner_prod_thresh , &
                te_proj_logical , te_proj_dir , te_proj_vec ,  &
                e_te, e_te_fountain, i_te, rr_te, ii_te, neigh_te, o_te, t_te )
                                                      !te as an output

  real(wp), intent(in)                :: rr(:,:)
  integer , intent(in)                :: ee(:,:) , neigh_m(:,:)
  real(wp), intent(in)                :: inner_prod_thresh
  logical   , intent(in)              :: te_proj_logical
  character(len=*), intent(in)        :: te_proj_dir
  real(wp)  , intent(in)              :: te_proj_vec(:)
  !> Actual arrays 
  integer , allocatable , intent(out) :: e_te(:,:), i_te(:,:), ii_te(:,:), e_te_fountain(:,:)
  integer , allocatable , intent(out) :: neigh_te(:,:), o_te(:,:)
  real(wp), allocatable , intent(out) :: rr_te(:,:), t_te(:,:)
  !> TMP arrays 
  integer , allocatable               :: e_te_tmp(:,:), i_te_tmp(:,:), ii_te_tmp(:,:)
  real(wp), allocatable               :: rr_te_tmp(:,:)
  integer , allocatable               :: i_el_nodes_tmp(:,:)
  real(wp), allocatable               :: nor(:,:), cen(:,:), area(:)
  integer                             :: nSides, nSides_b, nSides1, nSides2
  integer                             :: i_e, i_b, i_n, e1, e2
  integer                             :: n_el, n_p, ne_te, nn_te
  integer                             :: ind1, ind2
  integer                             :: i_node1, i_node2
  integer                             :: i_node1_max, i_node2_max
  real(wp)                            :: mi1, mi2
  real(wp)                            :: t_te_len
  integer                             :: t_te_nelem
  real(wp) , dimension(3)             :: side_dir
  ! TODO: read as an input of the component
  integer                             :: i1, i2
  character(len=*), parameter         :: this_sub_name = 'find_te_general'

  n_el = size(ee,2)
  n_p  = size(rr,2)

  ! allocate tmp structures --
  allocate( e_te_tmp(2,(n_el+1)/2)); e_te_tmp = 0
  allocate( i_el_nodes_tmp(2,n_p )); i_el_nodes_tmp = 0
  allocate( i_te_tmp      (2,n_p )); i_te_tmp = 0
  allocate(ii_te_tmp      (2,n_el)); ii_te_tmp = 0
  allocate(rr_te_tmp      (3,n_p )); rr_te_tmp = 0.0_wp

  !> Initialise counters 
  ne_te = 0 ; nn_te = 0

  !> Compute normals -> TODO: move out of this routine
  allocate(nor(3,n_el)); nor = 0.0_wp
  allocate(cen(3,n_el)); cen = 0.0_wp
  allocate(area( n_el)); area= 0.0_wp
  do i_e = 1 , n_el

    nSides = count( ee(:,i_e) .ne. 0 )

     nor(:,i_e) = 0.5_wp * cross( rr(:,ee(     3,i_e)) - rr(:,ee(1,i_e)) , &
                                  rr(:,ee(nSides,i_e)) - rr(:,ee(2,i_e))     )

    area( i_e) = norm2(nor(:,i_e))
    nor(:,i_e) = nor(:,i_e) / area(i_e)

    do i1 = 1 , nSides
      cen(:,i_e) = cen(:,i_e) + rr(:,ee(i1,i_e))
    end do
    cen(:,i_e) = cen(:,i_e) / real(nSides,wp)
  end do
  
  do i_e = 1 , n_el

    nSides = count( ee(:,i_e) .ne. 0 )

    do i_b = 1 , nSides

      if ( ( neigh_m(i_b,i_e) .ne. 0 ) .and. &
          ( all( e_te_tmp(1,1:ne_te) .ne. neigh_m(i_b,i_e) ) ) ) then
        !> Check normals
        !  1.
        if ( sum( nor(:,i_e) * nor(:,neigh_m(i_b,i_e)) ) .le. &
                                                      inner_prod_thresh ) then

          ne_te = ne_te + 1

          !> 2. ... other criteria to find te elements 

          !> Surface elements at the trailing edge 
          e_te_tmp(1,ne_te) = i_e
          e_te_tmp(2,ne_te) = neigh_m(i_b,i_e)

          !> Nodes on the trailing edge
          i_el_nodes_tmp(1,ne_te) = ee( i_b , i_e )
          i_el_nodes_tmp(2,ne_te) = ee( mod(i_b,nSides) + 1 , i_e )

          !> Find the corresponding node on the other side
          !  of the traling edge (for open te)
          ind1 = 1 ; mi1 = norm2( rr(:,i_el_nodes_tmp(1,ne_te)) -  &
                                                  rr(:,ee(1,neigh_m(i_b,i_e))) )
          ind2 = 1 ; mi2 = norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                                                  rr(:,ee(1,neigh_m(i_b,i_e))) )

          nSides_b = count( ee(:,neigh_m(i_b,i_e)) .ne. 0 )

          do i1 = 2 , nSides_b
            if (     norm2( rr(:,i_el_nodes_tmp(1,ne_te)) - &
                                rr(:,ee(i1,neigh_m(i_b,i_e)) ) ) .lt. mi1 ) then
              ind1 = i1
              mi1  = norm2( rr(:,i_el_nodes_tmp(1,ne_te)) - &
                                                rr(:,ee(i1,neigh_m(i_b,i_e)) ) )
            end if
            if (     norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                                rr(:,ee(i1,neigh_m(i_b,i_e)) ) ) .lt. mi2 ) then
              ind2 = i1
              mi2  = norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                                                rr(:,ee(i1,neigh_m(i_b,i_e)) ) )
            end if
          end do

          i_node1     = min(i_el_nodes_tmp(1,ne_te), ee(ind1,neigh_m(i_b,i_e)))
          i_node2     = min(i_el_nodes_tmp(2,ne_te), ee(ind2,neigh_m(i_b,i_e)))
          i_node1_max = max(i_el_nodes_tmp(1,ne_te), ee(ind1,neigh_m(i_b,i_e)))
          i_node2_max = max(i_el_nodes_tmp(2,ne_te), ee(ind2,neigh_m(i_b,i_e)))

          if ( all( i_te_tmp(1,1:nn_te) .ne. i_node2 ) ) then
            nn_te = nn_te + 1
            i_te_tmp(1,nn_te) = i_node2
            i_te_tmp(2,nn_te) = i_node2_max
            rr_te_tmp(:,nn_te) = 0.5_wp * ( &
                                  rr(:,i_node2) + rr(:,i_node2_max) )
            ii_te_tmp(1,ne_te) = nn_te
          else
            do i1 = 1 , nn_te   ! find ...
              if ( i_te_tmp(1,i1) .eq. i_node2 ) then
                ii_te_tmp(1,ne_te) = i1 ; exit
              end if
            end do
          end if
          if ( all( i_te_tmp(1,1:nn_te) .ne. i_node1 ) ) then
            nn_te = nn_te + 1
            i_te_tmp(1,nn_te) = i_node1
            i_te_tmp(2,nn_te) = i_node1_max
           rr_te_tmp(:,nn_te) = 0.5_wp * ( &
                                rr(:,i_node1) + rr(:,i_node1_max) )
            ii_te_tmp(2,ne_te) = nn_te
          else
            do i1 = 1 , nn_te   ! find ...
              if ( i_te_tmp(1,i1) .eq. i_node1 ) then
                ii_te_tmp(2,ne_te) = i1 ; exit
              end if
            end do
          end if
        end if
      end if
    end do
  end do

  ! From tmp to actual e_te array 
  if (allocated(e_te)) deallocate(e_te)
  allocate(e_te(2,ne_te)) ; e_te = e_te_tmp(:,1:ne_te)

  !> dummy fountain elements (TODO: currently only supported for parametric and pointwise elements) 
  if (allocated(e_te_fountain)) deallocate(e_te_fountain) 
  allocate(e_te_fountain(2,ne_te)) ; e_te_fountain = 0 
  
  ! From tmp to actual e_te array 
  if (allocated(i_te)) deallocate(i_te)
  allocate(i_te(2,nn_te)) ; i_te = i_te_tmp(:,1:nn_te)

  ! From tmp to actual r_te array 
  if (allocated(rr_te))  deallocate(rr_te)
  allocate(rr_te(3,nn_te)) ; rr_te = rr_te_tmp(:,1:nn_te)

  ! From tmp to actual r_te array 
  if (allocated(ii_te))  deallocate(ii_te)
  allocate(ii_te(2,ne_te)) ; ii_te = ii_te_tmp(:,1:ne_te)


  write(msg,'(A,I0,A,I0,A)') ' Trailing edge detection completed, found &
        &',ne_te,' elements at the trailing edge, with ',nn_te,' nodes &
        &on the trailing edge'
  call printout(msg)

  deallocate( e_te_tmp, i_te_tmp, rr_te_tmp, ii_te_tmp )

  !> Trailing edge elements connectivity 

  allocate(neigh_te(2,ne_te) ) ; neigh_te = 0
  allocate(    o_te(2,ne_te) ) ;     o_te = 0

  do e1 = 1 , ne_te
    nSides1 = 2
    do i1 = 1 , nSides1
      do e2 = e1+1 , ne_te
        nSides2 = 2
        do i2 = 1 , nSides2
          if ( ii_te(i1,e1) .eq. ii_te(i2,e2) ) then
            neigh_te(i1,e1) = e2
            neigh_te(i2,e2) = e1
            if ( i1 .ne. i2 ) then
              o_te(i1,e1) = 1
              o_te(i2,e2) = 1
            else
              o_te(i1,e1) = -1
              o_te(i2,e2) = -1
            end if
          end if
        end do
      end do
    end do
  end do

  ! Compute unit vector at the te nodes ----------------
  !  for the first implicit wake panel  ----------------
  allocate(t_te  (3,nn_te)) ; t_te = 0.0_wp
  do i_n = 1 , nn_te
    t_te_len = 0.0_wp
    t_te_nelem = 0
    do i_e = 1 , ne_te
      if ( any( i_n .eq. ii_te(:,i_e)) ) then
        t_te(:,i_n) =  t_te(:,i_n) + nor(:,e_te(1,i_e)) &
                                  + nor(:,e_te(2,i_e))
  ! **** mod 2018-07-04: avoid projection on v_rel direction.
  ! **** - avoid using hard-coded or params of the solver
  ! **** - for simulating fixed elements w/o free-stream velocity
        t_te_len = t_te_len + sqrt(area(e_te(1,i_e))) + sqrt(area(e_te(2,i_e)))
        t_te_nelem = t_te_nelem + 2
      end if
    end do

! tests: projection tests
!  ! 1. projection along the 'relative velocity'
!  t_te(:,i_n) = sum(t_te(:,i_n)*v_rel) * v_rel
!  ! 2. projection perpendicular to the side_slip direction
!  t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*side_dir) * side_dir

  ! normalisation
  t_te(:,i_n) = t_te(:,i_n) / norm2(t_te(:,i_n)) 

  t_te_len = 1.0_wp   
  t_te(:,i_n) = t_te(:,i_n) * t_te_len ! / dble(t_te_nelem)


    if ( trim(te_proj_dir) .eq. 'normal' ) then

      side_dir = te_proj_vec
      if ( te_proj_logical ) then
        t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*side_dir) * side_dir
      end if

    elseif ( trim(te_proj_dir) .eq. 'parallel' ) then

      side_dir = te_proj_vec
      if ( te_proj_logical ) then
        t_te(:,i_n) = sum(t_te(:,i_n)*side_dir) * side_dir
      end if

    end if
  end do


! WARNING: in mod_geo.f90 the inverse transformation was introduced in the
!  definition of the streamwise tangent unit vector at the te



end subroutine find_te_general



end subroutine build_te_general

!----------------------------------------------------------------------

!> Detects the trailing edge and generates the relevant structures,
!! for parametric elements
!!
!! Finds the trailing edge and fills relevant structures based on the
!! structure of parametric components
subroutine build_te_parametric ( ee , rr , y_fountain,  ElType , &
                npoints_chord_tot , nelems_span , &
                e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te, n_e_te_panel, e_te_fountain) !te as an output

  integer , intent(in)  :: ee(:,:)
  real(wp), intent(in)  :: rr(:,:)
  character,intent(in)  :: ElType
  real(wp), intent(in)  :: y_fountain 
  integer , intent(in)  :: npoints_chord_tot , nelems_span
  !> te structures
  integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
  integer , allocatable :: neigh_te(:,:) , o_te(:,:), ind(:) 
  real(wp), allocatable :: rr_te(:,:) , t_te(:,:), y_te_sorted(:), y_te_notsorted(:)
  integer               :: nelems_chord
  integer               :: i1
  integer               :: n_e_te_panel
  !> todo add as input 
  integer , allocatable :: e_te_fountain(:, :) 
  
  
  nelems_chord = npoints_chord_tot - 1
  n_e_te_panel = 0
  ! e_te: goes with the elements 
  allocate(e_te(2,nelems_span)) ; e_te = 0

  if ( ElType .eq. 'p' ) then
    e_te(1,:) = (/(   (i1  )*nelems_chord , i1=1,nelems_span)/)
    e_te(2,:) = (/( 1+(i1-1)*nelems_chord , i1=1,nelems_span)/)
    n_e_te_panel = size(e_te,2)
  else
    e_te(1,:) = (/(   (i1  )*nelems_chord , i1=1,nelems_span)/)
  end if

  !> i_te: goes with the nodes
  allocate(i_te(2,nelems_span+1)) ; i_te = 0
  if ( ElType .eq. 'p' ) then
    i_te(:,1) = (/ ee(2,1) , ee(3,nelems_chord)/)
    i_te(1,2:nelems_span+1) = (/ ( ee(1,1+(i1-1)*nelems_chord) ,  &
                                                            i1=1,nelems_span) /)
    i_te(2,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                            i1=1,nelems_span) /)
  else
    i_te(:,1) = (/ ee(3,nelems_chord) , ee(3,nelems_chord)/)
    i_te(1,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                            i1=1,nelems_span) /)
    i_te(2,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                            i1=1,nelems_span) /)
  end if
  
  !> rr_te 
  allocate(rr_te(3,nelems_span+1)) ; rr_te = 0.0_wp
  rr_te = 0.5_wp * ( rr(:,i_te(1,:)) + rr(:,i_te(2,:)) )
  
  ! Allocate id_fountain to the correct size
  allocate(e_te_fountain(2, nelems_span)); e_te_fountain = 0 
  !!> debug 
  !do i1 = 1, size(rr_te, 2)
  !  write(*,*) 'rr_te(i1) = ', rr_te(2, i1)
  !end do
  !> debug
  ! Populate id_fountain with the indices
  do i1 = 1, size(rr_te, 2)
    if ( rr_te(2, i1) .lt. y_fountain ) then !> along y axis
      if (i1 .gt. size(e_te, 2)) then 
        !> do nothing (only for mirrored cases) 
      else 
        if ( ElType .eq. 'p' ) then 
          e_te_fountain(1, i1) = e_te(1, i1)
          e_te_fountain(2, i1) = e_te(2, i1) 
        else
          e_te_fountain(1, i1) = e_te(1, i1)
        end if
      endif  
    end if
  end do 
    
  !> debug 
  !do i1 = 1, size(e_te_fountain, 2)
  !  write(*,*) 'e_te_fountain(i1) = ', e_te_fountain(1, i1)
  !end do
  
  !> ii_te 
  allocate(ii_te(2,nelems_span)) ; ii_te = 0
  if ( ElType .eq. 'p' ) then
    ii_te(1,:) = (/ ( i1 , i1 = 2,nelems_span+1 ) /)
    ii_te(2,:) = (/ ( i1 , i1 = 1,nelems_span   ) /)
  elseif ( ElType .eq. 'v' .or. ElType .eq. 'l') then
    ii_te(1,:) = (/ ( i1 , i1 = 2,nelems_span+1 ) /)
    ii_te(2,:) = (/ ( i1 , i1 = 1,nelems_span   ) /)
  end if

  !> neigh_te 
  allocate(neigh_te(2,nelems_span)) ; neigh_te = 0
  neigh_te(1,:) = (/ ( i1+1 , i1 = 1,nelems_span ) /)
  neigh_te(2,:) = (/ ( i1-1 , i1 = 1,nelems_span ) /)
  neigh_te(2,1) = 0
  neigh_te(1,nelems_span) = 0

  !> o_te 
  allocate(o_te(2,nelems_span)) ; o_te = 1
  o_te(2,1) = 0
  o_te(1,nelems_span) = 0


  !> t_te 
  allocate(t_te(3,nelems_span+1)) ; t_te = 0.0_wp
  if ( ElType .eq. 'p' ) then
    t_te(:,1) = 0.5_wp * ( rr(:,ee(3,e_te(1,1))) - rr(:,ee(2,e_te(1,1))) + &
                        rr(:,ee(2,e_te(2,1))) - rr(:,ee(3,e_te(2,1)))     )
    t_te(:,1) = t_te(:,1) / norm2(t_te(:,1))

    do i1 = 1 , nelems_span

      t_te(:,i1+1) = 0.5_wp*( rr(:,ee(4,e_te(1,i1))) - rr(:,ee(1,e_te(1,i1))) + &
                            rr(:,ee(1,e_te(2,i1))) - rr(:,ee(4,e_te(2,i1)))  )
      t_te(:,i1+1) = t_te(:,i1+1) / norm2(t_te(:,i1+1))
    end do
  else
    t_te(:,1) = 0.5_wp * (  rr(:,ee(3,e_te(1,1))) - rr(:,ee(2,e_te(1,1))) )
    t_te(:,1) = t_te(:,1) / norm2(t_te(:,1))
    do i1 = 1 , nelems_span
      t_te(:,i1+1) = 0.5_wp*( rr(:,ee(4,e_te(1,i1))) - rr(:,ee(1,e_te(1,i1))) )
      t_te(:,i1+1) = t_te(:,i1+1) / norm2(t_te(:,i1+1))
    end do

  end if

  !> clean up
  !deallocate(ind, y_te_sorted) 

end subroutine build_te_parametric

!----------------------------------------------------------------------

!> Updates lifting lines fields in case of symmetry
subroutine symmetry_update_ll_lists ( nelem_span_list , &
                theta_e, theta_p , chord_p , i_airfoil_e , normalised_coord_e )

  integer , allocatable , intent(inout) :: nelem_span_list(:)
  integer , allocatable , intent(inout) :: i_airfoil_e(:,:)
  real(wp), allocatable , intent(inout) :: normalised_coord_e(:,:)
  real(wp), allocatable , intent(inout) :: theta_e(:), theta_p(:), chord_p(:)

  integer , allocatable :: nelem_span_list_tmp(:)
  integer , allocatable :: i_airfoil_e_tmp(:,:)
  real(wp), allocatable :: normalised_coord_e_tmp(:,:)
  real(wp), allocatable :: theta_e_tmp(:), theta_p_tmp(:), chord_p_tmp(:)

  integer :: nelem_span_section
  integer :: nelems
  integer :: npts   ! nelem + 1

  integer :: i , siz

  ! Update dimensions ----------
  nelems = size(i_airfoil_e,2) * 2
  npts  = nelems + 1
  nelem_span_section = size(nelem_span_list) * 2

  !> Fill temporary arrays
  allocate(nelem_span_list_tmp( nelem_span_section ))
  siz = size(nelem_span_list)
  do i = 1 , siz
    nelem_span_list_tmp( siz+i   ) = nelem_span_list(i)
    nelem_span_list_tmp( siz-i+1 ) = nelem_span_list(i)
  end do

  allocate(i_airfoil_e_tmp( 2, nelems ))
  siz = size(i_airfoil_e,2)
  do i = 1 , siz
    i_airfoil_e_tmp( 1:2 , siz+i   ) = i_airfoil_e( 1:2    , i )
    i_airfoil_e_tmp( 1:2 , siz-i+1 ) = i_airfoil_e( 2:1:-1 , i )
  end do

  allocate(normalised_coord_e_tmp( 2, nelems ))
  siz = size(normalised_coord_e,2)
  do i = 1 , siz
    normalised_coord_e_tmp( 1:2 , siz+i   ) = normalised_coord_e( 1:2    , i )
    normalised_coord_e_tmp( 1:2 , siz-i+1 ) = 1.0_wp - normalised_coord_e( 2:1:-1 , i )
  end do

  allocate(theta_e_tmp( nelems ))
  siz = size(theta_e)
  do i = 1 , siz
    theta_e_tmp( siz+i   ) = theta_e( i )
    theta_e_tmp( siz-i+1 ) = theta_e( i )
  end do

  allocate(theta_p_tmp( npts ))
  allocate(chord_p_tmp( npts ))
  siz = size(theta_p)
  theta_p_tmp( siz ) = theta_p(1) ; chord_p_tmp( siz ) = chord_p(1)
  do i = 2 , siz
    theta_p_tmp( siz+i-1 ) = theta_p( i ) ; chord_p_tmp( siz+i-1 ) = chord_p( i )
    theta_p_tmp( siz-i+1 ) = theta_p( i ) ; chord_p_tmp( siz-i+1 ) = chord_p( i )
  end do

  ! Move_alloc to the original arrays -------------------------
  call move_alloc(    nelem_span_list_tmp ,     nelem_span_list )
  call move_alloc(        i_airfoil_e_tmp ,         i_airfoil_e )
  call move_alloc( normalised_coord_e_tmp ,  normalised_coord_e )
  call move_alloc(            theta_e_tmp ,             theta_e )
  call move_alloc(            theta_p_tmp ,             theta_p )
  call move_alloc(            chord_p_tmp ,             chord_p )

end subroutine symmetry_update_ll_lists

!----------------------------------------------------------------------

!> Updates lifting lines fields in case of mirroring
subroutine mirror_update_ll_lists ( nelem_span_list , &
          theta_e, theta_p, chord_p , i_airfoil_e , normalised_coord_e )

  integer , allocatable , intent(inout) :: nelem_span_list(:)
  integer , allocatable , intent(inout) :: i_airfoil_e(:,:)
  real(wp), allocatable , intent(inout) :: normalised_coord_e(:,:)
  real(wp), allocatable , intent(inout) :: theta_e(:), theta_p(:), chord_p(:)

  integer , allocatable :: nelem_span_list_tmp(:)
  integer , allocatable :: i_airfoil_e_tmp(:,:)
  real(wp), allocatable :: normalised_coord_e_tmp(:,:)
  real(wp), allocatable :: theta_e_tmp(:), theta_p_tmp(:), chord_p_tmp(:)

  integer :: nelem_span_section
  integer :: nelem
  integer :: npts   ! nelem + 1

  integer :: i , siz

  ! Update dimensions ----------
  nelem = sum(nelem_span_list)
  npts  = nelem + 1
  nelem_span_section = size(nelem_span_list)


  !> Fill temporary arrays
  allocate(nelem_span_list_tmp( nelem_span_section ))
  siz = size(nelem_span_list)
  do i = 1 , siz
    nelem_span_list_tmp( siz-i+1 ) = nelem_span_list(i)
  end do

  allocate(i_airfoil_e_tmp( 2, size(i_airfoil_e,2) ))
  siz = size(i_airfoil_e,2)
  do i = 1 , siz
    i_airfoil_e_tmp( 1:2 , siz-i+1 ) = i_airfoil_e( 2:1:-1 , i )
  end do

  allocate(normalised_coord_e_tmp( 2, size(normalised_coord_e,2) ))
  do i = 1 , siz
    normalised_coord_e_tmp( 1:2 , siz-i+1 ) = 1.0_wp - normalised_coord_e( 2:1:-1 , i )
  end do

  allocate(theta_e_tmp( size(theta_e) ))
  siz = size(theta_e)
  do i = 1 , siz
    theta_e_tmp( siz-i+1 ) = theta_e( i )
  end do

  allocate(theta_p_tmp( npts ))
  allocate(chord_p_tmp( npts ))
  siz = size(theta_p)
  theta_p_tmp( siz ) = theta_p(1) ; chord_p_tmp( siz ) = chord_p(1)
  do i = 2 , siz
    theta_p_tmp( siz-i+1 ) = theta_p( i ) ; chord_p_tmp( siz-i+1 ) = chord_p( i )
  end do

  ! Move_alloc to the original arrays -------------------------
  call move_alloc(    nelem_span_list_tmp ,     nelem_span_list )
  call move_alloc(        i_airfoil_e_tmp ,         i_airfoil_e )
  call move_alloc( normalised_coord_e_tmp ,  normalised_coord_e )
  call move_alloc(            theta_p_tmp ,             theta_p )
  call move_alloc(            chord_p_tmp ,             chord_p )
  call move_alloc(            theta_e_tmp ,             theta_e )

end subroutine mirror_update_ll_lists

!----------------------------------------------------------------------
!> Updates lifting lines fields in case of symmetry
subroutine symmetry_update_vl_lists ( nelem_span_list , &
                        i_airfoil_e , normalised_coord_e , thickness)

  integer , allocatable , intent(inout) :: nelem_span_list(:)
  integer , allocatable , intent(inout) :: i_airfoil_e(:,:)
  real(wp), allocatable , intent(inout) :: normalised_coord_e(:,:)
  real(wp), allocatable , intent(inout) :: thickness(:,:)

  integer , allocatable :: nelem_span_list_tmp(:)
  integer , allocatable :: i_airfoil_e_tmp(:,:)
  real(wp), allocatable :: normalised_coord_e_tmp(:,:)
  real(wp), allocatable :: thickness_tmp(:,:)
  
  integer :: nelem_span_section
  integer :: nelems
  integer :: npts   ! nelem + 1

  integer :: i , siz

  ! Update dimensions ----------
  nelems = size(i_airfoil_e,2) * 2
  npts  = nelems + 1
  nelem_span_section = size(nelem_span_list) * 2

  !> Fill temporary arrays
  allocate(nelem_span_list_tmp( nelem_span_section ))
  siz = size(nelem_span_list)
  do i = 1 , siz
    nelem_span_list_tmp( siz+i   ) = nelem_span_list(i)
    nelem_span_list_tmp( siz-i+1 ) = nelem_span_list(i)
  end do

  allocate(i_airfoil_e_tmp( 2, nelems ))
  siz = size(i_airfoil_e,2)
  do i = 1 , siz
    i_airfoil_e_tmp( 1:2 , siz+i   ) = i_airfoil_e( 1:2    , i )
    i_airfoil_e_tmp( 1:2 , siz-i+1 ) = i_airfoil_e( 2:1:-1 , i )
  end do

  allocate(normalised_coord_e_tmp( 2, nelems ))
  siz = size(normalised_coord_e,2)
  do i = 1 , siz
    normalised_coord_e_tmp( 1:2 , siz+i   ) = normalised_coord_e( 1:2    , i )
    normalised_coord_e_tmp( 1:2 , siz-i+1 ) = 1.0_wp - normalised_coord_e( 2:1:-1 , i )
  end do

  allocate(thickness_tmp( 2, nelems ))
  siz = size(thickness,2)
  do i = 1 , siz
    thickness_tmp( 1:2 , siz+i   ) = thickness( 1:2    , i )
    thickness_tmp( 1:2 , siz-i+1 ) = thickness( 2:1:-1 , i )
  end do

  ! Move_alloc to the original arrays -------------------------
  call move_alloc(    nelem_span_list_tmp ,     nelem_span_list )
  call move_alloc(        i_airfoil_e_tmp ,         i_airfoil_e )
  call move_alloc( normalised_coord_e_tmp ,  normalised_coord_e )
  call move_alloc( thickness_tmp ,  thickness )
  
  

end subroutine symmetry_update_vl_lists

!----------------------------------------------------------------------

!> Updates lifting lines fields in case of mirroring
subroutine mirror_update_vl_lists ( nelem_span_list , &
  i_airfoil_e , normalised_coord_e, thickness)

  integer , allocatable , intent(inout) :: nelem_span_list(:)
  integer , allocatable , intent(inout) :: i_airfoil_e(:,:)
  real(wp), allocatable , intent(inout) :: normalised_coord_e(:,:)
  real(wp), allocatable , intent(inout) :: thickness(:,:)

  integer , allocatable :: nelem_span_list_tmp(:)
  integer , allocatable :: i_airfoil_e_tmp(:,:)
  real(wp), allocatable :: normalised_coord_e_tmp(:,:)
  real(wp), allocatable :: thickness_tmp(:,:) 

  integer :: nelem_span_section
  integer :: nelem
  integer :: npts   ! nelem + 1

  integer :: i , siz

  ! Update dimensions ----------
  nelem = sum(nelem_span_list)
  npts  = nelem + 1
  nelem_span_section = size(nelem_span_list)


  !> Fill temporary arrays
  allocate(nelem_span_list_tmp( nelem_span_section ))
  siz = size(nelem_span_list)
  do i = 1 , siz
    nelem_span_list_tmp( siz-i+1 ) = nelem_span_list(i)
  end do

  allocate(i_airfoil_e_tmp( 2, size(i_airfoil_e,2) ))
  siz = size(i_airfoil_e,2)
  do i = 1 , siz
    i_airfoil_e_tmp( 1:2 , siz-i+1 ) = i_airfoil_e( 2:1:-1 , i )
  end do

  allocate(normalised_coord_e_tmp( 2, size(normalised_coord_e,2) ))
  do i = 1 , siz
    normalised_coord_e_tmp( 1:2 , siz-i+1 ) = 1.0_wp - normalised_coord_e( 2:1:-1 , i )
  end do

  allocate(thickness_tmp( 2, size(thickness,2) ))
  siz = size(thickness,2)
  do i = 1 , siz
    thickness_tmp( 1:2 , siz-i+1 ) = thickness( 2:1:-1 , i )
  end do
  
  ! Move_alloc to the original arrays -------------------------
  call move_alloc(    nelem_span_list_tmp ,     nelem_span_list )
  call move_alloc(        i_airfoil_e_tmp ,         i_airfoil_e )
  call move_alloc( normalised_coord_e_tmp ,  normalised_coord_e )
  call move_alloc( thickness_tmp ,  thickness )

end subroutine mirror_update_vl_lists


!> [X,Y] = cigar2D(L,R,RN,N)
!!
!! Generate the 2D curve describing half of the boom geometry. The boom is
!! composed of three parts a left nose, a central straight section and a
!! right nose. The geometry is symmetric and the nose is described by a NACA
!! 4-digit equation taken from the leading edge to the maximum thickness
!! coordinate (~30! of the chord).
!!
!!      left          straight        right
!!      nose          section          nose
!!
!!   |---RN---|----------L----------|---RN---|
!!
!! The overall length is 2*RN + L, and the thickness to chord ratio of the
!! corresponding NACA for digit airfoil is given by 2*st*R/RN. With st
!! roughly equal to 0.3.
!!
!! Inputs
!!  L    length of the central straight section
!!  R    radius of the cross section of the body of revolution
!!  RN   length of the nose fairing
!!
!! Outputs
!! X    longitudinal coordinates of the body of revolution geometry
!! Y    transversal coordinates of the body of revolution  geometry

subroutine cigar2D(L,R,RN,n,x,y)

  real(wp),               intent(in)  :: L, R, RN
  integer,                intent(in)  :: n
  real(wp), dimension(:), intent(out) :: x, y

  real(wp) :: tc, lhalf, cc, dphi, phi, s
  integer  :: i

  real(wp), parameter :: st = 0.29983_wp, &
  pi = acos(-1.0_wp)

  tc    = 2.0_wp*R/RN*st
  cc    = 2.0_wp*R/tc
  lhalf = ( st*cc + 0.5_wp*L )
  dphi = pi/real(n,wp)

  do i = 1,(n+1)
    phi = dphi*real(i-1,wp)

    x(i) = lhalf*cos(phi)

    s = ( lhalf-abs(x(i)) ) / cc

    if ( (s-st) .lt. 0.0_wp ) then
      y(i)= cc * 5.0_wp*tc*( 0.2969_wp*sqrt(s)-0.1260_wp*s  &
          -0.3516_wp*s**2 + 0.2843_wp*s**3 -0.1015_wp*s**4 )
    else
      y(i) = R
    endif
  enddo

end subroutine cigar2D

!----------------------------------------------------------------------

!> [RR,EE] = meshbyrev ( X, Y, NPHI )
!!
!! Generate a mesh by revolution around the x-axis of the input curve,
!! expressed by the input coordinates Y and Y.
!!
!! Inputs
!! X    longitudinal coordinates of the generating curve
!! Y    transversal coordinates of the generating curve, the first and the
!!      last point must be zero
!! NPHI number of points in the angular direction
!!
!! Outputs
!! RR   coordinates of the mesh nodes, size (3,nphi*nc+2)
!! EE   element to node connectivity, size (4,nphi*(n-1))

subroutine meshbyrev ( x, y, nphi, rr, ee, rr_virtual, ee_virtual )

  real(wp), dimension(:),   intent(in)  :: x, y
  integer,                  intent(in)  :: nphi
  real(wp), dimension(:,:), intent(out) :: rr, rr_virtual
  integer,  dimension(:,:), intent(out) :: ee, ee_virtual

  real(wp) :: dphi, phi
  integer  :: i, ip, j, e, n, nc

  real(wp), parameter :: pi = acos(-1.0_wp)
  character(len=*), parameter :: this_sub_name = 'meshbyrev'

  n    = size(x)
  nc   = n-2
  dphi = 2.0_wp*pi/dble(nphi)
  if ( any( y .lt. -1.0e-16_wp ) ) then
      call error(this_sub_name, this_mod_name, &
        'invalid mesh, input curve has negative points')
  endif

  ! Generation of the nodes

  e = 1
  do i = 1,nphi
    phi = dphi*dble(i-1)
    do j = 2,nc+1
      e = e + 1
      rr(1,e) = x(j)
      rr(2,e) = y(j)*sin(phi)
      rr(3,e) = y(j)*cos(phi)
    enddo
  enddo
  rr(1,1)          = x(1)
  rr(1,size(rr,2)) = x(n)

  if ( y(1) .gt. 1.0e-16_wp .or. y(n) .gt. 1.0e-16_wp ) then
      call warning(this_sub_name, this_mod_name, &
          'closing input curve, check the extrema')

      rr(2,1)          = 0.0_wp
      rr(2,size(rr,2)) = 0.0_wp
      rr(3,1)          = 0.0_wp
      rr(3,size(rr,2)) = 0.0_wp
  endif


  !> Generation of the elements

  e = 0
  do i = 1,nphi

    ip = i+1
    if ( ip > nphi ) ip = 1

    !> Left triangular elements
    e = e+1
    ee(2,e) = 1
    ee(1,e) = 2 + nc*(i-1)
    ee(3,e) = 2 + nc*(ip-1)

    !> Center of the body quadragles
    do j = 1,nc-1
      e = e+1
      ee(2,e) = 1 + nc*(i-1)  + j
      ee(1,e) = 2 + nc*(i-1)  + j
      ee(4,e) = 2 + nc*(ip-1) + j
      ee(3,e) = 1 + nc*(ip-1) + j
    enddo

    ! Right triangular elements
    e = e+1
    ee(2,e) = 1 + nc*i
    ee(1,e) = size(rr,2)
    ee(3,e) = 1 + nc*ip
  enddo
  ! Virtual Mesh
  rr_virtual = rr;
  ee_virtual = ee;
end subroutine meshbyrev

end module mod_build_geo
