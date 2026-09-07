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
!!          Davide Montagnani
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================

module mod_precice

use mod_param, only: &
  wp, pi, nl, max_char_len

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, printout

use mod_math, only: &
  cross, rotation_vector_combination, vec2mat

use mod_geometry, only: &
  t_geo, t_geo_component, t_tedge

use mod_aeroel, only: &
  t_pot_elem_p

use mod_virtual, only: &
  t_elem_virtual_p

use mod_wake, only: &
  t_wake, join_first_panels

use mod_liftlin, only: &
  t_liftlin, t_liftlin_p

use mod_hinges, only: &
  t_hinge

use mod_wind, only: &
  variable_wind

implicit none

private

public :: t_precice

!> Parameters
integer, parameter :: precice_mcl = max_char_len

!> PreCICE mesh -------------------------------------------------
type :: t_precice_mesh
  character(len=precice_mcl) :: mesh_name
  integer                    :: mesh_id
  integer , allocatable      :: node_ids(:)
  real(wp), allocatable      :: nodes(:,:)
  integer                    :: nnodes
  integer                    :: ndim
  integer , allocatable      :: elem_ids(:)
end type t_precice_mesh

!> PreCICE field ------------------------------------------------
type :: t_precice_field
  integer :: fid
  character(len=precice_mcl) :: fname
  character(len=precice_mcl) :: fio    ! 'read'/'write'
  character(len=precice_mcl) :: ftype  ! 'scalar'/'vector'
  real(wp), allocatable :: fdata(:,:)  ! (1,nnodes)/(nd,nnodes)
  real(wp), allocatable :: cdata(:,:)  ! data for storing/reloading fields
end type t_precice_field

!> PreCICE type -------------------------------------------------
type :: t_precice

  !> PreCICE configuration
  character(len=precice_mcl) :: config_file_name
  character(len=precice_mcl) :: solver_name
  character(len=precice_mcl) :: mesh_name     ! *** to do *** just one mesh?
  integer :: comm_rank, comm_size

  !> PreCICE mesh
  type(t_precice_mesh) :: mesh

  !> PreCICE fields
  integer :: n_fields
  type(t_precice_field), allocatable :: fields(:)

  !> PreCICE variables
  real(wp) :: dt_precice
  integer  :: is_ongoing

  contains

  procedure, pass(this) :: initialize
  procedure, pass(this) :: initialize_mesh
  procedure, pass(this) :: initialize_fields
  procedure, pass(this) :: update_force
  procedure, pass(this) :: update_force_coupled_hinge
  procedure, pass(this) :: update_elems
  procedure, pass(this) :: update_near_field_wake

end type t_precice

!> --------------------------------------------------------------
contains
!----------------------------------------------------------------
!>
subroutine initialize(this) ! first part of constructor creation for dust participant
  class(t_precice), intent(inout) :: this
  logical :: exists
  
  call printout(nl//'Using PreCICE')

  ! *** to do *** read %config_file_name as an input
  !> Default input for dust in a preCICE coupled simulation 
  this%config_file_name = sim_param%precice_config
  this%solver_name = 'dust' ! todo change solver to participant
  this%  mesh_name = 'dust_mesh'
  this%comm_rank = 0
  this%comm_size = 1
  
  inquire(file=trim(this%config_file_name), exist=exists)
  if (.not. exists) then
    call error('initialize', 'mod_precice', 'Precice configuration file "'&
              &//trim(this%config_file_name)//'" not found.'//nl//'If this is not a coupled &
              &simulation, please compile DUST using WITH_PRECICE=NO')
  end if
  
  !> Initialize PreCICE participant and mesh
  call precicef_create( this%solver_name, &
                        this%config_file_name, &
                        this%comm_rank, &
                        this%comm_size )

  this%mesh%mesh_name = this%mesh_name
  call precicef_get_mesh_dimensions( this%mesh%mesh_name, this%mesh%ndim )
  
end subroutine initialize

!----------------------------------------------------------------
!>
subroutine initialize_mesh( this, geo )
  class(t_precice), intent(inout) :: this
  type(t_geo)     , intent(in)    :: geo

  integer :: i_comp, n_comp
  integer :: dnnodes, nnodes, ih, n_hinges

  !> TODO: add component field, that describe if the component
  ! participates to the coupling. So far, all the components
  ! participate.
  logical :: coupling = .true.
  character(len=precice_mcl) :: coupling_type

  n_comp = size(geo%components)

  !> Count n. of nodes ==========================================
  nnodes = 0
  do i_comp = 1, n_comp
    coupling      = geo%components(i_comp)%coupling
    coupling_type = geo%components(i_comp)%coupling_type
    if ( coupling ) then
      if ( trim(coupling_type) .eq. 'rigid' ) then
        !> rigid coupling
        !> Increment number of nodes
        dnnodes = 1
        nnodes = nnodes + dnnodes
      elseif ( trim(coupling_type) .eq. 'rbf' ) then
        !> Set PreCICE nodes
        dnnodes = size(geo%components(i_comp)%rbf%nodes,2)
        !> Increment number of nodes
        nnodes = nnodes + dnnodes
        !> Hinges 
        n_hinges = size(geo%components(i_comp)%hinge)
        do ih = 1, n_hinges
          if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. &
                'coupling' ) then
            dnnodes = size(geo%components(i_comp)%hinge(ih)%i_points_precice)
            nnodes = nnodes + dnnodes
          end if
        end do
      else
        write(*,*) ' Error in initialize mesh. CouplingType= '// &
                    trim(coupling_type)//', while only rigid'&
                    ' CouplingType are implemented. Stop'; stop
      end if
    end if
  end do

  !> Allocate participant%mesh fields ===========================
  ! dust needs two grids:
  ! - node-centered for reading position from MBDyn
  ! - cell-centered for writing force to MBDyn
  allocate(this%mesh%node_ids(nnodes)); this%mesh%node_ids = 0
  allocate(this%mesh%nodes( this%mesh%ndim, nnodes ))
  nnodes = 0
  do i_comp = 1, n_comp
    coupling      = geo%components(i_comp)%coupling
    coupling_type = geo%components(i_comp)%coupling_type
    if ( coupling ) then
      if ( trim(coupling_type) .eq. 'rigid' ) then
        !> rigid coupling
        dnnodes = 1
        !> Here in the local reference frame ! ***to do***
        this%mesh%nodes(:,nnodes+1) = &
          geo%components(i_comp)%coupling_node
        nnodes = nnodes + dnnodes
      elseif ( trim(coupling_type) .eq. 'rbf' ) then
        !> rbf coupling
        dnnodes = size(geo%components(i_comp)%rbf%nodes,2)

        this%mesh%nodes(:,geo%components(i_comp)%i_points_precice) = &
          geo%components(i_comp)%rbf%nodes
        ! new, w/ hinges ---
        nnodes = nnodes + dnnodes
        !> Hinges ---
        n_hinges = size(geo%components(i_comp)%hinge)
        do ih = 1, n_hinges
          if ( trim(geo%components(i_comp)%hinge(ih)%input_type) .eq. &
                'coupling' ) then
            dnnodes = size(geo%components(i_comp)%hinge(ih)%i_points_precice)
            this%mesh%nodes(:,geo%components(i_comp)%hinge(ih)%i_points_precice) = &
              geo%components(i_comp)%hinge(ih)%nodes
            nnodes = nnodes + dnnodes
          end if
        end do
      end if
    end if
  end do


  !> Initialize mesh ============================================
  !> Nodes
  call precicef_set_vertices( this%mesh%mesh_name, nnodes,  &
                              this%mesh%nodes, this%mesh%node_ids )
  this%mesh%nnodes = nnodes

end subroutine initialize_mesh

!----------------------------------------------------------------
!>
subroutine initialize_fields( this )
  class(t_precice), intent(inout) :: this

  !> Expected fields
  integer, parameter :: n_max_fields = 6
  character(len=precice_mcl) :: field_list(n_max_fields)
  character(len=precice_mcl) :: type_list(n_max_fields)
  character(len=precice_mcl) :: io_list(n_max_fields)

  integer :: i, nf, hasdata

  !> Fields that can be exchanged through PreCICE. Add a field if needed
  field_list(1) = 'Position'       ; type_list(1) = 'vector'; io_list(1) = 'read'
  field_list(2) = 'Velocity'       ; type_list(2) = 'vector'; io_list(2) = 'read'
  field_list(3) = 'Rotation'       ; type_list(3) = 'vector'; io_list(3) = 'read'
  field_list(4) = 'AngularVelocity'; type_list(4) = 'vector'; io_list(4) = 'read'
  field_list(5) = 'Force'          ; type_list(5) = 'vector'; io_list(5) = 'write'
  field_list(6) = 'Moment'         ; type_list(6) = 'vector'; io_list(6) = 'write'

  !> Count n. of exchanged fields
  ! has_data is deprecated
  !nf = 0
  !do i = 1, n_max_fields
  !  call precicef_has_data( trim(field_list(i)), this%mesh%mesh_id, hasdata )
   ! if ( hasdata .eq. 1 ) nf = nf + 1
  !end do

  allocate( this%fields( 6 ) )
  nf = 0
  do i = 1, n_max_fields
    !call precicef_has_data( trim(field_list(i)), this%mesh%mesh_id, hasdata )
    !if ( hasdata .eq. 1 ) then
      nf = nf + 1
      this%fields(nf)%fname = trim(field_list(i))
      this%fields(nf)%fio   = trim(   io_list(i))
      this%fields(nf)%ftype = trim( type_list(i))
      !call precicef_get_data_id( trim(this%fields(nf)%fname), &
      !                          this%mesh%mesh_id, this%fields(nf)%fid )
      !> Allocate and initialize fields to zero *** to do *** where to write initial data?
      if ( trim(this%fields(nf)%ftype) .eq. 'scalar' ) then
        allocate( this%fields(nf)%fdata(1, this%mesh%nnodes) )
        this%fields(nf)%fdata = 0.0_wp
      else if ( trim(this%fields(nf)%ftype) .eq. 'vector' ) then
        allocate( this%fields(nf)%fdata(this%mesh%ndim, &
                                        this%mesh%nnodes))
        this%fields(nf)%fdata = 0.0_wp
      end if
    !end if
  end do

  !> Allocate cdata for storing/reloading data in implicit coupling
  ! and initialize it
  do i = 1, nf
    if ( trim(this%fields(i)%fio) .eq. 'write' ) then
      allocate( this%fields(i)%cdata( size(this%fields(i)%fdata,1), &
                                      size(this%fields(i)%fdata,2) ) )
      this%fields(i)%cdata = this%fields(i)%fdata
    end if
  end do

end subroutine initialize_fields

!----------------------------------------------------------------
!> Update force/moment fields
subroutine update_force( this, geo, elems )
  class(t_precice)  , intent(inout) :: this
  type(t_geo)       , intent(inout) :: geo
  type(t_pot_elem_p), intent(in)    :: elems(:)
  real(wp) :: n_rot(3), chord(3), chord_rot(3)
  real(wp) :: theta
  real(wp) :: eps = 1.0e-9_wp

  integer :: i, j, i_comp, iw, ip, ih

  integer :: j_for, j_mom, j_rot, j_pos


  do j = 1, size(this%fields)
    if ( trim(this%fields(j)%fname) .eq. 'Position') j_pos = j
    if ( trim(this%fields(j)%fname) .eq. 'Rotation') j_rot = j
    if ( trim(this%fields(j)%fname) .eq. 'Force'   ) j_for = j
    if ( trim(this%fields(j)%fname) .eq. 'Moment'  ) j_mom = j
  end do

  do i_comp = 1, size(geo%components)
    associate( comp => geo%components(i_comp) )

    if ( comp%coupling ) then

      !> rigid coupling --------------------------------------------------------
      if ( trim(comp%coupling_type) .eq. 'rigid' ) then
        !> rigid coupling. All the forces and moments are reduced to
        ! the coupling_node
        ip = comp%i_points_precice(1)

        !> Reset force and moment fields, to be filled by accumulation
        this%fields(j_for)%fdata(:, ip ) = 0.0_wp
        this%fields(j_mom)%fdata(:, ip ) = 0.0_wp
        do ih = 1, comp%n_hinges
          if ( trim(comp%hinge(ih)%input_type) .eq. 'coupling' ) then
            do i = 1, size( comp%hinge(ih)%i_points_precice )
              ip = comp%hinge(ih)%i_points_precice( i )
              this%fields(j_for)%fdata(:, ip) = 0.0_wp
              this%fields(j_mom)%fdata(:, ip) = 0.0_wp
            end do
          end if
        end do

        !> Forces
        do i = 1, size(comp%el)

          this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) + &
                                            comp%el(i)%dforce
        end do


        !> Rotation
        n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(1))
        theta = norm2( n_rot )

        if ( theta .lt. eps ) then
          n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
          theta = 0.0_wp
        else
          n_rot = n_rot / theta
        end if

        if ( comp%comp_el_type(1:1) .eq. 'l' ) then

          do i = 1, size(comp%el)

            !> vector between the center of the elements and the coupling node
            chord = 0.5_wp * (comp%c_ref_p(:,2*(i-1)+1) + &
                              comp%c_ref_p(:,2* i   +1) )
            chord_rot =  cos(theta) * chord + &
                         sin(theta) * cross( n_rot, chord ) + &
                       ( 1.0_wp - cos(theta) ) * sum( chord* n_rot ) * n_rot

            this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) + &
              comp%el(i)%dmom + &
              cross( chord_rot, comp%el(i)%dforce)

          end do

        else

          do i = 1, size(comp%el)

            !> vector between the center of the elements and the coupling node
            chord = comp%c_ref_c(:,i)
            chord_rot =  cos(theta) * chord + &
                         sin(theta) * cross( n_rot, chord ) + &
                       ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot

            this%fields(j_mom)%fdata(:,ip) =  this%fields(j_mom)%fdata(:,ip) + &
              comp%el(i)%dmom + &
              cross( chord_rot, comp%el(i)%dforce)

          end do

        end if

      !> rbf coupling ----------------------------------------------------------
      elseif ( trim(comp%coupling_type) .eq. 'rbf' ) then

        !> Forces and moments by accumulation, only those belonging to the comp
        do i = 1, size( comp%i_points_precice )
          ip = comp%i_points_precice( i )
          this%fields(j_for)%fdata(:,ip) = 0.0_wp
          this%fields(j_mom)%fdata(:,ip) = 0.0_wp
        end do
        do ih = 1, comp%n_hinges
          if ( trim(comp%hinge(ih)%input_type) .eq. 'coupling' ) then
            do i = 1, size( comp%hinge(ih)%i_points_precice )
              ip = comp%hinge(ih)%i_points_precice( i )
              this%fields(j_for)%fdata(:, ip) = 0.0_wp
              this%fields(j_mom)%fdata(:, ip) = 0.0_wp
            end do
          end if
        end do

        !> === Aerodynamic mesh to structural nodes ===
        !  w/o considering hinge rotations
        do i = 1, size(comp%el)
          do iw = 1, size(comp%rbf%cen%ind,1)

            ip = comp%i_points_precice( comp%rbf%cen%ind(iw,i) )

            !> Force
            this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) + &
                    comp%rbf%cen%wei(iw,i) * comp%el(i)%dforce

            !> Moments
            if ( comp%comp_el_type(1:1) .eq. 'l' ) then
              comp%el(i)%cen =  sum ( comp%el(i)%ver(:,1:2),2 ) / 2.0_wp
            endif

            this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) + &
                    comp%rbf%cen%wei(iw,i) *( ( comp%el(i)%dmom)  + &
                    cross( comp%el(i)%cen  - this%fields(j_pos)%fdata(:,ip) , &
                    comp%el(i)%dforce ) )

          end do
        end do

        ! -------------------------------------------------------------------------------
        !>  === Add hinge motion: START ===
        ! -------------------------------------------------------------------------------
        !> === Aerodynamic mesh to hinge nodes ===
        ! - reduce forces and moments to hinge nodes
        ! - update "structural" nodes
        do ih = 1, comp%n_hinges
          if ( trim(comp%hinge(ih)%input_type) .eq. 'coupling' ) then

            !> Update structural forcing, taking into account hinges
            call this%update_force_coupled_hinge( comp, comp%hinge(ih), j_for, j_mom )

          end if
        end do
        ! -------------------------------------------------------------------------------
        !>  === Add hinge motion: END ===
        ! -------------------------------------------------------------------------------
      end if
    end if

    end associate
  end do

end subroutine update_force

!----------------------------------------------------------------
!> Update structural forces, taking into account coupled hinges
subroutine update_force_coupled_hinge( this, comp, hinge, j_for, j_mom )
  class(t_precice)     , intent(inout) :: this
  type(t_geo_component), intent(inout) :: comp
  type(t_hinge)        , intent(inout) :: hinge
  integer              , intent(in)    :: j_for, j_mom

  integer :: i, j, n_nodes, nj
  integer :: a, h, b, ip
  real(wp) :: al_ah, w_ah, w_ab

  !> check input type
  if ( trim(comp%coupling_type) .ne. 'rbf' ) then
    write(*,*) ' Error in update_force_coupled_hinge( comp, hinge ). '
    write(*,*) ' So far, comp%coupling_type must be .eq. "rbf", but  '
    write(*,*) ' comp%coupling_type: ', trim(comp%coupling_type)
    write(*,*) ' Stop. '; stop
  end if

  !> Rot
  n_nodes = size( hinge%rot_cen%node_id )
  do i = 1, n_nodes

    a     = hinge%rot_cen%node_id(i)
    al_ah = hinge%rot_cen%span_wei(i)


    !> Update f_h, f_h += ...
    nj = size( hinge%rot_cen%ind, 1 )
    do j = 1, nj

      h    = hinge%rot_cen%ind(j,i)
      w_ah = hinge%rot_cen%wei(j,i)

      ip   = hinge%i_points_precice(h)


      !> Update f_h
      this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) &
                                     + al_ah * w_ah * comp%el(a)%dforce

      !> Update m_h
      this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) &
        + al_ah * w_ah * ( &
            comp%el(a)%dmom + &
            cross( comp%el(a)%cen - hinge%act%rr(:,h) , comp%el(a)%dforce ) )

    end do

    !> Update f_b, f_b -= ...
    nj = size( comp%rbf%cen%ind, 1 )
    do j = 1, nj

      b    = comp%rbf%cen%ind(j,a)
      w_ab = comp%rbf%cen%wei(j,a)

      ip   = comp%i_points_precice(b)


      !> Update f_b
      this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) &
                                     - al_ah * w_ab * comp%el(a)%dforce

      !> Update m_b
      this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) &
        - al_ah * w_ab * ( &
            comp%el(a)%dmom + &
            cross( comp%el(a)%cen - comp%rbf%rrb(:,b) , comp%el(a)%dforce ) )

    end do

  end do

  !> Blen
  n_nodes = size( hinge%blen_cen%node_id )
  do i = 1, n_nodes

    a     = hinge%blen_cen%node_id(i)
    al_ah = hinge%blen_cen%span_wei(i)

    !> Update f_h, f_h += ...
    nj = size( hinge%blen_cen%ind, 1 )
    do j = 1, nj

      h    = hinge%blen_cen%ind(j,i)
      w_ah = hinge%blen_cen%wei(j,i)

      ip   = hinge%i_points_precice(h)

      !> Update f_h
      this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) &
                                      + al_ah*w_ah*comp%el(a)%dforce

      !> Update m_h
      this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) &
        + al_ah * w_ah * ( &
            comp%el(a)%dmom + &
            cross( comp%el(a)%cen - hinge%act%rr(:,h) , comp%el(a)%dforce ) )

    end do

    !> Update f_b, f_b -= ...
    nj = size( comp%rbf%cen%ind, 1 )
    do j = 1, nj

      b    = comp%rbf%cen%ind(j,a)
      w_ab = comp%rbf%cen%wei(j,a)

      ip   = comp%i_points_precice(b)

      !> Update f_b
      this%fields(j_for)%fdata(:,ip) = this%fields(j_for)%fdata(:,ip) &
                                     - al_ah * w_ab * comp%el(a)%dforce

      !> Update m_b
      this%fields(j_mom)%fdata(:,ip) = this%fields(j_mom)%fdata(:,ip) &
        - al_ah * w_ab * ( &
            comp%el(a)%dmom + &
            cross( comp%el(a)%cen - comp%rbf%rrb(:,b) , comp%el(a)%dforce ) )

    end do

  end do

end subroutine update_force_coupled_hinge

!----------------------------------------------------------------
!> Update force/moment fields
subroutine update_elems( this, geo, elems, elems_virtual, te )
  class(t_precice)  , intent(inout)       :: this
  type(t_geo)       , intent(inout)       :: geo
  type(t_pot_elem_p), intent(inout)       :: elems(:)
  type(t_elem_virtual_p), intent(inout)   :: elems_virtual(:)
  type(t_tedge), optional, intent(inout)  :: te

  integer   :: i,j, i_comp, ip, iw, ih, ib, ii, it, il
  real(wp)  :: n_rot(3), chord(3), chord_rot(3), omega(3), pos(3), vel(3)
  real(wp)  :: r_drot(3), n_drot(3)
  real(wp)  :: theta
  real(wp)  :: Rot_mat(3,3)
  integer   :: j_pos, j_vel, j_rot, j_ome
  real(wp)  :: th1, yc, xq, yq, thp, xqp, yqp
  real(wp),  parameter :: eps = 1.0e-9_wp
  real(wp)             :: Rot(3,3), nx(3,3)
  
  ! Find rotation and angular velocity field id
  j_rot = 0; j_ome = 0
  do j = 1, size(this%fields)
    if ( trim(this%fields(j)%fname) .eq. 'Position' )       j_pos = j
    if ( trim(this%fields(j)%fname) .eq. 'Velocity' )       j_vel = j
    if ( trim(this%fields(j)%fname) .eq. 'Rotation' )       j_rot = j
    if ( trim(this%fields(j)%fname) .eq. 'AngularVelocity') j_ome = j
  end do

  if (present(te)) te%t_hinged = te%t
  !> Update elems
  ! *** to do *** build and exploit the connectivity preCICE-dust
  do i_comp = 1, size(geo%components)

    associate( comp => geo%components(i_comp) )

    if ( comp%coupling ) then

      !> rigid coupling --------------------------------------------------------
      if ( trim(comp%coupling_type) .eq. 'rigid' ) then

        !> === Coupling node ===
        !> Position
        pos = this%fields(j_pos)%fdata(:, comp%i_points_precice(1))

        !> Velocity
        vel = this%fields(j_vel)%fdata(:, comp%i_points_precice(1))

        !> Rotation
        n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(1))
        theta = norm2( n_rot )
        if ( theta .lt. eps ) then
          n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
          theta = 0.0_wp
        else
          n_rot = n_rot / theta
        end if

        !> Angular velocity of the point at the LE
        omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(1))

        !> === Grid nodes of the component ===
        ! Rigid motion of the component, defined by:
        ! - the motion of the coupling node,
        ! - the relative position of the component nodes and the
        !   coupling node.
        do i = 1, size(comp%c_ref_p,2) ! loop over comp points
          chord = comp%c_ref_p(:,i)
          chord_rot =  cos(theta) * chord + &
                       sin(theta) * cross( n_rot, chord ) + &
                     ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot

          !> Position and velocity of the nodes of the grid
          geo%points(    :, comp%i_points(i)) = pos + chord_rot
          geo%points_vel(:, comp%i_points(i)) = vel + cross( omega, chord_rot )

        end do

        !> === Control nodes of the elements ===
        ! *** to do *** avoid computing element quantities as the
        ! average value of node quantities
        do i = 1, size(comp%el)
          comp%el(i)%ub = 0.0_wp
          !> Compute the velocity of the element centre as the
          ! average value of the velocity of its nodes, by
          ! accumulation
          do j = 1, comp%el(i)%n_ver
            comp%el(i)%ub = comp%el(i)%ub + &
               1.0_wp / dble(comp%el(i)%n_ver) * &
                geo%points_vel(:, comp%el(i)%i_ver(j) )
          end do
          !> Velocity of the control point for LL components
          !> (exploit implicit connectivity of LL components)
          select type( el => comp%el(i) ); type is(t_liftlin)
            el%vel_ctr_pt = 0.5_wp * ( &
                  geo%points_vel(:, comp%i_points( 2*i-1 ) ) &
                + geo%points_vel(:, comp%i_points( 2*i+1 ) ) )
          end select
        end do

      !> rbf coupling ----------------------------------------------------------
      elseif ( trim(comp%coupling_type) .eq. 'rbf' ) then
        !> Nodes (rr)
        !> Reset, before accumulation, only nodes belonging to the component
        do i = 1, size(comp%i_points)
          ip = comp%i_points(i)
          geo%points    (:,ip) = 0.0_wp
          geo%points_vel(:,ip) = 0.0_wp
        end do
        
        !> Update surface quantities, as the weighted averages of the structure
        !  quantities, w/o considering rotations of the hinges
        do i = 1, size(comp%i_points)

          ip = comp%i_points(i)
          
          do iw = 1, size(comp%rbf%nod%ind,1)

            ! === Coupling Node ===
            !> Position of MBDyn nodes  
            pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%nod%ind(iw,i)))
            !write(*,*) 'pos', pos
            comp%rbf%rrb(:,comp%rbf%nod%ind(iw,i)) = pos
            !> Velocity
            vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%nod%ind(iw,i)))
            !> Rotation
            n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%nod%ind(iw,i)))
            comp%rbf%rrb_rot(:,comp%rbf%nod%ind(iw,i)) = n_rot
            
            theta = norm2( n_rot )
            if ( theta .lt. eps ) then
              n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
              theta = 0.0_wp
            else
              n_rot = n_rot / theta
            end if

            !> Angular velocity of the structural point
            omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%nod%ind(iw,i)))

            ! === Grid nodes of the components ===
            !> Reference difference
            chord = comp%loc_points(:,i) - comp%rbf%nodes(:,comp%rbf%nod%ind(iw,i))
            
            !> Rotated position difference
            chord_rot =  cos(theta) * chord + &
                         sin(theta) * cross( n_rot, chord ) + &
                       ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot
            
            !> Position
            geo%points(:, ip) = geo%points(:, ip) + &
                                comp%rbf%nod%wei(iw,i) * ( pos + chord_rot )
            !> Velocity
            geo%points_vel(:, ip) = geo%points_vel(:, ip) + &
                                    comp%rbf%nod%wei(iw,i) * ( vel + cross( omega, chord_rot ) )
          end do

        end do

        !> Points virtual 
        !> Reset, before accumulation, only nodes belonging to the component
        do i = 1, size(comp%i_points_virtual)
          ip = comp%i_points_virtual(i)
          geo%points_virtual(:,ip) = 0.0_wp
          geo%points_vel_virtual(:,ip) = 0.0_wp 
        end do
        
        !> Update surface quantities, as the weighted averages of the structure
        !  quantities, w/o considering rotations of the hinges
        do i = 1, size(comp%i_points_virtual)

          ip = comp%i_points_virtual(i)
          
          do iw = 1, size(comp%rbf%nod_virtual%ind,1)

            ! === Coupling Node ===
            !> Position of MBDyn nodes  
            pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%nod_virtual%ind(iw,i)))
            
            !> Velocity
            vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%nod_virtual%ind(iw,i)))
            !> Rotation
            n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%nod_virtual%ind(iw,i)))
            
            theta = norm2( n_rot )
            if ( theta .lt. eps ) then
              n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
              theta = 0.0_wp
            else
              n_rot = n_rot / theta
            end if

            !> Angular velocity of the structural point
            omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%nod_virtual%ind(iw,i)))

            ! === Grid nodes of the components ===
            !> Reference difference
            chord = comp%loc_points_virtual(:,i) - comp%rbf%nodes(:,comp%rbf%nod_virtual%ind(iw,i))
            
            !> Rotated position difference
            chord_rot =  cos(theta) * chord + &
                         sin(theta) * cross( n_rot, chord ) + &
                       ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot
            
            !> Position
            geo%points_virtual(:, ip) = geo%points_virtual(:, ip) + &
                                      comp%rbf%nod_virtual%wei(iw,i) * ( pos + chord_rot )
            !> Velocity
            geo%points_vel_virtual(:, ip) = geo%points_vel_virtual(:, ip) + &
                                    comp%rbf%nod_virtual%wei(iw,i) * ( vel + cross( omega, chord_rot ) )
          end do

        end do

        !> Center of panel 
        !> Reset, before accumulation, only nodes belonging to the component
        do i = 1, size(comp%el)
          comp%el(i)%cen = 0.0_wp
          comp%el(i)%ub = 0.0_wp
          comp%el(i)%ori = 0.0_wp
          comp%el(i)%R_cen = 0.0_wp
        end do

        !> Update surface quantities, as the weighted averages of the structure
        !  quantities, w/o considering rotations of the hinges
        
        do i = 1, size(comp%el) 
          
          do iw = 1, size(comp%rbf%cen%ind,1)

            ! === Coupling Node ===
            !> Position
            pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%cen%ind(iw,i)))
            !> Velocity
            vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%cen%ind(iw,i)))
            !> Rotation
            n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%cen%ind(iw,i)))

            theta = norm2( n_rot )
            if ( theta .lt. eps ) then
              n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
              theta = 0.0_wp
            else
              n_rot = n_rot / theta
            end if
            !> Angular velocity of the structural point
            omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%cen%ind(iw,i)))
            
            !> Reference difference
            chord = comp%loc_cen(:,i) - comp%rbf%nodes(:,comp%rbf%cen%ind(iw,i))

            !> Rotated position difference
            chord_rot = cos(theta) * chord + &
                        sin(theta) * cross(n_rot, chord ) + &
                        (1.0_wp - cos(theta)) * sum(chord*n_rot ) * n_rot

            !> Position
            comp%el(i)%cen = comp%el(i)%cen + &
                                comp%rbf%cen%wei(iw,i) * (pos + chord_rot)
            
            !> Orientation 
            comp%el(i)%ori = comp%el(i)%ori + &
                                comp%rbf%cen%wei(iw,i) * n_rot * theta  

                                
            !> Velocity
            comp%el(i)%ub = comp%el(i)%ub + &
                                comp%rbf%cen%wei(iw,i) * (vel + cross(omega, chord_rot))
            
          end do

            !> Orientation Matrix for surface panels 
            call vec2mat(comp%el(i)%ori, comp%el(i)%R_cen) 
            comp%el(i)%R_cen = matmul(comp%coupling_node_rot, comp%el(i)%R_cen)

        end do

        !> Center of panel virtual
        !> Reset, before accumulation, only nodes belonging to the component
        do i = 1, size(comp%el_virtual)
          comp%el_virtual(i)%cen = 0.0_wp
          comp%el_virtual(i)%ub = 0.0_wp
          comp%el_virtual(i)%ori = 0.0_wp
          comp%el_virtual(i)%R_cen = 0.0_wp
        end do

        !> Update surface quantities, as the weighted averages of the structure
        !  quantities, w/o considering rotations of the hinges
        
        do i = 1, size(comp%el_virtual) 
          
          do iw = 1, size(comp%rbf%cen_virtual%ind,1)

            ! === Coupling Node ===
            !> Position
            pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%cen_virtual%ind(iw,i)))
            !> Velocity
            vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%cen_virtual%ind(iw,i)))
            !> Rotation
            n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%cen_virtual%ind(iw,i)))

            theta = norm2( n_rot )
            if ( theta .lt. eps ) then
              n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
              theta = 0.0_wp
            else
              n_rot = n_rot / theta
            end if
            !> Angular velocity of the structural point
            omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%cen_virtual%ind(iw,i)))
            
            !> Reference difference
            chord = comp%loc_cen_virtual(:,i) - comp%rbf%nodes(:,comp%rbf%cen_virtual%ind(iw,i))

            !> Rotated position difference
            chord_rot = cos(theta) * chord + &
                        sin(theta) * cross(n_rot, chord ) + &
                        (1.0_wp - cos(theta)) * sum(chord*n_rot ) * n_rot

            !> Position
            comp%el_virtual(i)%cen = comp%el_virtual(i)%cen + &
                                comp%rbf%cen_virtual%wei(iw,i) * (pos + chord_rot)
            
            !> Orientation 
            comp%el_virtual(i)%ori = comp%el_virtual(i)%ori + &
                                comp%rbf%cen_virtual%wei(iw,i) * n_rot * theta  

            !> Velocity
            comp%el_virtual(i)%ub = comp%el_virtual(i)%ub + &
                                comp%rbf%cen_virtual%wei(iw,i) * (vel + cross(omega, chord_rot))
            
          end do

            !!> Orientation Matrix for surface panels 
            call vec2mat(comp%el_virtual(i)%ori, comp%el_virtual(i)%R_cen) 
            comp%el_virtual(i)%R_cen = matmul(comp%coupling_node_rot, comp%el_virtual(i)%R_cen)

        end do

        !> stripe ub 
        !> Reset, before accumulation, only nodes belonging to the component
        if (trim(comp%comp_el_type) .eq. 'v' .and. trim(comp%aero_correction) .eq. 'true') then 
          
          do i = 1, size(comp%stripe)
            comp%stripe(i)%ctr_pt = 0.0_wp
            comp%stripe(i)%ub  = 0.0_wp
          end do
          !> Update surface quantities, as the weighted averages of the structure
          !  quantities, w/o considering rotations of the hinges
          do i = 1, size(comp%stripe)
            do iw = 1, size(comp%rbf%ctr_pt%ind,1)
              ! === Coupling Node ===
              !> Position
              pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
              !> Velocity
              vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
              !> Rotation
              n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
              theta = norm2( n_rot )
            
              if ( theta .lt. eps ) then
                n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
                theta = 0.0_wp
              else
                n_rot = n_rot / theta
              end if

              !> Angular velocity of the structural point
              omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
              ! === Grid nodes of the components ===
              !> Reference difference
              chord = comp%stripe(i)%loc_ctr_pt - comp%rbf%nodes(:,comp%rbf%ctr_pt%ind(iw,i))
              
              !> Rotated position difference
              chord_rot =  cos(theta) * chord + &
                           sin(theta) * cross( n_rot, chord ) + &
                         ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot
  
              !> Position
              comp%stripe(i)%ctr_pt = comp%stripe(i)%ctr_pt + &
                                  comp%rbf%ctr_pt%wei(iw,i) * (pos + chord_rot)
              !> Velocity
              comp%stripe(i)%ub = comp%stripe(i)%ub + comp%rbf%ctr_pt%wei(iw,i) * (vel + cross(omega, chord_rot))
            end do
          end do

        elseif (trim(comp%comp_el_type) .eq. 'l') then

          do i = 1, size(comp%el)
            select type(el => comp%el(i))
              type is(t_liftlin) 
              el%ctr_pt = 0.0_wp
              el%vel_ctr_pt  = 0.0_wp
            end select
          end do
          !> Update surface quantities, as the weighted averages of the structure
          !  quantities, w/o considering rotations of the hinges
          do i = 1, size(comp%el)
            do iw = 1, size(comp%rbf%ctr_pt%ind,1)
              select type(el => comp%el(i))
                type is(t_liftlin) 
                ! === Coupling Node ===
                !> Position
                pos   = this%fields(j_pos)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
                !> Velocity
                vel   = this%fields(j_vel)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
                !> Rotation
                n_rot = this%fields(j_rot)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
                theta = norm2(n_rot)

                if ( theta .lt. eps ) then
                  n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
                  theta = 0.0_wp
                else
                  n_rot = n_rot / theta
                end if
                !> Angular velocity of the structural point
                omega = this%fields(j_ome)%fdata(:, comp%i_points_precice(comp%rbf%ctr_pt%ind(iw,i)))
                ! === Grid nodes of the components ===
                !> Reference difference
                chord = el%loc_ctr_pt - comp%rbf%nodes(:,comp%rbf%ctr_pt%ind(iw,i))

                !> Rotated position difference
                chord_rot =  cos(theta) * chord + &
                             sin(theta) * cross( n_rot, chord ) + &
                           ( 1.0_wp - cos(theta) ) * sum( chord*n_rot ) * n_rot
              
                !> Position
                el%ctr_pt = el%ctr_pt + comp%rbf%ctr_pt%wei(iw,i) * (pos + chord_rot)
                !> Velocity
                el%vel_ctr_pt = el%vel_ctr_pt + comp%rbf%ctr_pt%wei(iw,i) * (vel + cross(omega, chord_rot))

              end select
            end do
          end do

        endif 

        !!> *** to do *** move to a function: update_elems_coupled_hinges(), here in mod_precice
        !! -------------------------------------------------------------------------------
        !!>  === Add hinge motion: START ===
        ! -------------------------------------------------------------------------------

        do ih = 1, comp%n_hinges
          if ( trim(comp%hinge(ih)%input_type) .eq. 'coupling' ) then

            !> Update hinge nodes
            comp%hinge(ih)%act%rr = &
                this%fields(j_pos)%fdata(:,comp%hinge(ih)%i_points_precice)

            !> Update hinge node reference frame attached to the non-rotating structure
            comp%hinge(ih)%hin_rot = 0.0_wp
            
            do i = 1, size(comp%hinge(ih)%hin%ind,2)
              do j = 1, size(comp%hinge(ih)%hin%ind,1)
                n_rot = this%fields(j_rot)%fdata(:, &
                        comp%i_points_precice(comp%hinge(ih)%hin%ind(j,i)))
                comp%hinge(ih)%hin_rot(:,i) = comp%hinge(ih)%hin_rot(:,i) + &
                                                comp%hinge(ih)%hin%wei(j,i)*n_rot 
                
              end do

              !> Update h,v,n reference frame attached to the non-rotating hinge nodes
              ! From rotation vector to rotation matrix
              theta = norm2(comp%hinge(ih)%hin_rot(:,i))
              if ( theta .ne. 0.0_wp ) then

                Rot_mat(1,:) = ( 1.0_wp - cos(theta) ) * n_rot(1) * n_rot / theta**2.0_wp
                Rot_mat(2,:) = ( 1.0_wp - cos(theta) ) * n_rot(2) * n_rot / theta**2.0_wp
                Rot_mat(3,:) = ( 1.0_wp - cos(theta) ) * n_rot(3) * n_rot / theta**2.0_wp

                Rot_mat(1,:) = Rot_mat(1,:) + &
                  (/ cos(theta)*theta, -sin(theta)*n_rot(3), sin(theta)*n_rot(2)/)/theta
                Rot_mat(2,:) = Rot_mat(2,:) + &
                  (/ sin(theta)*n_rot(3), cos(theta)*theta   ,-sin(theta)*n_rot(1) /)/theta
                Rot_mat(3,:) = Rot_mat(3,:) + &
                  (/-sin(theta)*n_rot(2), sin(theta)*n_rot(1), cos(theta)*theta    /)/theta

                comp%hinge(ih)%act%h = matmul(Rot_mat, comp%hinge(ih)%ref%h )
  
                comp%hinge(ih)%act%v = matmul(Rot_mat, comp%hinge(ih)%ref%v )
                comp%hinge(ih)%act%n = matmul(Rot_mat, comp%hinge(ih)%ref%n )

              else
                comp%hinge(ih)%act%h = comp%hinge(ih)%ref%h
                comp%hinge(ih)%act%v = comp%hinge(ih)%ref%v
                comp%hinge(ih)%act%n = comp%hinge(ih)%ref%n

              end if

            end do

            !> 0. Update node position and velocity, before update with coupled hinge motion,
            !   r = ( 1 - alpha ) * r_beam + alpha * r_hinge =
            !     = r_beam + alpha * ( r_hinge - r_beam ) = r_bea m + alpha * dr
            !   v = ( 1 - alpha ) * v_beam + alpha * v_hinge =
            !     = v_beam + alpha * ( v_hinge - v_beam ) = v_beam + alpha * dv
            ! where:
            ! - <alpha> is the spanwise weight,
            ! - r,v_beam  the position and the velocity due to the motion of the structural
            !   part of the model,
            ! - r,v_hinge the position and the velocity due to the motion of the hinge nodes
            
            do i = 1, size( comp%hinge(ih)%rot%node_id )
              ip = comp%i_points( comp%hinge(ih)%rot%node_id(i) )
              geo%points(:,ip)     = geo%points(:,ip) * &
                                    ( 1.0_wp - comp%hinge(ih)%rot%span_wei(i))
              geo%points_vel(:,ip) = geo%points_vel(:,ip) * &
                                    ( 1.0_wp - comp%hinge(ih)%rot%span_wei(i))
            end do

            !> From motion of hinge nodes to surface motion
            ! ... see ~ geo/mod_hinges/hinge_deflection()
            !> 1. Update position
            do i = 1, comp%hinge(ih)%n_nodes ! hinge nodes

              !> === Hinge nodes ===
              !> Rotation vector and rotation matrix
              n_rot = this%fields(j_rot)%fdata(:,comp%hinge(ih)%i_points_precice(i))
              theta = norm2(n_rot)
              if ( theta .lt. eps ) then
                n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
                theta = 0.0_wp
              else
                n_rot = n_rot / theta
              end if

              nx(1,:) = (/    0.0_wp, -n_rot(3),  n_rot(2) /)
              nx(2,:) = (/  n_rot(3),    0.0_wp, -n_rot(1) /)
              nx(3,:) = (/ -n_rot(2),  n_rot(1),    0.0_wp /)

              Rot = reshape( (/1.0_wp, 0.0_wp, 0.0_wp, &
                              0.0_wp, 1.0_wp, 0.0_wp, &
                              0.0_wp, 0.0_wp, 1.0_wp /), (/3,3/) ) &
                  - sin(theta) * nx + ( 1.0_wp - cos(theta) ) * matmul(nx, nx)

              !> Evaluate hinge deflection, theta, from the relative position of the rotating
              ! and non-rotating ref.frames
              call rotation_vector_combination( &
                      this%fields(j_rot)%fdata(:, &
                      comp%hinge(ih)%i_points_precice( i ) ), &
                      -comp%hinge(ih)%hin_rot(:,i), r_drot, theta, n_drot )

              !> === Surface nodes ===
              !> 1.1. Update points: rigid rotation
              do ib = 1, size(comp%hinge(ih)%rot%n2h(i)%p2h)

                !> Reference difference
                ip = comp%hinge(ih)%rot%n2h(i)%p2h(ib)  ! Local numbering
                chord = comp%loc_points(:,ip) - comp%hinge(ih)%nodes(:,i)

                ip = comp%i_points(ip)  ! Local-to-global connectivity
                !> Position: absolute position (!)
                geo%points(:,ip) = geo%points(:,ip) + &
                     comp%hinge(ih)%rot%n2h(i)%s2h(ib) * &
                     comp%hinge(ih)%rot%n2h(i)%w2h(ib) * &
                    ( comp%hinge(ih) % act % rr(:,i) + matmul( transpose(Rot), chord ) )
              end do

              !> 1.2. Chordwise blending region
              do ib = 1, size(comp%hinge(ih)%blen%n2h(i)%p2h)

                ii = comp%hinge(ih)%blen%n2h(i)%p2h(ib)
                ip = comp%i_points(ii)  ! Local-to-global connectivity

                ! condition to identify if the hinge is moving up/down
                if (n_drot(2) .lt. 0) then
                  th1 = theta * comp%hinge(ih)%blen%n2h(i)%s2h(ib)
                else
                  th1 = -theta * comp%hinge(ih)%blen%n2h(i)%s2h(ib)
                endif

                if ( th1 .ne. 0.0_wp ) then
                  !> coordinate of the centre of the circle used for blending,
                  ! in the n-direction
                  yc = cos(th1)/sin(th1)*comp%hinge(ih)%offset*(1.0_wp + cos(th1)) + &
                                        comp%hinge(ih)%offset*sin(th1)
                  !> Some auxiliary quantities
                  xq = sum( ( geo%points(:,ip) - comp%hinge(ih)%act%rr(:,i) ) * &
                                                comp%hinge(ih)%act%v( :,i) )
                  yq = sum( ( geo%points(:,ip) - comp%hinge(ih)%act%rr(:,i) ) * &
                                                comp%hinge(ih)%act%n( :,i) )
                  thp = 0.5_wp * ( xq + comp%hinge(ih)%offset ) / comp%hinge(ih)%offset * th1
                  xqp = yc*sin(thp)   - comp%hinge(ih)%offset - yq*sin(thp) - xq
                  yqp = yc*(1.0_wp-cos(thp))                  + yq*cos(thp) - yq

                  !> Update coordinates
                  geo%points(:,ip) = geo%points(:,ip) + &
                             comp%hinge(ih)%blen%n2h(i)%s2h(ib) * &
                             comp%hinge(ih)%blen%n2h(i)%w2h(ib) * &
                           ( xqp * comp%hinge(ih)%act%v(:,i) + &
                             yqp * comp%hinge(ih)%act%n(:,i) )
                ! else do nothing
                end if

              end do

            end do
            
            if (present(te)) then 
              !> Update trailing edge direction
              do it = 1, size(te%t_hinged,2)
                do i = 1, comp%hinge(ih)%n_nodes ! hinge nodes

                  !> === Hinge nodes ===
                  !> Rotation vector and rotation matrix
                  n_rot = this%fields(j_rot)%fdata(:,comp%hinge(ih)%i_points_precice(i))
                  theta = norm2(n_rot)

                  if ( theta .lt. eps ) then
                  n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /); theta = 0.0_wp
                  else
                    n_rot = n_rot / theta
                  end if

                  nx(1,:) = (/    0.0_wp, -n_rot(3),  n_rot(2) /)
                  nx(2,:) = (/  n_rot(3),    0.0_wp, -n_rot(1) /)
                  nx(3,:) = (/ -n_rot(2),  n_rot(1),    0.0_wp /)
                
                  do ib = 1, size(comp%hinge(ih)%rot%n2h(i)%p2h)

                    ip = comp%hinge(ih)%rot%n2h(i)%p2h(ib)  ! Local numbering

                    
                    il = comp%i_points(ip)  ! Node of the hinge region        

                    th1 = theta * comp%hinge(ih)%rot%n2h(i)%w2h(ib) 
                    Rot = sin(th1) * nx + ( 1.0_wp - cos(th1) ) * matmul( nx, nx )
                  
                    if (te%i(1,it) .eq. il) then ! hinge node is also trailing edge node

                      te%is_hinged(it) = te%i(1,it)
                      te%t_hinged(:,it) = te%t_hinged(:,it) + comp%hinge(ih)%rot%n2h(i)%s2h(ib) &
                                        * matmul( Rot, te%t_hinged(:,it))   
                    
                      te%t_hinged(:,it) = te%t_hinged(:,it)/norm2(te%t_hinged(:,it))

                    end if

                  end do

                end do
              end do
            end if

            !> 2. Update velocity with weighted rigid motion, after the new position
            ! of the nodes has been evaluated
            do i = 1, comp%hinge(ih)%n_nodes ! hinge nodes

              !> === Hinge nodes ===
              !> Velocity and Angular Velocity of hinge nodes
              vel   = this%fields(j_vel)%fdata(:,comp%hinge(ih)%i_points_precice(i))
              omega = this%fields(j_ome)%fdata(:,comp%hinge(ih)%i_points_precice(i))

              !> === Surface nodes ===
              !> 2.1. Rigid rotation region
              do ib = 1, size(comp%hinge(ih)%rot%n2h(i)%p2h)

                !> Connectivity
                ip = comp%hinge(ih)%rot%n2h(i)%p2h(ib)  ! Local numbering
                il = comp%i_points(ip)                  ! Local-to-global connectivity
                !> Velocity: v = v_H + Omega_H x ( r - r_H )
                geo%points_vel(:,il) = geo%points_vel(:,il) + &
                     comp%hinge(ih)%rot%n2h(i)%s2h(ib) * &
                     comp%hinge(ih)%rot%n2h(i)%w2h(ib) * &
                    ( vel + cross( omega, &
                          geo%points(:,il) - comp%hinge(ih)%act%rr(:,i) ) )

              end do

              !!> 2.2. Chordwise blending region
              !do ib = 1, size(comp%hinge(ih)%blen%n2h(i)%p2h)
              !
              !  ! *** to do ***
              !
              !end do

            end do

          else
            write(*,*) ' comp%hinge(ih)%input_type must be equal to <coupling>. Stop '; stop
          end if
        end do

      else !
        call error('update_elems','mod_precice', &
                  ' Wrong CouplingType: '//trim(comp%coupling_type)// &
                  ' for component: '//trim(comp%comp_name)// &
                  '. So far, available CouplingType inputs are: rigid, rbf.'// &
                  ' Stop.'); stop
      end if

    end if ! if coupling
    end associate
  end do

end subroutine update_elems

!----------------------------------------------------------------
!> Update near field wake
subroutine update_near_field_wake( this, geo, wake, te )
  class(t_precice)  , intent(inout) :: this
  type(t_geo)       , intent(in)    :: geo
  type(t_wake)      , intent(inout) :: wake
  type(t_tedge)     , intent(inout) :: te

  real(wp) :: n_rot(3), dist(3), vel_te(3), wind(3)
  real(wp) :: theta
  real(wp) :: eps = 1.0e-9_wp

  integer :: icomp
  integer :: ip, ir, p1, p2, j, j_rot, iw, i

  ! Store old wake points to use in complete_wake
  if (wake%update_old_second_row) then
    wake%old_second_row = wake%pan_w_points(:,:,2)
    wake%update_old_second_row = .false.
  end if

  ! Find rotation and angular velocity field id
  j_rot = 0
  do j = 1, size(this%fields)
    if ( trim(this%fields(j)%fname) .eq. 'Rotation' )  j_rot = j
  end do

  !> First row of points of the panel wake
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%pan_gen_points(1,:)) + &
                                  geo%points(:,wake%pan_gen_points(2,:)))

  wake%pan_w_points(:,:,1) = wake%w_start_points

  !> Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1,wake%n_pan_points

    if ( geo%components( wake%pan_gen_icomp(ip) )%coupling ) then

      if ( trim(geo%components( wake%pan_gen_icomp(ip) )%coupling_type) &
              .eq. 'rigid' ) then

        !> rigid coupling
        !> Rotation
        n_rot = this%fields(j_rot)%fdata(:, &
                geo%components( wake%pan_gen_icomp(ip) )%i_points_precice(1) )
        theta = norm2( n_rot )
        if ( theta .lt. eps ) then;  n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /);  theta = 0.0_wp
        else                      ;  n_rot = n_rot / theta
        end if

        dist =  cos(theta) * wake%pan_gen_dir(:,ip) + &
                sin(theta) * cross( n_rot, wake%pan_gen_dir(:,ip) ) + &
              ( 1.0_wp - cos(theta) ) * sum( wake%pan_gen_dir(:,ip)*n_rot ) * n_rot

        vel_te = geo%points_vel(:, wake%pan_gen_icomp(ip))

        wind = variable_wind(geo%points(:,wake%pan_gen_icomp),sim_param%time)
        
        if ( norm2(wind-vel_te) .gt. sim_param%min_vel_at_te ) then
          
          if (sim_param%autoscale_te .and. te%id_lifting_line(ip) .ne. 1) then 
            wake%pan_gen_scaling(ip) = te%scaling(ip) / &
                                  norm2(wind-vel_te) 
          endif 

          wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
              dist*wake%pan_gen_scaling(ip)* &
              norm2(wind-vel_te)*sim_param%dt / norm2(dist) * &
              real(sim_param%ndt_update_wake,wp)
        else
          if (sim_param%autoscale_te .and. te%id_lifting_line(ip) .ne. 1) then 
            wake%pan_gen_scaling(ip) = te%scaling(ip) / &
                                  sim_param%min_vel_at_te 
          endif
          wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
              dist*wake%pan_gen_scaling(ip) * & ! next line may be commented
              sim_param%min_vel_at_te*sim_param%dt * &
              real(sim_param%ndt_update_wake,wp)
        end if

      elseif (trim(geo%components( wake%pan_gen_icomp(ip) )%coupling_type) &
              .eq. 'rbf' ) then

        !> rbf coupling (general coupling)
        !> Rotation
        icomp = wake%pan_gen_icomp(ip)
        associate( comp => geo%components(icomp) )
          ! trailing edge points ID in component index 
          iw = wake%pan_gen_points(1,ip) - minval(comp%i_points) + 1          
          ! find if the trailing edge node belongs to a hinge
          ! if so, it has already been rotated in update_geometry
          ! and we can use t_hinged
          if (any(te%is_hinged .eq. wake%pan_gen_points(1,ip))) then
            ! iw is the node index, but we need to find where that node is 
            ! in the trailing edge array
            ! (workaround for findloc)
            do i=1, size(te%is_hinged)
              if (te%is_hinged(i) .eq. wake%pan_gen_points(1,ip)) then              
                dist = te%t_hinged(:,i)
              end if
            end do
          else
            ! *** to do ***
            ! So far, t_te inherits orientation ONLY from the closest coupling node,
            ! without interpolation with rbf coefficients
            n_rot = this%fields(j_rot)%fdata(:, &
                    comp%i_points_precice( comp%rbf%nod%ind(1,iw) ) )
            theta = norm2( n_rot )

            if ( theta .lt. eps ) then
              n_rot = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
              theta = 0.0_wp
            else
              n_rot = n_rot / theta
            end if

               dist =  cos(theta) * wake%pan_gen_dir(:,ip) + &
                       sin(theta) * cross( n_rot, wake%pan_gen_dir(:,ip) ) + &
                     ( 1.0_wp - cos(theta) ) * sum( wake%pan_gen_dir(:,ip)*n_rot ) * n_rot
          end if
          ! trailing edge velocity (global numbering)
          vel_te = geo%points_vel(:, wake%pan_gen_points(1,ip))

          wind = variable_wind(geo%points(:,wake%pan_gen_icomp),sim_param%time)
        
          if ( norm2(wind-vel_te) .gt. sim_param%min_vel_at_te ) then
            if (sim_param%autoscale_te .and. te%id_lifting_line(ip) .ne. 1) then 
              wake%pan_gen_scaling(ip) = te%scaling(ip) / &
                                    norm2(wind-vel_te) 
            endif
            wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                                        dist*wake%pan_gen_scaling(ip)* &
                                        norm2(wind-vel_te)*sim_param%dt / norm2(dist) * &
                                        real(sim_param%ndt_update_wake,wp)
          
          else
            if (sim_param%autoscale_te .and. te%id_lifting_line(ip) .ne. 1) then
              wake%pan_gen_scaling(ip) = te%scaling(ip) / &
                                      sim_param%min_vel_at_te
            endif
            wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                dist*wake%pan_gen_scaling(ip) * & ! next line may be commented
                sim_param%min_vel_at_te*sim_param%dt * &
                real(sim_param%ndt_update_wake,wp)
          
          end if
        
        end associate

      else
        !> other coupling
        ! not implemented, so far
      end if

    end if

  enddo

  !==> Check if the panels need to be joined
  if(sim_param%join_te) then

    wake%joined_tes(:,:,2:wake%nmax_pan+1) = &
          wake%joined_tes(:,:,1:wake%nmax_pan)
    wake%joined_tes(:,:,1) = 0
    call join_first_panels(wake,sim_param%join_te_factor)
  endif  

  ! Calculate geometrical quantities of first 2 rows
  do ip = 1,wake%n_pan_stripes
    do ir = 1, wake%pan_wake_len
      p1 = wake%i_start_points(1,ip)
      p2 = wake%i_start_points(2,ip)
      call wake%wake_panels(ip,ir)%calc_geo_data( &
        reshape((/wake%pan_w_points(:,p1,ir),   wake%pan_w_points(:,p2,ir), &
                  wake%pan_w_points(:,p2,ir+1), wake%pan_w_points(:,p1,ir+1)/),&
                                                                    (/3,4/)))
    enddo
  enddo

end subroutine update_near_field_wake

end module mod_precice
