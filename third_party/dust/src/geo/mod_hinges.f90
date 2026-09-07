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

!> Module for introducing hinges and rotating parts in a component
module mod_hinges

use mod_param, only: &
  wp, max_char_len, pi

use mod_math, only: &
  cross, sort_vector_real

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, getsuboption, &
  countoption, finalizeparameters

#if USE_PRECICE
use mod_precice_rbf, only: &
  t_precice_rbf
#endif

implicit none

public :: t_hinge, t_hinge_input, build_hinges, hinge_input_parser, &
          initialize_hinge_config

private

! ---------------------------------------------------------------
!> Hinge input data type
type :: t_hinge_input
  character(len=max_char_len) :: tag
  character(len=max_char_len) :: nodes_input
  real(wp) :: node1(3), node2(3)
  integer  :: n_nodes
  real(wp), allocatable :: rr(:,:)
  character(len=max_char_len) :: node_file
  real(wp) :: ref_dir(3)
  real(wp) :: offset
  real(wp) :: span_blending
  character(len=max_char_len) :: rotation_input
  real(wp) :: rotation_amplitude
  real(wp) :: rotation_omega
  real(wp) :: rotation_phase
  real(wp) :: rotation_amplitude_init
  real(wp) :: rotation_cosine_cycl
  real(wp) :: rotation_initial_time
  integer, allocatable :: coupling_nodes(:)
  real(wp) :: le1(2) 
  real(wp) :: te1(2) 
  real(wp) :: le2(2)
  real(wp) :: te2(2)
  real(wp) :: chord1, chord2
  real(wp) :: csi1, csi2
  real(wp) :: merge_tol
  logical  :: adaptive_mesh
end type t_hinge_input

! ---------------------------------------------------------------
!> Hinge node to surface node connectivity. The different number of surface
! nodes per hinge nodes requires an array of obj, containing indices and
! weights
type :: t_n2h_conn
  !> Indices of the surface nodes
  integer, allocatable :: p2h(:)
  !> Weights
  real(wp),allocatable :: w2h(:)
  !> Spanwise weights
  real(wp),allocatable :: s2h(:)
end type t_n2h_conn

! ---------------------------------------------------------------
!> Hinge connectivity, meant for rigid rotation and bleding regions
type :: t_hinge_conn

  !> Local index of the nodes (or cell center) performing the desired motion
  integer, allocatable :: node_id(:)
  !> Surface node performing motion vs. hinge node connectivity
  integer , allocatable :: ind(:,:)
  !> Surface node performing motion vs. hinge node connectivity, weights
  ! for the weighted average of the motion
  real(wp), allocatable :: wei(:,:)
  !> Surface node performing motion vs. hinge node connectivity, weights
  ! for the weighted average of the motion in the spanwise direction
  real(wp), allocatable :: span_wei(:)
  !> Node to hinge connectivity array of objs
  type(t_n2h_conn), allocatable :: n2h(:)

end type t_hinge_conn

! ---------------------------------------------------------------
!> Hinge node configurations: to be used in defining reference
! and actual configurations
type :: t_hinge_config

  !> Position of the first and last point of the hinge in the local
  ! reference frame of the component
  real(wp) :: rr0(3), rr1(3)
  !> Local coordinates of the hinge nodes
  real(wp), allocatable :: rr(:,:)

  !> Unit vectors of the hinge node reference frames
  ! User input in local reference frame, for reference configuration: h, v
  !    h,
  !    n = cross(v,h) , n = n/norm2(n)
  !    v = cross(h,n)
  !> h: rotation axis
  real(wp), allocatable :: h(:,:)
  !> v: zero direction (user inputs are overwritten, in order to
  ! build ortonormal reference frameshinge%rot
  real(wp), allocatable :: v(:,:)
  !> n: normal direction
  real(wp), allocatable :: n(:,:)

  ! !> Spanwise weight to avoid irregular behaviors
  ! real(wp) :: span_wei

end type t_hinge_config

! ---------------------------------------------------------------
!> Hinge type
type :: t_hinge

  !> Number of points for transferring the motion from the hinge to the
  ! surface, through a weighted average
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 2 (or 1?)
  integer :: n_wei = 2

  !> Order of the norm used for computing distance-based weights
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 1 (or 2?)
  real(wp) :: w_order = 1

  !> Type: parametric, from_file
  character(len=max_char_len) :: nodes_input

  !> N.nodes of the hinge
  integer :: n_nodes

  !> Offset for avoiding irregular behavior
  real(wp) :: offset

  !> Spanwise blending to avoid irregular behaviors
  real(wp) :: span_blending

  !> Vector identifying zero-rotation direction
  real(wp) :: ref_dir(3)

  !> Type: constant, function, coupling
  character(len=max_char_len) :: input_type

  !> function:const, :sin, :cos amplitude, angular velocity and phase
  real(wp) :: f_ampl
  real(wp) :: f_omega
  real(wp) :: f_phase
  real(wp) :: f_ampl_init
  real(wp) :: f_cosine_cycl
  real(wp) :: f_initial_time

  !> Coupling nodes ( for Hinge_Rotation_input = coupling )
  ! replicating geo%components(i_comp)%i_points_precice(:)
  !             geo%components(i_comp)%rbf%nodes(:,:)
  integer , allocatable :: i_points_precice(:)
  real(wp), allocatable :: nodes(:,:)
  integer , allocatable :: i_coupling_nodes(:)

  !> Array of the rotation angle (read as an input, or from coupling nodes)
  real(wp), allocatable :: theta(:)
  real(wp), allocatable :: theta_old(:)

  !> Reference configuration
  type(t_hinge_config) :: ref
  !> Actual configuration
  type(t_hinge_config) :: act
  !> Actual configuration of the nodes attached to the non-rotating structure
  ! (defined and used only for coupled hinges, to determine the rotation of the
  !  rotating surface, w.r.t. the structure, and evaluate the motion of the
  !  surface in the blending region)
  type(t_hinge_config) :: fix

  !> Connectivity for the rigid rotation motion of a ``control surface'',
  ! both for grid nodes and cell centers
  type(t_hinge_conn) :: rot , rot_cen
  !> Connectivity of the blending region to avoid irregular behavior
  ! with a ``control surface'' large rotation (blending region extends
  ! from -offset to offset qualitatively in the ref_dir of the hinge)
  type(t_hinge_conn) :: blen, blen_cen
  !> Connectivity between hinge nodes and non-rotating structural nodes,
  ! read from mbdyn through preCICE
  type(t_hinge_conn) :: hin

  character(len=max_char_len) ::  tag 
  !> Interpolated rotation vector of the hinge nodes, attached to the non-
  ! rotating part of the component
  real(wp), allocatable :: hin_rot(:,:)


  contains

  procedure, pass(this) :: build_connectivity
  procedure, pass(this) :: build_connectivity_cen
  procedure, pass(this) :: build_connectivity_hin
  procedure, pass(this) :: init_theta
  procedure, pass(this) :: from_reference_to_actual_config ! empty (useless?)
  procedure, pass(this) :: update_hinge_nodes
  procedure, pass(this) :: update_theta
  procedure, pass(this) :: hinge_deflection

end type t_hinge

! ---------------------------------------------------------------

contains

! ---------------------------------------------------------------
!> Build hinge connectivity from component nodes rr, expressed in the
! local reference frame, and the hinge nodes.
subroutine build_connectivity(this, loc_points, coupling_node_rot)
  class(t_hinge), intent(inout) :: this
  real(wp),       intent(in)    :: loc_points(:,:)
  real(wp),       intent(in)    :: coupling_node_rot(3,3)

  real(wp) :: hinge_width
  integer  :: nb, nh, ib, ih, iw
  real(wp), allocatable :: rrb(:,:), rrh(:,:), rrb_wei(:,:)
  real(wp) :: Rot(3,3) = 0.0_wp

  real(wp) :: span_wei, tol 
  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::              ind_v(:)

  integer , allocatable :: rot_node_id(:) , ble_node_id(:)
  integer , allocatable :: rot_ind(:,:)   , ble_ind(:,:)
  real(wp), allocatable :: rot_wei(:,:)   , ble_wei(:,:)
  real(wp), allocatable :: rot_span_wei(:), ble_span_wei(:)

  integer , allocatable :: rot_p2h(:,:)   
  real(wp), allocatable :: rot_w2h(:,:)   
  integer , allocatable :: rot_i2h(:)     
  real(wp), allocatable :: rot_s2h(:,:)
  integer , allocatable :: ble_p2h(:,:)   
  real(wp), allocatable :: ble_w2h(:,:)   
  integer , allocatable :: ble_i2h(:)     
  real(wp), allocatable :: ble_s2h(:,:)

  integer :: nrot, nble
  real(wp) :: wei_hinge, x_hinge, z_hinge


  allocate(wei_v(this%n_wei))
  allocate(ind_v(this%n_wei))

  !> N. of surface and hinge nodes
  nb = size(loc_points,2)
  nh = this%n_nodes

  !> Coordinates in the hinge reference frame
  ! Rotation matrix, build with the local ortonormal ref.frame of
  ! the first hinge node
  Rot(1,:) = this%ref%v(:,1)
  Rot(2,:) = this%ref%h(:,1)
  Rot(3,:) = this%ref%n(:,1)

  allocate( rrb(3,nb) );  allocate( rrh(3,nh) )
  allocate( rrb_wei(3,nb) );
  !rotate the hinge nodes (from coupling_nodes file) to recover the value
  !of the component.in input
  do ih = 1, nh 
    this%ref%rr(:,ih) = matmul( (coupling_node_rot),this%ref%rr(:,ih))
  end do 
  !rotate the dust mesh points in dust reference
  do ib = 1, nb
    rrb(:,ib) =  matmul( (coupling_node_rot), loc_points(:,ib))
  end do

  do ib = 1, nh
    rrh(:,ib) =  this%ref%rr(:,ib) - this%ref%rr(:,1)
  end do

  ! hinge width, measured in the hinge direction
  hinge_width = rrh(2,nh) - rrh(2,1)

  !> Compute connectivity and weights
  ! Allocate auxiliary node_id(:), ind(:,:), wei(:,:) arrays
  allocate(rot_node_id(        nb));  rot_node_id = 0
  allocate(ble_node_id(        nb));  ble_node_id = 0
  allocate(rot_ind(this%n_wei, nb));      rot_ind = 0
  allocate(ble_ind(this%n_wei, nb));      ble_ind = 0
  allocate(rot_wei(this%n_wei, nb));      rot_wei = 0.0_wp
  allocate(ble_wei(this%n_wei, nb));      ble_wei = 0.0_wp
  allocate(rot_span_wei(       nb)); rot_span_wei = 0.0_wp
  allocate(ble_span_wei(       nb)); ble_span_wei = 0.0_wp

  ! diff_all and dist_all auxiliaary arrays
  allocate(dist_all(   nh)); dist_all = 0.0_wp

  ! auxiliary matrix for hinge to surf connectivity
  ! ***** to do ***** improve the implementation
  allocate(rot_p2h(nh,nb)); rot_p2h = 0      ;  allocate(ble_p2h(nh,nb)); ble_p2h = 0
  allocate(rot_w2h(nh,nb)); rot_w2h = 0.0_wp ;  allocate(ble_w2h(nh,nb)); ble_w2h = 0.0_wp
  allocate(rot_i2h(nh   )); rot_i2h = 0      ;  allocate(ble_i2h(nh   )); ble_i2h = 0
  allocate(rot_s2h(nh,nb)); rot_s2h = 0.0_wp ;  allocate(ble_s2h(nh,nb)); ble_s2h = 0.0_wp

  nrot = 0; nble = 0
  ! rrb_wei are the rrb re-rotate in the same reference of rrh in order to evaluate the weight
  ! while rrb are used to evaluate the flap region in the wind axis
  do ib = 1, nb
    rrb_wei(:,ib) = matmul( coupling_node_rot, &
                            (loc_points(:,ib) - &
                    matmul(transpose(coupling_node_rot) ,this%ref%rr(:,1)) ))
  enddo
  tol = abs(this%ref%rr(2,nh) - this%ref%rr(2,1))/1000.0_wp  ! tolerance for the hinge width (relative)
  ! Loop over all the surface points
  do ib = 1, nb
    !> wind axis conditions
    if ((rrb(2,ib) .ge. (this%ref%rr(2,1) - tol)) .and. (rrb(2,ib) .le. (this%ref%rr(2,nh)+tol))) then
      
      wei_hinge = (rrb(2,ib) - this%ref%rr(2,1)) / (this%ref%rr(2,nh)- this%ref%rr(2,1))
      x_hinge = this%ref%rr(1,1) + wei_hinge*(this%ref%rr(1,nh)- this%ref%rr(1,1))
      z_hinge = this%ref%rr(3,1) + wei_hinge*(this%ref%rr(3,nh)- this%ref%rr(3,1))

      if (rrb(1,ib) .gt. (this%offset + x_hinge)) then

        nrot = nrot + 1
        rot_node_id(nrot) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb_wei(2,ib) - rrh(2,ih) )
        end do

        !> Weights in chordwise direction
        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

        wei_v = 1.0_wp / max( wei_v, 1e-9_wp )**this%w_order
        wei_v = wei_v / sum(wei_v)

        !> Weights in spanwise direction
        if ( rrb(2,ib) .lt. (this%ref%rr(2,1) + this%span_blending) ) then
          span_wei = 1.0_wp + rrb_wei(2,ib) / this%span_blending
        elseif( rrb(2,ib) .lt. (this%ref%rr(2,nh) - this%span_blending)  ) then
          span_wei = 1.0_wp
        else
          span_wei = 1.0_wp - ( rrb_wei(2,ib) - hinge_width ) / this%span_blending
        endif

        rot_wei(   :,nrot) = wei_v
        rot_ind(   :,nrot) = ind_v
        rot_span_wei(nrot) = span_wei

        ! *****
        do iw = 1, this%n_wei
          rot_i2h(ind_v(iw)) = rot_i2h(ind_v(iw)) + 1
          rot_p2h(ind_v(iw),   rot_i2h(ind_v(iw))) = ib
          rot_w2h(ind_v(iw),   rot_i2h(ind_v(iw))) = wei_v(iw)
          rot_s2h(ind_v(iw),   rot_i2h(ind_v(iw))) = span_wei
        end do
        ! *****

      elseif (rrb(1,ib) .lt. (this%offset + x_hinge) .and. rrb(1,ib) .gt. (x_hinge - this%offset))  then ! blending region

        nble = nble + 1
        ble_node_id(nble) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb_wei(2,ib) - rrh(2,ih) )
        end do

        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

        !> Weights in chordwise direction
        wei_v = 1.0_wp / max( wei_v, 1e-9_wp ) **this%w_order
        wei_v = wei_v / sum(wei_v)

        !> Weights in spanwise direction
        if ( rrb(2,ib) .lt. (this%ref%rr(2,1) + this%span_blending) ) then
          span_wei = 1.0_wp + rrb_wei(2,ib) / this%span_blending
        elseif( rrb(2,ib) .lt. (this%ref%rr(2,nh) - this%span_blending)  ) then
          span_wei = 1.0_wp
        else
          span_wei = 1.0_wp - ( rrb_wei(2,ib) - hinge_width ) / this%span_blending
        endif

        ble_wei(   :,nble) = wei_v
        ble_ind(   :,nble) = ind_v
        ble_span_wei(nble) = span_wei

        ! *****
        do iw = 1, this%n_wei
          ble_i2h(ind_v(iw)) = ble_i2h(ind_v(iw)) + 1
          ble_p2h(ind_v(iw),   ble_i2h(ind_v(iw))) = ib
          ble_w2h(ind_v(iw),   ble_i2h(ind_v(iw))) = wei_v(iw)
          ble_s2h(ind_v(iw),   ble_i2h(ind_v(iw))) = span_wei
        end do
        ! *****
      else
      end if

    else
    end if

  end do

  !> Fill hinge object, with the connectivity and weight arrays
  allocate(this%rot %node_id(        nrot)); this%rot %node_id = rot_node_id( 1:nrot)
  allocate(this%rot %ind(this%n_wei, nrot)); this%rot %ind     = rot_ind(  :, 1:nrot)
  allocate(this%rot %wei(this%n_wei, nrot)); this%rot %wei     = rot_wei(  :, 1:nrot)
  allocate(this%rot %span_wei(       nrot)); this%rot %span_wei= rot_span_wei(1:nrot)
  allocate(this%blen%node_id(        nble)); this%blen%node_id = ble_node_id( 1:nble)
  allocate(this%blen%ind(this%n_wei, nble)); this%blen%ind     = ble_ind(  :, 1:nble)
  allocate(this%blen%wei(this%n_wei, nble)); this%blen%wei     = ble_wei(  :, 1:nble)
  allocate(this%blen%span_wei(       nble)); this%blen%span_wei= ble_span_wei(1:nble)
  allocate(this%rot %n2h(nh))
  allocate(this%blen%n2h(nh))
  do ih = 1, nh
    allocate(this%rot %n2h(ih)%p2h( rot_i2h(ih) )) ; this%rot %n2h(ih)%p2h = rot_p2h( ih, 1:rot_i2h(ih) )
    allocate(this%rot %n2h(ih)%w2h( rot_i2h(ih) )) ; this%rot %n2h(ih)%w2h = rot_w2h( ih, 1:rot_i2h(ih) )
    allocate(this%rot %n2h(ih)%s2h( rot_i2h(ih) )) ; this%rot %n2h(ih)%s2h = rot_s2h( ih, 1:rot_i2h(ih) )
    allocate(this%blen%n2h(ih)%p2h( ble_i2h(ih) )) ; this%blen%n2h(ih)%p2h = ble_p2h( ih, 1:ble_i2h(ih) )
    allocate(this%blen%n2h(ih)%w2h( ble_i2h(ih) )) ; this%blen%n2h(ih)%w2h = ble_w2h( ih, 1:ble_i2h(ih) )
    allocate(this%blen%n2h(ih)%s2h( ble_i2h(ih) )) ; this%blen%n2h(ih)%s2h = ble_s2h( ih, 1:ble_i2h(ih) )
  end do

  !> Explicit deallocations
  deallocate(rrb, rrh, dist_all, wei_v, ind_v)
  deallocate(rot_node_id, rot_ind, rot_wei, rot_span_wei)
  deallocate(ble_node_id, ble_ind, ble_wei, ble_span_wei)
  deallocate(rot_i2h, rot_p2h, rot_w2h, rot_s2h)
  deallocate(ble_i2h, ble_p2h, ble_w2h, ble_s2h)

end subroutine build_connectivity

! ---------------------------------------------------------------
!> Build connectivity betweeen hinge nodes and cell centers of the
! aerodynamic grid
subroutine build_connectivity_cen(this, rr, ee, coupling_node_rot)
  class(t_hinge), intent(inout) :: this
  real(wp),       intent(in)    :: rr(:,:) ! loc_points_in
  integer ,       intent(in)    :: ee(:,:)
  real(wp),       intent(in)    :: coupling_node_rot(3,3)

  !> Cell centers (to be evaluated), whose connectivity with the
  ! hinge nodes is created in this subroutine
  real(wp), allocatable :: loc_points(:,:)

  real(wp) :: hinge_width
  integer  :: n_nodes
  integer  :: nb, nh, ib, ih, iw, i
  real(wp), allocatable :: rrb(:,:), rrh(:,:), rrb_wei(:,:)
  real(wp) :: Rot(3,3) = 0.0_wp

  real(wp) :: span_wei, tol
  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::              ind_v(:)

  integer , allocatable :: rot_node_id(:) , ble_node_id(:)
  integer , allocatable :: rot_ind(:,:)   , ble_ind(:,:)
  real(wp), allocatable :: rot_wei(:,:)   , ble_wei(:,:)
  real(wp), allocatable :: rot_span_wei(:), ble_span_wei(:)

  integer , allocatable :: rot_p2h(:,:)   ! *** to do *** improve the actual inefficient
  real(wp), allocatable :: rot_w2h(:,:)   ! and memory intensive implementation
  integer , allocatable :: rot_i2h(:)     ! Look for ***** in this routine
  real(wp), allocatable :: rot_s2h(:,:)
  integer , allocatable :: ble_p2h(:,:)   ! *** to do *** improve the actual inefficient
  real(wp), allocatable :: ble_w2h(:,:)   ! and memory intensive implementation
  integer , allocatable :: ble_i2h(:)     ! Look for ***** in this routine
  real(wp), allocatable :: ble_s2h(:,:)

  integer :: nrot, nble
  real(wp) :: wei_hinge, x_hinge, z_hinge

  allocate(wei_v(this%n_wei))
  allocate(ind_v(this%n_wei))

  !> N. of surface elements (cell centers) and hinge nodes
  nb = size(ee,2)
  nh = this%n_nodes

  !> Cell centers, as the average of the element nodes
  allocate(loc_points(3,nb)); loc_points = 0.0_wp
  do ib = 1, nb

    n_nodes = 0

    do i = 1, size(ee,1)
      if ( ee(i,ib) .ne. 0 ) then
        n_nodes = n_nodes + 1
        loc_points(:,ib) = loc_points(:,ib) + rr(:,ee(i,ib))
      end if
    end do

    loc_points(:,ib) = loc_points(:,ib) / real(n_nodes,wp)

  end do
  
  !> Coordinates in the hinge reference frame
  ! Rotation matrix, build with the local ortonormal ref.frame of
  ! the first hinge node

  !for not coupled hinge
  Rot(1,:) = this%ref%v(:,1)
  Rot(2,:) = this%ref%h(:,1)
  Rot(3,:) = this%ref%n(:,1)

  allocate( rrb(3,nb) );  allocate( rrh(3,nh) )
  allocate( rrb_wei(3,nb) );

  !rotate the hinge nodes (from coupling_nodes file) to recover the value
  !of the component.in input
  ! NOTE: this%ref%rr are already rotated in build_connectivity

  do ib = 1, nb
    rrb(:,ib) =  matmul( (coupling_node_rot), loc_points(:,ib))
  end do

  ! hinge width, measured in the hinge direction (span)
  hinge_width = rrh(2,nh) - rrh(2,1)  !dh

  !> Compute connectivity and weights
  ! Allocate auxiliary node_id(:), ind(:,:), wei(:,:) arrays
  allocate(rot_node_id(        nb));  rot_node_id = 0
  allocate(ble_node_id(        nb));  ble_node_id = 0
  allocate(rot_ind(this%n_wei, nb));      rot_ind = 0
  allocate(ble_ind(this%n_wei, nb));      ble_ind = 0
  allocate(rot_wei(this%n_wei, nb));      rot_wei = 0.0_wp
  allocate(ble_wei(this%n_wei, nb));      ble_wei = 0.0_wp
  allocate(rot_span_wei(       nb)); rot_span_wei = 0.0_wp
  allocate(ble_span_wei(       nb)); ble_span_wei = 0.0_wp

  ! diff_all and dist_all auxiliaary arrays
  allocate(dist_all(   nh)); dist_all = 0.0_wp

  ! auxiliary matrix for hinge to surf connectivity
  ! ***** to do ***** improve the implementation
  allocate(rot_p2h(nh,nb)); rot_p2h = 0      ;  allocate(ble_p2h(nh,nb)); ble_p2h = 0
  allocate(rot_w2h(nh,nb)); rot_w2h = 0.0_wp ;  allocate(ble_w2h(nh,nb)); ble_w2h = 0.0_wp
  allocate(rot_i2h(nh   )); rot_i2h = 0      ;  allocate(ble_i2h(nh   )); ble_i2h = 0
  allocate(rot_s2h(nh,nb)); rot_s2h = 0.0_wp ;  allocate(ble_s2h(nh,nb)); ble_s2h = 0.0_wp

  nrot = 0; nble = 0

  do ib = 1, nb
    rrb_wei(:,ib) = matmul(coupling_node_rot, (loc_points(:,ib) - matmul(transpose(coupling_node_rot) ,this%ref%rr(:,1))))
  enddo

  tol = abs(this%ref%rr(2,nh) - this%ref%rr(2,1))/1000.0_wp  !tolerance for the hinge width (relative)

  ! Loop over all the surface points
  do ib = 1, nb

    if ((rrb(2,ib) .gt. this%ref%rr(2,1) - tol) .and. (rrb(2,ib) .lt. this%ref%rr(2,nh) + tol)) then
      
      wei_hinge = (rrb(2,ib) - this%ref%rr(2,1)) / (this%ref%rr(2,nh)- this%ref%rr(2,1))
      x_hinge = this%ref%rr(1,1) + wei_hinge*(this%ref%rr(1,nh)- this%ref%rr(1,1))
      z_hinge = this%ref%rr(3,1) + wei_hinge*(this%ref%rr(3,nh)- this%ref%rr(3,1))

      if (rrb(1,ib) .gt. (this%offset + x_hinge)) then

        nrot = nrot + 1
        rot_node_id(nrot) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb_wei(2,ib) - rrh(2,ih) )
        end do

        !> Weights in chordwise direction
        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

        wei_v = 1.0_wp / max( wei_v, 1e-9_wp )**this%w_order
        wei_v = wei_v / sum(wei_v)

        !> Weights in spanwise direction
        if ( rrb(2,ib) .lt. (this%ref%rr(2,1) + this%span_blending) ) then
          span_wei = 1.0_wp + rrb_wei(2,ib) / this%span_blending
        elseif( rrb(2,ib) .lt. (this%ref%rr(2,nh) - this%span_blending)  ) then
          span_wei = 1.0_wp
        else
          span_wei = 1.0_wp - ( rrb_wei(2,ib) - hinge_width ) / this%span_blending
        endif

        rot_wei(   :,nrot) = wei_v
        rot_ind(   :,nrot) = ind_v
        rot_span_wei(nrot) = span_wei

        do iw = 1, this%n_wei
          rot_i2h(ind_v(iw)) = rot_i2h(ind_v(iw)) + 1
          rot_p2h(ind_v(iw),   rot_i2h(ind_v(iw))) = ib
          rot_w2h(ind_v(iw),   rot_i2h(ind_v(iw))) = wei_v(iw)
          rot_s2h(ind_v(iw),   rot_i2h(ind_v(iw))) = span_wei
        end do

      elseif (rrb(1,ib) .lt. (this%offset + x_hinge) .and. rrb(1,ib) .gt. (x_hinge - this%offset))  then ! blending region

        nble = nble + 1
        ble_node_id(nble) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb_wei(2,ib) - rrh(2,ih) )
        end do

        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

        !> Weights in chordwise direction
        wei_v = 1.0_wp / max( wei_v, 1e-9_wp ) **this%w_order
        wei_v = wei_v / sum(wei_v)

        !> Weights in spanwise direction (is this needed?)
        if ( rrb(2,ib) .lt. (this%ref%rr(2,1) + this%span_blending) ) then
          span_wei = 1.0_wp + rrb_wei(2,ib) / this%span_blending
        elseif( rrb(2,ib) .lt. (this%ref%rr(2,nh) - this%span_blending)  ) then
          span_wei = 1.0_wp
        else
          span_wei = 1.0_wp - ( rrb_wei(2,ib) - hinge_width ) / this%span_blending
        endif

        ble_wei(   :,nble) = wei_v
        ble_ind(   :,nble) = ind_v
        ble_span_wei(nble) = span_wei

        ! *****
        do iw = 1, this%n_wei
          ble_i2h(ind_v(iw)) = ble_i2h(ind_v(iw)) + 1
          ble_p2h(ind_v(iw),   ble_i2h(ind_v(iw))) = ib
          ble_w2h(ind_v(iw),   ble_i2h(ind_v(iw))) = wei_v(iw)
          ble_s2h(ind_v(iw),   ble_i2h(ind_v(iw))) = span_wei
        end do
        ! *****
      else
      end if

    else
    end if

  end do

  !> Fill hinge object, with the connectivity and weight arrays
  allocate(this% rot_cen%node_id(        nrot)); this% rot_cen%node_id = rot_node_id( 1:nrot)
  allocate(this% rot_cen%ind(this%n_wei, nrot)); this% rot_cen%ind     = rot_ind(  :, 1:nrot)
  allocate(this% rot_cen%wei(this%n_wei, nrot)); this% rot_cen%wei     = rot_wei(  :, 1:nrot)
  allocate(this% rot_cen%span_wei(       nrot)); this% rot_cen%span_wei= rot_span_wei(1:nrot)
  allocate(this%blen_cen%node_id(        nble)); this%blen_cen%node_id = ble_node_id( 1:nble)
  allocate(this%blen_cen%ind(this%n_wei, nble)); this%blen_cen%ind     = ble_ind(  :, 1:nble)
  allocate(this%blen_cen%wei(this%n_wei, nble)); this%blen_cen%wei     = ble_wei(  :, 1:nble)
  allocate(this%blen_cen%span_wei(       nble)); this%blen_cen%span_wei= ble_span_wei(1:nble)
  allocate(this% rot_cen%n2h(nh))
  allocate(this%blen_cen%n2h(nh))
  do ih = 1, nh
    allocate(this% rot_cen%n2h(ih)%p2h ( rot_i2h(ih) )) ; &
            this% rot_cen%n2h(ih)%p2h = rot_p2h( ih, 1:rot_i2h(ih) )
    allocate(this% rot_cen%n2h(ih)%w2h ( rot_i2h(ih) )) ; &
            this% rot_cen%n2h(ih)%w2h = rot_w2h( ih, 1:rot_i2h(ih) )
    allocate(this% rot_cen%n2h(ih)%s2h ( rot_i2h(ih) )) ; &
            this% rot_cen%n2h(ih)%s2h = rot_s2h( ih, 1:rot_i2h(ih) )
    allocate(this%blen_cen%n2h(ih)%p2h ( ble_i2h(ih) )) ; &
            this%blen_cen%n2h(ih)%p2h = ble_p2h( ih, 1:ble_i2h(ih) )
    allocate(this%blen_cen%n2h(ih)%w2h ( ble_i2h(ih) )) ; &
            this%blen_cen%n2h(ih)%w2h = ble_w2h( ih, 1:ble_i2h(ih) )
    allocate(this%blen_cen%n2h(ih)%s2h ( ble_i2h(ih) )) ; &
            this%blen_cen%n2h(ih)%s2h = ble_s2h( ih, 1:ble_i2h(ih) )
  end do

  !> Explicit deallocations
  deallocate(rrb, rrh, dist_all, wei_v, ind_v)
  deallocate(rot_node_id, rot_ind, rot_wei, rot_span_wei)
  deallocate(ble_node_id, ble_ind, ble_wei, ble_span_wei)
  deallocate(rot_i2h, rot_p2h, rot_w2h, rot_s2h) 
  deallocate(ble_i2h, ble_p2h, ble_w2h, ble_s2h) 


end subroutine build_connectivity_cen

! ---------------------------------------------------------------
!> Build connectivity betweeen structural nodes and hinge nodes
! for evaluating the reference frames attached to the non-rotating
! structure and the deflection of the rotating surfaces, theta
! (needed for the blending regions)
subroutine build_connectivity_hin(this, rr_t, ind_h )
  class(t_hinge), intent(inout) :: this
  real(wp),       intent(in)    ::  rr_t(:,:) ! precice nodes of the component
  integer ,       intent(in)    :: ind_h(:)   ! local id of precice hinge nodes

  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::    ind_b(:), ind_v(:)
  real(wp), allocatable :: rr_h(:,:), rr_b(:,:)
  integer :: n_b, n_h, n_t
  integer :: i_b, i_h, i_t

  !integer, allocatable :: ind(:)
  integer :: ind

  allocate(wei_v(this%n_wei))
  allocate(ind_v(this%n_wei))

  !> Find hinge and structural nodes in rr_t array, collecting all the nodes
  n_t = size(rr_t,2)
  n_h = size(ind_h);  allocate(rr_h(3,n_h)) ;  rr_h = 0.0_wp
  n_b = n_t - n_h  ;  allocate(rr_b(3,n_b)) ;  rr_b = 0.0_wp
  allocate(ind_b(n_b)) ;  ind_b = -333
  i_h = 0; i_b = 0


  do i_t = 1, n_t
    if ( any( ind_h .eq. i_t ) ) then
      i_h = i_h + 1
      !ind = ind_h( findloc(ind_h, value=i_t) )
      !workaround for older compilers
      do ind = 1,size(ind_h)
        if (ind_h(ind) .eq. i_t) exit
      enddo

      !rr_h(:,i_h) = rr_t(:, ind(1) )
      rr_h(:,i_h) = rr_t(:, i_t)
    else
      i_b = i_b + 1
      rr_b(:,i_b) = rr_t(:,i_t)
      ind_b(i_b) = i_b
    end if
  end do


  !> Allocate and fill hinge%hin object
  allocate(this%hin%node_id(       n_h)); this%hin%node_id = -333 ! useless
  allocate(this%hin%ind(this%n_wei,n_h))
  allocate(this%hin%wei(this%n_wei,n_h))
  allocate(this%hin%span_wei(0)) ! useless
  allocate(this%hin%n2h(0))      ! useless
  ! diff_all and dist_all auxiliaary arrays
  allocate(dist_all(n_b)); dist_all = 0.0_wp

  ! Loop over all the surface points
  do i_h = 1, n_h

    do i_b = 1, n_b
      dist_all(i_b) = norm2( rr_b(:,i_b) - rr_h(:,i_h) )
    end do

    !> Weights in chordwise direction
    call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )
    wei_v = 1.0_wp / max( wei_v, 1e-9_wp ) **this%w_order
    wei_v = wei_v / sum(wei_v)

    ! this%hin%node_id(i_h) = i_h
    this%hin%ind(:,i_h) = ind_b(ind_v)
    this%hin%wei(:,i_h) = wei_v

  end do

  !> Allocate t_hinge%hin_rot
  allocate(this%hin_rot(3,n_h)); this%hin_rot = -333.3_wp

end subroutine build_connectivity_hin

! ---------------------------------------------------------------
!> Compute actual configuration of the hinge nodes
! *** to do *** unpredictable (or better, wrong) behavior when
! restarting a simulation
subroutine init_theta(this, t)
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    :: t

  if ( t .ne. 0.0_wp ) then
    write(*,*) ' Error in t_hinge % init_theta: t .ne. 0.0_wp. &
                &This argument is meant for future restart capabilities. &
                &So far, must be passed to the function equal to 0.0_wp. Stop'
    stop
  end if

  call this % update_theta( t )
  
  if ( trim(this%input_type) .eq. 'function:const' ) then
    this%theta_old = this%theta
  else
    this%theta_old = 0.0_wp
  end if

end subroutine init_theta

! ---------------------------------------------------------------
!> Compute actual configuration of the hinge nodes
subroutine from_reference_to_actual_config(this)
  class(t_hinge), intent(inout) :: this

  ! *** to do ***
  ! still usefull? anything else to do, that is missing in
  ! update_hinge_nodes() subroutine?


end subroutine from_reference_to_actual_config

! ---------------------------------------------------------------
!> Update hinge nodes, for non-coupled components
subroutine update_hinge_nodes( this, R, of )
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    ::  R(:,:)
  real(wp)      , intent(in)    :: of(:)

  if ( trim(this%input_type) .ne. 'coupling' ) then !> Prescribed motion
    !> Actual configuration: node position
    this % act % rr = matmul( R, this % ref % rr )
    this % act % rr(1,:) = this % act % rr(1,:) + of(1)
    this % act % rr(2,:) = this % act % rr(2,:) + of(2)
    this % act % rr(3,:) = this % act % rr(3,:) + of(3)

    !> Actual configuration: node orientation
    this % act % h = matmul( R, this % ref % h )
    this % act % v = matmul( R, this % ref % v )
    this % act % n = matmul( R, this % ref % n )

  else !> Coupling with dynamics solver

    ! see precice/mod_precice.f90/t_precice % update_elems(), l.1100 approx

  end if


end subroutine update_hinge_nodes

! ---------------------------------------------------------------
!> Update hinge nodes, for non-coupled components
subroutine update_theta( this, t )
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    :: t

  !> Update theta_old
  if (t .gt. 0) then 
    this % theta_old = this % theta
  end if

  if ( trim(this%input_type) .eq. 'function:const' ) then
    this%theta = this%f_ampl

  elseif ( trim(this%input_type) .eq. 'function:sin' ) then
    if ( t .le. this%f_initial_time ) then
      this%theta = this%f_ampl_init
    else
      this%theta = this%f_ampl_init + &
                   this%f_ampl * sin( this%f_omega * (t - this%f_initial_time) - this%f_phase )
    end if

  elseif ( trim(this%input_type) .eq. 'function:cos' ) then
    if ( t .le. this%f_initial_time ) then
      this%theta = this%f_ampl_init
    else
      this%theta = this%f_ampl_init + &
                   this%f_ampl * cos( this%f_omega * (t - this%f_initial_time) - this%f_phase )
    end if

  elseif ( trim(this%input_type) .eq. 'function:cosine_driver' ) then
    !> before the start of the driver
    if ( t .le. this%f_initial_time ) then
      this%theta = this%f_ampl_init
    !> driver
    elseif ((t .gt. this%f_initial_time) .and. (t .le. this%f_initial_time + 2.0_wp*pi/this%f_omega*this%f_cosine_cycl)) then
        this%theta = this%f_ampl_init + &
                    this%f_ampl * (1.0_wp - cos( this%f_omega * (t - this%f_initial_time)))
    !> final value
    elseif ( t .gt. this%f_initial_time +  2.0_wp*pi/this%f_omega*this%f_cosine_cycl) then
      !> check if this%f_cosine_cycl is a multiple of 0.5 (half cycle(s))
      if (mod(this%f_cosine_cycl, 0.5) .eq. 0.0 .and. modulo(this%f_cosine_cycl, 1.0) .ne. 0.0) then
        this%theta = 2.0_wp * this%f_ampl + this%f_ampl_init
      !> check if this%f_cosine_cycl is a multiple of 1 (full cycle(s))
      elseif (mod(this%f_cosine_cycl, 1.0_wp) .eq. 0.0  ) then
        this%theta = this%f_ampl_init
      !> otherwise, it takes the last calculated value
      else
        this%theta = this%theta_old
      endif
    endif

  !elseif ( trim(this%input_type) .eq. 'from_file' ) then
    ! TODO
    
  else
    this%theta_old = 0.0_wp
    this%theta     = 0.0_wp

  end if

end subroutine update_theta

! ---------------------------------------------------------------
!> Update the coordinates rr of the points of the surface, after
! hinge deflection
subroutine hinge_deflection(i_points, this,  rr, t, te_i, te_t, postpro )
  class(t_hinge),     intent(inout)   :: this
  real(wp), optional, intent(inout)   :: te_t(:,:)
  integer,  optional, intent(in)      :: te_i(:,:)
  real(wp),           intent(inout)   :: rr(:,:)
  real(wp),           intent(in)      :: t
  logical,  optional, intent(in)      :: postpro
  integer,            intent(in)      :: i_points(:)
  
  logical                             :: local_postpro
  logical, parameter                  :: default_postpro = .false.
  real(wp), allocatable               :: rr_in(:,:)
  real(wp)                            :: th, th1, thp, yc, xq, yq, xqp, yqp
  real(wp)                            :: nx(3,3), Rot_I(3,3), eye(3,3) 
  integer                             :: nrot, nble, ib, ih, ii, it

  eye(1,:) = (/1.0_wp, 0.0_wp, 0.0_wp/)
  eye(2,:) = (/0.0_wp, 1.0_wp, 0.0_wp/)
  eye(3,:) = (/0.0_wp, 0.0_wp, 1.0_wp/)

  if ( trim(this%input_type) .ne. 'coupling' ) then
    ! Old routine for hinges with prescribed motion

    !> Use the same routine for solver (incremental rotation) and
    ! postpro (non-incremental rotation)
    if ( present(postpro) ) then;  local_postpro = postpro
    else                        ;  local_postpro = default_postpro
    end if

    allocate(rr_in(size(rr,1),size(rr,2)))
    rr_in = rr

    !> n.nodes in the rigid-rotation and in the blending regions
    nrot = size(this%rot %node_id)
    nble = size(this%blen%node_id)

    do ih = 1, this%n_nodes

      th =   this % theta(ih) * pi/180.0_wp

      if ( th .ne. 0.0_wp ) then ! (equality check on real?)

        ! Rotation matrix
        nx(1,:) = (/            0.0_wp, -this%act%h(3,ih),  this%act%h(2,ih) /)
        nx(2,:) = (/  this%act%h(3,ih),            0.0_wp, -this%act%h(1,ih) /)
        nx(3,:) = (/ -this%act%h(2,ih),  this%act%h(1,ih),            0.0_wp /)

        !> Rigid rotation
        do ib = 1, size(this%rot%n2h(ih)%p2h)
          ii = this%rot%n2h(ih)%p2h(ib)
          th1 = th * this%rot%n2h(ih)%s2h(ib)
          
          Rot_I = sin(th1) * nx + ( 1.0_wp - cos(th1) ) * matmul( nx, nx )       
          rr(:,ii) = rr(:,ii) + &
                      this%rot%n2h(ih)%w2h(ib) * &
                      matmul( Rot_I, rr_in(:,ii)-this%act%rr(:,ih) )

        end do
        
        !> Blending region
        do ib = 1, size(this%blen%n2h(ih)%p2h)

          ii = this%blen%n2h(ih)%p2h(ib)
          th1 = -th * this%blen%n2h(ih)%s2h(ib)
          !> coordinate of the centre of the circle used for blending,
          ! in the n-direction
          yc = cos(th1)/sin(th1) * this%offset * ( 1.0_wp + cos(th1) ) + &
                                   this%offset * sin(th1)
          !> Some auxiliary quantities
          xq = sum( ( rr_in(:,ii) - this%act%rr(:,ih) ) * this%act%v(:,ih) )
          yq = sum( ( rr_in(:,ii) - this%act%rr(:,ih) ) * this%act%n(:,ih) )
          thp = 0.5_wp * ( xq + this%offset ) / this%offset * th1
          xqp = yc*sin(thp)          - this%offset - yq*sin(thp) - xq
          yqp = yc*(1.0_wp-cos(thp))               + yq*cos(thp) - yq

          !> Update coordinates
          rr(:,ii) = rr(:,ii) + &
                     this%blen%n2h(ih)%w2h(ib) * &
                   ( xqp * this%act%v(:,ih) + yqp * this%act%n(:,ih) )

        end do
      
      end if

    end do

    if (present(te_t)) then     
      ! Rotate trailing edge direction in the rigid-rotation region
      do it = 1, size(te_t,2)
        
        do ih = 1, this%n_nodes
          th =   (this%theta(ih) - this%theta_old(ih)) * pi/180.0_wp  !
          ! Rotation matrix
            nx(1,:) = (/            0.0_wp, -this%act%h(3,ih),  this%act%h(2,ih) /)
            nx(2,:) = (/  this%act%h(3,ih),            0.0_wp, -this%act%h(1,ih) /)
            nx(3,:) = (/ -this%act%h(2,ih),  this%act%h(1,ih),            0.0_wp /)

            do ib = 1, size(this%rot%n2h(ih)%p2h)


              ii = this%rot%n2h(ih)%p2h(ib)
              ii = i_points(ii)
              th1 = th * this%rot%n2h(ih)%s2h(ib)
              
              Rot_I = eye + sin(th1) * nx + ( 1.0_wp - cos(th1) ) * matmul( nx, nx )
              
              if (te_i(1 , it) .eq. ii) then ! hinge node is also trailing edge node

                te_t(:,it) = this%rot%n2h(ih)%s2h(ib) * matmul( Rot_I, te_t(:,it))
              
              end if
            
            end do
            
        end do
      end do

    end if  

    deallocate(rr_in)

  else  !> Coupled hinge

    ! see precice/mod_precice.f90/t_precice % update_elems()

  end if

end subroutine hinge_deflection

! ---------------------------------------------------------------
!> Build hinge configuration
! Used to build reference configuration in load_component(), and
! initialize actual configuration, as the reference configuration
subroutine initialize_hinge_config( h_config, hinge )
  type(t_hinge_config), intent(inout) :: h_config
  type(t_hinge)       , intent(inout) :: hinge

  real(wp) :: hv(3), vv(3), nv(3)
  integer :: i

  !> rr0, rr1 (useless?)
  h_config%rr0 = hinge%ref%rr(:,1)
  h_config%rr1 = hinge%ref%rr(:,hinge%n_nodes)

  hv = h_config%rr1 - h_config%rr0
  hv = hv / norm2(hv)
  vv = hinge%ref_dir / norm2(hinge%ref_dir)

  nv = cross( vv, hv )
  
  if (norm2(nv) .le. 1e-16_wp) then ! workaround for debug compiling 
    nv(3) = 1e-16_wp
  endif

  nv = nv / norm2(nv)
  vv = cross( hv, nv )
  vv = vv / norm2(vv)  ! useless normalization

  !> h, v, n
  allocate( h_config%h(3, hinge%n_nodes) )
  allocate( h_config%v(3, hinge%n_nodes) )
  allocate( h_config%n(3, hinge%n_nodes) )

  do i = 1, hinge%n_nodes

    h_config%h(:,i) = hv;  h_config%v(:,i) = vv;  h_config%n(:,i) = nv

  end do


end subroutine initialize_hinge_config


! ---------------------------------------------------------------
!> build hinges_input, used in dust_pre for generating geometry input file
! for the solver
subroutine build_hinges( geo_prs, n_hinges, hinges )
  type(t_parse), intent(inout) :: geo_prs
  integer      , intent(in)    :: n_hinges
  type(t_hinge_input), allocatable, intent(inout) :: hinges(:)

  type(t_parse), pointer :: hinge_prs, fun_prs, file_prs, coupling_prs
  character(len=max_char_len) :: hinge_node_subset
  integer :: i, j, id_1, id_2

  if ( allocated(hinges) ) deallocate(hinges)
  allocate( hinges(n_hinges) )

  do i = 1, n_hinges

    !> De-associate, before reading next hinge
    hinge_prs => null()
    fun_prs => null(); file_prs => null(); coupling_prs => null()
    !> Open hinge sub-parser
    call getsuboption(geo_prs, 'Hinge', hinge_prs)

    hinges(i)%tag = getstr(hinge_prs, 'hinge_tag')
    hinges(i)%nodes_input = getstr(hinge_prs, 'hinge_nodes_input')

    if ( trim(hinges(i) % nodes_input) .eq. 'parametric' ) then
      hinges(i)%node1 = getrealarray(hinge_prs, 'node_1', 3)
      hinges(i)%node2 = getrealarray(hinge_prs, 'node_2', 3)
      hinges(i)%le1 = 0.0_wp
      hinges(i)%le2 = 0.0_wp
      hinges(i)%te1 = 0.0_wp
      hinges(i)%te2 = 0.0_wp
      hinges(i)%chord1 = 0.0_wp
      hinges(i)%chord2 = 0.0_wp
      hinges(i)%csi1 = 0.0_wp
      hinges(i)%csi2 = 0.0_wp

      hinges(i)%n_nodes = getint(hinge_prs, 'n_nodes')
      hinges(i)%node_file = 'hinge with parametric input. If you read this &
                            &string, something probabily went wrong.'
    allocate( hinges(i)%rr( 3, hinges(i)%n_nodes ) )
      do j = 1, hinges(i)%n_nodes
          hinges(i)%rr(:,j) = hinges(i)%node1 + &
                             (hinges(i)%node2 - hinges(i)%node1 ) * &
                                dble(j-1)/dble(hinges(i)%n_nodes-1)
      end do

    elseif ( trim(hinges(i) % nodes_input) .eq. 'from_file' ) then
      ! *** to do *** fill dummy node1, node2, n_nodes fields
      hinges(i) % node_file = getstr(hinge_prs,'Node_File')
      call read_hinge_nodes( hinges(i)%node_file, &
                              hinges(i)%n_nodes  , &
                              hinges(i)%rr )
    else
      ! *** to do *** use error message handling
      write(*,*) ' Wrong Hinge_Nodes_Input input. Stop '; stop
    end if

    !> Hinge reference direction (zero rotation direction) and offset to
    ! avoid irregular behavior during hinge rotation
    hinges(i)%ref_dir = getrealarray(hinge_prs, 'hinge_ref_dir', 3)
    hinges(i)%offset  = getreal(hinge_prs, 'hinge_offset')
    hinges(i)%span_blending = getreal(hinge_prs, 'hinge_spanwise_blending')
    hinges(i)%merge_tol = getreal(hinge_prs,'hinge_merge_tol')
    hinges(i)%adaptive_mesh = getlogical(hinge_prs,'hinge_adaptive_mesh')
    !> Hinge input: function, amplitude
    !> Overwrite 'constant' rotation_input with 'function:const'
    hinges(i) % rotation_input = getstr(hinge_prs,'hinge_rotation_input')
    !> Overwrite old 'constant' input
    if ( trim( hinges(i)%rotation_input ) .eq. 'constant' ) then
      hinges(i)%rotation_input = 'function:const'
    end if
    if ( ( trim(hinges(i)%rotation_input) .ne. 'function:const' ) .and. &
          ( trim(hinges(i)%rotation_input) .ne. 'function:sin'   ) .and. &
          ( trim(hinges(i)%rotation_input) .ne. 'function:cos'   ) .and. &
          ( trim(hinges(i)%rotation_input) .ne. 'function:cosine_driver') .and. &
          ( trim(hinges(i)%rotation_input) .ne. 'from_file'      ) .and. &
          ( trim(hinges(i)%rotation_input) .ne. 'coupling'       ) ) then
            write(*,*) ' Error in t_hinge%build_hinge(): rotation_input = &
          & $ '//trim(hinges(i)%rotation_input)// &
                    &'not known.'; stop
    else
      if ( ( trim(hinges(i)%rotation_input) .eq. 'function:const' ) .or. &
            ( trim(hinges(i)%rotation_input) .eq. 'function:sin'   ) .or. &
            ( trim(hinges(i)%rotation_input) .eq. 'function:cos'   ) .or. &
            ( trim(hinges(i)%rotation_input) .eq. 'function:cosine_driver')) then

        call getsuboption(hinge_prs, 'Hinge_Rotation_Function', fun_prs)
        hinges(i) % rotation_amplitude = getreal(fun_prs,'amplitude')
        hinges(i) % rotation_omega     = getreal(fun_prs,'omega')
        hinges(i) % rotation_phase     = getreal(fun_prs,'phase')
        hinges(i) % rotation_amplitude_init     = getreal(fun_prs,'amplitude_initial')
        hinges(i) % rotation_cosine_cycl     = getreal(fun_prs,'number_of_cycles')
        hinges(i) % rotation_initial_time     = getreal(fun_prs,'initial_time')
      !elseif ( trim(hinges(i)%rotation_input) .eq. 'from_file' ) then
      !  write(*,*) ' Error in t_hinge%build_hinge(): rotation_input = from_file &
      !              &not implemented yet.'; stop

      elseif ( trim(hinges(i)%rotation_input) .eq. 'coupling' ) then

        !> Read coupling options
        call getsuboption(hinge_prs, 'Hinge_Rotation_Coupling', coupling_prs)
        hinge_node_subset = getstr(coupling_prs,'coupling_node_subset')

        if ( trim(hinge_node_subset) .eq. 'range' ) then
          id_1 = getint(coupling_prs,'coupling_node_first')
          id_2 = getint(coupling_prs,'coupling_node_last' )
           ! *** to do *** add some checks on node numbering ???

          allocate( hinges(i) % coupling_nodes( id_2-id_1+1 ) )
          do j = id_1, id_2; hinges(i)%coupling_nodes(j-id_1+1) = j; end do

        elseif ( trim(hinge_node_subset) .eq. 'from_file' ) then
           ! *** to do ***
          write(*,*) ' Coupling_Node_Subset = from_file, not implemented yet. Stop'
          stop
        else
          write(*,*) ' Coupling_Node_Subset must be = "range" or "from_file", but it &
                      &is = '// trim(hinge_node_subset)//'. Stop.'; stop
        end if

        !> Set "default" values of the function: inputs
        hinges(i) % rotation_amplitude = 1.0_wp
        hinges(i) % rotation_omega     = 0.0_wp
        hinges(i) % rotation_phase     = 0.0_wp
        hinges(i) % rotation_amplitude_init  = 0.0_wp
        hinges(i) % rotation_cosine_cycl     = 0.5_wp
        hinges(i) % rotation_initial_time     = 0.0_wp

      end if
    end if

  end do


end subroutine build_hinges

! ---------------------------------------------------------------
!> Read ascii file of hinge node coordinates, w/ input type from_file
!  Three-column ascii file is expected
subroutine read_hinge_nodes( filen, n_nodes, rr )
  character(len=max_char_len), intent(in)  :: filen
  integer                    , intent(out) :: n_nodes
  real(wp), allocatable      , intent(out) :: rr(:,:)

  integer :: i, io
  integer :: fid = 21

  if ( allocated(rr) ) deallocate(rr)

  n_nodes = 0;  io = 0

  !> Preliminary read for determining the n. of lines
  open(unit=fid, file=trim(filen))
  do while ( io .eq. 0 )
    read(fid,*,iostat=io) ! dummy
    n_nodes = n_nodes + 1
  end do
  close(fid)

  n_nodes = n_nodes - 1 ! *** to do *** check

  !> Allocate and fill rr array
  allocate(rr(3,n_nodes))
  open(unit=fid, file=trim(filen))
  do i = 1, n_nodes
    read(fid,*) rr(:,i)
  end do
  close(fid)


end subroutine read_hinge_nodes

! ---------------------------------------------------------------
!> Hinge input parser, called in mod_build_geo.f90 by dust_pre preprocessor
subroutine hinge_input_parser( geo_prs, hinge_prs, &
                              fun_prs, file_prs, coupling_prs )
  type(t_parse),          intent(inout) :: geo_prs
  type(t_parse), pointer, intent(inout) :: hinge_prs
  type(t_parse), pointer, intent(inout) :: fun_prs, file_prs, coupling_prs

  call geo_prs%CreateIntOption('n_hinges', &
              'N. of hinges and rotating parts (e.g. aileron) of the component', &
              '0') ! default: no hinges -> n_hinges = 0

  call geo_prs%CreateSubOption('hinge', 'Parser for hinge input', &
                                hinge_prs, multiple=.true. )

  call hinge_prs%CreateStringOption('hinge_tag','Name of the hinge')
  call hinge_prs%CreateStringOption('hinge_nodes_input', &
                              'Type of hinge nodes input: parametric or from_file.')
  call hinge_prs%CreateIntOption('n_nodes','N.hinge nodes')
  call hinge_prs%CreateRealArrayOption('node_1', &
      'First node of the hinge. Components in the local ref.frame of the component')
  call hinge_prs%CreateRealArrayOption('node_2', &
      'Second node of the hinge. Components in the local ref.frame of the component')
  call hinge_prs%CreateRealArrayOption('hinge_ref_dir', &
      'Reference direction of the hinges, indicating zero-deflection direction in &
      &the local ref.frame of the component')
  call hinge_prs%CreateRealOption('hinge_offset','Offset in the Ref_Dir needed for &
      &avoiding irregular behavior of the surface for large deflections')
  call hinge_prs%CreateRealOption('hinge_spanwise_blending', &
      'Blending in the spanwise direction needed for &
      &avoiding irregular behavior of the surface for large deflections','0.0')
  call hinge_prs%CreateRealOption('hinge_merge_tol', & 
      'Tolerance for adaptive hinge mesh in chord percentage', '0.01')
  call hinge_prs%CreateLogicalOption('hinge_adaptive_mesh', & 
      'Enable for adaptive hinge mesh', 'F')
  !> Hinge_Rotation_Input
  call hinge_prs%CreateStringOption('hinge_rotation_input', &
      'Input type of the rotation: function, from_file, coupling')
  !> Hinge_Rotation_Input = function:...
  call hinge_prs%CreateSubOption('hinge_rotation_function', &
                'Parser for hinge input w/ simple functions', fun_prs )
  call fun_prs%CreateRealOption('amplitude', &
      'Amplitude of the rotation, for constant, function:const, :sin, &
      &:cos Rotation_Input')
  call fun_prs%CreateRealOption('omega', &
      'Angular velocity of the rotation, for constant, function:const, &
      &:sin, :cos Rotation_Input', '0.0')
  call fun_prs%CreateRealOption('phase', &
      'Phase of the rotation, for constant, function:const, :sin, :cos &
      &Rotation_Input', '0.0')
  call fun_prs%CreateRealOption('amplitude_initial', &
      'Initial Amplitude of the rotation, for function:sin, &
      &:cos, cosine_driver Rotation_Input', '0.0')
  call fun_prs%CreateRealOption('number_of_cycles', &
      'Number of cycles for function:cosine_driver Rotation_Input', '0.5')
  call fun_prs%CreateRealOption('initial_time', &
      'Rotation Initial Time for function::sin, :cos &
      &Rotation_Input', '0.0')      

  !> Hinge_Rotation_Input = from_file
  call hinge_prs%CreateSubOption('hinge_rotation_file', &
              'Parser for hinge input from file', file_prs )
  call file_prs%CreateStringOption('file_name', &
      'Name of the file containing the input of the hinge rotation')
  !> Hinge_Rotation_Input = coupling
  call hinge_prs%CreateSubOption('hinge_rotation_coupling', &
              'Parser for hinge input from coupling', coupling_prs )
  call coupling_prs%CreateStringOption('coupling_node_subset', &
      'Optional. Define a subset of structural nodes to evaluate &
      &coupling: "range" or "from_file"')
  call coupling_prs%CreateIntOption('coupling_node_first','If node subset &
      &is defined through "range" input: first id of the nodes')
  call coupling_prs%CreateIntOption('coupling_node_last','If node subset &
      &is defined through "range" input: last id of the nodes')
  call coupling_prs%CreateStringOption('coupling_node_filename', &
      'File collecting the IDs of the coupling nodes for hinge coupling')

  ! *** to do ***
  ! add all the fields required for all the input types


end subroutine hinge_input_parser

! ---------------------------------------------------------------

end module mod_hinges
