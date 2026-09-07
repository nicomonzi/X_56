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

!> Module containing the upper level definitions of the aerodynamic
!! elements


module mod_aeroel

use mod_linsys_vars, only: t_linsys

use mod_param, only: &
  wp

implicit none

public :: c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
          t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

private

!----------------------------------------------------------------------

!> Pointers to different levels of classes
type :: t_elem_p
  class(c_elem), pointer :: p
end type

type :: t_pot_elem_p
  class(c_pot_elem), pointer :: p
end type

type :: t_vort_elem_p
  class(c_vort_elem), pointer :: p
end type

type :: t_impl_elem_p
  class(c_impl_elem), pointer :: p
end type

type :: t_expl_elem_p
  class(c_expl_elem), pointer :: p
end type

!----------------------------------------------------------------------

!> Abstract type defining a generic aerodynamic element
!!
!! It needs to be extended into more specific elements
type, abstract :: c_elem

  !> Magnitude of the singularity, pointer to the actual location
  real(wp), pointer :: mag => null()

  !> Center of the element
  real(wp) :: cen(3)

  contains

  procedure(i_compute_vel) , deferred, pass(this) :: compute_vel
  procedure(i_compute_grad), deferred, pass(this) :: compute_grad

end type c_elem

!----------------------------------------------------------------------

!> Class (not abstract, but to be further expanded) to contain potential
!! elements
type, abstract, extends(c_elem) :: c_pot_elem

  !> id of the component to which it belongs
  integer :: comp_id

  !> Id of the element, TODO consider removing this
  integer :: id

  !> Number of vertexes
  integer :: n_ver

  !> Vertexes coordinates
  real(wp), allocatable :: ver(:,:)

  !> Id of the vertexes
  integer, allocatable :: i_ver(:)

  !> Element area
  real(wp) :: area

  !> Element normal vector
  real(wp) :: nor(3)

  !> Element normal vector (at previous dt, to compute dn_dt)
  real(wp) :: nor_old(3)

  !> Element normal vector (at previous dt, to compute dn_dt)
  real(wp) :: dn_dt(3)

  !> Element tangential vectors
  real(wp) :: tang(3,2) ! tangent unit vectors as in PANAIR

  !> Vector of each edge
  real(wp), allocatable :: edge_vec(:,:)

  !> Length of each edge
  real(wp), allocatable :: edge_len(:)

  !> Unit vector of each edge
  real(wp), allocatable :: edge_uni(:,:)

  !> Body velocity at the centre
  real(wp) :: ub(3)

  !> Orientation vector of the center 
  real(wp) :: ori(3)

  !> Orientation matrix of the center 
  real(wp) :: R_cen(3,3)

  !> Vorticity induced velocity at the centre
  real(wp) :: uvort(3)

  !> Is the element moving during simulation?
  logical  :: moving

  !> Element neighbours global index (in the elements vector)
  type(t_pot_elem_p), allocatable :: neigh(:)

  !> Average pressure on the element
  real(wp)  :: pres

  !> Elementary force acting on the element (components in the base ref.sys.)
  real(wp)  :: dforce(3)

  !> Elementary force acting on the element (components in the base ref.sys.)
  real(wp)  :: dmom(3)
  integer   :: n_c, n_s
  real(wp)  :: dy
  !TODO: these three are used only by vortlatt and liftlin
  ! consider moving them there, but then change the implementation of
  ! create_strip_connectivity
  !> Previous element in a stripe
  type(t_pot_elem_p)  :: stripe_1
  
  !> Element indices in the component%strip_elem array
  type(t_elem_p), allocatable :: stripe_elem(:)
  
  !> Hinge motion
  ! Initialization to zero *** to do: restart??? ***
  !> Delta position, due to hinge motion of the element center
  real(wp) :: dcen_h(3)     ! = 0.0_wp
  real(wp) :: dcen_h_old(3) ! = 0.0_wp
  !> Delta velocity, due to hinge motion of the element center, evaluated
  ! with finite difference method: dvel_h = ( dcen_h - dcen_h_old ) / dt
  real(wp) :: dvel_h(3) ! = 0.0_wp
  real(wp) :: loc_ctr_pt(3)
  contains

  procedure(i_compute_pot)   , deferred, pass(this) :: compute_pot
  procedure(i_compute_linear_pot)   , deferred, pass(this) :: compute_linear_pot
  procedure(i_compute_psi)   , deferred, pass(this) :: compute_psi
  procedure(i_compute_pres)  , deferred, pass(this) :: compute_pres
  procedure(i_compute_dforce), deferred, pass(this) :: compute_dforce
  procedure(i_calc_geo_data) , deferred, pass(this) :: calc_geo_data

end type c_pot_elem

!----------------------------------------------------------------------

!> Class to contain vortical elements
type, abstract, extends(c_elem) :: c_vort_elem

end type c_vort_elem

!----------------------------------------------------------------------

!> Class to contain the implicit elements (Vortex lattices and surface panels)
type, abstract, extends(c_pot_elem) :: c_impl_elem

  real(wp)          :: didou_dt

  contains

  procedure(i_build_row)           , deferred, pass(this) :: build_row
  procedure(i_build_row_static)    , deferred, pass(this) :: build_row_static
  procedure(i_add_wake)            , deferred, pass(this) :: add_wake
  procedure(i_add_expl)            , deferred, pass(this) :: add_expl
  procedure(i_get_vort_vel)        , deferred, pass(this) :: get_vort_vel
  procedure(i_get_bernoulli_source), deferred, pass(this) :: get_bernoulli_source

end type c_impl_elem

!----------------------------------------------------------------------

!> Class to contain the explicit elements
type, abstract, extends(c_pot_elem) :: c_expl_elem

end type c_expl_elem

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!> Build a complete row of the linear system, including the right hand
!! side
!!
!! Since the static part of the sistem should already be initialized,
!! the function updates just the coefficients from ista to iend,
!! which should be the moving elements ordered at the end of the
!! elements array
!!
!! The static part od the rhs should already be initialized, and the
!! moving contribution is just added to the static one.
abstract interface
  subroutine i_build_row(this, elems, linsys, ie, ista, iend)
    import                              :: wp
    import                              :: c_impl_elem
    import                              :: t_impl_elem_p
    import                              :: t_linsys
    implicit none
    class(c_impl_elem), intent(inout)   :: this
    type(t_impl_elem_p), intent(in)     :: elems(:)
    type(t_linsys), intent(inout)       :: linsys
    integer, intent(in)                 :: ie
    integer, intent(in)                 :: ista, iend
  end subroutine
end interface

!----------------------------------------------------------------------

!> Build the static part of a row of the linear system
!!
!! It is necessary to calculate the aerodynamic influence coefficients as
!! well as the rhs contribution of the static elements just once at the
!! beginning of the simulation.
!!
!! For this reason this subroutine is called to calculate the static aic
!! and the static contribution to the rhs. The static contribution to the
!! rhs is stored into a 3 components * number of elements array so that
!! at each timestep the actual rhs contribution can be calculated
!! multiplying the 3 components for the velocity vector, which in principle
!! could vary.
!!
!! The function operates only on the components from ista to iend, which
!! should be the static ones ordered at the beginning of the array
abstract interface
  subroutine i_build_row_static(this, elems, expl_elems, linsys, ie, ista, iend)
    import                              :: wp
    import                              :: c_impl_elem
    import                              :: t_impl_elem_p
    import                              :: t_expl_elem_p
    import                              :: t_linsys
    implicit none
    class(c_impl_elem), intent(inout)   :: this
    type(t_impl_elem_p), intent(in)     :: elems(:)
    type(t_expl_elem_p), intent(in)     :: expl_elems(:)
    type(t_linsys), intent(inout)       :: linsys
    integer, intent(in)                 :: ie
    integer, intent(in)                 :: ista, iend
  end subroutine
end interface

!----------------------------------------------------------------------

!> Add the contribution of the wake to the linear system.
!!
!! The contribution is divided into an implicit part, due to the first row
!! of wake panels, and an explicit contribution due to all the other wake
!! elements.
!! The implicit part affects the linear system matrix A while the explicit
!! contributions affect only the right hand side b.
!!
!! The vector impl_wake_ind tells which are the elements affected by the
!! implicit wake panels.
!! The parameters ista and iend are employed to separate the update of the
!! stationary part of the matrix from the dynamic part.
!! The stationary part is updated once at the beginning of the simulation
!! while the moving one is updated each timestep
abstract interface
  subroutine i_add_wake(this, wake_elems, impl_wake_ind, linsys, ie, ista, iend)
    import                              :: wp
    import                              :: c_impl_elem
    import                              :: t_pot_elem_p
    import                              :: t_linsys
    implicit none
    class(c_impl_elem), intent(inout)   :: this
    type(t_pot_elem_p), intent(in)      :: wake_elems(:)
    integer, intent(in)                 :: impl_wake_ind(:,:)
    type(t_linsys), intent(inout)       :: linsys
    integer, intent(in)                 :: ie
    integer, intent(in)                 :: ista
    integer, intent(in)                 :: iend
  end subroutine
end interface

!----------------------------------------------------------------------


!> Add the contribution of explicit elements to the linear system.
!!
!! Since the actuator disks are completely explicit their contribution does
!! not affect the system matrix, but only the system right han side.
!! The contribution affecting the stationary part was already calculated
!! and just retrieved, while the contribution due to the moving components
!! is calculated and added
abstract interface
  subroutine i_add_expl(this, expl_elems, linsys, ie, ista, iend)
    import                              :: wp
    import                              :: c_impl_elem
    import                              :: t_expl_elem_p
    import                              :: t_linsys
    implicit none
    class(c_impl_elem), intent(inout)   :: this
    type(t_expl_elem_p), intent(in)     :: expl_elems(:)
    type(t_linsys), intent(inout)       :: linsys
    integer, intent(in)                 :: ie
    integer, intent(in)                 :: ista
    integer, intent(in)                 :: iend
  end subroutine
end interface

!----------------------------------------------------------------------

!> Get the velocity generated by vortical elements on implicit elements.
!!
abstract interface
  subroutine i_get_vort_vel(this, vort_elems)
    import                              :: c_impl_elem
    import                              :: t_vort_elem_p
    import                              :: wp
    implicit none
    class(c_impl_elem), intent(inout)   :: this
    type(t_vort_elem_p), intent(in)     :: vort_elems(:)
  end subroutine
end interface

!----------------------------------------------------------------------

!> Compute the potential induced by an aerodinamic element in a certain
!! position
!!
!! The structure of the subroutine is already intended to be used to calculate
!! the contribution to the linear system matrix and to the right hand side
!! for an equation for the potential
abstract interface
  subroutine i_compute_pot(this, A, b, pos,i,j)
    import                            :: c_pot_elem, wp, t_linsys
    implicit none
    class(c_pot_elem), intent(inout)  :: this
    real(wp), intent(out)             :: A
    real(wp), intent(out)             :: b
    real(wp), intent(in)              :: pos(:)
    integer , intent(in)              :: i, j
  end subroutine
end interface


!> Compute the linear potential induced by an aerodinamic element in a certain
!! position
!!
!! The structure of the subroutine is already intended to be used to calculate
!! the contribution to the linear system matrix and to the right hand side
!! for an equation for the potential
abstract interface
  subroutine i_compute_linear_pot(this, TL, TR, pos,i,j)
    import                            :: c_pot_elem, wp, t_linsys
    implicit none
    class(c_pot_elem), intent(inout)  :: this
    real(wp), intent(out)             :: TL
    real(wp), intent(out)             :: TR
    real(wp), intent(in)              :: pos(:)
    integer , intent(in)              :: i, j
  end subroutine
end interface
!----------------------------------------------------------------------

!> Compute the velocity induced by an aerodynamic element in a certain
!! position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
abstract interface
  subroutine i_compute_vel(this, pos, vel)
    import                    :: c_elem, wp
    implicit none
    class(c_elem), intent(in) :: this
    real(wp), intent(in)      :: pos(:)
    real(wp), intent(out)     :: vel(3)
  end subroutine

end interface

!----------------------------------------------------------------------

!> Compute the velocity induced by an aerodynamic element in a certain
!! position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
abstract interface
  subroutine i_compute_grad(this, pos, grad)
    import                    :: c_elem , wp
    implicit none
    class(c_elem), intent(in) :: this
    real(wp), intent(in)      :: pos(:)
    real(wp), intent(out)     :: grad(3,3)
  end subroutine
end interface

!----------------------------------------------------------------------

!> Compute the psi induced by an aerodinamic element in a certain
! position
!
! The psi is the velocity induced by the singularities multiplied by
! 4*pi.
! The structure of the subroutine is already intended to be used to calculate
! the contribution to the linear system matrix and to the right hand side
! for an equation for the velocity.
abstract interface
  subroutine i_compute_psi(this, A, b, pos, nor,i,j)
    import                            :: c_pot_elem, 
    import                            :: wp 
    import                            :: t_linsys
    implicit none
    class(c_pot_elem), intent(inout)  :: this
    real(wp), intent(out)             :: A
    real(wp), intent(out)             :: b
    real(wp), intent(in)              :: pos(:)
    real(wp), intent(in)              :: nor(:)
    integer , intent(in)              :: i, j
  end subroutine
end interface

!----------------------------------------------------------------------

!> Compute an approximation of the pressure acting on the acutal element
abstract interface
  subroutine i_compute_pres (this, R_g)
    import                            :: c_pot_elem , t_elem_p , wp
    implicit none
    class(c_pot_elem), intent(inout)  :: this
    real(wp)         , intent(in)     :: R_g(3,3)
  end subroutine
end interface

!----------------------------------------------------------------------

!> Compute the elementary force acting on the actual element
abstract interface
  subroutine i_compute_dforce (this)
    import                           :: c_pot_elem , t_elem_p , wp
    implicit none
    class(c_pot_elem), intent(inout) :: this
  end subroutine
end interface

!----------------------------------------------------------------------

!> Compute the geometrical quantities of the elemsnts
abstract interface
  subroutine i_calc_geo_data (this,vert)
    import                           :: c_pot_elem, t_elem_p, wp
    implicit none
    class(c_pot_elem), intent(inout) :: this
    real(wp), intent(in)             :: vert(:,:)
  end subroutine
end interface

!----------------------------------------------------------------------

!> Get the bernoulli source for the pressure equation
!
! This interface is used to avoid the use of select types to get the
! bernoulli source field which is defined only on surface panels
abstract interface
  function i_get_bernoulli_source(this) result(source)
    import                            :: c_impl_elem , wp
    implicit none
    class(c_impl_elem), intent(inout) :: this
    real(wp)                          :: source
  end function
end interface

!----------------------------------------------------------------------

end module mod_aeroel
