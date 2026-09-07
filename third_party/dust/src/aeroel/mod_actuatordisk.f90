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

!> Module containing the specific subroutines for the actuator disks
!! elements
module mod_actuatordisk

use mod_param, only: &
  wp, pi

use mod_handling, only: &
  error

use mod_sim_param, only: &
  sim_param

use mod_math, only: &
  cross

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet  , &
  gradient_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_actdisk, update_actdisk


!----------------------------------------------------------------------

type, extends(c_expl_elem) :: t_actdisk
  real(wp), allocatable :: traction
  real(wp), allocatable :: radius
contains

  procedure, pass(this) :: compute_pot      => compute_pot_actdisk
  procedure, pass(this) :: compute_linear_pot  => compute_linear_pot_actdisk
  procedure, pass(this) :: compute_vel      => compute_vel_actdisk
  procedure, pass(this) :: compute_grad     => compute_grad_actdisk
  procedure, pass(this) :: compute_psi      => compute_psi_actdisk
  procedure, pass(this) :: compute_pres     => compute_pres_actdisk
  procedure, pass(this) :: compute_dforce   => compute_dforce_actdisk
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_actdisk
end type

character(len=*), parameter :: this_mod_name='mod_actuatordisk'

integer :: it=0

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Compute the potential due to an actuator disk
!!
!! this subroutine employs doublets  to calculate
!! the AIC of an actuator disk on a surface panel, adding the contribution
!! to an equation for the potential.
subroutine compute_pot_actdisk (this, A, b, pos,i,j)
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b
 real(wp), intent(in) :: pos(:)
 integer , intent(in) :: i,j

 real(wp) :: dou

  if ( i .ne. j ) then
    call potential_calc_doublet(this, dou, pos)
  else
!   AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

  A = -dou

  b=0.0_wp

end subroutine compute_pot_actdisk

subroutine compute_linear_pot_actdisk (this, TL, TR, pos,i,j)
  class(t_actdisk), intent(inout) :: this
  real(wp), intent(out) :: TL
  real(wp), intent(out) :: TR
  real(wp), intent(in) :: pos(:)
  integer , intent(in) :: i,j
  

end subroutine compute_linear_pot_actdisk

!----------------------------------------------------------------------

!> Compute the velocity due to an actuator disk
!!
!! This subroutine employs doublets basic subroutines to calculate
!! the AIC coefficients of an actuator disk  to a vortex ring, adding
!! the contribution to an equation for the velocity
subroutine compute_psi_actdisk (this, A, b, pos, nor, i, j )
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: nor(:)
 integer , intent(in) :: i , j

 real(wp) :: vdou(3)

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)


  !  b = ... (from boundary conditions)
  !TODO: consider moving this outside
  if ( i .eq. j ) then
    b =  4.0_wp*pi
  else
    b = 0.0_wp
  end if

end subroutine compute_psi_actdisk

!----------------------------------------------------------------------

!> Compute the velocity induced by an actuator disk in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_actdisk (this, pos, vel)
 class(t_actdisk), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)


  vel = vdou*this%mag


end subroutine compute_vel_actdisk

!----------------------------------------------------------------------

!> Compute the velocity induced by an actuator disk in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_grad_actdisk (this, pos, grad)
 class(t_actdisk), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(out) :: grad(3,3)

 real(wp) :: grad_dou(3,3)


  ! doublet ---
  call gradient_calc_doublet(this, grad_dou, pos)


  grad = grad_dou*this%mag


end subroutine compute_grad_actdisk

!----------------------------------------------------------------------
subroutine compute_cp_actdisk (this, elems)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in) :: elems(:)

 character(len=*), parameter      :: this_sub_name='compute_cp_actdisk'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine compute_cp_actdisk

!----------------------------------------------------------------------

!> The computation of the pressure in the actuator disk is not meant to
!! happen, loads are retrieved from the tables
subroutine compute_pres_actdisk (this, R_g)
 class(t_actdisk) , intent(inout) :: this
 real(wp)         , intent(in)    :: R_g(3,3)
 !type(t_elem_p), intent(in) :: elems(:)

 character(len=*), parameter      :: this_sub_name='compute_pres_actdisk'

! Add this simple routine in order to easily include AD elements in postpro
 this%pres   = this%traction  / this%area

! call error(this_sub_name, this_mod_name, 'This was not supposed to &
! &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_pres_actdisk

!----------------------------------------------------------------------

subroutine compute_dforce_actdisk (this)
 class(t_actdisk), intent(inout) :: this
 !type(t_elem_p), intent(in) :: elems(:)

 character(len=*), parameter      :: this_sub_name='compute_dforce_actdisk'

! Add this simple routine in order to easily include AD elements in postpro
 this%dforce = this%traction * this%nor
!write(*,*) ' debug. this%dforce : ' , this%dforce

! call error(this_sub_name, this_mod_name, 'This was not supposed to &
! &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_dforce_actdisk

!----------------------------------------------------------------------

subroutine update_actdisk(elems_ad, linsys)
 type(t_expl_elem_p), intent(inout) :: elems_ad(:)
 type(t_linsys), intent(inout) :: linsys

 integer :: ie

 do ie=1,size(elems_ad)
   select type(el => elems_ad(ie)%p)
   type is(t_actdisk)

     el%mag = -el%traction*sim_param%dt/(sim_param%rho_inf*pi*el%radius**2)

   end select
 enddo


end subroutine update_actdisk

!----------------------------------------------------------------------

!TODO: is this really necessary?
subroutine solve_actdisk(elems_ll, elems_tot, elems_wake, airfoil_data)
 type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_elem_p), intent(in)    :: elems_tot(:)
 type(t_elem_p), intent(in)    :: elems_wake(:)
 type(t_aero_tab),  intent(in) :: airfoil_data(:)

end subroutine solve_actdisk

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an actuator disk
!!
!! The subroutine calculates all the relevant geometrical quantities of an
!! actuator disk
subroutine calc_geo_data_actdisk(this, vert)
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
 real(wp):: nor(3), tanl(3)
 integer :: nxt

  this%ver = vert
  nsides = this%n_ver

  ! center, for the lifting line is the mid-point
  this%cen =  sum ( this%ver,2 ) / real(nsides,wp)

  this%area = 0.0_wp; this%nor = 0.0_wp
  do is = 1, nsides
    nxt = 1+mod(is,nsides)
    nor = cross(this%ver(:,is) - this%cen,&
                this%ver(:,nxt) - this%cen )
    this%area = this%area + 0.5_wp * norm2(nor)
    this%nor = this%nor + nor/norm2(nor)
  enddo
    this%nor = this%nor/real(nsides,wp)

  ! local tangent unit vector: aligned with first node, normal to n
  tanl = (this%ver(:,1)-this%cen)-&
          sum((this%ver(:,1)-this%cen)*this%nor)*this%nor

  this%tang(:,1) = tanl / norm2(tanl)
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )

  ! vector connecting two consecutive vertices:
  do is = 1 , nsides
    nxt = 1+mod(is,nsides)
    this%edge_vec(:,is) = this%ver(:,nxt) - this%ver(:,is)
  end do

  ! edge: edge_len(:)
  do is = 1 , nsides
    this%edge_len(is) = norm2(this%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nsides
    this%edge_uni(:,is) = this%edge_vec(:,is) / this%edge_len(is)
  end do

  !TODO: is it necessary to initialize it here?
  this%dforce = 0.0_wp
  this%dmom   = 0.0_wp

end subroutine calc_geo_data_actdisk
!----------------------------------------------------------------------

end module mod_actuatordisk
