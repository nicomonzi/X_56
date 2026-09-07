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

!> Module containing the specific subroutines for a single vortex line
module mod_vortline

use mod_param, only: &
  wp, pi, max_char_len, prev_tri, next_tri, prev_qua, next_qua

use mod_handling, only: &
  error, printout

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  sim_param

use mod_math, only: &
  cross

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_vortline, initialize_vortline


!----------------------------------------------------------------------

type, extends(c_vort_elem) :: t_vortline
  real(wp) :: ver(3,2)
  real(wp) :: ver_vel(3,2)
  real(wp) :: edge_vec(3)
  real(wp) :: edge_uni(3)
  real(wp) :: edge_len
contains

  procedure, pass(this) :: compute_vel      => compute_vel_vortline
  procedure, pass(this) :: compute_grad     => compute_grad_vortline
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_vortline

end type

character(len=*), parameter :: this_mod_name='mod_vortline'

real(wp) :: r_Rankine
real(wp) :: r_cutoff

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize vortex line
subroutine initialize_vortline()

  r_Rankine = sim_param%RankineRad
  r_cutoff  = sim_param%CutoffRad

end subroutine initialize_vortline

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex line in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortline (this, pos, vel)
  class(t_vortline), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: vel(3)

  real(wp) :: vdou(3)
  real(wp) :: av(3) , hv(3)
  real(wp) :: ai    , hi
  real(wp) :: R1 , R2
  real(wp) :: r_Ran

  !TODO: add far field approximations
  !radius_v = pos-this%cen
  !radius   = norm2(radius_v)
  if(associated(this%mag)) then

    ! use this%ver instead of its projection this%verp
    av = pos-this%ver(:,1)
    ai = sum(av*this%edge_uni(:))
    R1 = norm2(av)
    R2 = norm2(pos-this%ver(:,2))
    hv = av - ai*this%edge_uni(:)
    hi = norm2(hv)
    if ( hi .gt. this%edge_len*r_Rankine ) then
      ! (a/r+(s-a)/r)/h
      vdou = ( (this%edge_len-ai)/r2 + ai/r1 )/(hi**2.0_wp) * &
                      cross(this%edge_uni(:),hv)
    else
!     ! (a/r+(s-a)/r)* h/r_Rankine^2.
      if ((R1 .gt. this%edge_len*r_cutoff) .and. &! avoid singularity
          (R2 .gt. this%edge_len*r_cutoff)) then
        r_Ran = r_Rankine * this%edge_len
        vdou = ((this%edge_len-ai)/R2 + ai/R1)/(r_Ran**2.0_wp)* &
                         cross(this%edge_uni(:),hv)
     else
      vdou = 0.0_wp
      end if
    end if
    vel = vdou*this%mag

  else
    vel = 0.0_wp
  endif

end subroutine compute_vel_vortline

!----------------------------------------------------------------------

subroutine compute_grad_vortline(this, pos, grad)
  class(t_vortline), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: grad(3,3)

  integer  :: i1 , i2
  real(wp) :: R1(3) , R2(3) , a1(3) , a2(3) , l(3) , a
  real(wp) :: R1v(3,1) , R2v(3,1) , a1v(3,1) , a2v(3,1) , lv(3,1)
  real(wp) :: lx(3,3) , aa1(3,3) , aa2(3,3) , ax1(3,3) , ax2(3,3) , al1(3,3) , al2(3,3)
  real(wp) :: del , a2del2

  !TODO: add far field approximations
  !radius_v = pos-this%cen
  !radius   = norm2(radius_v)
  del = sim_param % RankineRad
  !del = sim_param % VortexRad

  if(associated(this%mag)) then

    i1 = 1 ;  i2 = 2

    l = this%edge_uni(:)
    lv(:,1) = l

    R1 = pos-this%ver(:,i1) ;   a1 = cross( l , R1 )
    R2 = pos-this%ver(:,i2) ;   a2 = cross( l , R2 )
    a = norm2(a1)  ! = norm(a2)
    a2del2 = a**2.0_wp + del**2.0_wp

    R1v(:,1) = R1 ;  a1v(:,1) = a1
    R2v(:,1) = R2 ;  a2v(:,1) = a2

    lx(:,1) = (/  0.0_wp ,  l(3)   , -l(2)   /)
    lx(:,2) = (/ -l(3)   ,  0.0_wp ,  l(1)   /)
    lx(:,3) = (/  l(2)   , -l(1)   ,  0.0_wp /)

    aa1 = matmul( a1v , transpose(a1v) )
    ax1 = matmul( a1v , transpose(R1v) )
    al1 = matmul( a1v , transpose( lv) )
    aa2 = matmul( a2v , transpose(a2v) )
    ax2 = matmul( a2v , transpose(R2v) )
    al2 = matmul( a2v , transpose( lv) )

    grad = &
     + 1.0_wp / ( a2del2**1.5_wp * norm2(R1) ) * &
       ( (   a * lx &
           + ( 1.0_wp/a - 3.0_wp*a/a2del2 ) * matmul(aa1,lx) ) * sum(l*R1) &
       + a * ( al1 - ax1 * sum( l * R1 ) / norm2(R1)**2.0_wp ) )  &
     - 1.0_wp / ( a2del2**1.5_wp * norm2(R2) ) * &
       ( (   a * lx &
           + ( 1.0_wp/a - 3.0_wp*a/a2del2 ) * matmul(aa2,lx) ) * sum(l*R2) &
       + a * ( al2 - ax2 * sum( l * R2 ) / norm2(R2)**2.0_wp ) )

    grad = grad*this%mag

  else
    grad = 0.0_wp
  endif


end subroutine compute_grad_vortline

!----------------------------------------------------------------------

subroutine calc_geo_data_vortline(this, vert)
 class(t_vortline), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)


  this%ver = vert


  ! center, for lines  is the mid-point
  this%cen =  sum ( this%ver(:,1:2),2 ) / 2.0_wp

  !! vector connecting the two consecutive vertices:
  this%edge_vec(:) = this%ver(:,2) - this%ver(:,1)

  this%edge_len = norm2(this%edge_vec(:))

  this%edge_uni(:) = this%edge_vec(:) / this%edge_len

end subroutine calc_geo_data_vortline
!----------------------------------------------------------------------

end module mod_vortline
