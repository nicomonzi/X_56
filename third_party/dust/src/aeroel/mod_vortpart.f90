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
module mod_vortpart

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

public :: t_vortpart, t_vortpart_p, initialize_vortpart


!----------------------------------------------------------------------

type, extends(c_vort_elem) :: t_vortpart
  !> Orientation of the vorticity vector
  real(wp) :: dir(3)
  real(wp) :: vel(3)
  real(wp) :: cen_prev(3)
  real(wp) :: dir_prev(3)
  real(wp) :: mag_prev(3)
  real(wp) :: vel_prev(3) 
  real(wp) :: vel_old(3)
  real(wp) :: stretch(3)
  real(wp) :: stretch_old(3)
  real(wp) :: stretch_prev(3) 
  real(wp) :: stretch_alone(3) 
  logical  :: free=.true.
  real(wp) :: turbvisc
  real(wp) :: rotu(3)
  real(wp) :: r_Vortex
  real(wp) :: r_Vortex_prev
  real(wp) :: r_cutoff
  real(wp) :: vol
  real(wp) :: parent_id
  logical  :: initial_layer=.false.
contains

  procedure, pass(this) :: compute_vel        => compute_vel_vortpart
  procedure, pass(this) :: compute_grad       => compute_grad_vortpart
  procedure, pass(this) :: compute_stretch    => compute_stretch_vortpart
  procedure, pass(this) :: compute_rotu       => compute_rotu_vortpart
  procedure, pass(this) :: compute_diffusion  => compute_diffusion_vortpart
  procedure, pass(this) :: calc_geo_data      => calc_geo_data_vortpart
  procedure, pass(this) :: compute_quantities => compute_quantities_vortpart 

end type

type :: t_vortpart_p

  type(t_vortpart), pointer :: p

end type

character(len=*), parameter :: this_mod_name='mod_vortpart'



!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize vortex line
subroutine initialize_vortpart(this)
  class(t_vortpart), intent(inout) :: this
  
  this%r_Vortex = sim_param%VortexRad
  this%r_cutoff  = sim_param%CutoffRad

end subroutine initialize_vortpart

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex particle in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortpart (this, pos, vel)
  class(t_vortpart), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: vel(3)

  real(wp) :: dist(3), distn, c, d


  dist = pos-this%cen

  !generic
  distn = norm2(dist)
  call kernel_coeffs(this%r_Vortex, distn, c, d) 

  vel = c * cross(dist, this%dir)*this%mag

end subroutine compute_vel_vortpart

!----------------------------------------------------------------------

subroutine compute_grad_vortpart(this, pos, grad)
  class(t_vortpart), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: grad(3,3)

  grad = 0.0_wp

end subroutine compute_grad_vortpart

!----------------------------------------------------------------------

!> Compute the vortex stretching induced by a vortex particle
!! in a prescribed position with a prescribed vorticity (i.e. another particle)
!!
!! WARNING: the calculated term, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_stretch_vortpart (this, pos, alpha, r_Vortex_p, stretch)
  class(t_vortpart), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: alpha(3)
  real(wp), intent(out) :: stretch(3)
  real(wp), intent(in)          :: r_Vortex_p ! vortex rad of the particle p (induced on)
  real(wp) :: dist(3), distn, vecprod(3), c, d !, Sr, dSr
  real(wp) :: r_ave
  !TODO: add far field approximations

  dist = pos-this%cen
  distn = norm2(dist)
  r_ave = (this%r_Vortex + r_Vortex_p)*0.5_wp

  !transpose, generic
  vecprod = cross(alpha, this%dir*this%mag)
  call kernel_coeffs(this%r_Vortex, distn, c, d)
  stretch = -( c * vecprod + d * dist * sum(dist*vecprod) )

end subroutine compute_stretch_vortpart

!----------------------------------------------------------------------

!> Compute the vorticity created by a particle in a prescibed position
!!
!! WARNING: the calculated term, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_rotu_vortpart (this, pos, alpha, r_Vortex_p, rotu)
  class(t_vortpart), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: alpha(3)
  real(wp), intent(out) :: rotu(3)
  real(wp), intent(in)          :: r_Vortex_p ! vortex rad of the particle p (induced on)

  real(wp) :: dist(3), distn, c, d
  real(wp) :: r_ave

  dist = pos-this%cen
  distn = norm2(dist)
  r_ave = (this%r_Vortex + r_Vortex_p)*0.5_wp

  !generic
  call kernel_coeffs(this%r_Vortex, distn, c, d)
  rotu = -2.0_wp * c * this%dir*this%mag + d * cross(dist, cross(dist, this%dir*this%mag))

end subroutine compute_rotu_vortpart

!----------------------------------------------------------------------
!> Compute kernel derivatives coefficients
subroutine kernel_coeffs(r_vort, rr, c, d)
  real(wp), intent(in)  :: r_vort ! vortex radius 
  real(wp), intent(in)  :: rr ! distance between the two particles
  real(wp), intent(out) :: c, d ! coefficients

  real(wp) :: distn, r, r_squared, r_vort_squared

  r = rr 
#if(DUST_KERNEL==1)
  !> Rosenhead
  r_squared = r**2.0_wp
  r_vort_squared = r_vort**2.0_wp
  distn = sqrt(r_squared + r_vort_squared) 
  
  c = -1.0_wp/distn**3.0_wp
  d = 3.0_wp/distn**5.0_wp
#elif(DUST_KERNEL==2) 
  !> HOA (Winckelmans)
  r_squared = r**2.0_wp
  r_vort_squared = r_vort**2.0_wp
  distn = sqrt(r_squared + r_vort_squared)

  c = -(r_squared + 2.5_wp*r_vort_squared)/distn**5.0_wp 
  d = 3.0_wp*(r_squared + 3.5_wp*r_vort_squared)/distn**7.0_wp 
#endif 
  
  !> Other kernels to implement in FMM 
  !elseif (trim(sim_param%kernel) .eq. 'gaussian') then 
  !  if(r.gt.1e-13_wp) then
  !    c = -erf(r/(sqrt(2.0_wp)*r_vort))/r**3.0_wp + &
  !         2.0_wp/(r**2.0_wp * sqrt(2.0_wp*pi) * r_vort) * exp(-r**2.0_wp/(2.0_wp * r_vort**2.0_wp))
  !    d = 3.0_wp/r**5.0_wp * erf(r/(sqrt(2.0_wp)*r_vort)) + exp(-r**2.0_wp/(2.0_wp * r_vort**2.0_wp)) * &
  !        (-6.0_wp/(r**4.0_wp*r_vort*sqrt(2.0_wp*pi)) - 2.0_wp/(r**2.0_wp*r_vort**3.0_wp*sqrt(2.0_wp*pi)))
  !  else
  !    c=0.0_wp 
  !    d=0.0_wp 
  !  endif 
  !elseif (trim(sim_param%kernel) .eq. 'rankine')  then
  !  if (r .ge. r_vort) then
  !    c = -1.0_wp/r**3.0_wp 
  !    d = 3.0_wp/r**5.0_wp
  !  else
  !    c = -1.0_wp/r_vort**3
  !    d = 0.0_wp
  !  endif
  
end subroutine


!> compute quantities of particle (when stretching, divfilter are enabled) 
subroutine compute_quantities_vortpart (this, pos, alpha, r_Vortex_p, vel, stretch, rotu) 
  class(t_vortpart), intent(in) :: this
  real(wp), intent(in)  :: pos(:)
  real(wp), intent(in)  :: alpha(3)
  real(wp), intent(in)  :: r_Vortex_p ! vortex rad of the particle p (induced on)
  
  real(wp), intent(out) :: vel(3)
  real(wp), intent(out) :: stretch(3)
  real(wp), intent(out) :: rotu(3) 
  real(wp) :: dist(3), distn, c, d, r_ave, vecprod(3) 

  !> calc kernel coeffs (compute only once) 
  dist = pos-this%cen
  distn = norm2(dist) 
  r_ave = (this%r_Vortex + r_Vortex_p)*0.5_wp
  call kernel_coeffs(this%r_Vortex, distn, c, d)  

  !> compute velocity
  vel = c * cross(dist, this%dir)*this%mag 
  !> compute stretching 
  vecprod = cross(alpha, this%dir*this%mag) 
  stretch = -( c * vecprod + d * dist * sum(dist*vecprod) ) 
  !> compute rotu 
  rotu = -2.0_wp * c * this%dir*this%mag + d * cross(dist, cross(dist, this%dir*this%mag)) 

end subroutine compute_quantities_vortpart
!----------------------------------------------------------------------
!> Compute the vorticity diffusion induced by a vortex particle
!! in a prescribed position with a prescribed vorticity (i.e. another particle)
!!
subroutine compute_diffusion_vortpart (this, pos, alpha, r_Vortex_p, volp, diff)
  class(t_vortpart), intent(in) :: this ! particle q (inducing)
  real(wp), intent(in)          :: pos(:)
  real(wp), intent(in)          :: alpha(3)
  real(wp), intent(in)          :: r_Vortex_p ! vortex rad of the particle p (induced on)
  real(wp), intent(out)         :: diff(3)
  real(wp)                      :: dist(3), distn
  real(wp), intent(in)          :: volp
  real(wp)                      :: volq
  logical                       :: condition 
  dist = pos-this%cen
  distn = norm2(dist)

  volq = this%vol
  
  !> check singularities (vorticities on the pair of particles are zero) 
  !if (norm2(volp*this%dir*this%mag) .le. 1e-10_wp .and. norm2(volq*alpha) .le. 1e-10_wp ) then 
  !  diff = 0.0_wp 
  !  write(*,*) 'No diffusion due to vol*mag'
  !else 
  !  diff = 1.0_wp/(this%r_Vortex**2.0_wp)*(volp*this%dir*this%mag - volq*alpha)*etaeps(distn,this%r_Vortex) 
  !endif 
  diff = 1.0_wp/(this%r_Vortex**2.0_wp)*(volp*this%dir*this%mag - volq*alpha)*etaeps(distn,this%r_Vortex) 
  !> debug 
  !write(*,*) 'condition', condition
  !write(*,*) 'diff', diff
  !write(*,*) 'volp', volp
  !write(*,*) 'volq', volq
  !write(*,*) 'alpha', alpha
  !write(*,*) 'this%dir*this%mag', this%dir*this%mag
  !write(*,*) 'etaeps', etaeps(distn,this%r_Vortex) 
  !write(*,*) 'distn', distn
  !write(*,*) 'this%r_Vortex', this%r_Vortex
  !write(*,*) 'dist', dist 
  !write(*,*) 'pos', pos   

end subroutine compute_diffusion_vortpart

!----------------------------------------------------------------------

function etaeps(dist, eps) result(eta)
  real(wp), intent(in) :: dist
  real(wp), intent(in) :: eps

  real(wp)             :: eta 
  !> (eq 182 theory manual) 
#if(DUST_KERNEL==1)
  !> low order algebraic
  eta = 15.0_wp/(4.0_wp*pi) / ((dist/eps)**2.0_wp+1.0_wp)**(7.0_wp/2.0_wp)
#elif(DUST_KERNEL==2) 
  !> high order algebraic 
  eta = 105.0_wp/(8.0_wp*pi) / ((dist/eps)**2.0_wp+1.0_wp)**(9.0_wp/2.0_wp)
#endif

  eta = eta/(eps**3.0_wp)

end function etaeps
!----------------------------------------------------------------------
subroutine calc_geo_data_vortpart(this, vert)
  class(t_vortpart), intent(inout) :: this
  real(wp), intent(in) :: vert(:)

  ! center, it is the only coordinate available
  this%cen = vert

end subroutine calc_geo_data_vortpart

!----------------------------------------------------------------------

end module mod_vortpart
