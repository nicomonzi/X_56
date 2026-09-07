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


!> Module containing the specific subroutines for the vortex lattice
!! type of aerodynamic elements
module mod_vortlatt

use mod_handling, only: &
  internal_error

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet  , &
  gradient_calc_doublet  , &
  linear_potential_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  sim_param

use mod_param, only: &
  wp, pi, one_4pi, max_char_len, prev_tri, next_tri, prev_qua, next_qua

use mod_math, only: &
  cross, dot, linear_interp 

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

use mod_vortpart, only: &
  t_vortpart_p


use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: t_vortlatt


!----------------------------------------------------------------------

!> Planar aerodynamic elements with a surface distribution of doublets
! 
!  Aerodynamic elements characterized by a uniform surface distribution
!  of doublets, employed to model a planar surface.
!  Thanks to the equivalence of a surface distribution of doublets and a
!  closed vortex ring it can be also viewed as ring vortex laying on the
!  edges of the element.
type, extends(c_impl_elem) :: t_vortlatt

  real(wp) :: vel_ctr_pt(3)
  

  !TODO: consider applying the correct element pointer here
  contains

  procedure, pass(this) :: build_row                => build_row_vortlatt
  procedure, pass(this) :: build_row_static         => build_row_static_vortlatt
  procedure, pass(this) :: add_wake                 => add_wake_vortlatt
  procedure, pass(this) :: add_expl                 => add_expl_vortlatt
  procedure, pass(this) :: compute_pot              => compute_pot_vortlatt
  procedure, pass(this) :: compute_vel              => compute_vel_vortlatt
  procedure, pass(this) :: compute_grad             => compute_grad_vortlatt
  procedure, pass(this) :: compute_psi              => compute_psi_vortlatt
  !> Dummy for intel workaround
  procedure, pass(this) :: compute_linear_pot       => compute_linear_pot_vortlatt
  procedure, pass(this) :: compute_pres             => compute_pres_dummy
  procedure, pass(this) :: compute_pres_vortlatt    => compute_pres_vortlatt
  procedure, pass(this) :: compute_dforce           => compute_dforce_dummy
  procedure, pass(this) :: calc_geo_data            => calc_geo_data_vortlatt
  procedure, pass(this) :: get_vort_vel             => get_vort_vel_vortlatt
  procedure, pass(this) :: get_bernoulli_source     => get_bernoulli_source_vortlatt
  procedure, pass(this) :: get_vel_ctr_pt           => get_vel_ctr_pt_vortlatt
  procedure, pass(this) :: compute_dforce_jukowski  => compute_dforce_jukowski_vortlatt

end type

character(len=*), parameter :: this_mod_name='mod_vortlatt'
public 
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Build a row of the linear system for a vortex ring
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static
!! subroutine.
subroutine build_row_vortlatt(this, elems, linsys, ie, ista, iend)
  class(t_vortlatt), intent(inout) :: this
  type(t_impl_elem_p), intent(in)  :: elems(:)
  type(t_linsys), intent(inout)    :: linsys
  integer, intent(in)              :: ie
  integer, intent(in)              :: ista, iend

  integer :: j1
  real(wp) :: b1, wind(3)

  !Not moving components, in the rhs contribution there is no body velocity
  !linsys%b(ie) = sum(linsys%b_static(:,ie) * (-uinf))
  
  linsys%b(ie) = 0.0_wp

    do j1 = 1,ista-1

      wind = variable_wind(elems(j1)%p%cen, sim_param%time)
      linsys%b(ie) = linsys%b(ie) +  &
          linsys%b_static(ie,j1) *sum(elems(j1)%p%nor*(-wind-elems(j1)%p%uvort))
    
    enddo

  ! ista and iend will be the end of the unknowns vector, containing
  ! the moving elements
  do j1 = ista , iend

    call elems(j1)%p%compute_psi( linsys%A(ie,j1), b1,  &
                                  this%cen, this%nor, ie, j1 )

    if (ie .eq. j1) then
      
      !> Diagonal, we are certainly employing vortrin, enforce the b.c. on ie
      wind = variable_wind(elems(ie)%p%cen, sim_param%time)
      linsys%b(ie) = linsys%b(ie) + &
                b1*sum(elems(ie)%p%nor*(elems(ie)%p%ub-wind-elems(j1)%p%uvort))
    
    else

      !> off-diagonal: if it is a vortrin b1 is zero, if it is a surfpan
      !  enforce the boundary condition on it (j1)
      wind = variable_wind(elems(j1)%p%cen, sim_param%time)
      linsys%b(ie) = linsys%b(ie) + &
                b1*sum(elems(j1)%p%nor*(elems(j1)%p%ub-wind-elems(j1)%p%uvort))

    endif

  end do

end subroutine build_row_vortlatt

!----------------------------------------------------------------------

!> Build a static row of the linear system for a vortex ring
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_vortlatt(this, elems, expl_elems, linsys, &
                                  ie, ista, iend)
  class(t_vortlatt), intent(inout) :: this
  type(t_impl_elem_p), intent(in)       :: elems(:)
  type(t_expl_elem_p), intent(in)       :: expl_elems(:)
  type(t_linsys), intent(inout)    :: linsys
  integer, intent(in)              :: ie
  integer, intent(in)              :: ista, iend

  integer :: j1
  real(wp) :: b1

  linsys%b(ie) = 0.0_wp
  linsys%b_static(:,ie) = 0.0_wp

  !> Cycle just all the static elements, ista and iend will be the beginning of
  !  the result vector. Then save the rhs in b_static
  do j1 = ista , iend

    call elems(j1)%p%compute_psi( linsys%A(ie,j1), b1,  &
                                  this%cen, this%nor, ie, j1 )

    linsys%b_static(ie,j1) = b1

  end do

  !> Now build the static contribution from the lifting line elements
  do j1 = 1,linsys%nstatic_expl

    call expl_elems(j1)%p%compute_psi( linsys%L_static(ie,j1), b1,  &
                                  this%cen, this%nor,  1, 2 )

  enddo
  !The rest of the dynamic part will be completed during the first
  ! iteration of the assempling

end subroutine build_row_static_vortlatt

!----------------------------------------------------------------------

!> Add the contribution of the lifting lines to one equation for a vortex ring
!!
!! The rhs of the equation for a vortex ring is updated  adding the
!! the contribution of velocity due to the lifting lines
subroutine add_expl_vortlatt(this, expl_elems, linsys, &
                          ie, ista, iend)
  class(t_vortlatt), intent(inout) :: this
  type(t_expl_elem_p), intent(in)       :: expl_elems(:)
  type(t_linsys), intent(inout)    :: linsys
  integer, intent(in)              :: ie
  integer, intent(in)             :: ista
  integer, intent(in)             :: iend

  integer :: j1
  real(wp) :: a, b

  !> Static part: take what was already computed
  do  j1 = 1, ista-1
    linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*expl_elems(j1)%p%mag
  enddo

  !> Add the explicit vortex panel wake contribution to the rhs
  do j1 = ista, iend

    call expl_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )


    linsys%b(ie) = linsys%b(ie) - a*expl_elems(j1)%p%mag

  end do

end subroutine add_expl_vortlatt

!----------------------------------------------------------------------

!> Add the contribution of the wake to one equation for a vortex ring
!!
!! The rhs of the equation for a surface panel is updated  adding the
!! the contribution of velocity due to the wake
subroutine add_wake_vortlatt(this, wake_elems, impl_wake_ind, linsys, &
                          ie, ista, iend)
  class(t_vortlatt), intent(inout)  :: this
  type(t_pot_elem_p), intent(in)    :: wake_elems(:)
  integer, intent(in)               :: impl_wake_ind(:,:)
  type(t_linsys), intent(inout)     :: linsys
  integer, intent(in)               :: ie
  integer, intent(in)               :: ista
  integer, intent(in)               :: iend

  integer :: j1, ind1, ind2
  real(wp) :: a, b
  integer :: n_impl

  !Count the number of implicit wake contributions
  n_impl = size(impl_wake_ind,2)

  !Add the contribution of the implicit wake panels to the linear system
  !Implicitly we assume that the first set of wake panels are the implicit
  !ones since are at the beginning of the list
  do j1 = 1, n_impl
    ind1 = impl_wake_ind(1,j1) 
    ind2 = impl_wake_ind(2,j1)
    if ((ind1.ge.ista .and. ind1.le.iend)) then
      call wake_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )
      linsys%A(ie,ind1) = linsys%A(ie,ind1) + a
      if ( ind2 .ne. 0 ) then 
        linsys%A(ie,ind2) = linsys%A(ie,ind2) - a
      endif 
    endif
  end do

  ! Add the explicit vortex panel wake contribution to the rhs
  do j1 = n_impl+1 , size(wake_elems)
    call wake_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )
    linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%mag
  end do

end subroutine add_wake_vortlatt

!----------------------------------------------------------------------

!> Compute the potential due to a vortex ring
!!
!! this subroutine employs doublets basic subroutines to calculate
!! the AIC of a vortex ring on a surface panel. The contribution to its rhs
!! is zero since there are no sources (and no b.c. enforcing)
subroutine compute_pot_vortlatt(this, A, b, pos,i,j)
  class(t_vortlatt), intent(inout) :: this
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

  b = 0.0_wp

end subroutine compute_pot_vortlatt

!> Compute the linear potential due to a vortex ring
!! this subroutine employs doublets basic subroutines to calculate
!! the AIC of a vortex ring on a surface panel. The contribution to its rhs
!! is zero since there are no sources (and no b.c. enforcing)
subroutine compute_linear_pot_vortlatt(this, TL, TR, pos,i,j)
  class(t_vortlatt), intent(inout) :: this
  real(wp), intent(out) :: TL
  real(wp), intent(out) :: TR
  real(wp), intent(in) :: pos(:)
  integer , intent(in) :: i,j

  call linear_potential_calc_doublet(this, TL, TR, pos) 

end subroutine compute_linear_pot_vortlatt

!----------------------------------------------------------------------

!> Compute the velocity due to a vortex ring
!!
!! This subroutine employs basic doublet subroutines to calculate the AIC of
!! a vortex ring on a vortex ring. The contribution to the rhs is the boundary
!! condition and it is enforced only if the influencing and influenced
!! vortex rings are the same
subroutine compute_psi_vortlatt(this, A, b, pos, nor, i, j )
  class(t_vortlatt), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: nor(:)
  integer , intent(in) :: i , j

  real(wp) :: vdou(3)

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)

  ! b = ... (from boundary conditions)
  !TODO: consider moving this outside
  if ( i .eq. j ) then

    b =  4.0_wp*pi

  else

    b = 0.0_wp

  end if

end subroutine compute_psi_vortlatt

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex ring in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! doublets
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortlatt(this, pos, vel)
  class(t_vortlatt), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: vel(3)

  real(wp) :: vdou(3)

  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  vel = vdou*this%mag

end subroutine compute_vel_vortlatt

!----------------------------------------------------------------------

!> Compute the gradient induced by a vortex ring in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! doublets
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_grad_vortlatt(this, pos, grad)
class(t_vortlatt), intent(in) :: this
real(wp), intent(in) :: pos(:)
real(wp), intent(out) :: grad(3,3)

real(wp) :: grad_dou(3,3)

! doublet ---
call gradient_calc_doublet(this, grad_dou, pos)

grad = grad_dou*this%mag

end subroutine compute_grad_vortlatt


!----------------------------------------------------------------------

!> Compute an approximate value of the mean DELTA pressure on the actual
!! element
!!
!! pres = " DELTA pressure = ( pressure lower - pressure upper ) "
!!  s.t. vec{dforce} = pres * vec{n}  ( since vec{n} = vec{n_upper} )
!!
!! see compute_dforce_vortlatt
subroutine compute_pres_vortlatt(this) 
  class(t_vortlatt), intent(inout) :: this
  integer  :: i_stripe
  real(wp) :: wind(3)

  this%pres = 0.0_wp
  i_stripe = size(this%stripe_elem)

  if ( i_stripe .gt. 1 ) then

    wind = variable_wind(this%cen, sim_param%time)
    this%pres = - sim_param%rho_inf * &
          ( norm2(wind + this%uvort - this%ub) * this%dy / this%area * &
          ( this%mag - this%stripe_elem(i_stripe-1)%p%mag ) + &
              this%didou_dt ) 

  else

    wind = variable_wind(this%cen, sim_param%time)
    this%pres = - sim_param%rho_inf * &
          ( norm2(wind + this%uvort - this%ub) * this%dy / this%area * &
                this%mag + &
              this%didou_dt ) 
  end if

end subroutine compute_pres_vortlatt

!> Dummy version doing nothing
subroutine compute_pres_dummy(this, R_g)
  class(t_vortlatt), intent(inout) :: this
  real(wp)         , intent(in)    :: R_g(3,3)

end subroutine compute_pres_dummy


!> Compute the elementary force on the on the actual element
subroutine compute_dforce_vortlatt(this)
class(t_vortlatt), intent(inout) :: this

end subroutine compute_dforce_vortlatt

!> Dummy version doing nothing
subroutine compute_dforce_dummy(this)
class(t_vortlatt), intent(inout) :: this

end subroutine compute_dforce_dummy

!> Compute the elementary force on the on the actual element
!  using Kutta-Jukowski theorem as implemented in AVL:
!  the velocity
!
subroutine compute_dforce_jukowski_vortlatt(this, correction)

  class(t_vortlatt), intent(inout) :: this
  logical, intent(in)              :: correction 
  real(wp)                         :: gam(3), wind(3)

  integer                          :: i_stripe
  real(wp)                         :: mach
  
  ! Prandt -- Glauert correction for compressibility effect
  wind = variable_wind(this%cen, sim_param%time)
  
  mach = abs(sum(this%vel_ctr_pt*this%tang(:,2)) / sim_param%a_inf)
  
  ! === Steady contribution (KJ) ===
  gam = cross ( this%vel_ctr_pt, this%edge_vec(:,1) )
  i_stripe = size(this%stripe_elem)
  if (correction) then 
    !> apply the prantdl glauert correction (for non corrected vl)
    if ( i_stripe .gt. 1 ) then
      this%dforce = sim_param%rho_inf * gam &
                  * ( this%mag - this%stripe_elem(i_stripe-1)%p%mag )  / sqrt(1 - mach**2)
    else
      this%dforce = sim_param%rho_inf * gam * this%mag / sqrt(1 - mach**2)
    end if
  else 
    if ( i_stripe .gt. 1 ) then
      this%dforce = sim_param%rho_inf * gam &
                  * ( this%mag - this%stripe_elem(i_stripe-1)%p%mag )
    else
      this%dforce = sim_param%rho_inf * gam * this%mag 
    end if
  endif   
  
  ! === Unsteady contribution ===
  this%dforce = this%dforce &
              - sim_param%rho_inf * this%area * ( &
                this%didou_dt * this%nor)


end subroutine compute_dforce_jukowski_vortlatt




!----------------------------------------------------------------------
!> Compute an approximate value of the induced velocity at 1/4 of chord
! of an element. This routine re-computes the contributions of the
! potential elements only:
! - body panels
! - wake panels
! ------------------------------------------------------------------- !
! This routine uses the value in the centre of the panels of:         !
! - the free-stream and the body velocity                             !
! - the rotational part of the velocity, collected in this%uvort      !
! ------------------------------------------------------------------- !
subroutine get_vel_ctr_pt_vortlatt(this, elems, wake_elems, vort_elems)
  class(t_vortlatt),    intent(inout)   :: this
  type(t_pot_elem_p),   intent(in)      :: elems(:)
  type(t_pot_elem_p),   intent(in)      :: wake_elems(:)
  type(t_vort_elem_p),  intent(in)      :: vort_elems(:)

  real(wp) :: v(3), x0(3), wind(3)
  integer :: j

  !> Initialisation to zero
  this%vel_ctr_pt = 0.0_wp

  !> Control point at 1/4-fraction of the chord
  x0 = this%cen + (this%edge_vec(:,4) - this%edge_vec(:,2))/4.0_wp

  !=== Compute the velocity from all the elements ===
  do j = 1,size(wake_elems)  ! wake panels
    call wake_elems(j)%p%compute_vel(x0, v)
    this%vel_ctr_pt = this%vel_ctr_pt + v
  enddo

  do j = 1,size(vort_elems) ! wake vorticity elements 
    call vort_elems(j)%p%compute_vel(x0, v)
    this%vel_ctr_pt = this%vel_ctr_pt + v
  enddo

  do j = 1,size(elems) ! body elements
    call elems(j)%p%compute_vel(x0, v)
    this%vel_ctr_pt = this%vel_ctr_pt + v
  enddo

  wind = variable_wind(this%cen, sim_param%time)
  this%vel_ctr_pt = this%vel_ctr_pt/(4.0_wp * pi) &
                    + wind - this%ub

end subroutine get_vel_ctr_pt_vortlatt

!!!----------------------------------------------------------------------
!> Calculate the geometrical quantities of a vortex lattice
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! vortex lattice panel
subroutine calc_geo_data_vortlatt(this, vert)
  class(t_vortlatt), intent(inout) :: this
  real(wp), intent(in) :: vert(:,:)

  integer :: nsides, is
  real(wp):: nor(3), tanl(3)

  this%ver = vert
  nsides = this%n_ver
  this%cen =  sum ( this%ver,2 ) / real(nsides,wp)

  ! unit normal and area
  if ( nsides .eq. 4 ) then

    nor = cross( this%ver(:,3) - this%ver(:,1) , &
                this%ver(:,4) - this%ver(:,2)     )

  else if ( nSides .eq. 3 ) then

    nor = cross( this%ver(:,3) - this%ver(:,2) , &
                this%ver(:,1) - this%ver(:,2)     )

  end if
  
  this%area = 0.5_wp * norm2(nor)
  ! avoid numerical singularities when compiled in debug mode
  if (norm2(nor) .lt. 1e-16_wp) then 
    nor(3) = 1e-16_wp
  end if

  this%nor = nor / norm2(nor)
  
  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

  this%tang(:,1) = tanl / norm2(tanl)
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )

  ! vector connecting two consecutive this%verices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then

    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_tri(is)) - this%ver(:,is)
    end do

  else if ( nsides .eq. 4 ) then

    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_qua(is)) - this%ver(:,is)
    end do

  end if

  ! edge: edge_len(:)
  do is = 1 , nsides
    this%edge_len(is) = norm2(this%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    this%edge_uni(:,is) = this%edge_vec(:,is) / (this%edge_len(is))
  end do

  !TODO: is it necessary to initialize it here?
  this%dforce = 0.0_wp
  this%dmom   = 0.0_wp

end subroutine calc_geo_data_vortlatt

!> Calc geo data of the stripe: equivalent to calc_geo_data_liftlin 


!----------------------------------------------------------------------

subroutine get_vort_vel_vortlatt(this, vort_elems)
  class(t_vortlatt), intent(inout)  :: this
  type(t_vort_elem_p), intent(in)    :: vort_elems(:)

  integer :: iv
  real(wp) :: vel(3)

  !> Initialize to zero, before accumuulation
  this%uvort = 0.0_wp

  do iv=1,size(vort_elems)

    call vort_elems(iv)%p%compute_vel(this%cen, vel)
    this%uvort = this%uvort + vel*one_4pi

  enddo

end subroutine

!----------------------------------------------------------------------

function get_bernoulli_source_vortlatt(this) result(source)
  class(t_vortlatt), intent(inout) :: this
  real(wp) :: source

  character(len=*), parameter :: this_sub_name = 'get_bernoulli_source_vortlatt'

  !this is just a dummy, should never be used
  source = 0.0_wp
  !and for this reason here is an internal error
  call internal_error(this_sub_name, this_mod_name, &
                      'getting bernoulli source from a vortex lattice')

end function

!----------------------------------------------------------------------

end module mod_vortlatt
