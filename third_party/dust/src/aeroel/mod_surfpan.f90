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
!!          Alberto Savino
!!          Alessandro Cocco
!!=========================================================================


!> Module to treat surface doublet + source panels

module mod_surfpan

use mod_param, only: &
  wp , nl, &
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi, one_4pi, max_char_len

use mod_handling, only: &
  error, warning, printout

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet  , &
  gradient_calc_doublet  , & 
  linear_potential_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_math, only: &
  cross , compute_qr

use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: t_surfpan, initialize_surfpan

private

!----------------------------------------------------------------------

!> Surface panel (Morino kind) with uniform distribution of doublets
!! and sources
!!
!! This type inherits most of its members from the general \ref c_elem class
!! however implements its own subroutines to calculate the coefficients of the
!! linear system and the
type, extends(c_impl_elem) :: t_surfpan

  real(wp), allocatable :: pot_vel_stencil(:,:)
  real(wp), allocatable :: cosTi(:) , sinTi(:)
  real(wp), allocatable :: verp(:,:)
  real(wp)              :: surf_vel(3)
  !> stencil for the Constrained Hermite Taylor Series
  ! Least Square computation of derivatives
  real(wp), allocatable ::   chtls_stencil(:,:)

  !> surface quantities for Bernoulli integral equation
  real(wp) :: dUn_dt
  real(wp) :: bernoulli_source
  real(wp), pointer :: pres_sol

  !> boundary layer and flow separation
  real(wp) :: h_bl         ! height of the surface (boundary??) layer
  real(wp) :: al_free      ! ratio: shed vorticity / total vorticity
  real(wp) :: surf_vort(3) ! free vorticity
  real(wp) :: free_vort(3) ! free vorticity

contains

  procedure, pass(this) ::  build_row               => build_row_surfpan
  procedure, pass(this) ::  build_row_static        => build_row_static_surfpan
  procedure, pass(this) ::  add_wake                => add_wake_surfpan
  procedure, pass(this) ::  add_expl                => add_expl_surfpan
  procedure, pass(this) ::  compute_pot             => compute_pot_surfpan
  procedure, pass(this) ::  compute_linear_pot      => compute_linear_pot_surfpan
  procedure, pass(this) ::  compute_vel             => compute_vel_surfpan
  procedure, pass(this) ::  compute_grad            => compute_grad_surfpan
  procedure, pass(this) ::  compute_psi             => compute_psi_surfpan
  procedure, pass(this) ::  compute_pres            => compute_pres_surfpan
  procedure, pass(this) ::  compute_dforce          => compute_dforce_surfpan
  procedure, pass(this) ::  calc_geo_data           => calc_geo_data_surfpan
  procedure, pass(this) ::  get_vort_vel            => get_vort_vel_surfpan

  procedure, pass(this) ::  create_local_velocity_stencil => &
                            create_local_velocity_stencil_surfpan
  procedure, pass(this) ::  create_chtls_stencil => &
                            create_chtls_stencil_surfpan

  procedure, pass(this) ::  get_bernoulli_source => get_bernoulli_source_surfpan

end type

real(wp) :: ff_ratio

character(len=max_char_len) :: msg(3)

character(len=*), parameter :: this_mod_name = 'mod_surfpan'
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Subroutine to populate the module variables from input (TODO move to a separate module)
!!
subroutine initialize_surfpan()

  ff_ratio = sim_param%FarFieldRatioSource

end subroutine initialize_surfpan

!----------------------------------------------------------------------
!> compute AIC of panel <this>, on the control point <pos>
!! Compute I_ik : velocity potential in <pos>_i induced by a (-4*pi)
!!   -intensity surface source  on the panel <this>_k.         ^---------
!!
!! I_ik = int_{S_k} { 1 / |r_i - r| }
!!
!! The relation with unitary surface doublet D_ik is:
!!   S_ik = -(1/4pi) * I_ik
!!
subroutine potential_calc_sou_surfpan(this, sou, dou, pos)
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out) :: sou
  real(wp), intent(in)  :: dou
  real(wp), intent(in) :: pos(:)

  real(wp) :: radius

  real(wp), dimension(3) :: e3
  real(wp) :: zQ , souLog , vi , R1 , R2
  real(wp), dimension(3) :: Qp
  integer :: indm1 , indp1
  real(wp) :: den

  real(wp), parameter :: eps_sou  = 1.0e-6_wp
  real(wp), parameter :: ff_ratio = 10.0_wp

  integer :: i1
  character(len=*), parameter :: this_sub_name = 'potential_calc_sou_surfpan'

  radius = norm2(pos-this%cen)

  if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then
    
    ! far-field approximation (1)
    sou = this%area / radius

    if ( sou .ne. sou ) then
      call error(this_sub_name, this_mod_name, 'Divide by zero in sources &
        &far field approximation')
    end if

  else

    ! initialisation
    sou = 0.0_wp

    ! unit normal
    e3 = this%nor

    ! Control point (Q)
    zQ = dot_product((pos-this%cen) , e3)
    Qp = pos - zQ * e3

    do i1 = 1 , this%n_ver

      if ( this%n_ver .eq. 3 ) then
        indm1 = prev_tri(i1)
        indp1 = next_tri(i1)
      else if ( this%n_ver .eq. 4 ) then
        indm1 = prev_qua(i1)
        indp1 = next_qua(i1)
      end if

      R1 = norm2( pos - this%verp(:,i1) )
      R2 = norm2( pos - this%verp(:,indp1) )

      den = R1+R2-this%edge_len(i1)
      if ( den .lt. 1e-6_wp ) then

        if(sim_param%debug_level .ge.5) then
          write(msg(1),'(A)') 'Too small denominator in &
          &source computation with point projection, using actual &
          &points instead.'
          write(msg(2),'(A,F12.6,F12.6,F12.6)') 'Computing sources on point: ',&
          pos(1),pos(2),pos(3)
          write(msg(3),'(A)')'This is most likely due to severely warped &
          &quadrilateral elements adjacent to small elements.'//nl//&
          &'      === CHECK MESH QUALITY! ==='

          call warning(this_sub_name, this_mod_name, msg)
        endif

        R1 = norm2( pos - this%ver(:,i1) )
        R2 = norm2( pos - this%ver(:,indp1) )
        den = R1+R2-this%edge_len(i1)

      end if

      souLog = log((R1+R2+this%edge_len(i1)) / (den))

      vi = - dot_product(cross(Qp-this%verp(:,i1), this%edge_vec(:,i1)), e3) / this%edge_len(i1)
      sou = sou + vi * souLog

    end do

    sou = sou - zQ * dou

  end if

end subroutine potential_calc_sou_surfpan

!----------------------------------------------------------------------

!> compute velocity AIC of panel <this>, on the control point <pos>
!! Compute the velocity induced in <pos>_i by the surface source <this>_k
!!  with intentsity intensity -4*pi. <------
!!
!! V^{sou}_ik = grad_{r_i} { int_{S_k} { 1 / |r_i - r| } }
!!
subroutine velocity_calc_sou_surfpan(this, vel, pos)
  class(t_surfpan), intent(in) :: this
  real(wp), intent(out) :: vel(3)
  real(wp), intent(in) :: pos(:)

  real(wp) :: phix , phiy , pdou
  real(wp) :: R1 , R2 , souLog

  real(wp), parameter :: ff_ratio = 10.0_wp
  real(wp) :: radius_v(3)
  real(wp) :: radius

  integer :: indm1 , indp1
  integer :: i1
  character(len=max_char_len) :: message
  character(len=*), parameter :: this_sub_name='velocity_calc_sou_surfpan'

  radius_v = pos-this%cen
  radius   = norm2(radius_v)

  if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1)

    phix = - this%area * sum( radius_v*this%tang(:,1) ) / radius**3.0_wp
    phiy = - this%area * sum( radius_v*this%tang(:,2) ) / radius**3.0_wp
    pdou = - this%area * sum( radius_v*this%nor       ) / radius**3.0_wp

  else

    phix = 0.0_wp
    phiy = 0.0_wp

    do i1 = 1 , this%n_ver

      if ( this%n_ver .eq. 3 ) then
        indm1 = prev_tri(i1)
        indp1 = next_tri(i1)
      else if ( this%n_ver .eq. 4 ) then
        indm1 = prev_qua(i1)
        indp1 = next_qua(i1)
      end if

      R1 = norm2(pos - this%verp(:,i1))
      R2 = norm2(pos - this%verp(:,indp1))

      if(sim_param%debug_level .ge.5) then
        if ( abs(R1+R2-this%edge_len(i1)) .lt. 1e-6_wp ) then
          call warning(this_sub_name, this_mod_name, &
            'too small denominator in calculation of velocity')
          write(message,*) 'R1, R2, this%edge_len, i1', R1, R2, this%edge_len(i1), i1
          call printout(message)
      end if
        if ( abs(this%edge_len(i1) ) .lt. 1e-6_wp ) then
          call warning(this_sub_name, this_mod_name, &
            'too small edge length in calculation of velocity')
          write(message,*) 'R1, R2, this%edge_len, i1', R1, R2, this%edge_len(i1), i1
          call printout(message)
        end if
        if ( abs(R1+R2+this%edge_len(i1)) .lt. 1e-6_wp ) then
          call warning(this_sub_name, this_mod_name, &
            'too small numerator in calculation of velocity')
          write(message,*) 'R1, R2, this%edge_len, i1', R1, R2, this%edge_len(i1), i1
          call printout(message)
        end if
      endif

      if ( R1+R2-this%edge_len(i1) .lt. 1e-12_wp ) then
        souLog = 0.0_wp
      else
        !> katz eqs 10.95 - 10.96 (the sign in the book is wrong)
        souLog = log( (R1+R2 + this%edge_len(i1)) / (R1+R2 - this%edge_len(i1)) )
      endif

      phix = phix + this%sinTi(i1) * souLog
      phiy = phiy - this%cosTi(i1) * souLog

    end do

    call potential_calc_doublet(this, pdou, pos)

    phix = - phix 
    phiy = - phiy 
    pdou = - pdou 

  end if

  ! vsou = (/ phix , phiy , pdou /) 
  !> Rotate the velocity vector from the local coordinate system to the global one 
  !> rot = (tang(:,:), nor) 
  vel(1) = this%tang(1,1)*phix + this%tang(1,2)*phiy + this%nor(1)* pdou
  vel(2) = this%tang(2,1)*phix + this%tang(2,2)*phiy + this%nor(2)* pdou
  vel(3) = this%tang(3,1)*phix + this%tang(3,2)*phiy + this%nor(3)* pdou

end subroutine velocity_calc_sou_surfpan

!----------------------------------------------------------------------
!> TODO: compute the gradient of the velocity induced by a source and
!        write this routine
subroutine gradient_calc_sou_surfpan(this, grad, pos)
  class(t_surfpan), intent(in) :: this
  real(wp), intent(in)         :: pos(:)
  real(wp), intent(out)        :: grad(3,3)

  real(wp)                     :: rot(3,3), grad_loc(3,3), dlog
  real(wp)                     :: csi, eta, zeta, v_dou(3) 
  real(wp)                     :: csi_1, eta_1, zeta_1
  real(wp)                     :: csi_2, eta_2, zeta_2 
  real(wp)                     :: phicsi_csi,  phieta_csi,  phizeta_csi 
  real(wp)                     :: phicsi_eta,  phieta_eta,  phizeta_eta
  real(wp)                     :: phicsi_zeta, phieta_zeta, phizeta_zeta
  real(wp)                     :: R1, R2 
  integer                      :: indm1, indp1 
  integer                      :: i1
  character(len=max_char_len)  :: message
  character(len=*), parameter  :: this_sub_name='velocity_calc_sou_surfpan'

  !> rotation matrix from local to global coordinates 
  rot(:,1) = this%tang(:,1)
  rot(:,2) = this%tang(:,2)
  rot(:,3) = this%nor(:)
  
  !> initialization
  grad_loc = 0.0_wp; grad = 0.0_wp; dlog = 0.0_wp
  !> structure of the gradient matrix (local coordinates)  
  phicsi_csi = 0.0_wp;  phieta_csi = 0.0_wp; phizeta_csi = 0.0_wp 
  phicsi_eta = 0.0_wp;  phieta_eta = 0.0_wp; phizeta_eta = 0.0_wp
  phicsi_zeta = 0.0_wp; phieta_zeta = 0.0_wp; phizeta_zeta = 0.0_wp

  !> TODO add far field approximation for the gradient of the source 
  !> loop over the vertices of the panel 
  do i1 = 1 , this%n_ver

    if (this%n_ver .eq. 3) then
      indm1 = prev_tri(i1)
      indp1 = next_tri(i1)
    else if (this%n_ver .eq. 4) then
      indm1 = prev_qua(i1)
      indp1 = next_qua(i1)
    end if

    R1 = norm2(pos - this%verp(:,i1))    !> R_{k}
    R2 = norm2(pos - this%verp(:,indp1)) !> R_{k+1} 

    !> R_{k} and R_{k+1} in local coordinates {csi, eta, zeta}
    csi = dot_product(rot(:,1), pos)
    eta = dot_product(rot(:,2), pos)
    zeta = dot_product(rot(:,3), pos)

    csi_1 = dot_product(rot(:,1), this%verp(:,i1))
    eta_1 = dot_product(rot(:,2), this%verp(:,i1))
    zeta_1 = dot_product(rot(:,3), this%verp(:,i1)) !> should be 0.0
    
    csi_2 = dot_product(rot(:,1), this%verp(:,indp1)) 
    eta_2 = dot_product(rot(:,2), this%verp(:,indp1))
    zeta_2 = dot_product(rot(:,3), this%verp(:,indp1)) !> should be 0.0
    
    !> log derivative first part  
    dlog = 2.0_wp*this%edge_len(i1)/((R1+R2)**2.0_wp - this%edge_len(i1)**2)

    !> grad(1,1) = dphicsi/dcsi
    phicsi_csi = phicsi_csi + this%sinTi(i1) * dlog * &
                  ((csi - csi_1)/R1 + (csi - csi_2)/R2) 
    !> grad(2,1) = dphicsi/deta
    phicsi_eta = phicsi_eta + this%sinTi(i1) * dlog * &
                  ((eta - eta_1)/R1 + (eta - eta_2)/R2) 
    !> grad(3,1) = dphicsi/dzeta
    phicsi_zeta = phicsi_zeta + this%sinTi(i1) * dlog * &
                  (zeta/R1 + zeta/R2) 
    !> grad(1,2) = dphicsi/dcsi
    phieta_csi = phieta_csi - this%cosTi(i1) * dlog * &
                  ((csi - csi_1)/R1 + (csi - csi_2)/R2) 
    !> grad(2,2) = dphicsi/deta
    phieta_eta = phieta_eta - this%cosTi(i1) * dlog * &
                  ((eta - eta_1)/R1 + (eta - eta_2)/R2) 
    !> grad(3,2) = dphicsi/dzeta
    phieta_zeta = phicsi_zeta - this%cosTi(i1) * dlog * &
                  (zeta/R1 + zeta/R2) 
    
  enddo

  !> the source is minus the integral of the doublet, therefore the zeta gradient components 
  !> is the induced velocity of a doublet in the local zeta direction   
  call velocity_calc_doublet(this, v_dou, pos) !> velocity induced by the doublet

  !> velocity induced by the doublet in local coordinates
  v_dou = -matmul(transpose(rot), v_dou) 

  phizeta_csi = v_dou(1)
  phizeta_eta = v_dou(2)
  phizeta_zeta = v_dou(3) 

  !> assembly of the gradient matrix 
  grad_loc(1,1) = phicsi_csi;  grad_loc(1,2) = phieta_csi;  grad_loc(1,3) = phizeta_csi   
  grad_loc(2,1) = phicsi_eta;  grad_loc(2,2) = phieta_eta;  grad_loc(2,3) = phizeta_eta  
  grad_loc(3,1) = phicsi_zeta; grad_loc(3,2) = phieta_zeta; grad_loc(3,3) = phizeta_zeta
  
  !> transformation of the gradient matrix from local to global coordinates 
  grad = matmul(transpose(rot), matmul(grad_loc, rot)) 

end subroutine gradient_calc_sou_surfpan

!----------------------------------------------------------------------

!> Build a row of the linear system for a surface panel
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static
!! subroutine.
subroutine build_row_surfpan(this, elems, linsys, ie, ista, iend)
  class(t_surfpan), intent(inout) :: this
  type(t_impl_elem_p), intent(in)      :: elems(:)
  type(t_linsys), intent(inout)   :: linsys
  integer, intent(in)             :: ie
  integer, intent(in)             :: ista, iend

  integer :: j1 , ipres
  real(wp) :: b1
  real(wp) :: wind(3)

  ! RHS for \phi equation --------------------------
  linsys%b(ie) = 0.0_wp

  !Components not moving, no body velocity in the boundary condition
  !linsys%b(ie) = sum(linsys%b_static(:,ie) * (-uinf))
  ! RHS for Bernoulli polynomial equation ----------
  ipres = linsys%idSurfPanG2L(ie)
  linsys%b_pres(ipres) = 0.0_wp

  ! Static part ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! + \phi equation ------------------------------
  do j1 = 1,ista-1
    wind = variable_wind(elems(j1)%p%cen, sim_param%time)
    linsys%b(ie) = linsys%b(ie) + &
          linsys%b_static(ie,j1) *sum(elems(j1)%p%nor*(-wind-elems(j1)%p%uvort))
  enddo

  ! + Bernoulli polynomial equation --------------
  do j1 = 1, min(ista-1 , size(linsys%b_static_pres,1))
      linsys%b_pres( ipres ) = linsys%b_pres( ipres ) + &
                               linsys%b_static_pres( ipres , j1 ) * &
                        elems(linsys%idSurfPan(j1))%p%get_bernoulli_source()
  end do

  ! Moving part ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ista and iend will be the end of the unknowns vector, containing
  ! the moving elements
  do j1 = ista , iend

    !Compute the AIC and the contribution of j1 to the rhs of ie
    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1, &
                                  this%cen, ie, j1 )

    ! + \phi equation ----------------------------
    !Add the contribution to the rhs with the
    wind = variable_wind(elems(j1)%p%cen, sim_param%time)
    linsys%b(ie) = linsys%b(ie) &
              + b1* sum(elems(j1)%p%nor*(elems(j1)%p%ub-wind-elems(j1)%p%uvort))

    ! + Bernoulli polynomial equation ------------
    !Add the contribution to the rhs

    select type( el => elems(j1)%p ) ; class is(t_surfpan)
      linsys%b_pres( ipres ) =  linsys%b_pres( ipres )    + &
                                b1* el%bernoulli_source

    end select

  end do

end subroutine build_row_surfpan

!----------------------------------------------------------------------

!> Build a static row of the linear system for a surface panel
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_surfpan(this, elems, expl_elems, linsys, &
                                    ie, ista, iend)
  class(t_surfpan), intent(inout) :: this
  type(t_impl_elem_p), intent(in)      :: elems(:)
  type(t_expl_elem_p), intent(in)      :: expl_elems(:)
  type(t_linsys), intent(inout)   :: linsys
  integer, intent(in)             :: ie
  integer, intent(in)             :: ista, iend

  integer :: j1
  real(wp) :: b1

  linsys%b(ie) = 0.0_wp

  !Cycle just all the static elements, ista and iend will be the beginning of
  !the result vector. Then save the rhs in b_static
  do j1 = ista , iend

    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1,  &
                                  this%cen, ie, j1 )

    linsys%b_static(ie,j1)  = b1

  end do

  !Now build the static contribution from the explicit elements
  do j1 = 1,linsys%nstatic_expl
    call expl_elems(j1)%p%compute_pot( linsys%L_static(ie,j1), b1,  &
                                  this%cen, 1, 2 )
  enddo
  !The rest of the dynamic part will be completed during the first
  ! iteration of the assembling

end subroutine build_row_static_surfpan

!----------------------------------------------------------------------

!> Add the contribution of the wake to one equation for a surface panel
!!
!! The rhs of the equation for a surface panel is updated  adding the
!! the contribution of potential due to the wake
subroutine add_wake_surfpan(this, wake_elems, impl_wake_ind, linsys, &
                            ie, ista, iend)
  class(t_surfpan), intent(inout) :: this
  type(t_pot_elem_p), intent(in)  :: wake_elems(:)
  integer, intent(in)             :: impl_wake_ind(:,:)
  type(t_linsys), intent(inout)   :: linsys
  integer, intent(in)             :: ie
  integer, intent(in)             :: ista
  integer, intent(in)             :: iend

  integer :: j1, ind1, ind2
  real(wp) :: a, b
  integer :: n_impl
  real (wp) :: TR, TL 

  !Count the number of implicit wake contributions
  n_impl = size(impl_wake_ind,2)

  !Add the contribution of the implicit wake panels to the linear system
  !Implicitly we assume that the first set of wake panels are the implicit
  !ones since are at the beginning of the list

  if (sim_param%kutta_correction) then
    do j1 = 1 , n_impl
      ind1 = impl_wake_ind(1,j1)
      ind2 = impl_wake_ind(2,j1)
      if ((ind1.ge.ista .and. ind1.le.iend) .and. &
          (ind2.ge.ista .and. ind2.le.iend)) then
        !todo: find a more elegant solution to avoid i=j
        call wake_elems(j1)%p%compute_linear_pot(TL, TR, this%cen, 1, 2 ) 
        linsys%TL(ie, j1) = TL 
        linsys%TR(ie, j1) = TR  
        linsys%A(ie,ind1) = linsys%A(ie,ind1) + TL
        linsys%A(ie,ind2) = linsys%A(ie,ind2) - TL
      endif
    end do
  else !> old formulation (with constant potential)
    do j1 = 1 , n_impl
      ind1 = impl_wake_ind(1,j1)
      ind2 = impl_wake_ind(2,j1)
      if ((ind1.ge.ista .and. ind1.le.iend) .and. &
          (ind2.ge.ista .and. ind2.le.iend)) then
        !todo: find a more elegant solution to avoid i=j
        call wake_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
        linsys%A(ie,ind1) = linsys%A(ie,ind1) + a
        linsys%A(ie,ind2) = linsys%A(ie,ind2) - a
      endif
    end do
  endif

  ! Add the explicit vortex panel wake contribution to the rhs
  do j1 = n_impl+1 , size(wake_elems)
    !todo: find a more elegant solution to avoid i=j
    call wake_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
    linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%mag
  end do

end subroutine add_wake_surfpan

!----------------------------------------------------------------------

!> Add the contribution of the lifing lines to one equation for a surface panel
!!
!! The rhs of the equation for a surface panel is updated  adding the
!! the contribution of potential due to the lifting lines
subroutine add_expl_surfpan(this, expl_elems, linsys, &
                            ie,ista, iend)
  class(t_surfpan), intent(inout) :: this
  type(t_expl_elem_p), intent(in) :: expl_elems(:)
  type(t_linsys), intent(inout)   :: linsys
  integer, intent(in)             :: ie
  integer, intent(in)             :: ista
  integer, intent(in)             :: iend

  integer                         :: j1
  real(wp)                        :: a, b

  !Static part: take what was already computed
  do  j1 = 1, ista-1
    linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*expl_elems(j1)%p%mag
  enddo

  !Dynamic part: compute the things now
  do j1 = ista , iend
    !todo: find a more elegant solution to avoid i=j
    call expl_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
    linsys%b(ie) = linsys%b(ie) - a*expl_elems(j1)%p%mag
  end do

end subroutine add_expl_surfpan

!----------------------------------------------------------------------

!> Compute the potential due to a surface panel
!!
!! this subroutine employs doublets and sources basic subroutines to calculate
!! the AIC of a suface panel on another surface panel, and the contribution
!! to its rhs
subroutine compute_pot_surfpan(this, A, b, pos , i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out)           :: A
  real(wp), intent(out)           :: b
  real(wp), intent(in)            :: pos(:)
  integer , intent(in)            :: i, j

  real(wp)                        :: dou, sou

  if ( i .ne. j ) then
    call potential_calc_doublet(this, dou, pos)
  else
  ! AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

  A = -dou

  call potential_calc_sou_surfpan(this, sou, dou, pos)

  b =  sou

end subroutine compute_pot_surfpan 

!> Compute the potential due to a surface panel
!!
!! this subroutine employs doublets and sources basic subroutines to calculate
!! the AIC of a suface panel on another surface panel, and the contribution
!! to its rhs
subroutine compute_linear_pot_surfpan(this, TL, TR, pos , i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out)           :: TL
  real(wp), intent(out)           :: TR
  real(wp), intent(in)            :: pos(:)
  integer , intent(in)            :: i, j

  real(wp)                        :: dou
  
  if ( i .ne. j ) then
    call linear_potential_calc_doublet(this, TL, TR, pos)
  else
  ! AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

end subroutine compute_linear_pot_surfpan 

!----------------------------------------------------------------------

!> Compute the velocity due to a surface panel
!!
!! This subroutine employs doublets and sources basic subroutines to calculate
!! the AIC coefficients of a surface panel to a vortex ring and the
!! contribution to its rhs
subroutine compute_psi_surfpan(this, A, b, pos, nor, i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out)           :: A
  real(wp), intent(out)           :: b
  real(wp), intent(in)            :: pos(:)
  real(wp), intent(in)            :: nor(:)
  integer , intent(in)            :: i, j

  real(wp)                        :: vdou(3), vsou(3)

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)

  call velocity_calc_sou_surfpan(this, vsou, pos)

  ! b = (sources from doublets)
  b =   sum(-vsou * nor )

end subroutine compute_psi_surfpan

!----------------------------------------------------------------------

!> Compute the velocity induced by a surface panel in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! both the doublets and the sources
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_surfpan(this, pos, vel )
  class(t_surfpan), intent(in)  :: this
  real(wp), intent(in)          :: pos(:)
  real(wp), intent(out)         :: vel(3)

  real(wp)                      :: vdou(3), vsou(3), wind(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  ! source ----
  call velocity_calc_sou_surfpan(this, vsou, pos)
  
  wind = variable_wind(this%cen, sim_param%time)
  vel = vdou*this%mag - vsou*(sum(this%nor*(this%ub-wind-this%uvort)))

end subroutine compute_vel_surfpan

!----------------------------------------------------------------------

!> Compute the gradient of the induced velocity in the position pos
!!
!! The velocity in the position is calculated considering the influece of
!! both the doublets and the sources
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_grad_surfpan(this, pos, grad )
  class(t_surfpan), intent(in) :: this
  real(wp), intent(in)         :: pos(:)
  real(wp), intent(out)        :: grad(3,3)

  real(wp)                     :: grad_dou(3,3), grad_sou(3,3)
  real(wp)                     :: wind(3)

  ! doublet ---
  call gradient_calc_doublet(this, grad_dou, pos)

  ! source ----
  call gradient_calc_sou_surfpan(this, grad_sou, pos)

  wind = variable_wind(this%cen,sim_param%time)
  grad = grad_dou*this%mag &
        - grad_sou*( sum(this%nor*(this%ub-wind-this%uvort)) )

end subroutine compute_grad_surfpan

!----------------------------------------------------------------------

!> Compute an approximate value of the mean pressure on the actual element
! TODO: move the velocity update outside this routine, in a dedicated routine
subroutine compute_pres_surfpan(this, R_g)
  class(t_surfpan) , intent(inout) :: this
  real(wp)         , intent(in)    :: R_g(3,3)

  real(wp) :: vel_phi(3)
  real(wp) :: f(5)    ! <- max n_ver of a surfpan = 4 ; +1 for the constraint eqn

  integer :: i_e , n_neigh
  real(wp) :: mach
  real(wp) :: wind(3)

! This routine contains the velocity update as well. TODO, move to a dedicated routine
! Method implemented for surface velocity computation:
! stencil relying on a CHTLS* method
! *CHTLS: Constrained Hermite Taylor series Least Square method Constrained Hermite TLS for Mesh-free Derivative
! Estimation Near and On Boundaries  https://digitalcommons.calpoly.edu/cgi/viewcontent.cgi?article=1095&context=aero_fac

  n_neigh = 0 ; f = 0.0_wp
  do i_e = 1 , this%n_ver
    if ( associated(this%neigh(i_e)%p) ) then
      n_neigh = n_neigh + 1
      f(n_neigh) = - ( this%neigh(i_e)%p%mag - this%mag )
    end if
  end do
  wind = variable_wind(this%cen, sim_param%time)
  f(n_neigh+1) = sum(this%nor * (-wind - this%uvort + this%ub) )

  vel_phi = matmul( this%chtls_stencil , f(1:n_neigh+1) )

  vel_phi = matmul( (R_g) , vel_phi )
  
  this%surf_vel = wind + vel_phi + this%uvort

  !> pressure -------------------------------------------------
  ! unsteady problems  : P =   P_inf 
  !                          + 0.5*rho_inf*V_inf^2
  !                          +(- 0.5*rho_inf*V^2 
  !                          + rho * ub.u_phi
  !                          - rho_inf*dphi/dt)      
  !                                     ^   
  !!        with: idou = -phi         __|
  ! Prandt -- Glauert correction for compressibility effect
  ! it applies to the dynamic pressure terms only
  ! so that only the *relative* static pressure is affected
  mach = abs(norm2(wind - this%ub) / sim_param%a_inf)
  
  this%pres =   sim_param%P_inf &
              +(0.5_wp * sim_param%rho_inf * norm2(sim_param%u_inf)**2.0_wp &
              - 0.5_wp * sim_param%rho_inf * norm2(  this%surf_vel)**2.0_wp  &
              + sim_param%rho_inf * sum(this%ub*(vel_phi+this%uvort)) &
              + sim_param%rho_inf * this%didou_dt)/sqrt(1 - mach**2)

  ! compute dforce
  call compute_dforce_surfpan(this)

end subroutine compute_pres_surfpan

!----------------------------------------------------------------------

!> Compute the elementary force on the on the actual element
!!  computing the eleent force from the pressure
subroutine compute_dforce_surfpan(this)
  class(t_surfpan), intent(inout) :: this

  this%dforce = - (this%pres - sim_param%P_inf) * this%area * this%nor

end subroutine compute_dforce_surfpan

!----------------------------------------------------------------------

!> Compute the coefficients (pot_vel_stencil, for surfpan elements) for
!!  computing the velocity from the velocity potential (phi = -mu).
!! On-body analysis for 3dPanels (surfpan). This coefficients are constant
!!  in the local frame, associated with the component. In order to obtain
!!  the components of the velcoity in the base frame, the global rotation
!!  matrix is needed.
subroutine create_local_velocity_stencil_surfpan ( this , R_g )
  class(t_surfpan)  , intent(inout) :: this
  real(wp)          , intent(in)    :: R_g(3,3)

  real(wp) :: bubble_surf
  integer  :: i_v

  if ( .not. allocated(this%pot_vel_stencil) ) then
    allocate(this%pot_vel_stencil(3,this%n_ver) )
  end if

  ! method #2: w/o averaging on neighbouring elements ; 0.5 factor added
  bubble_surf = this%area

  do i_v = 1 , this%n_ver

    this%pot_vel_stencil(:,i_v) = &
     0.5_wp * cross( this%edge_vec(:,i_v) , this%nor )

  end do

  this%pot_vel_stencil = this%pot_vel_stencil / bubble_surf

  ! chtls stencil need to be defined in the local ref.sys.
  ! it is built in the global ref.sys. at the first timestep
  !  -> rotation is needed
  this%pot_vel_stencil = matmul( transpose(R_g) , this%pot_vel_stencil )

end subroutine create_local_velocity_stencil_surfpan

!----------------------------------------------------------------------
!> create stencil in the "local" reference system. The components of the
!  velocity vector in the "global" ref.sys. are obtained with the rotation
!  matrix of the "local" reference system
subroutine create_chtls_stencil_surfpan( this , R_g )
  class(t_surfpan)  , intent(inout) :: this
  real(wp)          , intent(in)    :: R_g(3,3)

  real(wp), allocatable             :: A(:,:), B(:,:), W(:,:), V(:,:)
  real(wp)                          :: dx(3)
  real(wp), allocatable             :: C(:,:), CQ(:,:)
  real(wp), allocatable             :: Cls_tilde(:,:), iCls_tilde(:,:), chtls_tmp(:,:)
  real(wp)                          :: det_cls
  real(wp)                          :: r1
  real(wp), allocatable             :: Q(:,:), R(:,:)

  integer                           :: n_neigh
  integer                           :: i_n , i_nn

  ! # of neighbouring elements
  n_neigh = 0
  do i_n = 1 , this%n_ver
    if ( associated(this%neigh(i_n)%p) ) then ! the neighbouring element exists
      n_neigh = n_neigh + 1
    end if
  end do

  ! arrays A (differences), B (constraints), W (weights) -----------
  allocate(A(n_neigh,3) , W(n_neigh+1,n_neigh+1)); W = 0.0_wp
  allocate(B(1,3));                                B(1,:) = this%nor
  
  i_n = 0
  do i_nn = 1 , this%n_ver
    if ( associated(this%neigh(i_nn)%p) ) then ! the neighbouring element exists
      i_n = i_n + 1
      dx = this%neigh(i_nn)%p%cen - this%cen
      A(i_n, : ) = dx
      W(i_n,i_n) = 1.0_wp/norm2(dx) * &
                    max(0.0_wp, abs(sum(this%nor*this%neigh(i_nn)%p%nor)))
    end if
  end do
  W(n_neigh+1,n_neigh+1) = sum( W ) / real(n_neigh,wp)

  allocate( V(3,1) ) ; V(:,1) = B(1,:)

  call compute_qr( V , Q , R )

  allocate(         C(n_neigh+1,3) ) ; C(1:n_neigh,:) = A ; C(n_neigh+1,:) = B(1,:)
  allocate(        CQ(n_neigh+1,3) ) ; CQ = matmul( C , Q )
  allocate( Cls_tilde(        2,2) ) ; allocate( iCls_tilde(2,2) )
  Cls_tilde = matmul( transpose(CQ(:,2:3)) , matmul( W , CQ(:,2:3) ) )


  ! inverse Cls ----
  det_cls = Cls_tilde(1,1) * Cls_tilde(2,2) - Cls_tilde(1,2) * Cls_tilde(2,1)
  if (det_cls .eq. 0.0_wp) then
    iCls_tilde = 0.0_wp
  else
    iCls_tilde(1,1) =  Cls_tilde(2,2) / det_cls ; iCls_tilde(1,2) = -Cls_tilde(1,2) / det_cls
    iCls_tilde(2,1) = -Cls_tilde(2,1) / det_cls ; iCls_tilde(2,2) =  Cls_tilde(1,1) / det_cls
  endif 
  if ( .not. allocated(this%chtls_stencil) ) then
    allocate(this%chtls_stencil( 3 , n_neigh + 1 ) ) ; this%chtls_stencil = 0.0_wp
  end if

  r1 = R(1,1)
  allocate(chtls_tmp( 3 , n_neigh + 1 )) ; chtls_tmp = 0.0_wp
  chtls_tmp(  1,n_neigh + 1 ) = 1.0_wp / R(1,1)
  chtls_tmp(2:3,1 : n_neigh ) = matmul( iCls_tilde , &
        matmul(transpose(CQ(1:n_neigh,2:3)), W(1:n_neigh,1:n_neigh) ) )
  chtls_tmp(2:3, n_neigh+1:n_neigh+1 ) = matmul( iCls_tilde , &
        matmul(transpose(CQ( n_neigh+1:n_neigh+1,2:3)), W(n_neigh+1:n_neigh+1, n_neigh+1:n_neigh+1) ) &
      - matmul( matmul( transpose(CQ(1:n_neigh+1,2:3)), W(1:n_neigh+1,1:n_neigh+1) ) , &
                CQ(1:n_neigh+1,1:1) / r1 )  )

  this%chtls_stencil = matmul( Q , chtls_tmp )

  ! chtls stencil need to be defined in the local ref.sys.
  ! it is built in the global ref.sys. at the first timestep
  !  -> rotation is needed
  this%chtls_stencil = matmul( transpose(R_g) , this%chtls_stencil )


  deallocate(C,CQ,Cls_tilde,iCls_tilde,chtls_tmp)
  deallocate(A,B,V,W,Q,R)


end subroutine create_chtls_stencil_surfpan



!----------------------------------------------------------------------

!> Calculate the geometrical quantities of a surface panel
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! surface panel
subroutine calc_geo_data_surfpan(this,vert)
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(in) :: vert(:,:)

  integer :: nsides, is
  real(wp):: nor(3), tanl(3)

  this%ver = vert
  nsides = this%n_ver

  ! center
  this%cen =  sum ( this%ver,2 ) / real(nsides,wp)

  ! unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross(this%ver(:,3) - this%ver(:,1) , &
                this%ver(:,4) - this%ver(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross(this%ver(:,3) - this%ver(:,2) , &
                this%ver(:,1) - this%ver(:,2)     )
  end if

  this%area = 0.5_wp * norm2(nor)
  this%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

  this%tang(:,1) = tanl / norm2(tanl)
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )

  ! vector connecting two consecutive vertices:
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
    this%edge_uni(:,is) = this%edge_vec(:,is) / this%edge_len(is)
  end do

  !surface panels own fields
  do is = 1 , nsides
    ! cosTi , sinTi
    this%cosTi(is) = sum( this%edge_uni(:,is) * this%tang(:,1) )
    this%sinTi(is) = sum( this%edge_uni(:,is) * this%tang(:,2) )
    ! projection of the vertices on the mean plane
    this%verp(:,is) = this%ver(:,is) - this%nor * &
                    sum( (this%ver(:,is) - this%cen ) * this%nor )
  end do

  !TODO: is it necessary to initialize it here?
  this%dforce = 0.0_wp
  this%dmom   = 0.0_wp

end subroutine calc_geo_data_surfpan

!----------------------------------------------------------------------

!> Calculate the vorticity induced velocity from vortical elements
subroutine get_vort_vel_surfpan(this, vort_elems)
  class(t_surfpan), intent(inout)   :: this
  type(t_vort_elem_p), intent(in)   :: vort_elems(:)

  integer :: iv
  real(wp) :: vel(3)

  !> Initialize to zero, before accumulation
  this%uvort = 0.0_wp

  do iv=1,size(vort_elems)
    call vort_elems(iv)%p%compute_vel(this%cen, vel)
    this%uvort = this%uvort + vel*one_4pi
  enddo

end subroutine

!----------------------------------------------------------------------

function get_bernoulli_source_surfpan(this) result(source)
  class(t_surfpan), intent(inout) :: this
  real(wp) :: source

  source = this%bernoulli_source

end function

!----------------------------------------------------------------------
end module mod_surfpan
