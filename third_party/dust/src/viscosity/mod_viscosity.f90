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
!!          Matteo Tugnoli
!!=========================================================================

module mod_viscosity

use mod_param, only: &
  wp , &
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi, max_char_len

use mod_math, only: &
  cross

use mod_sim_param, only: &
  sim_param

use mod_geometry, only: &
  t_geo , t_tedge

use mod_aeroel, only: &
  t_pot_elem_p , t_impl_elem_p

use mod_surfpan, only: &
  t_surfpan

!----------------------------------------------------------------------

implicit none

public :: viscosity_effects

private

character(len=*), parameter :: this_mod_name = 'mod_viscosity'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
! subroutine viscosity_effects( geo , sim_param )
! Viscous and separations effects for panel elements only
! ( A loop over components is needed, because elems contains both panels
!   and vortex lattice elements )
! 1. preliminary compuation of surface quantities:
!    - surf_vel (already computed in dust.f90, elems%p%compute_pres
!      around l.536;
!    - approximation of the wall bounded vorticity.
! 2. computation of the velocity normal to the wall to release vorticity
!    in the domain.
! 3. computation of the amount of vorticity to be released in the domain
! 4. create new particles with the associated velocity to be evolved NOT HERE
!    but instead where the whole wake is evolved
!
subroutine viscosity_effects( geo , elems , te )
 type(t_geo)         , intent(inout) , target :: geo
 type(t_impl_elem_p) , intent(inout)          :: elems(:)
 type(t_tedge)       , intent(in)             :: te

 real(wp) :: vc , vd
 integer :: i_comp , n_comp , i_elem , n_elem
 integer :: i_e , ie_te , ne_te
 integer :: ie1 , ie2

 real(wp) , parameter :: c = 1.136_wp   ! Ojima-Kamemoto(2000)
 real(wp) :: c2_2

 real(wp) :: OmV(3) , OmV_free(3) , OmV_bound(3)
 real(wp) :: al_free , al_bound
 real(wp) :: edge_te(3)

 ! preliminary "fixed" parameters
 real(wp) , parameter :: h = 0.025_wp          ! h = 0.100_wp
 real(wp) , parameter :: tol_velSep = 0.0_wp  ! 0.00_wp for the cylinder ! 0.10_wp , 0.05_wp for the airfoil  ! 0.5_wp
 real(wp) , parameter :: vc_ratio   = 0.50_wp  ! 1.00_wp for the cylinder ! 0.50_wp for the airfoil

! airfoil ------------------
!   h = 0.025_wp   ! 0.05_wp   ! 0.1_wp   ! 0.025_wp
!   tol_velSep = 0.20_wp
! cylinder, R = 1.0 --------
!   h = 0.200_wp
!   tol_velSep = 0.00_wp

!write(*,*)
!write(*,*) ' ---- Viscosity effects ----------- '
!write(*,*)


 ! some parameters -------
 c2_2 = 0.5_wp * c**2  ! coefficient c^2 needed for viscous diffusion (Ojima & Kamemoto)

 n_comp = size( geo%components )
! debug -----
!write(*,*) ' n. components : ' , n_comp
! debug -----

 do i_comp = 1 , n_comp     ! ***** loop #1 over components *****

! debug -----
!  write(*,*) ' component id. : ' , i_comp , &
!    ' , comp_name : ' , trim(geo%components(i_comp)%comp_name) , &
!    ' , comp_el_type: ' , trim(geo%components(i_comp)%comp_el_type)
! debug -----

   ! flow separation allowed only for surfpan elements -----
   if ( trim( geo%components(i_comp)%comp_el_type ) .eq. 'p' ) then

     n_elem = size( geo%components(i_comp)%el )
!    write(*,*) '   n_elem : ' , n_elem

     do i_elem = 1 , n_elem     ! ***** loop #2 over elements   *****

       ! Initialisation -------------------------------------------
       ! volume integral contribution of vorticity
       OmV = 0.0_wp ; OmV_free = 0.0_wp ; OmV_bound = 0.0_wp
       al_free = 0.0_wp

       ! normal "convective" velocity
       vc = 0.0_wp

       select type( el => geo%components(i_comp)%el(i_elem) ) ; type is (t_surfpan)

       ! height of the surface layer (boundary layer??)
       el%h_bl = h

       do i_e = 1 , el % n_ver     ! ***** loop #3 over neighbors  *****

         if ( associated(el%neigh(i_e)%p) ) then
           select type( el_neigh => el%neigh(i_e)%p  ) ; type is (t_surfpan)
           vc = vc - &
                el%h_bl * sum( &
                   cross( el % edge_vec(:,i_e) , el % nor ) * &
                 ( el % surf_vel + el_neigh%surf_vel ) ) * 0.5_wp
           OmV = OmV + &
                   0.5_wp * el % edge_vec(:,i_e) * ( &
                                            el % mag - el_neigh % mag )
           end select
         else
           vc = vc - &
                el%h_bl * sum( &
                   cross( el % edge_vec(:,i_e) , el % nor ) * &
                 ( el % surf_vel ) )
           OmV = OmV + &
                   el % edge_vec(:,i_e) * el % mag
         end if

       end do     ! ***** loop #3 over neighbors  *****

       vc = ( vc - sum( el%ub * el%nor ) ) / el%area

! debug -----
!      write(*,*) ' el. i_elem : ' , i_elem , ' ,   vc :' , vc
! debug -----

       ! normal "diffusive" velocity
       vd = c2_2 * sim_param%mu_inf / el%h_bl     ! 0.0_wp

       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! criterion for flow separation
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      if ( vc + vd .gt. tol_velSep ) then ! release free vorticity and
!                                          !   update the intensity of bounded vorticity
       if ( vc + vd .gt. tol_velSep * norm2(el%surf_vel) ) then ! release free vorticity and
                                                                !   update the intensity of bounded vorticity

         al_bound = el%h_bl / ( ( vc + vd ) * sim_param%dt + el%h_bl )
         al_free  = 1.0_wp - al_bound

       end if
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! criterion for flow separation
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! to be corrected for te elements
       el%al_free = al_free
       el%surf_vort = OmV
       el%free_vort = al_free * OmV

       end select ! <- type(el) = t_surfpan

     end do     ! ***** loop #2 over elements   *****

   end if

!  write(*,*)

 end do     ! ***** loop #1 over components *****

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<
 ! -------------------------
 ! te correction
 ! -------------------------
 ne_te = size(te%e,2)
 do ie_te = 1 , ne_te
   if(associated(te%e(2,ie_te)%p)) then
   ie1 = te%e(1,ie_te)%p%id
   ie2 = te%e(2,ie_te)%p%id

   select type( el1 => elems(ie1)%p ) ; type is (t_surfpan)
   select type( el2 => elems(ie2)%p ) ; type is (t_surfpan)

     edge_te = 0.5_wp * ( geo%points( :,te%i( 1,te%ii(2,ie_te) ) ) + &
                          geo%points( :,te%i( 2,te%ii(2,ie_te) ) ) - &
                          geo%points( :,te%i( 1,te%ii(1,ie_te) ) ) - &
                          geo%points( :,te%i( 2,te%ii(1,ie_te) ) )  )
     ! correction of panel ie1 ----
     el1%surf_vort = el1%surf_vort - el1%mag * edge_te
     el1%free_vort = el1%al_free * el1%surf_vort

     ! correction of panel ie2 ----
     el2%surf_vort = el2%surf_vort + el2%mag * edge_te
!        el2%mag * ( te%rr(:,te%ii(1,ie_te)) - te%rr(:,te%ii(1,ie_te)) )
     el2%free_vort = el2%al_free * el2%surf_vort

   end select
   end select
   endif

 end do
 ! -------------------------
 ! te correction
 ! -------------------------
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<

end subroutine viscosity_effects

!----------------------------------------------------------------------

end module mod_viscosity
