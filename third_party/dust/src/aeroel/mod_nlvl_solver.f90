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
!!          Alberto Savino
!!          Alessandro Cocco
!!=========================================================================

!> Module collecting the nonlinear VL correction loop logic
module mod_nlvl_solver

  use mod_param, only: &
    wp, pi, max_char_len, extended_char_len

  use mod_sim_param, only: &
    sim_param

  use mod_handling, only: &
    warning, printout, dust_time, t_realtime

  use mod_geometry, only: &
    t_geo

  use mod_wake, only: &
    t_wake

  use mod_aeroel, only: &
    t_impl_elem_p, t_pot_elem_p

  use mod_linsys_vars, only: &
    t_linsys

  use mod_linsys, only: &
    solve_linsys, dump_linsys

  use mod_c81, only: &
    t_aero_tab

  use mod_vortlatt, only: &
    t_vortlatt

  use mod_math, only: &
    dot

  implicit none

  private
  public :: solve_nlvl_correction

contains

  subroutine solve_nlvl_correction(geo, wake, elems, elems_non_corr, elems_corr, linsys, airfoil_data, &
                                   res_old, sel, it, time_2_debug_out, basename_debug)
    type(t_geo),        intent(inout) :: geo
    type(t_wake),       intent(inout) :: wake
    type(t_impl_elem_p),intent(inout) :: elems(:)
    type(t_pot_elem_p), intent(inout) :: elems_non_corr(:)
    type(t_pot_elem_p), intent(inout) :: elems_corr(:)
    type(t_linsys),     intent(inout) :: linsys
    type(t_aero_tab),   intent(in)    :: airfoil_data(:)
    real(wp),           intent(in)    :: res_old(:)
    integer,            intent(in)    :: sel
    integer,            intent(in)    :: it
    logical,            intent(in)    :: time_2_debug_out
    character(len=*),   intent(in)    :: basename_debug

    integer                         :: i_el, i_c, i_s, i_p, i_c2, i_s2
    integer                         :: it_vl, it_stall
    real(wp)                        :: tol, diff, max_diff
    real(wp)                        :: d_cd(3), vel(3), v(3), a_v, area_stripe, dforce_stripe(3)
    real(wp)                        :: nor(3), tang_cen(3), u_v, q_inf
    real(wp), allocatable           :: residual_vl(:), residual_vl_old(:), residual_vl_delta(:), gamma_tmp(:)
    real(wp)                        :: rel_aitken, den_aitken
    real(t_realtime)                :: t0, t1
    character(len=max_char_len)     :: frmt, frmt_vl
    character(len=extended_char_len):: message

    if (.not. sim_param%vl_correction) return

    tol = sim_param%vl_tol
    it_vl = 0
    it_stall = 0
    max_diff = tol + 1e-6_wp
    linsys%skip = .true.
    t0 = dust_time()

    !> Select time step to start the VL correction
    !  (avoid strange behaviour at the begining of simulation)
    if (it .gt. sim_param%vl_startstep) then

      !> allocate residual terms for Aitken acceleration
      allocate(residual_vl(size(linsys%b)));        residual_vl = 0.0_wp
      allocate(residual_vl_old(size(linsys%b)));    residual_vl_old = 0.0_wp
      allocate(residual_vl_delta(size(linsys%b)));  residual_vl_delta = 0.0_wp
      allocate(gamma_tmp(size(linsys%b)));          gamma_tmp = 0.0_wp
      rel_aitken = sim_param%vl_relax

      do while (max_diff .gt. tol .and. it_vl .lt. sim_param%vl_maxiter)

        max_diff = 0.0_wp
        residual_vl = 0.0_wp

        !> geo quantities must be updated for all components before computing anything
        do i_c = 1, size(geo%components)
          if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
              trim(geo%components(i_c)%aero_correction) .eq. 'true') then
            if (it_vl .eq. 0) then
              !$omp parallel do private(i_s)
              do i_s = 1, size(geo%components(i_c)%stripe)
                call geo%components(i_c)%stripe(i_s)%calc_geo_data(geo%components(i_c)%stripe(i_s)%ver)
              end do
              !$omp end parallel do
            endif
          endif
        end do

        do i_c = 1, size(geo%components)
          if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
              trim(geo%components(i_c)%aero_correction) .eq. 'true') then

            !> Freeze external wake-induced velocity during the inner nonlinear loop.
            if (it_vl .eq. 0) then
              !$omp parallel do private(i_s)
              do i_s = 1, size(geo%components(i_c)%stripe)
                call geo%components(i_c)%stripe(i_s)%get_vel_ctr_pt(elems_non_corr, (/ wake%pan_p, wake%rin_p /), wake%vort_p)
              end do
              !$omp end parallel do
            endif

            !$omp parallel do private(i_s, i_s2, i_c2, vel, v, diff) schedule(dynamic, 4) reduction(max:max_diff)
            do i_s = 1, size(geo%components(i_c)%stripe)
              vel = 0.0_wp
              do i_c2 = 1, size(geo%components)
                if (trim(geo%components(i_c2)%comp_el_type) .eq. 'v' .and. &
                    trim(geo%components(i_c2)%aero_correction) .eq. 'true') then
                  do i_s2 = 1, size(geo%components(i_c2)%stripe)
                    call geo%components(i_c2)%stripe(i_s2)%compute_vel_stripe(geo%components(i_c)%stripe(i_s)%cen, v)
                    vel = vel + v
                  end do
                endif
              end do

              geo%components(i_c)%stripe(i_s)%vel = vel
              call geo%components(i_c)%stripe(i_s)%correction_c81_vortlatt(airfoil_data, linsys, diff, residual_vl, it_vl, i_s)
              max_diff = max(diff, max_diff)
            end do
            !$omp end parallel do
          end if
        end do

        !> Debug output of the system
        if ((sim_param%debug_level .ge. 50) .and. time_2_debug_out) then
          write(frmt,'(I4.4)') it
          write(frmt_vl,'(I4.4)') it_vl
          call dump_linsys(linsys,  &
                           trim(basename_debug)//'A_'//trim(frmt)//'_it_'//trim(frmt_vl)//'.dat', &
                           trim(basename_debug)//'b_'//trim(frmt)//'_it_'//trim(frmt_vl)//'.dat')
        endif

        residual_vl_delta = residual_vl - residual_vl_old
        den_aitken = dot(residual_vl_delta, residual_vl_delta)

        if (sim_param%rel_aitken .and. it_vl .gt. 2 .and. den_aitken .gt. epsilon(1.0_wp)) then
          rel_aitken = -rel_aitken * dot(residual_vl_old, residual_vl_delta) / den_aitken
        else
          rel_aitken = sim_param%vl_relax
        endif

        linsys%b = linsys%b + rel_aitken * residual_vl
        residual_vl_old = residual_vl

        !> Solve the factorized system
        if (linsys%rank .gt. 0) then
          call solve_linsys(linsys)
        endif

        it_vl = it_vl + 1

        !> average intensity for stall condition
        if (sim_param%vl_ave) then
          if (it_vl .gt. sim_param%vl_maxiter - sim_param%vl_iter_ave) then
            it_stall = it_stall + 1
            gamma_tmp = gamma_tmp + linsys%res
          endif

          if (it_vl .eq. sim_param%vl_maxiter) then
            linsys%res = gamma_tmp / real(it_stall,wp)
          endif
        endif

        !> update unsteady term
        !$omp parallel do private(i_el)
        do i_el = 1 , sel
          elems(i_el)%p%didou_dt = (linsys%res(i_el) - res_old(i_el)) / sim_param%dt
        enddo
        !$omp end parallel do

        do i_el = 1, size(elems_non_corr)
          select type(el => elems_non_corr(i_el)%p)
            class is(t_vortlatt)
              !> compute dforce using AVL formula with prandtl glauert for non corrected vl
              call el%compute_dforce_jukowski(.true.)
          end select
        end do

        do i_el = 1, size(elems_corr)
          select type(el => elems_corr(i_el)%p)
            class is(t_vortlatt)
              !> compute dforce using AVL formula without prandtl glauert correction
              call el%compute_dforce_jukowski(.false.)
          end select
        end do

      end do

      deallocate(residual_vl, residual_vl_old, residual_vl_delta, gamma_tmp)

      if(it_vl .eq. sim_param%vl_maxiter) then
        call warning('dust','dust','max iteration reached for non linear vl:&
                    &increase VLmaxiter!')
        write(message,*) 'Last iteration error: ', max_diff
        call printout(message)
      endif

      !> Viscous and pressure drag correction
      do i_c = 1, size(geo%components)
        if (trim(geo%components(i_c)%comp_el_type) .eq. 'v' .and. &
            trim(geo%components(i_c)%aero_correction) .eq. 'true') then

          do i_s = 1, size(geo%components(i_c)%stripe)

            nor         = geo%components(i_c)%stripe(i_s)%nor
            tang_cen    = geo%components(i_c)%stripe(i_s)%tang_cen
            a_v         = geo%components(i_c)%stripe(i_s)%alpha*pi/180.0_wp
            area_stripe = geo%components(i_c)%stripe(i_s)%area
            u_v         = geo%components(i_c)%stripe(i_s)%vel_2d
            q_inf       = 0.5_wp*sim_param%rho_inf * u_v ** 2.0_wp * area_stripe

            d_cd = 0.5_wp * sim_param%rho_inf * u_v**2.0_wp * &
                   geo%components(i_c)%stripe(i_s)%cd * &
                   (tang_cen * cos(a_v) + nor * sin(a_v))

            dforce_stripe = 0.0_wp
            do i_p = 1, size(geo%components(i_c)%stripe(i_s)%panels)
              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce = &
                  geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce + &
                  d_cd * geo%components(i_c)%stripe(i_s)%panels(i_p)%p%area

              geo%components(i_c)%stripe(i_s)%panels(i_p)%p%pres = &
                  sum(geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce * &
                      geo%components(i_c)%stripe(i_s)%panels(i_p)%p%nor) / &
                  geo%components(i_c)%stripe(i_s)%panels(i_p)%p%area

              dforce_stripe = dforce_stripe + geo%components(i_c)%stripe(i_s)%panels(i_p)%p%dforce
            end do
          end do
        end if
      end do
    endif

    linsys%skip = .false.
    t1 = dust_time()
    if(sim_param%debug_level .ge. 1) then
      write(message,'(A,F9.3,A)') 'Solved nonlinear vortex lattice in: ' , t1 - t0,' s.'
      call printout(message)
    endif

  end subroutine solve_nlvl_correction

end module mod_nlvl_solver
