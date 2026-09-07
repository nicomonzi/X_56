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
!!          Federico Gentile
!!          Matteo Dall'Ora
!!          Alessandro Cocco
!!=========================================================================


!> Module to treat the whole wake
module mod_wake_test

use mod_param, only: &
  wp, nl, pi, one_4pi, max_char_len, eps

use mod_math, only: &
  cross, infinite_plate_spline, tessellate

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_aeroel, only: &
  c_elem, c_vort_elem, &
  t_elem_p, t_vort_elem_p

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_hdf5_io, only: &
  h5loc, &
  new_hdf5_file, &
  open_hdf5_file, &
  close_hdf5_file, &
  new_hdf5_group, &
  open_hdf5_group, &
  close_hdf5_group, &
  write_hdf5, &
  write_hdf5_attr, &
  read_hdf5, &
  read_hdf5_al, &
  check_dset_hdf5

use mod_octree_test, only: &
  t_octree, sort_particles, calculate_multipole, &
  apply_multipole, apply_multipole_optimized

use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: t_wake, initialize_wake, update_wake, &
          prepare_wake, complete_wake
private


!> Type containing wake panels information
type :: t_wake

  !! Particles data

  !> Maximum number of particles
  integer :: nmax_prt

  !> Actual number of particles
  integer :: n_prt

  !> Wake particles
  type(t_vortpart), allocatable :: wake_parts(:)

  !> Magnitude of particles vorticity
  real(wp), allocatable :: prt_ivort(:)

  !> Wake particles pointer
  type(t_vortpart_p), allocatable :: part_p(:)

  !> Bounding box
  real(wp) :: part_box_min(3), part_box_max(3)

  type(t_vort_elem_p), allocatable :: vort_p(:)

end type

!> Class to change methods from different wake implementations
type, abstract :: c_wake_mov
  contains
  procedure(i_get_vel), deferred, pass(this) :: get_vel
end type

abstract interface
  subroutine i_get_vel(this, wake, pos, vel)
    import                                :: c_wake_mov, wp, t_wake
    class(c_wake_mov)                     :: this
    type(t_wake), intent(in)              :: wake
    real(wp), intent(in)                  :: pos(3)
    real(wp), intent(out)                 :: vel(3)
  end subroutine
end interface

type, extends(c_wake_mov) :: t_free_wake
contains
  procedure, pass(this) :: get_vel => get_vel_free
end type

class(c_wake_mov), allocatable  :: wake_movement
character(len=max_char_len)     :: msg
real(t_realtime)                :: t1 , t0
character(len=*), parameter     :: this_mod_name='mod_wake_test'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the particle wake
subroutine initialize_wake(wake)
  type(t_wake), intent(out),target     :: wake
  integer                              :: ip


allocate(t_free_wake::wake_movement)

!Particles

wake%nmax_prt = sim_param%part_n0
allocate(wake%wake_parts(wake%nmax_prt))
allocate(wake%prt_ivort(wake%nmax_prt))

wake%n_prt = wake%nmax_prt
allocate(wake%part_p(wake%n_prt))


do ip = 1,wake%n_prt
  wake%wake_parts(ip)%mag      => wake%prt_ivort(ip)
  wake%wake_parts(ip)%cen      = sim_param%part_pos0(ip,:)
  wake%wake_parts(ip)%dir      = sim_param%part_vort0_dir(ip,:)
  wake%wake_parts(ip)%mag      = sim_param%part_vort0_mag(ip)
  wake%wake_parts(ip)%vol      = sim_param%part_vol(ip)
  wake%wake_parts(ip)%r_Vortex = sim_param%VortexRad
  wake%wake_parts(ip)%free     = .false.
  wake%part_p(ip)%p            => wake%wake_parts(ip)
enddo

wake%part_box_min = sim_param%particles_box_min
wake%part_box_max = sim_param%particles_box_max


allocate(wake%vort_p(wake%n_prt))
end subroutine initialize_wake

!----------------------------------------------------------------------


!----------------------------------------------------------------------

!> Prepare the wake before the timestep
!!
!! Mainly prepare all the structures for the octree
subroutine prepare_wake(wake, octree)
  type(t_wake), intent(inout), target   :: wake
  type(t_octree), intent(inout)         :: octree
  integer                               :: k, ip, ir, iw, ie, n_end_vort

  if (sim_param%use_fmm) then
    call sort_particles(wake%wake_parts, wake%n_prt, octree)
    call calculate_multipole(wake%part_p, octree)
  endif

  !==>Recreate structures and pointers, if particles are present
  if( wake%n_prt.gt.0 ) then

    !Recreate the pointer vector
    if(allocated(wake%part_p)) then 
      deallocate(wake%part_p)
    endif

    allocate(wake%part_p(wake%n_prt))

    
    k = 1
    do ip = 1, wake%n_prt
      do ir=k,wake%nmax_prt
        if(.not. wake%wake_parts(ir)%free) then
          k = ir+1
          wake%part_p(ip)%p => wake%wake_parts(ir)
          wake%vort_p(ip)%p => wake%wake_parts(ir)
          exit
        endif
      enddo
    enddo
  endif
end subroutine prepare_wake

!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels 
!  Brings them to the next time step
!  Only updates the "old" panels, ie not the first two, and existing particles
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll)
!! elements
subroutine update_wake(wake, octree)
  type(t_wake), intent(inout), target :: wake
  type(t_octree), intent(inout)       :: octree

  integer                             :: iw, ipan, ie, ip, np, iq
  integer                             :: id, ir
  real(wp)                            :: pos_p(3), vel_p(3)
  real(wp)                            :: str(3), stretch(3)
  real(wp)                            :: ru(3), rotu(3)
  real(wp)                            :: df(3), diff(3)
  real(wp)                            :: hcas_vel(3)
  real(wp), allocatable               :: point_old(:,:,:)
  real(wp), allocatable               :: points(:,:,:)
  logical                             :: increase_wake
  integer                             :: size_old
  character(len=*), parameter         :: this_sub_name='update_wake'

  !==>    Particles: evolve the position in time

  !calculate the velocities at the points
!$omp parallel do private(pos_p, vel_p, ip, iq,  stretch, diff, df, str, ru, rotu)
  do ip = 1, wake%n_prt
    wake%part_p(ip)%p%vel_old = wake%part_p(ip)%p%vel
    wake%part_p(ip)%p%stretch_old = wake%part_p(ip)%p%stretch
    wake%part_p(ip)%p%stretch = 0.0_wp
    wake%part_p(ip)%p%rotu = 0.0_wp

    !If not using the fast multipole, update particles position now
    if (.not.sim_param%use_fmm) then
      pos_p = wake%part_p(ip)%p%cen

      call wake_movement%get_vel(wake, pos_p, vel_p)

      wake%part_p(ip)%p%vel =  vel_p
      !if using vortex stretching, calculate it now
      if(sim_param%use_vs) then
        stretch = 0.0_wp
        rotu = 0.0_wp
        do iq = 1, wake%n_prt
          if (ip.ne.iq) then
            call wake%part_p(iq)%p%compute_stretch(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, str)
            ! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            stretch = stretch + str*one_4pi

            if(sim_param%use_divfilt) then
              call wake%part_p(iq)%p%compute_rotu(wake%part_p(ip)%p%cen, &
                    wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, wake%part_p(ip)%p%r_Vortex, ru)
              rotu = rotu + ru*one_4pi
            endif
            
          endif
        enddo
!        if (ip .eq. 1) then
!          write(*,*) 'mag = ', wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag
!          write(*,*) 'stretch = ', stretch
!        endif
        
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + stretch
        wake%part_p(ip)%p%stretch_alone = wake%part_p(ip)%p%stretch ! Used in reformulated

        if(sim_param%use_divfilt) then 
          wake%part_p(ip)%p%rotu = wake%part_p(ip)%p%rotu + rotu
        endif
      
      endif !use_vs

      !if using the vortex diffusion, calculate it now
      if(sim_param%use_vd) then
        diff = 0.0_wp

        do iq = 1, wake%n_prt

          if (ip .ne. iq) then
            call wake%part_p(iq)%p%compute_diffusion(wake%part_p(ip)%p%cen, &
                  wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, &
                  wake%part_p(ip)%p%r_Vortex, wake%part_p(ip)%p%vol, df)
            if ( (norm2(wake%part_p(ip)%p%mag*wake%part_p(ip)%p%dir + 2.0_wp*df*sim_param%nu_inf*sim_param%dt)) .le. &
                1.05_wp*max( wake%part_p(ip)%p%mag, wake%part_p(iq)%p%mag)) then
                diff = diff + 2.0_wp*df*sim_param%nu_inf    ! 21/12/2023 Added factor 2 (see Winckelmans)
            else
                write(*,*) 'Unphysical diffusion'
            endif
          endif

        enddo !iq
        wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch + diff
!        if (ip .eq. 1) then
!          write(*,*) 'diff = ', diff
!        endif
      endif !use_vd

    end if ! not use_fmm

  enddo
!$omp end parallel do

  if (sim_param%use_fmm) then
    t0 = dust_time()
    if (sim_param%use_vs .and. sim_param%use_divfilt .and. sim_param%use_vd) then 
      call apply_multipole_optimized(wake%part_p, octree)
    else 
      call apply_multipole(wake%part_p, octree)
    endif 
    t1 = dust_time()
    write(msg,'(A,F9.3,A)') 'Multipoles calculation: ' , t1 - t0,' s.'
    if(sim_param%debug_level.ge.3) call printout(msg)
    write(msg,'(A,I0)') 'Number of particles: ' , wake%n_prt
    if(sim_param%debug_level.ge.1) call printout(msg)
  endif


end subroutine update_wake

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!! Completes the updating to the next time step begun in update_wake
!! first and second row are updated to the next step and new particles are
!! created if necessary; they will appear at the save_date in the next time step
subroutine complete_wake(wake, octree)
  type(t_wake), target, intent(inout)   :: wake
  type(t_octree), intent(inout)         :: octree
  integer                               :: p1, p2
  integer                               :: ip, iw, id, is, nprev
  real(wp)                              :: dist(3) , vel_te(3), pos_p(3)
  real(wp)                              :: dir(3), partvec(3), ave, alpha_p(3), alpha_p_n
  real(wp)                              :: q_1(3), q_2(3), q_3(3)
  real(wp)                              :: alpha_q_1(3), alpha_q_2(3), alpha_q_3(3)
  real(wp)                              :: alpha_p_1(3), alpha_p_2(3), alpha_p_3(3)
  real(wp)                              :: alpha_p_1_mag, alpha_p_2_mag, alpha_p_3_mag
  real(wp)                              :: alpha_p_1_dir(3), alpha_p_2_dir(3), alpha_p_3_dir(3)
  real(wp)                              :: r_Vortex_q_1, r_Vortex_q_2, r_Vortex_q_3 
  real(wp)                              :: r_Vortex_p_1, r_Vortex_p_2, r_Vortex_p_3, r_Vortex 
  integer                               :: n_part, count_free
  real(wp)                              :: vel_in(3), vel_out(3), wind(3), filt_eta
  real(wp)                              :: sigma_dot
  real(wp), allocatable                 :: alpha_pedrizzetti(:,:)

  character(len=max_char_len)           :: msg
  character(len=*), parameter           :: this_sub_name='complete_wake'

!==> Particles: update the position and intensity in time, avoid penetration
!               and chech if remain into the boundaries
  n_part = wake%n_prt
  filt_eta = sim_param%alpha_divfilt/sim_param%dt
select case (sim_param%integrator)
  case('euler') ! Explicit Euler
!$omp parallel do schedule(dynamic,4) private(ip,pos_p,alpha_p,alpha_p_n,vel_in,vel_out, sigma_dot, r_Vortex)
  do ip = 1, n_part

    if(.not. wake%part_p(ip)%p%free) then 
      if( wake%part_p(ip)%p%mag .ge. sim_param%mag_threshold) then ! to avoid negative magnitudes (and too small)

        pos_p = wake%part_p(ip)%p%cen + wake%part_p(ip)%p%vel* &
                sim_param%dt*real(sim_param%ndt_update_wake,wp)

        if(all(pos_p .ge. wake%part_box_min) .and. &
            all(pos_p .le. wake%part_box_max)) then
          wake%part_p(ip)%p%cen = pos_p

          if(sim_param%use_vs .or. sim_param%use_vd) then

            ! add reformulated contribution (Alvarez rVPM 2023)
            sigma_dot = 0.0_wp
            if(sim_param%use_reformulated) then
              wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch & 
                                          - (sim_param%g + sim_param%f)/(1.0_wp/3.0_wp + sim_param%f) & 
                                          * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir) * wake%part_p(ip)%p%dir

              sigma_dot = - (sim_param%g + sim_param%f)/(1.0_wp + 3.0_wp*sim_param%f) &
                          * wake%part_p(ip)%p%r_Vortex/wake%part_p(ip)%p%mag & 
                          * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir)
            endif

            !add divergence filtering (Pedrizzetti Relaxation)
            if(sim_param%use_divfilt .and. norm2(wake%part_p(ip)%p%rotu) .ge. 1.0e-9_wp) then
              wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch - &
                filt_eta/real(sim_param%ndt_update_wake,wp)*( wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag - &
                wake%part_p(ip)%p%rotu*wake%part_p(ip)%p%mag/norm2(wake%part_p(ip)%p%rotu))
            endif

            !Explicit Euler
            alpha_p = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag + &
                            wake%part_p(ip)%p%stretch* &
                            sim_param%dt*real(sim_param%ndt_update_wake,wp)
            alpha_p_n = norm2(alpha_p)

            if(sim_param%use_reformulated) then
              !r_Vortex update
              r_Vortex = wake%part_p(ip)%p%r_Vortex &
                        + sigma_dot * sim_param%dt*real(sim_param%ndt_update_wake,wp)
            else 
              r_Vortex = wake%part_p(ip)%p%r_Vortex
            endif
            
            !Magnitude check to avoid division by zero and negative magnitudes
            if(alpha_p_n .ge. sim_param%mag_threshold .and. r_Vortex .ge. sim_param%mag_threshold) then
              wake%part_p(ip)%p%dir = alpha_p/alpha_p_n
              wake%part_p(ip)%p%mag = alpha_p_n 
              if(sim_param%use_reformulated) then
                wake%part_p(ip)%p%r_Vortex = r_Vortex
              endif
            else
              wake%part_p(ip)%p%free = .true.
!$omp atomic update
              wake%n_prt = wake%n_prt -1
!$omp end atomic
            endif !magnitude check
          endif !use_vs 
        else !part_box
          wake%part_p(ip)%p%free = .true.
!$omp atomic update
          wake%n_prt = wake%n_prt -1
!$omp end atomic
        endif !part box  
      else !magnitude
        wake%part_p(ip)%p%free = .true.
!$omp atomic update
        wake%n_prt = wake%n_prt -1
!$omp end atomic       
      endif !magnitude
    endif ! not free
  enddo
!$omp end parallel do

  case('low_storage') ! Low storage Runge-Kutta 
    !> 1st stage
    count_free = 0
    if(sim_param%use_divfilt) then 
      allocate(alpha_pedrizzetti(n_part,3)); alpha_pedrizzetti = 0.0_wp
    endif
!$omp parallel do schedule(dynamic,4) private(ip, q_1, alpha_q_1, alpha_p_1, sigma_dot, r_Vortex_q_1, r_Vortex_p_1)
    do ip = 1, n_part
      if ( .not. wake%part_p(ip)%p%free) then
        if( wake%part_p(ip)%p%mag .ge. sim_param%mag_threshold) then ! to avoid negative magnitudes (and too small)

          if(sim_param%use_divfilt .and. norm2(wake%part_p(ip)%p%rotu) .ge. 1.0e-9_wp) then 
            ! In this case (low_storage) the contribution to alpha is calculated
            ! directly, instead of acting on stretch. Added to alpha_p_3 after 3rd stage
            alpha_pedrizzetti(ip,:) = (1.0_wp-sim_param%alpha_divfilt)*wake%part_p(ip)%p%mag*wake%part_p(ip)%p%dir &
                                    + sim_param%alpha_divfilt*wake%part_p(ip)%p%mag &
                                    * wake%part_p(ip)%p%rotu/norm2(wake%part_p(ip)%p%rotu) &
                                    - wake%part_p(ip)%p%mag*wake%part_p(ip)%p%dir ! Only the increment is considered
          endif

          q_1 = wake%part_p(ip)%p%vel*sim_param%dt*real(sim_param%ndt_update_wake,wp)
          wake%part_p(ip)%p%cen = wake%part_p(ip)%p%cen + 1.0_wp/3.0_wp*q_1 
          
          sigma_dot = 0.0_wp
          if(sim_param%use_reformulated) then
            wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch & 
                                        - (sim_param%g + sim_param%f)/(1.0_wp/3.0_wp + sim_param%f) & 
                                        * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir) * wake%part_p(ip)%p%dir
            sigma_dot = - (sim_param%g + sim_param%f)/(1.0_wp + 3.0_wp*sim_param%f) &
                        * wake%part_p(ip)%p%r_Vortex/wake%part_p(ip)%p%mag & 
                        * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir)
          endif        

          alpha_q_1 = wake%part_p(ip)%p%stretch*sim_param%dt*real(sim_param%ndt_update_wake,wp) 
          alpha_p_1 = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag + 1.0_wp/3.0_wp*alpha_q_1  
          if(sim_param%use_reformulated) then
            !> r_Vortex update
            r_Vortex_q_1 = sigma_dot * sim_param%dt*real(sim_param%ndt_update_wake,wp)
            r_Vortex_p_1 = wake%part_p(ip)%p%r_Vortex & 
                                        + 1.0_wp/3.0_wp*r_Vortex_q_1
          else 
            r_Vortex_p_1 = wake%part_p(ip)%p%r_Vortex
          endif
          if(norm2(alpha_p_1) .ge. sim_param%mag_threshold .and. r_Vortex_p_1 .ge. sim_param%mag_threshold) then 
            wake%part_p(ip)%p%mag = norm2(alpha_p_1) !> mag
            wake%part_p(ip)%p%dir = alpha_p_1/(wake%part_p(ip)%p%mag) !> direction 

            if(sim_param%use_reformulated) then
              wake%part_p(ip)%p%r_Vortex = r_Vortex_p_1 
            endif
            !> assign old values
            wake%part_p(ip)%p%cen_prev = wake%part_p(ip)%p%cen
            wake%part_p(ip)%p%dir_prev = wake%part_p(ip)%p%dir
            wake%part_p(ip)%p%mag_prev = wake%part_p(ip)%p%mag
            wake%part_p(ip)%p%vel_prev = q_1 !> velocity*dt
            wake%part_p(ip)%p%stretch_prev = alpha_q_1 !> stretch*dt
            wake%part_p(ip)%p%r_Vortex_prev = r_Vortex_q_1 !> sigma_dot*dt
          else !mag check
            wake%part_p(ip)%p%free = .true.
!$omp atomic update
            count_free = count_free + 1
!$omp end atomic   
          endif !mag check
        else !magnitude
          wake%part_p(ip)%p%free = .true.
!$omp atomic update
          count_free = count_free + 1
!$omp end atomic     
        endif !magnitude
      endif !not free
    enddo
!$omp end parallel do
     
    !> 2nd stage 
    if (sim_param%use_vs .and. sim_param%use_divfilt .and. sim_param%use_vd) then 
      call apply_multipole_optimized(wake%part_p, octree)
    else 
      call apply_multipole(wake%part_p, octree)
    endif 
!$omp parallel do schedule(dynamic,4) private(ip, q_2, alpha_q_2, alpha_p_2, sigma_dot, r_Vortex_q_2, r_Vortex_p_2)                        
    do ip = 1, n_part
      if ( .not. wake%part_p(ip)%p%free) then 
        if( wake%part_p(ip)%p%mag .ge. sim_param%mag_threshold) then ! to avoid negative magnitudes (and too small) 
          q_2 = wake%part_p(ip)%p%vel*sim_param%dt*real(sim_param%ndt_update_wake,wp) - &
                5.0_wp/9.0_wp*wake%part_p(ip)%p%vel_prev  
          wake%part_p(ip)%p%cen = wake%part_p(ip)%p%cen_prev + 15.0_wp/16.0_wp*q_2 

          sigma_dot = 0.0_wp
          if(sim_param%use_reformulated) then
            wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch & 
                                        - (sim_param%g + sim_param%f)/(1.0_wp/3.0_wp + sim_param%f) & 
                                        * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir) * wake%part_p(ip)%p%dir

            sigma_dot = - (sim_param%g + sim_param%f)/(1.0_wp + 3.0_wp*sim_param%f) &
                        * wake%part_p(ip)%p%r_Vortex/wake%part_p(ip)%p%mag & 
                        * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir)
          endif
          alpha_q_2 = wake%part_p(ip)%p%stretch*sim_param%dt*real(sim_param%ndt_update_wake,wp) - &
                      5.0_wp/9.0_wp*wake%part_p(ip)%p%stretch_prev 
          alpha_p_2 = wake%part_p(ip)%p%dir_prev*wake%part_p(ip)%p%mag_prev + 15.0_wp/16.0_wp*alpha_q_2 
          !> r_Vortex update
          if(sim_param%use_reformulated) then
            r_Vortex_q_2 = sigma_dot * sim_param%dt*real(sim_param%ndt_update_wake,wp) - &
                            5.0_wp/9.0_wp*wake%part_p(ip)%p%r_Vortex_prev
            r_Vortex_p_2 = wake%part_p(ip)%p%r_Vortex + 15.0_wp/16.0_wp*r_Vortex_q_2
          else
            r_Vortex_p_2 = wake%part_p(ip)%p%r_Vortex
          endif
          if(norm2(alpha_p_2) .ge. sim_param%mag_threshold .and. r_Vortex_p_2 .ge. sim_param%mag_threshold) then 
            wake%part_p(ip)%p%mag = norm2(alpha_p_2)
            wake%part_p(ip)%p%dir = alpha_p_2/(wake%part_p(ip)%p%mag)  

            if(sim_param%use_reformulated) then
              wake%part_p(ip)%p%r_Vortex = r_Vortex_p_2
            endif
            !> assign old values
            wake%part_p(ip)%p%cen_prev = wake%part_p(ip)%p%cen
            wake%part_p(ip)%p%dir_prev = wake%part_p(ip)%p%dir
            wake%part_p(ip)%p%mag_prev = wake%part_p(ip)%p%mag
            wake%part_p(ip)%p%vel_prev = q_2 !> velocity*dt
            wake%part_p(ip)%p%stretch_prev = alpha_q_2 !> stretch*dt
            wake%part_p(ip)%p%r_Vortex_prev = r_Vortex_q_2 !> sigma_dot*dt
          else !mag check
            wake%part_p(ip)%p%free = .true.
!$omp atomic update
            count_free = count_free + 1
!$omp end atomic    
          endif !mag check
        else !magnitude
          wake%part_p(ip)%p%free = .true.
!$omp atomic update
          count_free = count_free + 1
!$omp end atomic    
        endif !magnitude
      endif !not free
    enddo
!$omp end parallel do

    !> 3rd stage 
    if (sim_param%use_vs .and. sim_param%use_divfilt .and. sim_param%use_vd) then 
      call apply_multipole_optimized(wake%part_p, octree)
    else 
      call apply_multipole(wake%part_p, octree)
    endif 
  
!$omp parallel do schedule(dynamic,4) private(ip, pos_p, vel_in,vel_out, q_3, &
!$omp& alpha_q_3, alpha_p_3_mag, alpha_p_3_dir, alpha_p_3, sigma_dot, r_Vortex_q_3, r_Vortex_p_3)         
    do ip = 1, n_part 
    
      if ( .not. wake%part_p(ip)%p%free) then
        if( wake%part_p(ip)%p%mag .ge. sim_param%mag_threshold) then ! to avoid negative magnitudes (and too small) 
          q_3 = wake%part_p(ip)%p%vel*sim_param%dt*real(sim_param%ndt_update_wake,wp) - &
                153.0_wp/128.0_wp*wake%part_p(ip)%p%vel_prev 
          pos_p = wake%part_p(ip)%p%cen_prev + 8.0_wp/15.0_wp*q_3 
          sigma_dot = 0.0_wp
          if(sim_param%use_reformulated) then
            wake%part_p(ip)%p%stretch = wake%part_p(ip)%p%stretch & 
                                      - (sim_param%g + sim_param%f)/(1.0_wp/3.0_wp + sim_param%f) & 
                                      * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir) * wake%part_p(ip)%p%dir

            sigma_dot = - (sim_param%g + sim_param%f)/(1.0_wp + 3.0_wp*sim_param%f) &
                      * wake%part_p(ip)%p%r_Vortex/wake%part_p(ip)%p%mag & 
                      * sum(wake%part_p(ip)%p%stretch_alone*wake%part_p(ip)%p%dir)
          endif
          alpha_q_3 = wake%part_p(ip)%p%stretch*sim_param%dt - 153.0_wp/128.0_wp*wake%part_p(ip)%p%stretch_prev  
          alpha_p_3 = wake%part_p(ip)%p%dir_prev*wake%part_p(ip)%p%mag_prev + 8.0_wp/15.0_wp*alpha_q_3

          !add Pedrizzetti contribution
          if(sim_param%use_divfilt) then 
            alpha_p_3 = alpha_p_3 + alpha_pedrizzetti(ip,:)
          endif
          if (sim_param%use_reformulated) then
            !r_Vortex update
            r_Vortex_q_3 = sigma_dot * sim_param%dt*real(sim_param%ndt_update_wake,wp) - &
                            153.0_wp/128.0_wp*wake%part_p(ip)%p%r_Vortex_prev
            r_Vortex_p_3 = wake%part_p(ip)%p%r_Vortex + 8.0_wp/15.0_wp*r_Vortex_q_3
          else 
            r_Vortex_p_3 = wake%part_p(ip)%p%r_Vortex
          endif
          if(norm2(alpha_p_3) .ge. sim_param%mag_threshold .and. r_Vortex_p_3 .ge. sim_param%mag_threshold) then
            alpha_p_3_mag = norm2(alpha_p_3)
            alpha_p_3_dir = alpha_p_3/(alpha_p_3_mag)
              
            if(all(pos_p .ge. wake%part_box_min) .and. &
                all(pos_p .le. wake%part_box_max)) then
              wake%part_p(ip)%p%cen = pos_p
              wake%part_p(ip)%p%dir = alpha_p_3_dir
              wake%part_p(ip)%p%mag = alpha_p_3_mag
              if (sim_param%use_reformulated) then
                wake%part_p(ip)%p%r_Vortex = r_Vortex_p_3
              endif
            else !part_box
              wake%part_p(ip)%p%free = .true.
!$omp atomic update
              count_free = count_free + 1
!$omp end atomic
            endif !part_box
          else !mag check
            wake%part_p(ip)%p%free = .true.
!$omp atomic update
            count_free = count_free + 1
!$omp end atomic
          endif !mag check
        else !magnitude
          wake%part_p(ip)%p%free = .true.
!$omp atomic update
          count_free = count_free + 1
!$omp end atomic
        endif !magnitude
      endif !not free
    enddo
!$omp end parallel do

wake%n_prt = wake%n_prt - count_free
if(sim_param%use_divfilt) then
  deallocate(alpha_pedrizzetti)
endif 
end select 

end subroutine complete_wake

!----------------------------------------------------------------------

subroutine compute_vel_from_all(wake, pos, vel)
  type(t_wake), intent(in)        :: wake
  real(wp), intent(in)            :: pos(3)
  real(wp), intent(out)           :: vel(3)

  integer                         :: ie
  real(wp)                        :: v(3)

  vel = 0.0_wp

  !calculate the influence of particles
  do ie=1,size(wake%part_p)
    call wake%part_p(ie)%p%compute_vel(pos, v)
    vel = vel + v*one_4pi
  enddo

end subroutine compute_vel_from_all

!----------------------------------------------------------------------

subroutine get_vel_free(this, wake, pos, vel)
  class(t_free_wake)                    :: this
  type(t_wake), intent(in)              :: wake
  real(wp), intent(in)                  :: pos(3)
  real(wp), intent(out)                 :: vel(3)

  call compute_vel_from_all(wake, pos, vel)

  !vel = vel + sim_param%u_inf
  vel = vel + variable_wind(pos, sim_param%time) 

end subroutine get_vel_free

!---------------------------------------------------------------------

end module mod_wake_test
