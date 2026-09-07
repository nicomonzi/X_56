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


!> Module to treat several parts of the pressure integral equation

module mod_pressure_equation


use mod_param, only: &
  wp, nl, max_char_len, pi

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_math, only: &
  cross, vec2mat

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_aeroel, only: &
  t_impl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_linsys_vars, only: t_linsys

use mod_geometry, only: &
  t_geo

use mod_wake, only: &
  t_wake

!----------------------------------------------------------------------

implicit none

public :: dump_linsys_pres, press_normvel_der, initialize_pressure_sys, &
          assemble_pressure_sys, solve_pressure_sys

private

character(len=*), parameter :: this_mod_name = 'mod_pressure_equation'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Initialize the linear system of the pressure
!!
!! DISCONTINUED: operations performed inside initialize linsys
subroutine  initialize_pressure_sys(linsys, geo, elems)
 type(t_linsys), intent(inout), target :: linsys
 type(t_geo), intent(in) :: geo
 type(t_impl_elem_p), intent(inout) :: elems(:)
 integer :: ie

  !Set the number of surface panels
  linsys%n_sp = geo%nSurfPan
  linsys%nstatic_sp = geo%nstatic_SurfPan
  linsys%nmoving_sp = geo%nSurfPan - geo%nstatic_SurfPan

  !Allocate matrix and rhs for pressure integral equation
  allocate( linsys%idSurfPan(geo%nSurfPan) )
  linsys%idSurfPan    = geo%idSurfPan
  allocate( linsys%idSurfPanG2L(geo%nSurfPan) )
  linsys%idSurfPanG2L = geo%idSurfPanG2L
  allocate( linsys%A_pres(geo%nSurfPan,geo%nSurfPan) )
  linsys%A_pres = 0.0_wp
  allocate( linsys%b_pres(geo%nSurfPan) )
  linsys%b_pres = 0.0_wp
  allocate( linsys%res_pres(geo%nSurfPan) )
  linsys%res_pres = 0.0_wp
  allocate( linsys%b_static_pres(geo%nstatic_SurfPan,geo%nstatic_SurfPan) )
  linsys%b_static_pres = linsys%b_static( &
         geo%idSurfPan(1:geo%nstatic_SurfPan), &
         geo%idSurfPan(1:geo%nstatic_SurfPan) )

  do ie = 1 , geo%nSurfpan
    select type( el => elems(geo%idSurfPan(ie))%p ) ; class is(t_surfpan)
      el%pres_sol => linsys%res_pres(ie)
    end select
  end do

end subroutine initialize_pressure_sys


!----------------------------------------------------------------------

!> Assemble the pressure equation system
subroutine assemble_pressure_sys(linsys, geo, elems, wake)
 type(t_linsys), intent(inout) :: linsys
 type(t_geo), intent(in) :: geo
 type(t_impl_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(in) ::  wake

 integer :: ie, ip, iw, p1, p2, inext, is
 integer :: ipp(4) , iww(4), ntot
 real(wp) :: elcen(3), dist(3), dist2(3)
 real(wp) :: Pinf, rhoinf

  ! Free-stream conditions
  Pinf   = sim_param%P_inf
  rhoinf = sim_param%rho_inf
  ntot = linsys%rank


  ! Pressure integral equation +++++++++++++++++++++++++++++++++++++++++

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! Assemble the RHS of the linear system for the Bernoulli polynomial     !
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! RHS = B_infty +                                 (a) far field
  !       + \oint_{Sb} { G ( du/dt - ... ) } +      (b) time-dependent
  !       +  \int_{V } { DG . \omega x U } +        (c) rotational
  !       +  \int_{Sb} { viscous terms }            (d) viscous

  ! RHS assembled as the sum of:
  ! (b) time-dependent source contribution: computed and assembled in the
  !     %build_row routines above
  ! do nothing here
  ! (c) rotational effects (up to now, ignoring ring wakes)
!$omp parallel do private(ie, iw, dist, elcen) schedule(dynamic,32)
  do ie = 1 , geo%nSurfpan
    elcen = elems( geo%idSurfpan(ie) )%p%cen
    ! (c.1) particles ( part_p )
    if ( allocated( wake%part_p ) ) then
      do iw = 1 , size( wake%part_p )
        if  ( .not. ( wake%part_p(iw)%p%free ) ) then

          dist = wake%part_p(iw)%p%cen - elcen
          linsys%b_pres(ie) = linsys%b_pres(ie) + &
            sum(dist * cross( wake%part_p(iw)%p%dir , wake%part_p(iw)%p%vel ) ) * &
                wake%part_p(iw)%p%mag / ( (sqrt(sum(dist**2)+wake%part_p(iw)%p%r_Vortex**2))**3.0_wp )

        end if
      end do
    end if
  enddo
!$omp end parallel do

!$omp parallel do private(ie, iw, dist, dist2, p1, p2, ipp, iww,ip, is, inext, elcen) schedule(dynamic, 64)
  do ie = 1 , geo%nSurfpan
    elcen = elems( geo%idSurfpan(ie) )%p%cen
    ! (c.2) line elements ( end_vorts )
    do iw = 1 , wake%n_pan_stripes
      if( associated( wake%end_vorts(iw)%mag ) ) then
        dist  = wake%end_vorts(iw)%ver(:,1) - elcen
        dist2 = wake%end_vorts(iw)%ver(:,2) - elcen
        linsys%b_pres(ie) = linsys%b_pres(ie) - &
          0.5_wp * wake%end_vorts(iw)%mag * sum( wake%end_vorts(iw)%edge_vec * &
              ( cross(dist , wake%end_vorts(iw)%ver_vel(:,1) ) /(norm2(dist )**3.0_wp) + &
                cross(dist2, wake%end_vorts(iw)%ver_vel(:,2) ) /(norm2(dist2)**3.0_wp) ) )
      end if
    end do
    ! (c.3) constant surface doublets = vortex rings ( wake_panels )
    do ip = 1 , wake%pan_wake_len
      do iw = 1 , wake%n_pan_stripes

        p1 = wake%i_start_points(1,iw)
        p2 = wake%i_start_points(2,iw)

        ipp = (/ ip , ip , ip+1, ip+1 /)
        iww = (/ p1 , p2 , p2  , p1   /)

        do is = 1 , wake%wake_panels(iw,ip)%n_ver ! do is = 1 , 4

          inext = mod(is,wake%wake_panels(iw,ip)%n_ver) + 1
          dist  = wake%wake_panels(iw,ip)%ver(:,is   ) - elcen
          dist2 = wake%wake_panels(iw,ip)%ver(:,inext) - elcen

          linsys%b_pres(ie) = linsys%b_pres(ie) - &
            0.5_wp * wake%wake_panels(iw,ip)%mag * sum( wake%wake_panels(iw,ip)%edge_vec(:,is) * &
                ( cross(dist , wake%pan_w_vel(:,iww(is   ),ipp(is   )) ) /(norm2(dist )**3.0_wp) + &
                  cross(dist2, wake%pan_w_vel(:,iww(inext),ipp(inext)) ) /(norm2(dist2)**3.0_wp) ) )

        end do
      end do
    end do
    ! (c.4) rings from actuator disks
    ! TODO: ...

    ! (d) viscous effects
    ! TODO: to be implemented

  end do
!$omp end parallel do

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! END Assemble the RHS of the linear system for the Bernoulli polynomial !
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

end subroutine assemble_pressure_sys

!----------------------------------------------------------------------

!> Solve the linear system for the pressure
!!
!! The linear system for the pressure with a similar procedure to the one
!! used for the main linear system: the stored LU static decomposition is
!! completed by the decomposition of the dynamic part and the factorized
!! system is solved
subroutine solve_pressure_sys(linsys)
  type(t_linsys), intent(inout) :: linsys

  integer              :: INFO
  character(len=max_char_len) :: msg
  character(len=*), parameter :: this_sub_name = 'solve_linsys_pressure'

  ! Operations on the side band matrices: done only if the system is
  ! mixed static/dynamic and those matrices exists
  if (linsys%nstatic_sp .gt. 0 .and. linsys%nmoving_sp .gt.0) then

    !=>Create the upper-diagonal block Usd
    !Swap in place Asd to get PssAsd
#if (DUST_PRECISION==1)
    call slaswp(linsys%nmoving_sp, &
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp,1,linsys%nstatic_sp, &
          linsys%P_pres(1:linsys%nstatic_sp),1)
#elif (DUST_PRECISION==2)
    call dlaswp(linsys%nmoving_sp, &
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp,1,linsys%nstatic_sp, &
          linsys%P_pres(1:linsys%nstatic_sp),1)
#endif /*DUST_PRECISION*/

    !Solve Lss Usd = Pss Asd to get Usd and put it in place of Asd
#if (DUST_PRECISION==1)
    call strsm('L','L','N','U',linsys%nstatic_sp,linsys%nmoving_sp,1.0d+0,   &
          linsys%A_pres(1:linsys%nstatic_sp,1:linsys%nstatic_sp), &
          linsys%nstatic_sp, &
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp)
#elif (DUST_PRECISION==2)
    call dtrsm('L','L','N','U',linsys%nstatic_sp,linsys%nmoving_sp,1.0d+0,   &
          linsys%A_pres(1:linsys%nstatic_sp,1:linsys%nstatic_sp), &
          linsys%nstatic_sp, &
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp)
#endif /*DUST_PRECISION*/

    !==>Solve the lower-diagoal block Lds
    !Solve Pdd-1 Lds Uss = Ads for Pdd-1 Lds and put it in place of Ads
#if (DUST_PRECISION==1)
    call strsm('R','U','N','N',linsys%nmoving_sp,linsys%nstatic_sp,1.0d+0,   &
          linsys%A_pres(1:linsys%nstatic_sp,1:linsys%nstatic_sp), &
          linsys%nstatic_sp,  &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp),&
          linsys%nmoving_sp)
#elif (DUST_PRECISION==2)
    call dtrsm('R','U','N','N',linsys%nmoving_sp,linsys%nstatic_sp,1.0d+0,   &
          linsys%A_pres(1:linsys%nstatic_sp,1:linsys%nstatic_sp), &
          linsys%nstatic_sp,  &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp),&
          linsys%nmoving_sp)
#endif /*DUST_PRECISION*/

    !==>Modify the dynamic square block
    !Modify the square block from Add to Add - Pdd-1Lds Usd
#if (DUST_PRECISION==1)
    call sgemm('N','N',linsys%nmoving_sp,linsys%nmoving_sp,linsys%nstatic_sp, &
          -1.0d+0, &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp),&
          linsys%nmoving_sp,&
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp,1.0d+0,&
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nmoving_sp)
#elif (DUST_PRECISION==2)
    call dgemm('N','N',linsys%nmoving_sp,linsys%nmoving_sp,linsys%nstatic_sp, &
          -1.0d+0, &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp),&
          linsys%nmoving_sp,&
          linsys%A_pres(1:linsys%nstatic_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nstatic_sp,1.0d+0,&
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,linsys%nstatic_sp+1:linsys%n_sp),&
          linsys%nmoving_sp)
#endif /*DUST_PRECISION*/

  endif

  ! If the system has a dynamic part, factorize such part
  if (linsys%nmoving_sp .gt. 0) then

    !==>Factorize and put in place the square dynamic block
#if (DUST_PRECISION==1)
    call sgetrf(linsys%nmoving_sp,linsys%nmoving_sp, &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,linsys%nstatic_sp+1:linsys%n_sp), &
          linsys%nmoving_sp,linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp),info)
#elif (DUST_PRECISION==2)
    call dgetrf(linsys%nmoving_sp,linsys%nmoving_sp, &
          linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,linsys%nstatic_sp+1:linsys%n_sp), &
          linsys%nmoving_sp,linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp),info)
#endif /*DUST_PRECISION*/
    if ( info .ne. 0 ) then
      write(msg,*) 'error while factorizing the dynamic  block of the &
                    &pressure linear system, Lapack DGETRF error code ', info
      call error(this_sub_name, this_mod_name, trim(msg))
    end if

  endif

  ! If the system is mixed finish the operation on the band blocks
  if (linsys%nstatic .gt. 0 .and. linsys%nmoving .gt.0) then
  !==> Permute the lower mixed bloc
#if (DUST_PRECISION==1)
  call slaswp(linsys%nstatic_sp, &
              linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp), &
              linsys%nmoving_sp,1,linsys%nmoving_sp,&
              linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp),1)
#elif (DUST_PRECISION==2)
  call dlaswp(linsys%nstatic_sp, &
              linsys%A_pres(linsys%nstatic_sp+1:linsys%n_sp,1:linsys%nstatic_sp), &
              linsys%nmoving_sp,1,linsys%nmoving_sp,&
              linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp),1)
#endif /*DUST_PRECISION*/
  endif

  !==> Fix the lower part of the permutation matrix to make it global
  linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp) = &
  linsys%P_pres(linsys%nstatic_sp+1:linsys%n_sp) + linsys%nstatic_sp

  !==> Solve the factorized system
  linsys%res_pres = linsys%b_pres
#if (DUST_PRECISION==1)
  call sgetrs('N',linsys%n_sp,1,linsys%A_pres,linsys%n_sp,linsys%P_pres, &
              linsys%res_pres, linsys%n_sp,info)
#elif (DUST_PRECISION==2)
  call dgetrs('N',linsys%n_sp,1,linsys%A_pres,linsys%n_sp,linsys%P_pres, &
              linsys%res_pres, linsys%n_sp,info)
#endif /*DUST_PRECISION*/
  if ( info .ne. 0 ) then
    write(msg,*) 'error while solving the factorized pressure system &
      &Lapack DGETRS error code ', info
    call error(this_sub_name, this_mod_name, trim(msg))
  end if

end subroutine solve_pressure_sys

!----------------------------------------------------------------------

!> Compute the normal velocity derivative
!!
!! compute the time derivative of the normal component of the velocity on
!!  surfpan to be used in the source rhs of the Bernoulli integral equation.
!! surf_vel_SurfPan_old should be saved at the end of the time step
subroutine press_normvel_der(geo, elems, surf_vel_SurfPan_old)
  type(t_geo), intent(in) :: geo
  type(t_impl_elem_p), intent(inout) :: elems(:)
  real(wp)             :: R_cen(3,3) 
  real(wp), intent(in) :: surf_vel_SurfPan_old(:,:)

  integer :: i_el, i_e
  real(wp) :: GradS_Un(3), DivS_U
  character(len=*), parameter :: this_sub_name = 'press_normvel_der'

  do i_el = 1 , geo%nSurfPan

    select type ( el => elems(geo%idSurfPan(i_el))%p ) ; class is ( t_surfpan )

#if USE_PRECICE
      if (geo%components(elems(i_el)%p%comp_id)%coupling) then
        R_cen = el%R_cen
      else
        R_cen = geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g
      endif 
#else 
      R_cen = geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g
#endif 

      el%dUn_dt = sum(el%nor * ( el%ub - &
                      surf_vel_SurfPan_old( i_el , : ) ) ) / sim_param%dt

      ! Compute GradS_Un
      GradS_Un = 0.0_wp
      do i_e = 1 , el%n_ver
        if ( associated(el%neigh(i_e)%p) ) then !  .and. &
          select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            GradS_Un = GradS_Un + &
              matmul(R_cen, &
                el%pot_vel_stencil(:,i_e) * ( &
                  sum(el%nor* (el_neigh%surf_vel - el%surf_vel) ) ) )
          end select
        else

            GradS_Un = GradS_Un + &
              matmul(R_cen, &
               el%pot_vel_stencil(:,i_e) * ( &
                  sum(el%nor* ( - 2.0_wp * el%surf_vel) ) ) )

        end if
      end do

      ! Compute DivS_U
      DivS_U = 0.0_wp
      do i_e = 1 , el%n_ver
        if ( associated(el%neigh(i_e)%p) ) then !  .and. &
          select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            
            DivS_U = DivS_U + &
                      sum( &
                      matmul(R_cen,   &
                                                            el%pot_vel_stencil(:,i_e) ) * &
                            ( el_neigh%surf_vel - el%surf_vel )   )
          end select
        else
          DivS_U = DivS_U + &
            sum( &
              matmul(R_cen,   &
                                                         el%pot_vel_stencil(:,i_e) ) * &
                 ( - 2.0_wp * el%surf_vel ) )

        end if
      end do

      ! Compute "source intensity" of Bernoulli equations
      el%bernoulli_source = + el%dUn_dt &                    !    n . DU/Dt
                            - sum( GradS_Un * ( el%ub ))   & !  - GradS_Un . el%ub
                            + DivS_U * sum(el%ub*el%nor)     !  + Un * Div_S U

    end select
  end do

end subroutine press_normvel_der


!> Dump the matrix and rhs of the pressure linear system to file
!!
!! TODO: consider moving these functionalities to the i/o modules
subroutine dump_linsys_pres(linsys , filen_A , filen_b )
  type(t_linsys), intent(in) :: linsys
  character(len=*) , intent(in) :: filen_A , filen_b

  integer :: fid
  integer :: i1


  fid = 23
  open(unit=fid, file=trim(adjustl(filen_A)) )
  do i1 = 1 , size(linsys%A_pres,1)
    write(fid,*) linsys%A_pres(i1,:)
  end do
  close(fid)

  fid = 24
  open(unit=fid, file=trim(adjustl(filen_b)) )
  do i1 = 1 , size(linsys%b_pres,1)
    write(fid,*) linsys%b_pres(i1)
  end do
  close(fid)

end subroutine dump_linsys_pres

!----------------------------------------------------------------------

end module mod_pressure_equation