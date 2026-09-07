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


!> Module to handle the input/output of the solution of DUST in hdf5 format,
!! meant principally to exchange data with the post processor and restart
module mod_dust_io

use mod_param, only: &
  wp, max_char_len, nl

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_liftlin, only: &
  t_liftlin

use mod_virtual, only: &
  c_elem_virtual, t_elem_virtual_p

use mod_reference, only: &
  t_ref

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

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_wake, only: &
  t_wake

use mod_stringtools, only: &
  stricmp

use mod_version, only: &
  git_sha1, version


use mod_parse, only: &
  t_parse, getstr, getint, getintarray, getreal, getrealarray, getlogical, countoption, &
  getsuboption,&
  finalizeparameters, t_link, check_opt_consistency


!----------------------------------------------------------------------

implicit none

public :: save_status, load_solution

private


character(len=*), parameter :: this_mod_name = 'mod_dust_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine save_status(geo, wake,  it, time, run_id)
  type(t_geo), intent(in)           :: geo
  type(t_wake), intent(in)          :: wake
  integer, intent(in)               :: it
  real(wp), intent(in)              :: time
  integer, intent(in)               :: run_id(10)

  integer(h5loc)                    :: floc, gloc1, gloc2, gloc3, gloc4, ploc
  character(len=max_char_len)       :: comp_name, sit
  integer                           :: icomp, ncomp
  integer                           :: iref, nref
  character(len=max_char_len)       :: ref_name
  integer                           :: ie, ne
  real(wp), allocatable             :: vort(:), cp(:) , pres(:), didou_dt(:)
  real(wp), allocatable             :: dforce(:,:), dmom(:,:), surf_vel(:,:)
  real(wp), allocatable             :: turbvisc(:), v_rad(:), vp_parentid(:)
  real(wp), allocatable             :: points_w(:,:,:), cent(:,:,:) , vel_w(:,:,:)
  real(wp), allocatable             :: vort_v(:,:)
  integer,  allocatable             :: conn_pe(:)
  real(wp), allocatable             :: ori(:,:)
  !LL data
  real(wp), allocatable             :: alpha(:), vel_2d(:), vel_outplane(:)
  real(wp), allocatable             :: up_x(:), up_y(:), up_z(:)
  real(wp), allocatable             :: aero_coeff(:,:)
  real(wp), allocatable             :: Gamma_old(:), Gamma_old_old(:), dGamma_dt(:)
  !VL corrected data
  real(wp), allocatable             :: alpha_vl(:), vel_2d_vl(:), vel_outplane_vl(:)
  real(wp), allocatable             :: up_x_vl(:), up_y_vl(:), up_z_vl(:)
  real(wp), allocatable             :: alpha_isolated_vl(:), vel_2d_isolated_vl(:)
  real(wp), allocatable             :: vel_outplane_isolated_vl(:), aero_coeff_vl(:,:)
  integer                           :: ir, id, ip, np, ih
  ! Hinges
  character(len=2)                  :: hinge_id_str
  character(len=5)                  :: aero_table 

  !> create the output file
  write(sit,'(I4.4)') it
  call new_hdf5_file(trim(sim_param%basename)//'_res_'//trim(sit)//'.h5', &
                      floc)
  call write_hdf5_attr(git_sha1, 'git_sha1', floc)
  call write_hdf5_attr(version, 'version', floc)
  call sim_param%save_param(floc)
  call write_hdf5(time,'time',floc)

  call new_hdf5_group(floc, 'Parameters', ploc)
  call write_hdf5(sim_param%u_inf,'u_inf', ploc)
  call write_hdf5(sim_param%P_inf,'P_inf', ploc)
  call write_hdf5(sim_param%rho_inf,'rho_inf', ploc)
  call write_hdf5(sim_param%mu_inf,'mu_inf', ploc)
  call write_hdf5(sim_param%a_inf,'a_inf', ploc)
  
  
  call close_hdf5_group(ploc)

  ! 1) %%%%% Component solution:
  ! mimic the structure of the components inside the geometry input file

  call new_hdf5_group(floc, 'Components', gloc1)
  ncomp = size(geo%components)
  call write_hdf5(ncomp, 'NComponents', gloc1)
  !Cycle the components
  do icomp = 1, ncomp
    !write(comp_name,'(A,I3.3)')'Comp',icomp
    write(comp_name,'(A,I3.3)')'Comp',geo%components(icomp)%comp_id
    call new_hdf5_group(gloc1, trim(comp_name), gloc2)

    call write_hdf5(trim(geo%components(icomp)%comp_name),'CompName',gloc2)
    call write_hdf5(trim(geo%components(icomp)%ref_tag),'RefTag',gloc2)
    call write_hdf5(geo%components(icomp)%ref_id,'RefId',gloc2)

    call new_hdf5_group(gloc2, 'Solution', gloc3)

    ne = size(geo%components(icomp)%el)
    allocate(vort(ne), cp(ne) , pres(ne) , dforce(3,ne), dmom(3,ne))
    do ie = 1,ne
      vort(ie) = geo%components(icomp)%el(ie)%mag
      pres(ie) = geo%components(icomp)%el(ie)%pres
      dforce(:,ie) = geo%components(icomp)%el(ie)%dforce
      dmom(:,ie) = geo%components(icomp)%el(ie)%dmom
    enddo
    call write_hdf5(vort,'Vort',gloc3)
    call write_hdf5(pres,'Pres',gloc3)
    call write_hdf5(dforce,'dF',gloc3)
    call write_hdf5(dmom,'dMom',gloc3)
    deallocate(vort, cp, pres, dforce, dmom)

    !> Output the surface velocity
    if ( trim( geo%components(icomp)%comp_el_type ) .eq. 'p' ) then
      allocate(surf_vel(3,ne))
      do ie = 1,ne
        select type( el => geo%components(icomp)%el(ie) ) ; type is (t_surfpan)
          surf_vel(:,ie) = el%surf_vel
        end select
      end do
      call write_hdf5(surf_vel,'surf_vel',gloc3)
      deallocate(surf_vel)
    endif

    !> Output the lifting lines data
    if ( trim( geo%components(icomp)%comp_el_type ) .eq. 'l' ) then
      allocate(alpha(ne), vel_2d(ne), vel_outplane(ne))
      allocate(up_x(ne), up_y(ne), up_z(ne))
      allocate(aero_coeff(3,ne))
      allocate(Gamma_old(ne), Gamma_old_old(ne), dGamma_dt(ne))
      
      do ie = 1,ne
        select type( el => geo%components(icomp)%el(ie) ) ; type is (t_liftlin)
          alpha(ie) = el%alpha
          vel_2d(ie) = el%vel_2d
          up_x(ie) = el%up_x
          up_y(ie) = el%up_y
          up_z(ie) = el%up_z
          vel_outplane(ie) = el%vel_outplane
          aero_coeff(:,ie) = el%aero_coeff
          Gamma_old(ie) = el%Gamma_old
          Gamma_old_old(ie) = el%Gamma_old_old
          dGamma_dt(ie) = el%dGamma_dt
        end select
      end do    
      call write_hdf5(alpha,                  'alpha',                  gloc3)
      call write_hdf5(vel_2d,                 'vel_2d',                 gloc3)
      call write_hdf5(up_x,                   'up_x',                   gloc3)
      call write_hdf5(up_y,                   'up_y',                   gloc3)
      call write_hdf5(up_z,                   'up_z',                   gloc3)
      call write_hdf5(vel_outplane,           'vel_outplane',           gloc3)
      call write_hdf5(transpose(aero_coeff),  'aero_coeff',             gloc3)
      call write_hdf5(Gamma_old,              'Gamma_old',              gloc3)
      call write_hdf5(Gamma_old_old,          'Gamma_old_old',          gloc3)
      call write_hdf5(dGamma_dt,              'dGamma_dt',              gloc3)
      
      deallocate(alpha, vel_2d, vel_outplane, up_x, up_y, up_z)
      deallocate(aero_coeff)
      deallocate(Gamma_old, Gamma_old_old, dGamma_dt)
      
    !> Output the vl data
    elseif (trim(geo%components(icomp)%comp_el_type) .eq. 'v') then
      ne = size(geo%components(icomp)%el)
      allocate(didou_dt(ne))
      do ie = 1,ne
        select type( el => geo%components(icomp)%el(ie) ) ; type is (t_vortlatt)
          didou_dt(ie) = el%didou_dt          
        end select
      end do
      call write_hdf5(didou_dt,'didou_dt',gloc3)
      deallocate(didou_dt)
      
      ! if it is corrected
      if(trim(geo%components(icomp)%aero_correction) .eq. 'true') then 
      
        ne = size(geo%components(icomp)%stripe)
        allocate(alpha_vl(ne), vel_2d_vl(ne), vel_outplane_vl(ne))
        allocate(up_x_vl(ne), up_y_vl(ne), up_z_vl(ne))
        
        allocate(alpha_isolated_vl(ne), vel_2d_isolated_vl(ne), &
                  vel_outplane_isolated_vl(ne))
        allocate(aero_coeff_vl(3,ne))
        aero_table = geo%components(icomp)%aero_correction
        do ie = 1, ne
          alpha_vl(ie) = geo%components(icomp)%stripe(ie)%alpha
          vel_2d_vl(ie) = geo%components(icomp)%stripe(ie)%vel_2d
          up_x_vl(ie) = geo%components(icomp)%stripe(ie)%up_x
          up_y_vl(ie) = geo%components(icomp)%stripe(ie)%up_y
          up_z_vl(ie) = geo%components(icomp)%stripe(ie)%up_z
          vel_outplane_vl(ie) = geo%components(icomp)%stripe(ie)%vel_outplane
          alpha_isolated_vl(ie) = geo%components(icomp)%stripe(ie)%alpha_isolated
          vel_2d_isolated_vl(ie) = geo%components(icomp)%stripe(ie)%vel_2d_isolated
          vel_outplane_isolated_vl(ie) = geo%components(icomp)%stripe(ie)%vel_outplane_isolated
          aero_coeff_vl(:,ie) = geo%components(icomp)%stripe(ie)%aero_coeff
        end do     
        call write_hdf5(alpha_vl,                 'alpha_vl',                 gloc3)
        call write_hdf5(vel_2d_vl,                'vel_2d_vl',                gloc3)
        call write_hdf5(up_x_vl,                  'up_x_vl',                  gloc3)
        call write_hdf5(up_y_vl,                  'up_y_vl',                  gloc3)
        call write_hdf5(up_z_vl,                  'up_z_vl',                  gloc3)
        call write_hdf5(vel_outplane_vl,          'vel_outplane_vl',          gloc3)
        call write_hdf5(alpha_isolated_vl,        'alpha_isolated_vl',        gloc3)
        call write_hdf5(vel_2d_isolated_vl,       'vel_2d_isolated_vl',       gloc3)
        call write_hdf5(vel_outplane_isolated_vl, 'vel_outplane_isolated_vl', gloc3)
        call write_hdf5(transpose(aero_coeff_vl), 'aero_coeff_vl',            gloc3)
        call write_hdf5(aero_table,               'aero_table',               gloc3)

        deallocate(alpha_vl, vel_2d_vl, vel_outplane_vl, up_x_vl, up_y_vl, up_z_vl)
        deallocate(alpha_isolated_vl, vel_2d_isolated_vl, vel_outplane_isolated_vl)
        deallocate(aero_coeff_vl)

      endif
    !> Output the surfpan data  
    elseif (trim(geo%components(icomp)%comp_el_type) .eq. 'p') then
      ne = size(geo%components(icomp)%el)
      allocate(didou_dt(ne))
      do ie = 1,ne
        select type( el => geo%components(icomp)%el(ie) ) ; type is (t_surfpan)
          didou_dt(ie) = el%didou_dt
        end select
      end do
      call write_hdf5(didou_dt,'didou_dt',gloc3)
      deallocate(didou_dt)   
    endif
    call close_hdf5_group(gloc3)

    !> Hinges: node position, node orientation(v,h,n), theta, theta_old
    call new_hdf5_group(gloc2, 'Hinges', gloc3)
    call write_hdf5(geo%components(icomp)%n_hinges, 'n_hinges', gloc3)
    do ih = 1, geo%components(icomp)%n_hinges
      !> Open hinge group
      write(hinge_id_str,'(I2.2)') ih
      call new_hdf5_group(gloc3, 'Hinge_'//hinge_id_str, gloc4)
      associate( comp => geo%components(icomp) )
        call write_hdf5(comp%hinge(ih)%act%rr, 'act_rr', gloc4)
        call write_hdf5(comp%hinge(ih)%act%v , 'act_v ', gloc4)
        call write_hdf5(comp%hinge(ih)%act%h , 'act_h ', gloc4)
        call write_hdf5(comp%hinge(ih)%act%n , 'act_n ', gloc4)
        call write_hdf5(comp%hinge(ih)%theta     , 'theta'     , gloc4 )
        call write_hdf5(comp%hinge(ih)%theta_old , 'theta_old' , gloc4 )
      end associate
      call close_hdf5_group(gloc4)
    end do
    call close_hdf5_group(gloc3)

    if(sim_param%output_detailed_geo) then
      call new_hdf5_group(gloc2, 'Geometry', gloc3)
      call write_hdf5(geo%points(:,geo%components(icomp)%i_points), &
                      'rr',gloc3)
      call write_hdf5(geo%points_virtual(:,geo%components(icomp)%i_points_virtual), &
                      'rr_virtual',gloc3)
      call close_hdf5_group(gloc3)
    endif

#if USE_PRECICE
    if ( geo%components(icomp)%coupling ) then
      call new_hdf5_group(gloc2, 'Geometry', gloc3)
      call write_hdf5(geo%points(:,geo%components(icomp)%i_points), &
                      'rr', gloc3)
      call write_hdf5(geo%points_virtual(:,geo%components(icomp)%i_points_virtual), &
                      'rr_virtual',gloc3)   
                                         
      !> get orientation for sectional loads 
      ne = size(geo%components(icomp)%el)
      allocate(ori(ne, 3))
      do ie = 1,ne
        ori(ie, :) = geo%components(icomp)%el(ie)%ori
      end do
      call write_hdf5(ori, 'ori', gloc3)
      deallocate(ori)
      call close_hdf5_group(gloc3)
    end if
#endif

    call close_hdf5_group(gloc2)
  enddo
  call close_hdf5_group(gloc1)

  ! 2) %%%% Wake:
  ! just print the whole points, the solution and the starting row
  ! connectivity to build connectivity after
  !=== Panels ===
  call new_hdf5_group(floc, 'PanelWake', gloc1)
  call write_hdf5(wake%pan_w_points(:,:,1:wake%pan_wake_len+1),'WakePoints',gloc1 )
  call write_hdf5(wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1),'WakeVels'  ,gloc1 ) ! <<<< restart with Bernoulli integral equation
  call write_hdf5(wake%i_start_points,'StartPoints',gloc1)
  call write_hdf5(wake%pan_idou(:,1:wake%pan_wake_len),'WakeVort',gloc1)

  call close_hdf5_group(gloc1)

  !=== Rings ===
  call new_hdf5_group(floc, 'RingWake', gloc1)
  allocate(points_w(3,wake%np_row,wake%rin_wake_len))
  allocate(conn_pe(wake%np_row))
  allocate(cent(3,wake%ndisks,wake%rin_wake_len))
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    do ir = 1,wake%rin_wake_len
      points_w(:,ip:ip+np-1,ir) = wake%wake_rings(id,ir)%ver(:,:)
      cent(:,id,ir) = wake%wake_rings(id,ir)%cen(:)
    enddo
    conn_pe(ip:ip+np-1) = id
    ip = ip+np
  enddo
  call write_hdf5(points_w,'WakePoints',gloc1)
  call write_hdf5(conn_pe,'Conn_pe',gloc1)
  call write_hdf5(cent,'WakeCenters',gloc1)
  call write_hdf5(wake%rin_idou(:,1:wake%rin_wake_len),'WakeVort',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, conn_pe, cent)

  !=== Particles ===
  call new_hdf5_group(floc, 'ParticleWake', gloc1)
  allocate(points_w(3,wake%n_prt,1))
  allocate(vort_v(3,wake%n_prt))
  allocate(turbvisc(wake%n_prt))
  allocate(vel_w(3,wake%n_prt,1))
  allocate(v_rad(wake%n_prt))
  allocate(vp_parentid(wake%n_prt))

  do ip = 1, wake%n_prt
    points_w(:,ip,1) = wake%part_p(ip)%p%cen
    vel_w(:,ip,1) = wake%part_p(ip)%p%vel
    vort_v(:,ip) = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag
    turbvisc(ip) = wake%part_p(ip)%p%turbvisc
    v_rad(ip) = wake%part_p(ip)%p%r_Vortex
    vp_parentid(ip) = wake%part_p(ip)%p%parent_id
  enddo
  call write_hdf5(points_w(:,:,1),'WakePoints',gloc1)
  call write_hdf5(   vel_w(:,:,1),'WakeVels'  ,gloc1)
  call write_hdf5( turbvisc,'turbvisc'  ,gloc1)
  call write_hdf5( v_rad,'VortexRad'  ,gloc1)
  call write_hdf5( vp_parentid,'ParentId'  ,gloc1)
  call write_hdf5(vort_v,'WakeVort',gloc1)
  call write_hdf5(wake%last_pan_idou,'LastPanIdou',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, vort_v, vel_w, turbvisc, v_rad, vp_parentid)

  ! 3) %%%% References
  call new_hdf5_group(floc, 'References', gloc1)
  nref = size(geo%refs)
  call write_hdf5(nref, 'NReferences', gloc1)
  do iref = 0, nref-1
    write(ref_name,'(A,I3.3)')'Ref',iref
    call new_hdf5_group(gloc1, trim(ref_name), gloc2)

    call write_hdf5(geo%refs(iref)%tag, 'Tag', gloc2)
    call write_hdf5(geo%refs(iref)%of_g, 'Offset',gloc2)
    call write_hdf5(geo%refs(iref)%R_g, 'R',gloc2)
    call write_hdf5(geo%refs(iref)%f_g, 'Vel',gloc2)
    call write_hdf5(geo%refs(iref)%G_g, 'RotVel',gloc2)
    call write_hdf5(geo%refs(iref)%relative_pol_pos, 'RelativePolPos',gloc2)
    call write_hdf5(geo%refs(iref)%relative_rot_pos, 'RelativeRotPos',gloc2)

    call close_hdf5_group(gloc2)
  enddo
  call close_hdf5_group(gloc1)
    
  call close_hdf5_file(floc)
  
end subroutine save_status

!----------------------------------------------------------------------
!> Load the solution from a pre-existent solution file
!!
!! Loads just the bodies solution, which is stored into the components,
!! the loading of the wakes is left to other subroutines
subroutine load_solution(filename,comps,refs)
  character(len=max_char_len), intent(in) :: filename
  type(t_geo_component), intent(inout)    :: comps(:)
  type(t_ref), intent(in)                 :: refs(0:)

  integer(h5loc)                          :: floc, gloc1, gloc2, gloc3
  integer                                 :: ncomp, icomp, icomp2
  integer                                 :: ne, ie
  character(len=max_char_len)             :: comp_name_read, comp_id
  real(wp), allocatable                   :: idou(:), pres(:), dF(:,:), didou_dt(:)
  real(wp), allocatable                   :: Gamma_old(:), Gamma_old_old(:), dGamma_dt(:)
  character(len=*), parameter             :: this_sub_name = 'load_solution'

  call open_hdf5_file(filename, floc)

  !TODO: check the run id with the geo file
  call open_hdf5_group(floc, 'Components', gloc1)

  !Check if the number of components is consistent
  call read_hdf5(ncomp, 'NComponents', gloc1)
  if(ncomp .ne. size(comps)) call error(this_sub_name, this_mod_name, &
  'Mismatching number of components between geometry and solution files')

  !Cycle on all the components on the file and then on all the local
  !components: they come from different files and can be mixed up in the
  !order
  ! TODO check if this means that the elem can switch
  do icomp = 1, ncomp
    write(comp_id,'(A,I3.3)')'Comp',icomp
    call open_hdf5_group(gloc1, trim(comp_id), gloc2)
    
    call read_hdf5(comp_name_read,'CompName',gloc2)
    do icomp2 = 1,size(comps)
      if (stricmp(comp_name_read, comps(icomp2)%comp_name)) then
        !We found the correct local component, read all
        call check_ref(gloc2, floc, refs(comps(icomp2)%ref_id))
        call open_hdf5_group(gloc2, 'Solution', gloc3)
        call read_hdf5_al(idou,'Vort',gloc3)
        call read_hdf5_al(pres,'Pres',gloc3)
        call read_hdf5_al(dF,'dF',gloc3)
        call close_hdf5_group(gloc3)
        ne = size(comps(icomp2)%el)
        do ie =1,ne
          comps(icomp2)%el(ie)%mag   = idou(ie)
          comps(icomp2)%el(ie)%pres   = pres(ie)
          comps(icomp2)%el(ie)%dforce = dF(:,ie)
        enddo
        deallocate(idou, pres, dF)
        ! lifting line specific fields
        if (comps(icomp2)%comp_el_type .eq. 'l') then
          call open_hdf5_group(gloc2, 'Solution', gloc3)
          call read_hdf5_al(Gamma_old,'Gamma_old',gloc3)
          call read_hdf5_al(Gamma_old_old,'Gamma_old_old',gloc3)
          call read_hdf5_al(dGamma_dt,'dGamma_dt',gloc3)
          call close_hdf5_group(gloc3)
          do ie =1,ne
            select type( el => comps(icomp2)%el(ie) ) ; type is (t_liftlin)
              el%Gamma_old   = Gamma_old(ie)
              el%Gamma_old_old   = Gamma_old_old(ie)
              el%dGamma_dt = dGamma_dt(ie)
            end select
          enddo
          deallocate(Gamma_old, Gamma_old_old, dGamma_dt)
        ! vl specific fields
        else if(comps(icomp2)%comp_el_type .eq. 'v') then
          call open_hdf5_group(gloc2, 'Solution', gloc3)
          call read_hdf5_al(didou_dt,'didou_dt',gloc3)
          do ie =1,ne
            select type( el => comps(icomp2)%el(ie) ) ; type is (t_vortlatt)
              el%didou_dt = didou_dt(ie)
            end select
          enddo
          deallocate(didou_dt)
        ! panel elements specific fields
        else if(comps(icomp2)%comp_el_type .eq. 'p') then
          call open_hdf5_group(gloc2, 'Solution', gloc3)
          call read_hdf5_al(didou_dt,'didou_dt',gloc3)
          do ie =1,ne
            select type( el => comps(icomp2)%el(ie) ) ; type is (t_surfpan)
              el%didou_dt = didou_dt(ie)
            end select
          enddo
          deallocate(didou_dt)
        end if
        exit
      endif

    enddo

    !If the component was not found there is a problem
    if(icomp2 .gt. size(comps)) call error(this_sub_name, this_mod_name, &
    'Component loaded from geometry file not found in result file')

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)
  call close_hdf5_file(floc)

end subroutine load_solution

!----------------------------------------------------------------------

!> Chech that the component that we are about to load will be placed, in
!! the actual reference frames, in the same position in which it was saved
!!
!! This is done comparing offsets and rotations between the actual references
!! and the one saved
subroutine check_ref(gloc, floc, ref)
  integer(h5loc), intent(in)  :: gloc
  integer(h5loc), intent(in)  :: floc
  type(t_ref), intent(in)     :: ref

  integer                     :: iref
  integer(h5loc)              :: refs_gloc, ref_loc
  character(len=max_char_len) :: ref_title
  real(wp)                    :: R(3,3), of(3)
  character(len=*), parameter :: this_sub_name='check_ref'

  call read_hdf5(iref, 'RefId', gloc)
  call open_hdf5_group(floc, 'References', refs_gloc)
  write(ref_title,'(A,I3.3)')'Ref',iref
  call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)

  call read_hdf5(R,'R',ref_loc)
  call read_hdf5(of,'Offset',ref_loc)

  if ((.not. all(real(R) .eq. real(ref%R_g))) .or. &
      (.not. all(real(of) .eq. real(ref%of_g))) ) then
    call warning(this_sub_name, this_mod_name, 'Mismatching initial position &
    & while loading reference '//trim(ref%tag)//' This may lead to &
    &unexpected behaviour')
  endif

  call close_hdf5_group(ref_loc)
  call close_hdf5_group(refs_gloc)

end subroutine check_ref

!----------------------------------------------------------------------

end module mod_dust_io
