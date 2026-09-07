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

module mod_dat_out

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len, ascii_real

use mod_handling, only: &
  error, warning, internal_error, new_file_unit

!---------------------------------------------------------------------
implicit none

public :: dat_out_probes_header, dat_out_loads_header, dat_out_hinge_header, dat_out_aa_header, &
          dat_out_sectional, dat_out_sectional_ll, dat_out_sectional_vl, dat_out_aa, &
          dat_out_chordwise 

private

character(len=*), parameter :: &
  this_mod_name = 'mod_dat_output'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine dat_out_loads_header ( fid , comps_meas , ref_sys , average)
  integer , intent(in)          :: fid
  character(len=*), intent(in)  :: comps_meas(:)
  character(len=*), intent(in)  :: ref_sys
  logical, intent(in)           :: average

  character(len=max_char_len)   :: istr
  integer :: n_comps , ic

  n_comps = size(comps_meas)

  write(istr,'(I0)') n_comps
  write(fid,*) '# Integral loads: N.components: ' , trim(istr)
  write(fid,*) '#                 Ref.sys     : ' , trim(ref_sys)
  
  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #                 Components  : '
  
  do ic = 1 , n_comps - 1
    write(fid,'(A)',advance='no') trim(comps_meas(ic))//' , '
  end do

  write(fid,'(A)') trim(comps_meas(n_comps))
  
  if (.not. average)  then
    write(fid,*) '#  t , Fx , Fy , Fz , Mx , My , Mz , ref_mat(9) , ref_off(3) '
  else
    write(fid,*) '# Fx_average , Fy_average , Fz_average ,&
                  & Mx_average , My_average , Mz_average ,&
                  & ref_mat(9) , ref_off(3) '
  endif

end subroutine dat_out_loads_header

subroutine dat_out_hinge_header ( fid , comps_meas , hinge_tag , average)
  integer , intent(in)          :: fid
  character(len=*), intent(in)  :: comps_meas
  character(len=*), intent(in)  :: hinge_tag
  logical, intent(in)           :: average

  write(fid,*) '# Hinge Moment: '
  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #               Components  : '
  write(fid,'(A)') trim(comps_meas)

  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #                    Hinge  : '
  write(fid,'(A)') trim(hinge_tag)

  if (.not. average)  then
    write(fid,*) '#  t , Fv , Fh , Fn , Mv , Mh , Mn , axis_mat(9) , node_hinge(3) '
  else
    write(fid,*) '# Fv_average , Fh_average , Fn_average ,&
                  & Mv_average , Mh_average , Mn_average ,&
                  & axis_mat(9) , node_hinge(3) '
  endif

end subroutine dat_out_hinge_header

!---------------------------------------------------------------------

subroutine dat_out_aa_header ( fid , num_el, t, p_inf, rho_inf, a_inf, mu_inf, u_inf)
  integer , intent(in) :: fid
  real(wp), intent(in) :: t
  real(wp), intent(in) :: p_inf, rho_inf, a_inf, mu_inf
  real(wp), intent(in) :: u_inf(3)
  integer,  intent(in) :: num_el

  write(fid,'(A)') '# Aeroacoustic Data'

  write(fid,'(A)') '# Element, Time, Pressure, Density, Sound speed, Dynamic viscosity, Flow velocity'
  write(fid,'(I0,8'//ascii_real//')') num_el, t, p_inf, rho_inf, a_inf, mu_inf, u_inf

  write(fid,'(A)') '# cx, cy, cz, nx, ny, nz, area, rho, pressure, rhoux, rhouy, rhouz, '// &
                   'svx, svy, svz, cx_u, cy_u, cz_u, cx_l, cy_l, cz_l, nx_u, ny_u, nz_u, '// &
                   'nx_l, ny_l, nz_l, area_u, area_l, svx_u, svy_u, svz_u, svx_l, svy_l, svz_l'
  
end subroutine dat_out_aa_header

!---------------------------------------------------------------------

subroutine dat_out_aa ( fid , cen, n, vel, area, pres, rho_inf, &
                        cen_u, cen_l, n_u, n_l, area_u, area_l, vel_u, vel_l)
  integer , intent(in) :: fid
  real(wp), intent(in) :: cen(3), n(3), vel(3), area, pres, rho_inf
  real(wp), intent(in), optional :: cen_u(3), cen_l(3), n_u(3), n_l(3)
  real(wp), intent(in), optional :: area_u, area_l
  real(wp), intent(in), optional :: vel_u(3), vel_l(3)

  real(wp) :: cen_u_loc(3), cen_l_loc(3), n_u_loc(3), n_l_loc(3)
  real(wp) :: vel_u_loc(3), vel_l_loc(3)
  real(wp) :: area_u_loc, area_l_loc

  cen_u_loc = 0.0_wp
  cen_l_loc = 0.0_wp
  n_u_loc = 0.0_wp
  n_l_loc = 0.0_wp
  vel_u_loc = 0.0_wp
  vel_l_loc = 0.0_wp
  area_u_loc = 0.0_wp
  area_l_loc = 0.0_wp

  if (present(cen_u)) cen_u_loc = cen_u
  if (present(cen_l)) cen_l_loc = cen_l
  if (present(n_u)) n_u_loc = n_u
  if (present(n_l)) n_l_loc = n_l
  if (present(vel_u)) vel_u_loc = vel_u
  if (present(vel_l)) vel_l_loc = vel_l
  if (present(area_u)) area_u_loc = area_u
  if (present(area_l)) area_l_loc = area_l

  write(fid,'(36'//ascii_real//')') cen, n, area, rho_inf, pres, vel*rho_inf, vel, &
                                    cen_u_loc, cen_l_loc, n_u_loc, n_l_loc, &
                                    area_u_loc, area_l_loc, vel_u_loc, vel_l_loc

end subroutine dat_out_aa

!---------------------------------------------------------------------

subroutine dat_out_probes_header ( fid , rr_probes , vars_str, nt )
  integer , intent(in)         :: fid
  real(wp), intent(in)         :: rr_probes(:,:)
  character(len=*), intent(in) :: vars_str
  integer                      :: nt

  character(len=max_char_len) :: istr
  integer :: n_probes , ic
  character(len=8) :: nnum

  n_probes = size(rr_probes,2)

  if ( size(rr_probes,1) .ne. 3 ) then
    call error(trim(this_mod_name),'','Wrong format of the rr_probes inputs.&
            & size(rr_probes,1) .ne. 3. Stop ')
  end if

  write(fid,*) '# N. of point probes:' , n_probes
  !> three-dimensional space
  do ic = 1 , 3
    write(nnum,'(I0)') size(rr_probes,2)
    write(fid,'('//trim(nnum)//ascii_real//')') rr_probes(ic,:)
  end do

  write(istr,'(I0)') n_probes
  write(fid,'(A,I0)') '# n_time: ', nt
  write(fid,'(A,A,A,A,A)') '#    t     '//trim(istr)//' ('//trim(vars_str)//' )'

end subroutine dat_out_probes_header

!---------------------------------------------------------------------

subroutine dat_out_sectional (basename, compname, y_cen, y_span, chord, time, &
                              sec_loads, ref_mat, off_mat, average )
  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in) :: y_cen(:)
  real(wp) , intent(in) :: y_span(:)
  real(wp) , intent(in) :: chord(:)
  real(wp) , intent(in) :: time(:)
  real(wp) , intent(in) :: sec_loads(:,:,:)
  real(wp) , intent(in) :: ref_mat(:,:)
  real(wp) , intent(in) :: off_mat(:,:)
  logical,   intent(in) :: average

  character(len=2) :: load_str(4)
  character(len=8) :: nnum
  character(len=max_char_len) :: filename
  integer :: it , nt , fid , i1

  load_str = (/ 'Fx' , 'Fy' , 'Fz' , 'Mo' /)

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(sec_loads,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,2) .ne. size(y_cen). Stop ')
  end if
  
  if ( size(sec_loads,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,1) .ne. size(time). Stop ')
  end if
  
  if ( size(sec_loads,3) .ne. 4 ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,3) .ne. 4. Stop ')
  end if

  ! Print out .dat files
  fid = 21
  do i1 = 1 , 4
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(i1))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(i1))//'.dat'
    endif

    
    open(unit=fid,file=trim(filename))
    ! Header -----------
    write(fid,'(A)') '# Sectional load '//trim(load_str(i1))//&
                    &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(sec_loads,2) , ' ; n_time : ' , nt , & 
                                          &'. Next lines: y_cen , y_span, chord'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span
    write(fid,'('//trim(nnum)//ascii_real//')') chord

    if(average) then
      write(fid,'(A)') '#sec(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') sec_loads(1,:,i1)
    else
      write(fid,'(A)') '# t , sec(n_sec) , ref_mat(9) , ref_off(3) '
      ! Dump data --------
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)+9+3
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), &
                          sec_loads(it,:,i1) , ref_mat(it,:) , off_mat(it,:)
      end do
    endif
    close(fid)

  end do


end subroutine dat_out_sectional

!---------------------------------------------------------------------

subroutine dat_out_sectional_ll (basename, compname, y_cen, y_span, chord, time, &
                                ll_sec, average )
  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in)         :: y_cen(:)
  real(wp) , intent(in)         :: y_span(:)
  real(wp) , intent(in)         :: chord(:)  
  real(wp) , intent(in)         :: time(:)
  real(wp) , intent(in)         :: ll_sec(:,:,:)
  logical,   intent(in)         :: average

  character(len=8)              :: nnum
  character(len=max_char_len)   :: filename
  integer                       :: it , nt , fid , ierr, il
  character(len=*), parameter   :: this_sub_name = 'dat_out_sectional_ll'
  character(len=21)             :: load_str(9)
  character(len=max_char_len)   :: description_str(9)

  load_str(1) = 'Cl' 
  load_str(2) = 'Cd' 
  load_str(3) = 'Cm'
  load_str(4) = 'alpha'
  load_str(5) = 'vel_2d'
  load_str(6) = 'up_x'
  load_str(7) = 'up_y'
  load_str(8) = 'up_z'
  load_str(9) = 'vel_outplane'
  
  description_str(1) = 'lift coefficient'
  description_str(2) = 'drag coefficient'
  description_str(3) = 'moment coefficient'
  description_str(4) = 'angle of attack'
  description_str(5) = 'in section plane velocity (module)'
  description_str(6) = 'in section plane velocity x component'
  description_str(7) = 'in section plane velocity y component'
  description_str(8) = 'in section plane velocity z component'
  description_str(9) = 'out of section (spanwise) velocity'

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(ll_sec,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & different length of nodes and solution ')
  end if
  if ( size(ll_sec,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,1) .ne. size(time). Stop ')
  end if

  !> Print out .dat files
  do il = 1 , size(load_str)
    call new_file_unit(fid, ierr)
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'.dat'
    endif
    open(unit=fid,file=trim(filename))
    !> Header 
    write(fid,'(A)') '# Sectional '//trim(description_str(il))//&
                    &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(y_cen) , ' ; n_time : ' , nt , & 
                    &'. Next lines: y_cen , y_span, chord'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span
    write(fid,'('//trim(nnum)//ascii_real//')') chord 
    if(average) then
      write(fid,'(A)') '# '//trim(load_str(il))//'(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') ll_sec(1,:,il)
    else
      write(fid,'(A)') '# t , '//trim(load_str(il))//'(n_sec) '
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), ll_sec(it,:,il)
      end do
    endif
    close(fid)
  enddo


end subroutine dat_out_sectional_ll


subroutine dat_out_sectional_vl (basename, compname, y_cen, y_span, chord, time, &
                                vl_sec, average )

  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in)         :: y_cen(:)
  real(wp) , intent(in)         :: y_span(:)
  real(wp) , intent(in)         :: chord(:)
  
  real(wp) , intent(in)         :: time(:)
  real(wp) , intent(in)         :: vl_sec(:,:,:)
  logical,   intent(in)         :: average

  character(len=8)              :: nnum
  character(len=max_char_len)   :: filename
  integer                       :: it , nt , fid , ierr, il
  character(len=*), parameter   :: this_sub_name = 'dat_out_sectional_vl'
  character(len=21)             :: load_str(12)
  character(len=max_char_len)   :: description_str(12)

  load_str(1) = 'Cl'; 
  load_str(2) = 'Cd'; 
  load_str(3) = 'Cm'
  load_str(4) = 'alpha'; 
  load_str(5) = 'alpha_isolated'
  load_str(6) = 'vel_2d'; 
  load_str(7) = 'up_x'
  load_str(8) = 'up_y'
  load_str(9) = 'up_z'
  load_str(10) = 'vel_2d_isolated'
  load_str(11) = 'vel_outplane'; 
  load_str(12) = 'vel_outplane_isolated'

  description_str(1) = 'lift coefficient'
  description_str(2) = 'drag coefficient'
  description_str(3) = 'moment coefficient'
  description_str(4) = 'angle of attack'
  description_str(5) = 'isolated angle of attack'
  description_str(6) = 'in section plane velocity (module)'
  description_str(7) = 'in section plane velocity x component'
  description_str(8) = 'in section plane velocity y component'
  description_str(9) = 'in section plane velocity z component'
  description_str(10) = 'in section plane velocity isolated'
  description_str(11) = 'out of section (spanwise) velocity'
  description_str(12) = 'out of section (spanwise) isolated velocity'

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(vl_sec,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
    & different length of nodes and solution ')
  end if
  if ( size(vl_sec,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
    & size(sec_loads,1) .ne. size(time). Stop ')
  end if

  ! Print out .dat files
  do il = 1 , size(load_str)
    call new_file_unit(fid, ierr)
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'.dat'
    endif
    open(unit=fid,file=trim(filename))
    ! Header -----------
    write(fid,'(A)') '# Sectional '//trim(description_str(il))// &
                                &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(y_cen) , ' ; n_time : ' , nt , &
                                &'. Next lines: y_cen , y_span, chord'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span
    write(fid,'('//trim(nnum)//ascii_real//')') chord
    
    if(average) then
      write(fid,'(A)') '# '//trim(load_str(il))//'(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') vl_sec(1,:,il)
    else
      write(fid,'(A)') '# t , '//trim(load_str(il))//'(n_sec) '
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), vl_sec(it,:,il)
      end do
    endif
    close(fid)
  enddo

end subroutine dat_out_sectional_vl

!---------------------------------------------------------------------
subroutine dat_out_chordwise (basename, compname, time, &
                              force_int, tang_int, nor_int, cen_int, pres_int, & 
                              cp_int, average, n_station, station, chord_length)
  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp), intent(in)          :: time(:)
  real(wp), intent(in)          :: force_int(:,:,:,:)
  real(wp), intent(in)          :: tang_int(:,:,:,:)
  real(wp), intent(in)          :: nor_int(:,:,:,:)
  real(wp), intent(in)          :: cen_int(:,:,:,:)
  real(wp), intent(in)          :: pres_int(:,:,:)
  real(wp), intent(in)          :: cp_int(:,:,:)
  logical,  intent(in)          :: average
  integer,  intent(in)          :: n_station 
  real(wp), intent(in)          :: station(:) 
  real(wp), intent(in)          :: chord_length(:)  
  
  character(len=8)              :: nnum
  character(len=max_char_len)   :: filename
  integer                       :: it, ista, icase, nt, fid 
  character(len=5)              :: load_str(10)
  character(len=max_char_len)   :: description_str(10)

  load_str(1) = 'Pres' 
  load_str(2) = 'Cp'
  load_str(3) = 'dFx'
  load_str(4) = 'dFz'
  load_str(5) = 'dNx'
  load_str(6) = 'dNz'
  load_str(7) = 'dTx'
  load_str(8) = 'dTz' 
  load_str(9) = 'x_cen'
  load_str(10) = 'z_cen'
  
  description_str(1) = 'Panel pressure'
  description_str(2) = 'Panel coefficient of pressure'
  description_str(3) = 'Panel force per unit length in chordwise direction'
  description_str(4) = 'Panel force per unit length in flapwise direction'
  description_str(5) = 'Panel local normal in chordwise direction'
  description_str(6) = 'Panel local normal in flapwise direction'
  description_str(7) = 'Panel local tangent in chordwise direction'
  description_str(8) = 'Panel local tangent in flapwise direction'
  description_str(9) = 'Panel center chordwise coordinate'
  description_str(10) = 'Panel center flapwise coordinate'
  

  nt = size(time)

  ! Print out .dat files
  fid = 21
  do icase = 1, 10
    do ista = 1, n_station 
      if(average) then
        write(filename,'(A,I0,A)') trim(basename)//'_', ista, '_'//trim(load_str(icase))//'_ave.dat'
      else
        write(filename,'(A,I0,A)') trim(basename)//'_', ista, '_'//trim(load_str(icase))//'.dat'
      endif

      open(unit=fid,file=trim(filename))
      ! Header -----------
      write(fid,'(A)') '# Chordwise load of component: '//trim(compname)
      write(fid,'(A,F5.3,A,F5.3)') '# spanwise_location: ', station(ista), ' ; chord_length: ',  &
                                    chord_length(ista)
      write(fid,'(A,I0,A,I0,A)') '# n_chord : ' ,size(cen_int,3) , ' ; n_time : ' , nt ,& 
                                & '. Next lines: x_chord , z_chord'
      write(nnum,'(I0)') size(cen_int,3)
      write(fid,'('//trim(nnum)//ascii_real//')') cen_int(1,ista,:,1)
      write(fid,'('//trim(nnum)//ascii_real//')') cen_int(1,ista,:,3)
    
      select case(trim(load_str(icase)))
        case('Pres')
          write(fid,'(A)') '# t, Pres'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1 + size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), pres_int(it,ista,:) 
          end do
        case('Cp')
          write(fid,'(A)') '# t, Cp'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1 + size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), cp_int(it,ista,:) 
          end do
        case('dFx')
          write(fid,'(A)') '# t, dFz'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), force_int(it,ista,:,1)            
          end do
        case('dFz')
          write(fid,'(A)') '# t, dFz'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), force_int(it,ista,:,3)            
          end do
        case('dNx')
          write(fid,'(A)') '# t, dNx'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), nor_int(it,ista,:,1)            
          end do
        case('dNz')
          write(fid,'(A)') '# t, dNz'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), nor_int(it,ista,:,3)            
          end do
        case('dTx')
            write(fid,'(A)') '# t, dTx'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), tang_int(it,ista,:,1)            
          end do
        case('dTz')
          write(fid,'(A)') '# t, dTz'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), tang_int(it,ista,:,3)            
          end do
        case('x_cen')
          write(fid,'(A)') '# t, x_cen'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), cen_int(it,ista,:,1)      
          end do
        case('z_cen')
          write(fid,'(A)') '# t, z_cen'
          ! Dump data --------
          do it = 1 , nt
            write(nnum,'(I0)') 1+size(cen_int,3)
            write(fid,'('//trim(nnum)//ascii_real//')') time(it), cen_int(it,ista,:,3)      
          end do
        end select
      close(fid)
    end do 
  end do


end subroutine dat_out_chordwise 

end module mod_dat_out
