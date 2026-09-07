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

!> Module containing the subroutines to perform flowfield visualizations
!! during postprocessing
module mod_post_flowfield_test

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi, one_4pi

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, printout

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, getintarray, getrealarray, &
  countoption

use mod_wake_test, only: &
  t_wake

use mod_hdf5_io, only: &
  h5loc, &
  open_hdf5_file, &
  close_hdf5_file, &
  open_hdf5_group, &
  close_hdf5_group, &
  read_hdf5

use mod_stringtools, only: &
  LowCase

use mod_tecplot_out, only: &
  tec_out_box

use mod_vtk_out_test, only: &
  vtr_write

use mod_post_load_test, only: &
  load_wake_post

use mod_wind, only: &
  variable_wind
  
use mod_aeroel, only: &
  c_elem, c_vort_elem,  &
  t_elem_p, t_vort_elem_p 


implicit none

public :: post_flowfield

private

character(len=*), parameter :: this_mod_name = 'mod_post_flowfield_test'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_flowfield( sbprms, basename, data_basename, an_name, ia, &
                            out_frmt, an_start, an_end, an_step, average )
                            
  type(t_parse), pointer                                   :: sbprms
  character(len=*) , intent(in)                            :: basename
  character(len=*) , intent(in)                            :: data_basename
  character(len=*) , intent(in)                            :: an_name
  integer          , intent(in)                            :: ia
  character(len=*) , intent(in)                            :: out_frmt
  integer , intent(in)                                     :: an_start , an_end , an_step
  logical, intent(in)                                      :: average

  integer , parameter                                      :: n_max_vars = 3 !vel,p,vort, ! TODO: 4 with cp
  character(len=max_char_len), allocatable                 :: var_names(:)
  integer , allocatable                                    :: vars_n(:)
  integer                                                  :: i_var , n_vars
  logical                                                  :: probe_vel , probe_p , probe_vort, probe_cp
  real(wp)                                                 :: u_inf(3)
  real(wp)                                                 :: P_inf , rho
  real(wp)                                                 :: vel_probe(3) = 0.0_wp 
  real(wp)                                                 :: vort_probe(3) = 0.0_wp
  real(wp)                                                 :: pres_probe
  real(wp)                                                 :: pot_probe, pot_probe_old
  real(wp)                                                 :: v(3) = 0.0_wp
  real(wp)                                                 :: phi = 0.0_wp

  integer(h5loc)                                           :: floc , ploc

  real(wp), allocatable                                    :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                    :: refs_G(:,:,:), refs_f(:,:)
  real(wp), allocatable                                    :: vort(:), cp(:)
  real(wp), allocatable                                    :: vort_old(:), cp_old(:) ! TODO needed?

  type(t_wake)                                             :: wake, wake_old
  type(t_elem_p), allocatable                              :: wake_elems(:), wake_elems_old(:)

  real(wp)                                                 :: t, t_old

  integer                                                  :: nxyz(3)
  real(wp)                                                 :: minxyz(3), maxxyz(3)
  real(wp), allocatable                                    :: xbox(:), ybox(:), zbox(:)
  real(wp)                                                 :: dxbox, dybox, dzbox
  real(wp), allocatable                                    :: box_vel(:,:), box_p(:), box_vort(:,:)
  integer                                                  :: ix , iy , iz
  real(wp), allocatable                                    :: vars(:,:)
  real(wp), allocatable                                    :: ave_vars(:,:)
  integer                                                  :: i_var_v , i_var_p , i_var_w

  integer                                                  :: ip , ipp, ic , ie , i1 , it, itave, nelems_comp
  
  logical                                                  :: ex
  real(wp)                                                 :: dummy 

  character(len=max_char_len)                              :: str_a, var_name
  character(len=max_char_len)                              :: filename
  character(len=*), parameter                              :: this_sub_name = 'post_flowfield'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' flowfield'//nl
  call printout(trim(msg))


  ! Read variables to save : velocity | pressure | vorticity
  probe_vel = .false. ; probe_p = .false. ; probe_vort = .false. ; probe_cp = .false.
  n_vars = countoption(sbprms,'variable')

  if ( n_vars .eq. 0 ) then ! default: velocity | pressure | vorticity
    probe_vel = .true. ; probe_p = .true. ; probe_vort = .true.
  else
    do i_var = 1 , n_vars
      var_name = getstr(sbprms,'variable')
      call LowCase(var_name)
      select case(trim(var_name))
        case ( 'velocity' ) ; probe_vel = .true.
        case ( 'pressure' ) ; probe_p   = .true.
        case ( 'vorticity') ; probe_vort= .true.
        case ( 'cp'       ) ; probe_cp = .true. ! TODO behaviour of probe cp with 'all'?
        case ( 'all') ; probe_vel = .true. ; probe_p   = .true. ; probe_vort= .true. ; probe_cp = .false.
        case default
          write(str_a,*) ia
          call error('dust_post','','Unknown Variable: '//trim(var_name)//&
                    ' for analysis n.'//trim(str_a)//'.'//nl//&
                    'Choose "velocity", "pressure", "vorticity".')
      end select
    end do
  end if
 

  ! Read box dimensions ...
  nxyz   = getintarray( sbprms,'n_xyz'  ,3)
  minxyz = getrealarray(sbprms,'min_xyz',3)
  maxxyz = getrealarray(sbprms,'max_xyz',3)

  ! ... allocate 'box' and deal with inconsistent input
  if ( nxyz(1) .gt. 1 ) then ! x-coord
    allocate( xbox(nxyz(1)) )
    dxbox = ( maxxyz(1) - minxyz(1) ) / real(nxyz(1) - 1, wp)
    xbox = (/ ( minxyz(1) + real(i1-1,wp) * dxbox , i1 = 1 , nxyz(1) )/)
  else
    allocate( xbox(1) )
    xbox(1) = minxyz(1)
  end if
  if ( nxyz(2) .gt. 1 ) then ! y-coord
    allocate( ybox(nxyz(2)) )
    dybox = ( maxxyz(2) - minxyz(2) ) / real(nxyz(2) - 1, wp)
    ybox = (/ ( minxyz(2) + real(i1-1, wp) * dybox , i1 = 1 , nxyz(2) )/)
  else
    allocate( ybox(1) )
    ybox(1) = minxyz(2)
  end if
  if ( nxyz(3) .gt. 1 ) then ! z-coord
    allocate( zbox(nxyz(3)) )
    dzbox = ( maxxyz(3) - minxyz(3) ) / real(nxyz(3) - 1, wp)
    zbox = (/ ( minxyz(3) + real(i1-1, wp) * dzbox , i1 = 1 , nxyz(3) )/)
  else
    allocate( zbox(1) )
    zbox(1) = minxyz(3)
  end if

  ! Allocate box_vel, box_p, box_vort, ...
  allocate(var_names(n_max_vars)) ; var_names = ' '
  allocate(vars_n   (n_max_vars)) ; vars_n = 0
  i_var = 0
  i_var_v = 0 ; i_var_p = 0 ; i_var_w = 0
  if ( probe_vel ) then
    allocate(box_vel (product(nxyz),3))
    i_var = i_var + 1
    var_names(i_var) = 'velocity'
    vars_n(i_var) = 3
    i_var_v = 3
  end if
  if ( probe_p   ) then
    allocate(box_p   (product(nxyz)  ))
    i_var = i_var + 1
    if (probe_cp) then
      var_names(i_var) = 'cp'
    else
      var_names(i_var) = 'pressure'
    end if  
    vars_n(i_var) = 1
    i_var_p = i_var_v + 1
  end if
  if ( probe_vort) then
    allocate(box_vort(product(nxyz),3))
    i_var = i_var + 1
    var_names(i_var) = 'vorticity'
    vars_n(i_var) = 3
    i_var_w = i_var_p + 3
  end if

  ! Allocate and fill vars array, for output +++++++++++++++++
  !  sum(vars_n): # of scalar fields to be plotted
  !  product(nxyz) # of points where the vars are plotted
  allocate(vars(sum(vars_n),product(nxyz))) ; vars = 0.0_wp
  if(average) then
    allocate(ave_vars(sum(vars_n),product(nxyz)))
    ave_vars = 0.0_wp
  endif

  write(*,'(A,I0,A,I0,A,I0)') nl//' it_start,it_end,an_step : ' , &
    an_start , ' , ' , an_end , ' , ' , an_step
  itave = 0
  do it = an_start, an_end, an_step ! Time history
    itave = itave + 1

    ! Show timing, since this analysis is quite slow
    write(*,'(A,I0,A,I0)') ' it : ' , itave , ' / ' , &
      ( an_end - an_start + 1 ) / an_step

    ! Open the result file ----------------------
    write(filename,'(A,I4.4,A)') trim(data_basename)// &
                                          '_res_',it,'.h5'
    call open_hdf5_file(trim(filename),floc)

    ! Load u_inf --------------------------------
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(u_inf,'u_inf',ploc)
    call read_hdf5(P_inf,'P_inf',ploc)
    call read_hdf5(rho,'rho_inf',ploc)
    call close_hdf5_group(ploc)

    sim_param%u_inf = u_inf
    sim_param%P_inf = P_inf
    sim_param%rho_inf = rho


    ! Load the wake -----------------------------
    call load_wake_post(floc, wake)
    call read_hdf5(t,'time',floc)
    call close_hdf5_file(floc)
    
  !> Compute fields to be plotted +++++++++++++++++++++++++++++
!$omp parallel do collapse(3) private(iz,iy,ix,vel_probe, pres_probe, vort_probe, &
                                      !$omp ic, ie, v, ipp, pot_probe, pot_probe_old, phi)
    ! Loop over the nodes of the box
    do iz = 1 , size(zbox)  ! z-coord
      do iy = 1 , size(ybox)  ! y-coord
        do ix = 1 , size(xbox)  ! x-coord

          ipp = ix + (iy-1)*size(xbox) + (iz-1)*size(xbox)*size(ybox)

          if ( probe_vel .or. probe_p ) then

            ! Compute velocity
            vel_probe = 0.0_wp ; pres_probe = 0.0_wp ; vort_probe = 0.0_wp ; pot_probe = 0.0_wp; pot_probe_old = 0.0_wp
           
            !> wake
            do ie = 1, size(wake%wake_parts)
              call wake%wake_parts(ie)%compute_vel((/ xbox(ix), ybox(iy), zbox(iz)/), v)
              vel_probe = vel_probe + v
            enddo

            !> + u_inf
            vel_probe = vel_probe*one_4pi + variable_wind((/ xbox(ix) , ybox(iy) , zbox(iz) /), t)
          end if

          if ( probe_vel ) then
            vars(1:3,ipp) = vel_probe
          end if
        end do  ! x-coord
      end do  ! y-coord
    end do  ! z-coord
!$omp end parallel do

    !> Output or average +++++++++++++++++++++++++++++++++++++++++++++++++
    if (average) then
      ave_vars = ave_vars*(real(itave-1,wp)/real(itave,wp)) + &
                      vars/real(itave,wp)

    else
      select case (trim(out_frmt))

      case ('vtk')
        write(filename,'(A,I4.4,A)') trim(basename)//'_'//&
                                      trim(an_name)//'_',it,'.vtr'
        call vtr_write ( filename , xbox , ybox , zbox , &
                        vars_n(1:i_var) , var_names(1:i_var) , &
                        vars )
      case('tecplot')
        write(filename,'(A,I4.4,A)') trim(basename)//'_'//&
                                trim(an_name)//'_',it,'.plt'
        i_var = 0
        deallocate(var_names)
        allocate(var_names(7))
        if(probe_vel) then
          var_names(i_var + 1) = 'ux'
          var_names(i_var + 2) = 'uy'
          var_names(i_var + 3) = 'uz'
          i_var = i_var + 3
        endif
        if(probe_p) then
          var_names(i_var + 1) = 'p'
          i_var = i_var + 1
        endif
        if(probe_vort) then
          var_names(i_var + 1) = 'omx'
          var_names(i_var + 2) = 'omy'
          var_names(i_var + 3) = 'omz'
          i_var = i_var + 3
        endif

        call tec_out_box(filename, t, xbox, ybox, zbox, &
                        vars, var_names(1:i_var))

      case default
        call error('dust_post','','Unknown format '//trim(out_frmt)//&
                    ' for flowfield output. Choose: vtk or tecplot.')
      end select
    endif !average or not
  end do  ! Time loop

  !output if average
  if (average) then
    select case (trim(out_frmt))    
    case ('vtk')
      write(filename,'(A)') trim(basename)//'_'//&
                            trim(an_name)//'_ave.vtr'
      call vtr_write (filename , xbox , ybox , zbox , &
                      vars_n(1:i_var) , var_names(1:i_var) , &
                      ave_vars )
    case('tecplot')
      write(filename,'(A)') trim(basename)//'_'//&
                              trim(an_name)//'_ave.plt'
      i_var = 0
      deallocate(var_names)
      allocate(var_names(7))
      if(probe_vel) then
        var_names(i_var + 1) = 'ux'
        var_names(i_var + 2) = 'uy'
        var_names(i_var + 3) = 'uz'
        i_var = i_var + 3
      endif
      if(probe_p) then
        var_names(i_var + 1) = 'p'
        i_var = i_var + 1
      endif
      if(probe_vort) then
        var_names(i_var + 1) = 'omx'
        var_names(i_var + 2) = 'omy'
        var_names(i_var + 3) = 'omz'
        i_var = i_var + 3
      endif
      call tec_out_box(filename, t, xbox, ybox, zbox, &
                      ave_vars, var_names(1:i_var))
    case default
      call error( 'dust_post','','Unknown format '//trim(out_frmt)//&
                  ' for flowfield output. Choose: vtk or tecplot.')
    end select
  
  endif

  if ( allocated(box_vel ) ) deallocate(box_vel )
  if ( allocated(box_p   ) ) deallocate(box_p   )
  if ( allocated(box_vort) ) deallocate(box_vort)
  !if ( allocated(box_cp  ) ) deallocate(box_cp  )
  deallocate(xbox,ybox,zbox)
  deallocate(var_names,vars_n)

  write(msg,'(A,I0,A)') nl//'++++++++++ Flowfield done'//nl
  call printout(trim(msg))

end subroutine post_flowfield

! ----------------------------------------------------------------------

end module mod_post_flowfield_test
