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

!> Module containing the subroutines to perform integral loads calculations
!! during postprocessing
module mod_post_integral

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi, ascii_real

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, &
  countoption

use mod_hdf5_io, only: &
  h5loc, &
  open_hdf5_file, &
  close_hdf5_file, &
  open_hdf5_group, &
  close_hdf5_group, &
  read_hdf5

use mod_stringtools, only: &
  LowCase, isInList, stricmp

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro

use mod_geometry, only: &
  t_geo, t_geo_component, destroy_elements

use mod_post_load, only: &
  load_refs, load_res , &
  check_if_components_exist

use mod_tecplot_out, only: &
  tec_out_loads

use mod_dat_out, only: &
  dat_out_loads_header, dat_out_hinge_header

use mod_math, only: &
  cross, vec2mat

implicit none

public :: post_integral, post_hinge_loads

private

character(len=*), parameter :: this_mod_name = 'mod_post_integral'
character(len=max_char_len) :: msg

contains

! ----------------------------------------------------------------------

subroutine post_integral( sbprms, basename, data_basename, an_name , ia , &
                          out_frmt, components_names, all_comp , &
                          an_start, an_end, an_step, average )
  type(t_parse), pointer                                   :: sbprms
  character(len=*), intent(in)                             :: basename
  character(len=*), intent(in)                             :: data_basename
  character(len=*), intent(in)                             :: an_name
  integer, intent(in)                                      :: ia
  character(len=*), intent(in)                             :: out_frmt
  character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
  logical, intent(in)                                      :: all_comp
  integer, intent(in)                                      :: an_start , an_end , an_step
  logical, intent(in)                                      :: average

  type(t_geo_component), allocatable                       :: comps(:)
  integer(h5loc)                                           :: floc , ploc
  real(wp), allocatable                                    :: points(:,:), points_virtual(:,:)
  integer                                                  :: nelem, nelem_virtual

  character(len=max_char_len)                              :: filename
  character(len=max_char_len), allocatable                 :: refs_tag(:)
  real(wp), allocatable                                    :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                    :: refs_G(:,:,:), refs_f(:,:)
  real(wp), allocatable                                    :: vort(:), cp(:)
  real(wp)                                                 :: R_cen(3,3)
  character(len=max_char_len)                              :: ref_tag
  integer                                                  :: ref_id
  real(wp)                                                 :: F_ref(3), F_bas(3), F_bas1(3)
  real(wp)                                                 :: M_ref(3), M_bas(3), ac(3)
  real(wp)                                                 :: F_ave(3), M_ave(3)
  real(wp), allocatable                                    :: force(:,:), moment(:,:)
  real(wp)                                                 :: P_inf , rho
  integer                                                  :: ic , it , ie , ierr , ires , fid_out , nstep
  real(wp), allocatable                                    :: time(:)
  real(wp)                                                 :: t

  character(len=*), parameter :: this_sub_name = 'post_integral'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' integral loads'//nl
  call printout(trim(msg))

  !debug
  if (all_comp) then
    call printout('  Analysing all components.')
  else
    call printout('  Analysing the following components:')
    do ic = 1 , size(components_names)
      write(msg,'(A,I0,A)') '   ',ic,') '//trim(components_names(ic))
      call printout(trim(msg))
    end do
  endif

  ! load the geo components just once just once
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  !TODO: here get the run id
  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
                                components_names,  all_comp)
  call close_hdf5_file(floc)

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  ! Reference system where the loads are projected:
  !  read the input and check if it exist
  ref_tag = getstr(sbprms,'reference_tag')

  write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_', an_start,'.h5'
  call open_hdf5_file(trim(filename),floc)
  call load_refs(floc, refs_R, refs_off, refs_G,refs_f,refs_tag)
  call close_hdf5_file(floc)

  ref_id = -333  ! check input
  do it = lbound(refs_tag,1) , ubound(refs_tag,1)
    if ( stricmp(refs_tag(it),  ref_tag) ) ref_id = it
  end do
  if ( ref_id .eq. -333 ) then
    call warning(this_sub_name, this_mod_name, 'Unknown reference system &
    &requested for the analysis, these are the valid reference frames: ')
    call printout('   Available references systems: ')
    do it = lbound(refs_tag,1) , ubound(refs_tag,1)
      write(msg,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
      call printout(trim(msg))
    end do
    call warning(this_sub_name, this_mod_name, 'Analysis skipped')
    return ! jump this analysis if the reference frame is not available
  end if

  nstep = (an_end-an_start)/an_step + 1
  ! Output format
  select case(trim(out_frmt))

    case('dat')
      ! Open output .dat file
      call new_file_unit(fid_out, ierr)
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.dat'
      open(unit=fid_out,file=trim(filename))
      call dat_out_loads_header( fid_out , components_names , ref_tag, average )

    case('tecplot')
      if (.not. average) then
        allocate(force(3,nstep), moment(3,nstep), time(nstep))
      else
        call error(this_sub_name, this_mod_name, 'Cannot output in tecplot format&
        & when averaging the results')
      endif

    case default
      call error('dust_post','','Unknown format '//trim(out_frmt)//&
                ' for integral_loads analysis. Choose dat or tecplot.')

  end select

  ires = 0
  if(average) then
    F_ave = 0.0_wp; M_ave = 0.0_wp
  endif

  do it=an_start, an_end, an_step ! Time loop
    ires = ires+1

    ! Open the file:
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'

    call open_hdf5_file(trim(filename),floc)

    ! Load u_inf --------------------------------
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(P_inf,'P_inf',ploc)
    call read_hdf5(rho,'rho_inf',ploc)
    call close_hdf5_group(ploc)

    ! Load the references and move the points ---
    call load_refs(floc,refs_R,refs_off)
    ! Move the points ---------------------------
    call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                filen = trim(filename) )
    ! Load the results --------------------------
    call load_res(floc, comps, vort, cp, t)

    call close_hdf5_file(floc)


    ! Initialise integral loads in the desired ref.frame
    F_ref = 0.0_wp ; M_ref = 0.0_wp

    ! Update the overall load with the comtribution from all the components
    do ic = 1 , size(comps)

      ! Initialise integral loads in the local ref.frame
      F_bas = 0.0_wp ; M_bas = 0.0_wp

      ! Loads from the ic-th component in the base ref.frame
      do ie = 1 , size(comps(ic)%el)

#if USE_PRECICE
        if (comps(ic)%coupling) then
          call vec2mat(comps(ic)%el(ie)%ori, R_cen)
          R_cen = matmul(comps(ic)%coupling_node_rot, R_cen)
          F_bas1 = comps(ic)%el(ie)%dforce 
          F_bas = F_bas + matmul(transpose(R_cen) , F_bas1 )
          if (trim(comps(ic)%comp_el_type) .eq. 'l') then 
              ac = sum ( comps(ic)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
              M_bas  = M_bas + cross( matmul(transpose(R_cen), ac), F_bas1) &
                            + comps(ic)%el(ie)%dmom 
            else 
              M_bas  = M_bas + cross( matmul(transpose(R_cen), &
                              comps(ic)%el(ie)%cen) , F_bas1 )      &
                              + comps(ic)%el(ie)%dmom 
          endif          

        else 

          F_bas1 = comps(ic)%el(ie)%dforce 
          F_bas = F_bas +  F_bas1

          if (trim(comps(ic)%comp_el_type) .eq. 'l') then
            ac = sum ( comps(ic)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
            M_bas = M_bas + cross( ac   &
                        - refs_off(:,ref_id) , F_bas1 )  &
                        + comps(ic)%el(ie)%dmom  
          else
            M_bas = M_bas + cross( comps(ic)%el(ie)%cen    &
                        - refs_off(:,ref_id) , F_bas1 )  &
                        + comps(ic)%el(ie)%dmom  
          endif
        endif 
#else   
        F_bas1 = comps(ic)%el(ie)%dforce 
        F_bas = F_bas +  F_bas1
        
        if (trim(comps(ic)%comp_el_type) .eq. 'l') then
          ac = sum ( comps(ic)%el(ie)%ver(:,1:2),2 ) / 2.0_wp
          M_bas = M_bas + cross( ac   &
                        - refs_off(:,ref_id) , F_bas1 )  &
                        + comps(ic)%el(ie)%dmom  
        else
          M_bas = M_bas + cross( comps(ic)%el(ie)%cen    &
                        - refs_off(:,ref_id) , F_bas1 )  &
                        + comps(ic)%el(ie)%dmom  
        endif
#endif
        
      end do !ie

      ! From the base ref.sys to the chosen ref.sys (offset and rotation)
#if USE_PRECICE
      if (comps(ic)%coupling) then 
        F_ref = F_ref + F_bas
        M_ref = M_ref + M_bas
      else
        F_ref = F_ref + matmul( &
            transpose( refs_R(:,:, ref_id) ) , F_bas )
        M_ref = M_ref + matmul( &
            transpose( refs_R(:,:, ref_id) ) , M_bas )
      endif 
#else
      F_ref = F_ref + matmul( &
            transpose( refs_R(:,:, ref_id) ) , F_bas )
      M_ref = M_ref + matmul( &
            transpose( refs_R(:,:, ref_id) ) , M_bas )
#endif
    end do !ic

    if(.not. average) then
      ! Update output dat file / update output arrays for tecplot
      select case(trim(out_frmt))

      case ('dat')
        write(fid_out, '('//ascii_real//')',advance='no') t
        write(fid_out,'(3'//ascii_real//')',advance='no') F_ref
        write(fid_out,'(3'//ascii_real//')',advance='no') M_ref
        write(fid_out,'(9'//ascii_real//')',advance='no') refs_R(:,:, ref_id)
        write(fid_out,'(3'//ascii_real//')',advance='no') refs_off(:, ref_id)
        write(fid_out,'(A)',advance='no') nl

      case('tecplot')
        time(ires) = t
        force(:,ires) = F_ref
        moment(:,ires) = M_ref

      end select
    else
      F_ave = F_ave*(real(ires-1,wp)/real(ires,wp)) + F_ref/real(ires,wp)
      M_ave = M_ave*(real(ires-1,wp)/real(ires,wp)) + M_ref/real(ires,wp)
    endif

  end do ! Time loop

  ! Close dat file / write tec file
  select case(trim(out_frmt))

  case('dat')
    if(.not. average) then
      close(fid_out)
    else
        write(fid_out,'(3'//ascii_real//')' ,advance='no') F_ave
        write(fid_out,'(3'//ascii_real//')' ,advance='no') M_ave
        write(fid_out,'(9'//ascii_real//')',advance='no') refs_R(:,:, ref_id)
        write(fid_out,'(3'//ascii_real//')',advance='no') refs_off(:, ref_id)
    endif

  case('tecplot')
    write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
    call tec_out_loads(filename, time, force, moment)
    deallocate(time, force, moment)

  end select

  !TODO: move deallocate(comps) outside this routine,
  !      because it is common to all the analyses
  call destroy_elements(comps)
  deallocate(comps,components_names)

  write(msg,'(A,I0,A)') nl//'++++++++++ Integral loads done'//nl
  call printout(trim(msg))

end subroutine post_integral

! ----------------------------------------------------------------------

subroutine post_hinge_loads( sbprms, basename, data_basename, an_name , ia , &
                            out_frmt, components_names, all_comp , hinge_tag, all_hinge, &
                            an_start, an_end, an_step, average )
  type(t_parse), pointer                                    :: sbprms
  character(len=*), intent(in)                              :: basename
  character(len=*), intent(in)                              :: data_basename
  character(len=*), intent(in)                              :: an_name
  integer, intent(in)                                       :: ia
  character(len=*), intent(in)                              :: out_frmt
  character(len=max_char_len), allocatable , intent(inout)  :: components_names(:)
  character(len=max_char_len), allocatable, intent(inout)   :: hinge_tag(:) 
  logical, intent(in)                                       :: all_comp, all_hinge 
  integer, intent(in)                                       :: an_start , an_end , an_step
  logical, intent(in)                                       :: average
  
  integer                                                   :: i_comp 
  type(t_geo_component), allocatable                        :: comps(:)
  integer(h5loc)                                            :: floc, ploc
  real(wp), allocatable                                     :: points(:,:), points_virtual(:,:)
  integer                                                   :: nelem, nelem_virtual, i_hinge
  character(len=max_char_len)                               :: filename
  character(len=max_char_len), allocatable                  :: refs_tag(:)
  real(wp), allocatable                                     :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                     :: refs_G(:,:,:), refs_f(:,:)
  real(wp), allocatable                                     :: vort(:), cp(:)
  character(len=max_char_len)                               :: ref_tag
  real(wp)                                                  :: F_ref(3), F_bas(3), F_bas1(3)
  real(wp)                                                  :: M_ref(3), M_bas(3)
  real(wp)                                                  :: M_ave(3), F_ave(3) 
  real(wp), allocatable                                     :: force(:,:), moment(:,:)
  real(wp)                                                  :: hinge_R(3,3) = 0.0_wp 
  real(wp)                                                  :: P_inf , rho
  integer                                                   :: ic , it , ie , ierr , ires , fid_out , nstep
  real(wp), allocatable                                     :: time(:)
  real(wp)                                                  :: t

  character(len=*), parameter                               :: this_sub_name = 'post_hinge_loads'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' hinge loads'//nl
  call printout(trim(msg))

  !debug
  if (all_comp) then
    call printout('  Analysing all components.')
  else
    call printout('  Analysing the following components:')
    do ic = 1 , size(components_names)
      write(msg,'(A,I0,A)') '   ',ic,') '//trim(components_names(ic))
      call printout(trim(msg))
    end do
  endif

  ! load the geo components just once just once
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  
  !TODO: here get the run id
  call load_components_postpro(comps, points, points_virtual, nelem, nelem_virtual, floc, &
        components_names,  all_comp)

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  ! Reference system where the loads are projected:
  !  read the input and check if it exist
  ref_tag = '0'

  write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_', an_start,'.h5'
  call open_hdf5_file(trim(filename),floc)
  call load_refs(floc,refs_R,refs_off,refs_G,refs_f,refs_tag)
  call close_hdf5_file(floc)

  ! Hinge tag where the hinge moment is calculated:
  !  read the input and check if it exist
  nstep = (an_end-an_start)/an_step + 1
  ! Output format
  do i_comp = 1, size(components_names)
    do i_hinge = 1, size(hinge_tag)
      select case(trim(out_frmt))

      case('dat')
        ! Open output .dat file
        call new_file_unit(fid_out, ierr)
        write(filename,'(A)') trim(basename)//'_'//trim(an_name)// & 
                              '_'//trim(components_names(i_comp))//'_'//& 
                              trim(hinge_tag(i_hinge))//'.dat'

        open(unit=fid_out,file=trim(filename))
        call dat_out_hinge_header( fid_out , components_names(i_comp) , hinge_tag(i_hinge), average )

      case('tecplot')
        if (.not. average) then
          allocate(force(3,nstep), moment(3,nstep), time(nstep))
        else
          call error(this_sub_name, this_mod_name, 'Cannot output in tecplot format&
          & when averaging the results')
        endif

      case default
        call error('dust_post','','Unknown format '//trim(out_frmt)//&
        ' for integral_loads analysis. Choose dat or tecplot.')

      end select
    enddo
  enddo 

  ires = 0
  if(average) then
    M_ave = 0.0_wp
    F_ave = 0.0_wp
  endif

  do it=an_start, an_end, an_step ! Time loop

    ires = ires + 1
    !> Open the file:
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'

    call open_hdf5_file(trim(filename),floc)

    ! Load u_inf --------------------------------
    call open_hdf5_group(floc,'Parameters',ploc)
    call read_hdf5(P_inf,'P_inf',ploc)
    call read_hdf5(rho,'rho_inf',ploc)
    call close_hdf5_group(ploc)

    ! Load the references and move the points ---
    ! Move the points ---------------------------
    call update_points_postpro(comps, points, points_virtual, refs_R, refs_off, &
                                filen = trim(filename) )
    ! Load the results --------------------------
    call load_res(floc, comps, vort, cp, t)
  
    call close_hdf5_file(floc)

    !> Initialise integral loads in the desired ref.frame
    F_ref = 0.0_wp ; M_ref = 0.0_wp

    !> Update the overall load with the comtribution from all the components
    do ic = 1 , size(comps)
      do i_hinge = 1, size(comps(ic)%hinge)
        if (stricmp(comps(ic)%hinge(i_hinge)%tag, hinge_tag(1)))then 
          !> Initialise integral loads in the local ref.frame
          F_bas = 0.0_wp ; M_bas = 0.0_wp
          !> Loads from the ic-th component in the base ref.frame for the rotated reagion 
          do ie = 1 , size(comps(ic)%hinge(i_hinge)%rot_cen%node_id, 1)
            F_bas1 =  comps(ic)%el(comps(ic)%hinge(i_hinge)%rot_cen%node_id(ie))%dforce 
            F_bas = F_bas + F_bas1
            M_bas = M_bas + cross( comps(ic)%el(comps(ic)%hinge(i_hinge)%rot_cen%node_id(ie))%cen &
                            - comps(ic)%hinge(i_hinge)%act%rr(:,1) , F_bas1 )  &
                            + comps(ic)%el(comps(ic)%hinge(i_hinge)%rot_cen%node_id(ie))%dmom  
          end do 

          !> Loads from the ic-th component in the base ref.frame for the blending reagion 
          do ie = 1 , size(comps(ic)%hinge(i_hinge)%blen_cen%node_id, 1)
            F_bas1 =  comps(ic)%el(comps(ic)%hinge(i_hinge)%blen_cen%node_id(ie))%dforce 
            F_bas = F_bas + F_bas1
            M_bas = M_bas + cross( comps(ic)%el(comps(ic)%hinge(i_hinge)%blen_cen%node_id(ie))%cen &
                            - comps(ic)%hinge(i_hinge)%act%rr(:,1) , F_bas1 )  &
                            + comps(ic)%el(comps(ic)%hinge(i_hinge)%blen_cen%node_id(ie))%dmom  
          end do 

          !> From the base ref.sys to the chosen ref.sys (offset and rotation)
          hinge_R(1,:) = comps(ic)%hinge(i_hinge)%act%v(:,1)
          hinge_R(2,:) = comps(ic)%hinge(i_hinge)%act%h(:,1)
          hinge_R(3,:) = comps(ic)%hinge(i_hinge)%act%n(:,1)
          F_ref = F_ref + matmul(hinge_R, F_bas )
          M_ref = M_ref + matmul(hinge_R, M_bas )

          if(.not. average) then
          ! Update output dat file / update output arrays for tecplot
            select case(trim(out_frmt))
            case ('dat')
              write(fid_out, '('//ascii_real//')',advance='no') t
              write(fid_out,'(3'//ascii_real//')',advance='no') F_ref
              write(fid_out,'(3'//ascii_real//')',advance='no') M_ref
              write(fid_out,'(9'//ascii_real//')',advance='no') hinge_R
              write(fid_out,'(3'//ascii_real//')',advance='no') comps(ic)%hinge(i_hinge)%act%rr(:,1) 
              write(fid_out,'(A)',advance='no') nl
            case('tecplot')
              time(ires) = t
              force(:,ires) = F_ref
              moment(:,ires) = M_ref
            end select
          else
            F_ave = F_ave*(real(ires-1,wp)/real(ires,wp)) + F_ref/real(ires,wp)
            M_ave = M_ave*(real(ires-1,wp)/real(ires,wp)) + M_ref/real(ires,wp)
          endif
        endif
      enddo 
    enddo 
  end do ! Time loop
  
  do ic = 1 , size(comps)
    do i_hinge = 1, size(comps(ic)%hinge)
      if (stricmp(comps(ic)%hinge(i_hinge)%tag, hinge_tag(1))) then 
        ! Close dat file / write tec file
        select case(trim(out_frmt))
        
        case('dat')
          if(.not. average) then
            close(fid_out)
          else
            write(fid_out,'(3'//ascii_real//')' ,advance='no') F_ave
            write(fid_out,'(3'//ascii_real//')' ,advance='no') M_ave
            write(fid_out,'(9'//ascii_real//')',advance='no') hinge_R
            write(fid_out,'(3'//ascii_real//')',advance='no') comps(ic)%hinge(i_hinge)%act%rr(:,1)
          endif
        
        case('tecplot')
          write(filename,'(A)')   trim(basename)//'_'//trim(an_name)// & 
                                  '_'//trim(components_names(i_comp))//'_'//& 
                                  trim(hinge_tag(i_hinge))//'.plt'
        
          call tec_out_loads(filename, time, force, moment)
          deallocate(time, force, moment)
        
        end select
      endif 
    enddo
  enddo
  call destroy_elements(comps)
  deallocate(comps, components_names, hinge_tag)

  write(msg,'(A,I0,A)') nl//'++++++++++ Hinge loads done'//nl
  call printout(trim(msg))


end subroutine post_hinge_loads
! ----------------------------------------------------------------------

end module mod_post_integral
