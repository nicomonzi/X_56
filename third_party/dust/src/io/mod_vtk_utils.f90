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

module mod_vtk_utils

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, internal_error, warning, info, printout, new_file_unit


!---------------------------------------------------------------------
implicit none

public :: t_output_var, vtk_isize, vtk_fsize, add_output_var, &
          copy_output_vars, clear_output_vars, vtk_print_piece_header, &
          vtk_print_piece_data

private

type :: t_output_var
  character(max_char_len) :: var_name
  logical :: vector
  logical :: skip
  real(wp), allocatable :: var(:,:)
end type t_output_var

interface add_output_var
  module procedure add_output_var_scal
  module procedure add_output_var_vec
end interface add_output_var

integer, parameter :: vtk_isize = 4
integer, parameter :: vtk_fsize = 4
character(len=*), parameter :: &
  this_mod_name = 'mod_vtk_utils'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine add_output_var_scal(output_var, var, var_name, skip)
  type(t_output_var), intent(inout) :: output_var
  real(wp), intent(in) :: var(:)
  character(len=*), intent(in) :: var_name
  logical, intent(in) :: skip

  output_var%var_name = trim(var_name)
  output_var%vector = .false.
  output_var%skip = skip

  if(.not. skip) then
    allocate(output_var%var(1,size(var)))
    output_var%var(1,:) = var
  endif

end subroutine add_output_var_scal

!---------------------------------------------------------------------

subroutine add_output_var_vec(output_var, var, var_name, skip)
  type(t_output_var), intent(inout) :: output_var
  real(wp), intent(in) :: var(:,:)
  character(len=*), intent(in) :: var_name
  logical, intent(in) :: skip

  character(len=*), parameter :: this_sub_name = 'add_output_var_vec'

  output_var%var_name = trim(var_name)
  output_var%vector = .true.
  output_var%skip = skip

  if(.not. skip) then
    allocate(output_var%var(size(var,1),size(var,2)))
    output_var%var = var
  endif

  if(size(var,1).ne.3) call internal_error(this_sub_name, this_mod_name, &
  'Printing vector variable with size different from 3')

end subroutine add_output_var_vec

!---------------------------------------------------------------------

subroutine copy_output_vars(output_vars, copy_vars, reset)
  type(t_output_var), intent(in)  :: output_vars(:)
  type(t_output_var), intent(out) :: copy_vars(:)
  logical, intent(in)             :: reset

  integer                         :: i
  character(len=*), parameter     :: this_sub_name = 'copy_output_vars'

  if(size(output_vars).ne.size(copy_vars)) call internal_error(this_sub_name,&
    this_mod_name, 'Copying output vectors of different lengths')

  do i=1,size(output_vars)
    copy_vars(i)%var_name = output_vars(i)%var_name
    copy_vars(i)%vector = output_vars(i)%vector
    copy_vars(i)%skip = output_vars(i)%skip

    allocate(copy_vars(i)%var(size(output_vars(i)%var,1), &
                              size(output_vars(i)%var,2)))
    if(reset) then
      copy_vars(i)%var = 0.0_wp
    else
      copy_vars(i)%var = output_vars(i)%var
    endif
  enddo

end subroutine copy_output_vars

!---------------------------------------------------------------------

subroutine clear_output_vars(output_vars)
  type(t_output_var), intent(out)  :: output_vars(:)

  integer                          :: i
  character(len=*), parameter      :: this_sub_name = 'clear_output_vars'

  do i=1,size(output_vars)
    output_vars(i)%var_name = ''
    output_vars(i)%vector = .false.

    if(allocated(output_vars(i)%var)) deallocate(output_vars(i)%var)
  enddo


end subroutine clear_output_vars

!---------------------------------------------------------------------

subroutine vtk_print_piece_header(fu, offset, npoints, ncells, nquad, ntria, &
                                  npart, out_vars)
  integer, intent(in)            :: fu
  integer, intent(inout)         :: offset
  integer, intent(in)            :: npoints
  integer, intent(in)            :: ncells
  integer, intent(in)            :: nquad
  integer, intent(in)            :: ntria
  integer, intent(in)            :: npart
  type(t_output_var), intent(in) :: out_vars(:)

  character(len=200) :: buffer
  integer :: i_v
  character(len=20) :: ostr, str1, str2
  character(len=1)   :: lf
  character(len=*), parameter :: this_sub_name = 'vtk_print_piece_header'

  lf = char(10) !line feed char
  write(str1,'(I0)') npoints
  write(str2,'(I0)') ncells
  buffer = '  <Piece NumberOfPoints="'//trim(str1)//'" NumberOfCells="'//&
                    trim(str2)//'">'//lf; write(fu) trim(buffer)
  !Points
  buffer =  '   <Points>'//lf; write(fu) trim(buffer)
  write(ostr,'(I0)') offset
  buffer='    <DataArray type="Float32" Name="Coordinates" &
              &NumberOfComponents="3" format="appended" &
              &offset="'//trim(ostr)//'"/>'//lf;write(fu) trim(Buffer)

  buffer =  '   </Points>'//lf; write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*3*npoints

  !Cells
  buffer =  '   <Cells>'//lf; write(fu) trim(buffer)

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="connectivity" &
            &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*3*ntria + vtk_isize*4*nquad + &
                                vtk_isize*1*npart

  write(ostr,'(I0)') offset
  buffer =  '    <DataArray type="Int32" Name="offsets" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*ncells

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="types" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*ncells

  buffer =  '   </Cells>'//lf; write(fu) trim(buffer)

  !Data
  if (size(out_vars) .gt. 0) then
    buffer =  '   <CellData Scalars="scalars">'//lf;
    write(fu) trim(buffer)
    do i_v = 1,size(out_vars)
      if(.not. out_vars(i_v)%skip) then
        if(ncells.ne.size(out_vars(i_v)%var,2)) then
          call internal_error(this_sub_name, this_mod_name, 'Number of output &
                              &variables different from number of cells')
        end if
      end if
      if(out_vars(i_v)%vector) then
        write(ostr,'(I0)') offset
        buffer = '    <DataArray type="Float32" Name="'//trim(out_vars(i_v)%var_name)//'" &
          & NumberOfComponents="3" Format="appended" offset="'//trim(ostr)//'"/>'//lf
        write(fu) trim(buffer)
        offset = offset + vtk_isize + vtk_fsize*3*ncells
      else
        write(ostr,'(I0)') offset
        buffer = '    <DataArray type="Float32" Name="'//trim(out_vars(i_v)%var_name)//'" &
          &Format="appended" offset="'//trim(ostr)//'"/>'//lf
        write(fu) trim(buffer)
        offset = offset + vtk_isize + vtk_fsize*ncells
      endif
    enddo
    buffer = '   </CellData>'//lf; write(fu) trim(buffer)
  endif

  buffer = '  </Piece>'//lf; write(fu) trim(buffer)

end subroutine vtk_print_piece_header

!----------------------------------------------------------------------

subroutine vtk_print_piece_data(fu, out_vars, nquad, ntria, &
                                npart, rr, ee)
  integer, intent(in)            :: fu
  type(t_output_var), intent(in) :: out_vars(:)
  integer, intent(in)            ::  nquad, ntria, npart
  real(wp), intent(in)           :: rr(:,:)
  integer, intent(in), optional  :: ee(:,:)

  integer                        :: nbytes, ne, i_shift, etype
  integer                        :: i, i_v
  logical                        :: particles
  character(len=*), parameter    :: this_sub_name='vtk_print_piece_data'

  if(nquad+ntria .gt. 0 .and. .not. present(ee)) call internal_error( &
  this_sub_name,this_mod_name,'Not passed connectivity to piece with elements')

  particles = .false.
  if(nquad+ntria .le. 0) particles = .true.
  ne = nquad+ntria+npart

  !> Points
  nbytes = vtk_fsize *  &
                size(rr,1)*size(rr,2)
  write(fu) nbytes

  do i=1,size(rr,2)
    write(fu) real(rr(:,i),vtk_fsize)
  enddo

  !> Connectivity
  if(particles) then
    nbytes = vtk_isize*ne; write(fu) nbytes
    do i=1,npart
      write(fu) i-1
    enddo
  else
    nbytes = vtk_isize*3*ntria+vtk_isize*4*nquad; write(fu) nbytes
    do i=1,size(ee,2)
      if (ee(4,i) .eq. 0) then
        write(fu) ee(1:3,i)-1
      else
        write(fu) ee(1:4,i)-1
      endif
    enddo
  endif

  !> Offset
  nbytes =  vtk_isize*ne; write(fu) nbytes
  if(particles) then
    do i=1,npart
      write(fu) i
    enddo
  else
    i_shift = 0
    do i=1,size(ee,2)
      if (ee(4,i) .eq. 0) then
        i_shift = i_shift + 3
      else
        i_shift = i_shift + 4
      endif
      write(fu) i_shift
    enddo
  endif

  !> Cell types
  nbytes = vtk_isize*ne; write(fu) nbytes
  if(particles) then
    do i=1,ne
      etype = 1
      write(fu) etype
    enddo
  else
    do i=1,size(ee,2)
      if (ee(4,i) .eq. 0) then
        etype = 5
      else
        etype = 9
      endif
      write(fu) etype
    enddo
  endif

  !> Variables
  do i_v = 1,size(out_vars)
    if(.not.out_vars(i_v)%vector) then
      nbytes =  vtk_fsize*ne; write(fu) nbytes
      if(.not.out_vars(i_v)%skip) then
        do i=1,size(out_vars(i_v)%var,2)
          write(fu) real(out_vars(i_v)%var(1,i), vtk_fsize)
        enddo
      else
          do i=1,ne
            write(fu) real(0.0_wp, vtk_fsize)
          enddo
      endif
    else
      nbytes =  3*vtk_fsize*ne; write(fu) nbytes
      if(.not.out_vars(i_v)%skip) then
        do i=1,size(out_vars(i_v)%var,2)
          write(fu) real(out_vars(i_v)%var(:,i), vtk_fsize)
        enddo
      else
          do i=1,ne
            write(fu) real((/0.0_wp, 0.0_wp, 0.0_wp/), vtk_fsize)
          enddo
      endif
    endif
  enddo

end subroutine vtk_print_piece_data

!----------------------------------------------------------------------

end module mod_vtk_utils

