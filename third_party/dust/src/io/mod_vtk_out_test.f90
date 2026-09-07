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

module mod_vtk_out_test

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, warning, info, printout, new_file_unit

use mod_vtk_utils, only: &
  t_output_var, vtk_isize, vtk_fsize, vtk_print_piece_header, &
  vtk_print_piece_data


!---------------------------------------------------------------------
implicit none

public :: vtk_out_viz , vtr_write

private

character(len=*), parameter :: &
  this_mod_name = 'mod_vtk_output'
!---------------------------------------------------------------------

contains

! subroutine vtr_write      <--- for RectuangularGrid files
! subroutine vtk_out_bin
! subroutine vtk_out_viz

!> RectilinearGrid: binary xml .vtr file

subroutine vtr_write ( filen , x , y , z , vars_n , vars_name , vars )
  character(len=*), intent(in) :: filen
  real(wp), intent(in) :: x(:) , y(:) , z(:)
  integer , intent(in) :: vars_n(:)
  character(len=*), intent(in) :: vars_name(:)
  real(wp), intent(in) :: vars(:,:)

  integer :: nx , ny , nz , n_points
  integer :: n_vars
  character(len=200) :: buffer
  character(len=20)  :: str1 , ostr , istr
  character(len=1)   :: lf
  integer :: offset , nbytes
  integer :: fid , ierr , i1 , i2 , irow

  character(len=*), parameter :: this_sub_name = 'wtr_write'


  lf = char(10) !line feed char

  nx = size(x) ; ny = size(y) ; nz = size(z)
  n_points = nx * ny * nz
  n_vars = size(vars_n)

  call new_file_unit(fid,ierr)

  open(fid,file=trim(filen), &
      status='replace',access='stream',iostat=ierr)

  ! Header
  buffer = '<?xml version="1.0"?>'//lf ; write(fid) trim(buffer)
  buffer = '<VTKFile type="RectilinearGrid" version="0.1"&
          & byte_order="LittleEndian">'//lf ; write(fid) trim(buffer)
  write(str1,'(I0,a,I0,a,I0,a,I0,a,I0,a,I0)') &
            0,' ',size(x)-1,' ',0,' ',size(y)-1, ' ',0,' ',size(z)-1
  buffer = ' <RectilinearGrid WholeExtent="'//trim(str1)//'">'//lf ; write(fid) trim(buffer)
  buffer = '  <Piece Extent="'//trim(str1)//'">'//lf ; write(fid) trim(buffer)
  ! Coordinates
  buffer = '   <Coordinates>'//lf ; write(fid) trim(buffer)
  offset = 0
  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Float32" Name="x" format="appended"&
          & offset="'//trim(ostr)//'"/>'//lf ; write(fid) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*size(x)
  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Float32" Name="y" format="appended"&
          & offset="'//trim(ostr)//'"/>'//lf ; write(fid) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*size(y)
  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Float32" Name="z" format="appended"&
          & offset="'//trim(ostr)//'"/>'//lf ; write(fid) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*size(z)
  buffer = '   </Coordinates>'//lf ; write(fid) trim(buffer)
  write(ostr,'(I0)') offset
  ! Data
  buffer = '   <PointData>'//lf ; write(fid) trim(buffer)
  do i1 = 1 , n_vars
    if ( vars_n(i1) .eq. 1 ) then ! Scalar variable

      buffer = '    <DataArray type="Float32" Name="'//trim(vars_name(i1))//'"&
              & format="appended" offset="'//trim(ostr)//'"/>'//lf ; write(fid) trim(buffer)
      offset = offset + vtk_isize + vtk_fsize*size(x)*size(y)*size(z)
      write(ostr,'(I0)') offset

    elseif ( vars_n(i1) .eq. 3 ) then ! Vector variable

      buffer = '    <DataArray type="Float32" Name="'//trim(vars_name(i1))//'"&
              & NumberOfComponents="3" format="appended" offset="'//trim(ostr)//'"/>'//lf
      write(fid) trim(buffer)
      offset = offset + vtk_isize + 3*vtk_fsize*size(x)*size(y)*size(z)
      write(ostr,'(I0)') offset

    !elseif ( vars_n(i1) .eq. 9 ) then ! Tensor Variable

    else
      write(istr,'(I0)') i1
      call error(this_sub_name, this_mod_name, 'Wrong dimension of the variable n.'//istr)
    end if
  end do
  buffer = '   </PointData>'//lf ; write(fid) trim(buffer)
  buffer = '  </Piece>'//lf ; write(fid) trim(buffer)
  buffer = ' </RectilinearGrid>'//lf; write(fid) trim(buffer)
  ! Appended data
  buffer = ' <AppendedData encoding="raw">'//lf ; write(fid) trim(buffer)
  buffer = '_'; write(fid) trim(buffer) !mark the beginning of the data
  ! Coordinates
  nbytes = vtk_fsize * nx
  write(fid) nbytes
  write(fid) real(x,vtk_fsize)
  nbytes = vtk_fsize * ny
  write(fid) nbytes
  write(fid) real(y,vtk_fsize)
  nbytes = vtk_fsize * nz
  write(fid) nbytes
  write(fid) real(z,vtk_fsize)
  ! Data
  irow = 1
  do i1 = 1 , n_vars
    if ( vars_n(i1) .eq. 1 ) then
      nbytes = vtk_fsize * n_points
      write(fid) nbytes
      write(fid) real( vars(irow,:) , vtk_fsize )
      irow = irow + 1
    elseif ( vars_n(i1) .eq. 3 ) then
      nbytes = vtk_fsize * 3 * n_points
      write(fid) nbytes
      do i2 = 1 , n_points
        write(fid) real( vars(irow:irow+2,i2) , vtk_fsize )
      end do
      irow = irow + 3
    !elseif ( vars_n(i1) .eq. 9 ) then ! Tensor Variable
    else
      write(istr,'(I0)') i1
      call error(this_sub_name, this_mod_name, 'Wrong dimension of the variable n.'//istr)
    end if
  end do
  buffer = ' </AppendedData>'//lf ; write(fid) trim(buffer)
  buffer = '</VTKFile>' ; write(fid) trim(buffer)

  close(fid)

end subroutine vtr_write

!---------------------------------------------------------------------

!> Output the processed data into a binary xml  .vtu file
!!
!! This is a more advanced version, supporting different variables
subroutine vtk_out_viz (out_filename, &
                        vp_rr, vp_vars)
  character(len=*), intent(in)              :: out_filename
  real(wp), intent(in), optional            :: vp_rr(:,:)
  type(t_output_var), intent(in), optional  :: vp_vars(:)

  integer :: fu, ierr
  integer :: npoints, ncells, ne
  integer :: offset
  integer :: slen
  character(len=200) :: buffer
  character(len=1)  :: lf

  integer :: ie, nquad, ntria, nquad_w, ntria_w, nvp
  integer :: npoints_w, nw

  lf = char(10) !line feed char


  nvp = size(vp_rr,2)

  call new_file_unit(fu,ierr)
  open(fu,file=trim(out_filename), &
        status='replace',access='stream',iostat=ierr)
  offset = 0
  buffer = '<?xml version="1.0"?>'//lf; write(fu) trim(buffer)
  buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
  write(fu) trim(buffer)
  buffer = ' <UnstructuredGrid>'//lf; write(fu) trim(buffer)
  !=== Vortex Particles
  call vtk_print_piece_header(fu,offset,size(vp_rr,2),size(vp_rr,2),0,0, &
                                size(vp_rr,2),vp_vars)
  buffer = ' </UnstructuredGrid>'//lf; write(fu) trim(buffer)

  !All the appended data
  buffer = '  <AppendedData encoding="raw">'//lf; write(fu) trim(buffer)
  buffer = '_'; write(fu) trim(buffer) !mark the beginning of the data

  !=  == Particles
  call vtk_print_piece_data(fu, vp_vars, 0, 0, nvp, vp_rr)

  buffer = '  </AppendedData>'//lf; write(fu) trim(buffer)
  buffer = '</VTKFile>'//lf; write(fu) trim(buffer)


  close(fu,iostat=ierr)

end subroutine vtk_out_viz

end module mod_vtk_out_test

