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
!!          Matteo Tugnoli
!!=========================================================================



!-----------------------------------------------------------------------
module mod_hdf5_io

!-----------------------------------------------------------------------
 use HDF5

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, warning, info, printout, new_file_unit

!-----------------------------------------------------------------------

 implicit none

!-----------------------------------------------------------------------

! Module interface

 public :: &
   h5loc, &
   initialize_hdf5, &
   destroy_hdf5, &
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
   read_hdf5_attr, &
   append_hdf5, &
   h5t_mem_float, &
   h5t_file_float, &
   check_dset_hdf5, &
   get_dset_dimensions_hdf5


 private

!-----------------------------------------------------------------------
 !module variables

 character(len=*), parameter :: &
   this_mod_name = 'mod_hdf5_io'


 integer, parameter :: h5loc = HID_T
 integer, parameter :: h5sz = HSIZE_T


 !hdf5 types
 integer(h5loc) :: h5t_file_float, h5t_file_int, h5t_file_char
 integer(h5loc) :: h5t_mem_float, h5t_mem_int, h5t_mem_char
!-----------------------------------------------------------------------
 !Interfaces

   interface write_hdf5
   module procedure write_string_hdf5, &
                     write_1d_string_hdf5, &
                     write_2d_string_hdf5, &
                     write_int_hdf5, &
                     write_1d_int_hdf5, &
                     write_2d_int_hdf5, &
                     write_3d_int_hdf5, &
                     write_real_hdf5, &
                     write_1d_real_hdf5, &
                     write_2d_real_hdf5, &
                     write_3d_real_hdf5
   end interface

 interface write_hdf5_attr
  module procedure write_1d_int_hdf5_attr, &
                   write_int_hdf5_attr, &
                   write_1d_real_hdf5_attr, &
                   write_real_hdf5_attr, &
                   write_logical_hdf5_attr, &
                   write_string_hdf5_attr
 end interface

 interface append_hdf5
  module procedure append_1d_real_hdf5, &
                   append_1d_int_hdf5
 end interface

 interface read_hdf5
  module procedure read_int_hdf5, &
                   read_1d_int_hdf5, &
                   read_2d_int_hdf5, &
                   read_3d_int_hdf5, &
                   read_real_hdf5, &
                   read_1d_real_hdf5, &
                   read_2d_real_hdf5, &
                   read_3d_real_hdf5, &
                   read_string_hdf5, &
                   read_1d_string_hdf5, &
                   read_2d_string_hdf5
 end interface read_hdf5

 interface read_hdf5_attr
  module procedure read_1d_int_hdf5_attr, &
                   read_int_hdf5_attr, &
                   read_1d_real_hdf5_attr, &
                   read_real_hdf5_attr, &
                   read_logical_hdf5_attr, &
                   read_string_hdf5_attr
 end interface

 interface read_hdf5_al
  module procedure read_1d_real_hdf5_al, &
                   read_2d_real_hdf5_al, &
                   read_3d_real_hdf5_al, &
                   read_1d_int_hdf5_al, &
                   read_2d_int_hdf5_al, &
                   read_3d_int_hdf5_al, &
                   read_1d_string_hdf5_al, &
                   read_2d_string_hdf5_al
 end interface read_hdf5_al

 interface check_dset_hdf5
   module procedure check_dset_hdf5
 end interface check_dset_hdf5


!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!===== Files/Groups Operations =====
!-----------------------------------------------------------------------

 subroutine initialize_hdf5()

  integer :: h5err
  logical :: equal
  character(len=*), parameter :: &
    this_sub_name = 'initialize_hdf5'


   ! Load Hdf5 library and create parallel access properties
   call h5open_f(h5err)

   ! Select the output types
   h5t_file_float = h5kind_to_type(wp, H5_REAL_KIND)
   h5t_file_int   = H5T_STD_I32LE
   h5t_file_char  = H5T_NATIVE_CHARACTER ! TODO: set a universal value

   h5t_mem_float = h5kind_to_type(wp, H5_REAL_KIND)
   h5t_mem_int   = H5T_NATIVE_INTEGER
   h5t_mem_char  = H5T_NATIVE_CHARACTER

   ! check the datatypes
   equal = .false.
   call h5Tequal_f(h5t_file_float, h5t_mem_float, equal, h5err)
   if (.not. equal) then
       call warning(this_mod_name, this_sub_name, 'Selected and '// &
       'native real types do not match, consider making another '//&
       'choiche for efficency')
   endif
   call h5Tequal_f(h5t_file_int, h5t_mem_int, equal, h5err)
   if (.not. equal) then
       call warning(this_mod_name, this_sub_name, 'Selected and '// &
       'native integer types do not match, consider making another '//&
       'choiche for efficency')
   endif
   call h5Tequal_f(h5t_file_char, h5t_mem_char, equal, h5err)
   if (.not. equal) then
       call warning(this_mod_name, this_sub_name, 'Selected and '// &
       'native character types do not match, consider making another'//&
       ' choiche for efficency')
   endif

 end subroutine initialize_hdf5

!-----------------------------------------------------------------------

 subroutine destroy_hdf5()
  integer :: h5err
  character(len=*), parameter :: &
    this_sub_name = 'destroy_hdf5'

   ! Unload Hdf5 library
   call h5close_f(h5err)
 end subroutine destroy_hdf5

!-----------------------------------------------------------------------

!> Create a new hdf5 file and return the hdf5 location
subroutine new_hdf5_file(filename, file_id)
 character(len=*), intent(in) :: filename
 integer(HID_T), intent(out)   :: file_id
 integer :: err_setting

 character(len=*), parameter :: &
    this_sub_name = 'new_hdf5_file'
 integer :: h5err
  
  err_setting = 0
  call h5eset_auto_f(err_setting, h5err)
  call h5Fcreate_f (filename, H5F_ACC_TRUNC_F, file_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                          'Problems creating hdf5 file '//filename)

  err_setting=1
  call h5eset_auto_f(err_setting, h5err)
  
end subroutine new_hdf5_file

!-----------------------------------------------------------------------

!> Open an existing hdf5 file and return the hdf5 location
subroutine open_hdf5_file(filename, file_id)
 character(len=*), intent(in) :: filename
 integer(HID_T), intent(out)   :: file_id
 
 character(len=*), parameter :: &
    this_sub_name = 'open_hdf5_file'
 integer :: h5err
 integer :: err_setting
 
  err_setting = 0
  call h5eset_auto_f(err_setting, h5err)
  call h5Fopen_f (filename, H5F_ACC_RDWR_F, file_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                            'Problems opening hdf5 file '//filename)
  err_setting=1
  call h5eset_auto_f(err_setting, h5err)
  
end subroutine open_hdf5_file

!-----------------------------------------------------------------------

!> Close an opened hdf5 file
subroutine close_hdf5_file(file_id)
 integer(HID_T) :: file_id

 character(len=*), parameter :: &
    this_sub_name = 'close_hdf5_file'
 integer :: h5err

  call h5Fclose_f (file_id, h5err)
end subroutine close_hdf5_file

!-----------------------------------------------------------------------

!> Create a new hdf5 group (folder) in an existing open hdf5 location
!! and return the group location
subroutine new_hdf5_group (loc_id, groupname, group_id)
 integer(HID_T), intent(in)   :: loc_id
 character(len=*), intent(in) :: groupname
 integer(HID_T), intent(out)  :: group_id

 character(len=*), parameter :: &
    this_sub_name = 'new_hdf5_group'
 integer :: h5err

  call h5Gcreate_f(loc_id, groupname, group_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                            'Problems creating hdf5 group '//groupname)

end subroutine new_hdf5_group

!-----------------------------------------------------------------------

!> Open an existing group (folder) from its location and return the
!! group location
subroutine open_hdf5_group (loc_id, groupname, group_id)
 integer(HID_T), intent(in)   :: loc_id
 character(len=*), intent(in) :: groupname
 integer(HID_T), intent(out)  :: group_id

 character(len=*), parameter :: &
    this_sub_name = 'open_hdf5_group'
 integer :: h5err

  call h5Gopen_f(loc_id, groupname, group_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                             'Problems opening hdf5 group '//groupname)

end subroutine open_hdf5_group

!-----------------------------------------------------------------------

!> Close an open hdf5 group location
subroutine close_hdf5_group (group_id)
 integer(HID_T), intent(in) :: group_id

 character(len=*), parameter :: &
    this_sub_name = 'close_hdf5_group'
 integer :: h5err

  call h5Gclose_f(group_id, h5err)

end subroutine close_hdf5_group

!-----------------------------------------------------------------------

!> function to check if a dataset is actually present inside an open
!! file/group, returns a logical
function check_dset_hdf5(dsetname, loc_id) result(dset_exists)
 character(len=*), intent(in)  :: dsetname
 integer(HID_T), intent(in)    :: loc_id
 logical          :: dset_exists


 character(len=*), parameter :: &
    this_sub_name = 'collective_append_1d_real_hdf5'
 integer :: h5err
 logical :: l_exists

  !check the presence of the requested dataset (actually checking the
  !link to it)
  call h5Lexists_f(loc_id, dsetname, l_exists, h5err)
  dset_exists = l_exists

end function check_dset_hdf5


!-----------------------------------------------------------------------

subroutine get_dset_dimensions_hdf5(dsetname, loc_id, dims)
 character(len=*), intent(in)  :: dsetname
 integer(HID_T), intent(in)    :: loc_id
 integer, allocatable, intent(out) :: dims(:)

 integer(HID_T) :: dset_id
 integer(HID_T) :: filespace_id
 integer :: h5err
 integer :: dset_rank
 integer(HSIZE_T), allocatable :: max_dims(:), act_dims(:)
 character(len=*), parameter :: this_sub_name = 'get_dset_dimensions_hdf5_hdf5'

  !open the dataset
  call h5Dopen_f(loc_id, dsetname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//dsetname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, dset_rank, h5err)
  allocate(act_dims(dset_rank), max_dims(dset_rank))
  call h5Sget_simple_extent_dims_f(filespace_id, act_dims, &
                                   max_dims, h5err)

  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

  allocate(dims(dset_rank))
  dims = int(act_dims)
  deallocate(act_dims, max_dims)
end subroutine get_dset_dimensions_hdf5

!-----------------------------------------------------------------------
!==== String Operations ====
!-----------------------------------------------------------------------

!> Write a single string.
!!
subroutine write_string_hdf5(outdata, outname, file_id)
 character(len=*), intent(in), target :: outdata
 character(len=*), intent(in)         :: outname
 integer(HID_T), intent(in)           :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_string_hdf5'
 integer :: h5err
 integer ::strlen
 logical :: l_exists

  strlen = len(outdata)
  out_size(1) = 1

  !create the string datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tcopy_f(h5t_file_char, filetype_id, h5err)

  if(strlen .gt. 0) then !create the string type and a scalar space

    call h5Tset_size_f(memtype_id, int(strlen,HSIZE_T),h5err)
    call h5Tset_size_f(filetype_id, int(strlen,HSIZE_T),h5err)

    !create a (single) scalar dataspace
    call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  else !create a 1 rank space of lenth 0

    out_size = 0
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  endif

  !check if the dataset on file exists
  call h5Lexists_f(file_id, trim(outname), l_exists, h5err)

  if (l_exists) then
    call h5gunlink_f(file_id, trim(outname), h5err)
  endif
  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), filetype_id, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, memtype_id, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
end subroutine write_string_hdf5

!-----------------------------------------------------------------------

!> Read a single string.
!!
subroutine read_string_hdf5(input, inputname, loc_id)
 character(len=*), intent(out)        :: input
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'read_string_hdf5'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  in_size(1) = 1
  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  !get the type
  call h5Dget_type_f(dset_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  !WARNING: no dimensions check for the moment
  call h5Dread_f(dset_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading dataset '//inputname)

  !release resources
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Dclose_f(dset_id, h5err)
end subroutine read_string_hdf5

!-----------------------------------------------------------------------
!> Write a 1D array of strings
!!
!! If the array is empty nothing is written (not an empty dataset as in
!! other cases), due to known issues to load empty character dataset by
!! hdf5
subroutine write_1d_string_hdf5(outdata, outname, file_id)
 character(len=*), intent(in)  :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_string_hdf5'
 integer :: h5err
 integer ::strlen

  strlen = len(outdata)
  out_size(1) = int(size(outdata, 1),h5sz)

  !write only if the array is not empty
  ! (empty arrays of strings create problems while loading)
  if (out_size(1) .gt. 0) then
    !create the string datatypes
    call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
    call h5Tset_size_f(memtype_id, int(strlen,HSIZE_T),h5err)
    call h5Tcopy_f(h5t_file_char, filetype_id, h5err)
    call h5Tset_size_f(filetype_id, int(strlen,HSIZE_T),h5err)

    !create a 1D scalar dataspace
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

    !create the dataset on the file
    call h5Dcreate_f(file_id, trim(outname), filetype_id, dspace_id, &
                     dset_id, h5err)

    !write
    call h5Dwrite_f(dset_id, memtype_id, outdata, out_size, h5err)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

    !release resources
    call h5Dclose_f(dset_id, h5err)
    call h5Sclose_f(dspace_id, h5err)
    call h5Tclose_f(memtype_id, h5err)
    call h5Tclose_f(filetype_id, h5err)
  endif
end subroutine write_1d_string_hdf5

!-----------------------------------------------------------------------

!> Read a string array.
!!
subroutine read_1d_string_hdf5(input, inputname, loc_id)
 character(len=*), intent(out)        :: input(:)
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 character(len=*), parameter :: &
    this_sub_name = 'read_string_hdf5'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get the type
  call h5Dget_type_f(dset_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if(in_size(1) .ne. int(size(input,1), h5sz)) &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)

  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  call h5Dread_f(dset_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading dataset '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Dclose_f(dset_id, h5err)
end subroutine read_1d_string_hdf5

!-----------------------------------------------------------------------

!> Read and allocate string array.
!!
subroutine read_1d_string_hdf5_al(input, inputname, loc_id)
 character(len=*), allocatable, intent(out) :: input(:)
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 character(len=*), parameter :: &
    this_sub_name = 'read_1d_string_hdf5_al'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get the type
  call h5Dget_type_f(dset_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)

  ! allocate memory data
  allocate(input(in_size(1)))
  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  call h5Dread_f(dset_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading dataset '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Dclose_f(dset_id, h5err)
end subroutine read_1d_string_hdf5_al

!-----------------------------------------------------------------------

!> Write a 2D array of strings
!!
!! If the array is empty nothing is written (not an empty dataset as in
!! other cases), due to known issues to load empty character dataset by
!! hdf5
subroutine write_2d_string_hdf5(outdata, outname, file_id)
 character(len=*), intent(in)  :: outdata(:,:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=2
 integer(HSIZE_T) :: out_size(2)
 character(len=*), parameter :: &
    this_sub_name = 'write_string_hdf5'
 integer :: h5err
 integer ::strlen

  strlen = len(outdata)
  out_size(1) = int(size(outdata, 1), h5sz)
  out_size(2) = int(size(outdata, 2), h5sz)

  !write only if the array is not empty
  ! (empty arrays of strings create problems while loading)
  if (out_size(1) .gt. 0) then
    !create the string datatypes
    call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
    call h5Tset_size_f(memtype_id, int(strlen,HSIZE_T),h5err)
    call h5Tcopy_f(h5t_file_char, filetype_id, h5err)
    call h5Tset_size_f(filetype_id, int(strlen,HSIZE_T),h5err)

    !create a 1D scalar dataspace
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

    !create the dataset on the file
    call h5Dcreate_f(file_id, trim(outname), filetype_id, dspace_id, &
                     dset_id, h5err)

    !write
    call h5Dwrite_f(dset_id, memtype_id, outdata, out_size, h5err)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

    !release resources
    call h5Dclose_f(dset_id, h5err)
    call h5Sclose_f(dspace_id, h5err)
    call h5Tclose_f(memtype_id, h5err)
    call h5Tclose_f(filetype_id, h5err)
  endif
end subroutine write_2d_string_hdf5

!-----------------------------------------------------------------------

!> Read a string 2d array.
!!
subroutine read_2d_string_hdf5(input, inputname, loc_id)
 character(len=*), intent(out)        :: input(:,:)
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 character(len=*), parameter :: &
    this_sub_name = 'read_string_hdf5'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get the type
  call h5Dget_type_f(dset_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 2 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if((in_size(1) .ne. int(size(input,1), h5sz)) .or. &
     (in_size(2) .ne. int(size(input,2), h5sz))) &
       call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)

  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  call h5Dread_f(dset_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading dataset '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Dclose_f(dset_id, h5err)
end subroutine read_2d_string_hdf5

!-----------------------------------------------------------------------

!> Read a string 2d array.
!!
subroutine read_2d_string_hdf5_al(input, inputname, loc_id)
 character(len=*), allocatable, intent(out) :: input(:,:)
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 character(len=*), parameter :: &
    this_sub_name = 'read_string_hdf5'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get the type
  call h5Dget_type_f(dset_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 2 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1),in_size(2)))
  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  call h5Dread_f(dset_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading dataset '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Dclose_f(dset_id, h5err)
end subroutine read_2d_string_hdf5_al

!-----------------------------------------------------------------------
! ==== Write Reals ====
!-----------------------------------------------------------------------

!> Writes a single real value
!!
subroutine write_real_hdf5(outdata, outname, file_id)
 real(wp), intent(in)          :: outdata
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_real_hdf5'
 integer :: h5err

  !create a (single) scalar dataspace
  out_size(1) = 1
  call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_float, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_real_hdf5

!-----------------------------------------------------------------------

!> Write a rank 1 array of reals
!!
subroutine write_1d_real_hdf5(outdata, outname, file_id)
 real(wp), intent(in)          :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_1d_real_hdf5'
 integer :: h5err

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1),h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_float, dspace_id, &
                   dset_id, h5err)
  !write
  call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_1d_real_hdf5

!-----------------------------------------------------------------------

!> Write a rank 2 array of reals
!!
subroutine write_2d_real_hdf5(outdata, outname, file_id)
 real(wp), intent(in)          :: outdata(:,:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=2
 integer(HSIZE_T) :: out_size(2)
 character(len=*), parameter :: &
    this_sub_name = 'write_2d_real_hdf5'
 integer :: h5err

  !create a 2D dataspace
  out_size(1) = int(size(outdata,1),h5sz)
  out_size(2) = int(size(outdata,2),h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_float, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_2d_real_hdf5


!-----------------------------------------------------------------------

!> Write a rank 3 array of reals
!!
subroutine write_3d_real_hdf5(outdata, outname, file_id)
 real(wp), intent(in)          :: outdata(:,:,:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=3
 integer(HSIZE_T) :: out_size(3)
 character(len=*), parameter :: &
    this_sub_name = 'write_3d_real_hdf5'
 integer :: h5err

  !create a 3D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  out_size(2) = int(size(outdata,2), h5sz)
  out_size(3) = int(size(outdata,3), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_float, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_3d_real_hdf5

!-----------------------------------------------------------------------
! ==== Append Reals ====
!-----------------------------------------------------------------------

!> Subroutine to append a 1D real array to an existing dataset, to
!! create a time series in a single file. If the dataset does not exist
!! (first time used) the dataset is created.
subroutine append_1d_real_hdf5(outdata, outname, file_id)
 real(wp), intent(in)          :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer(HID_T) :: prp_id
 integer, parameter :: rank=2
 integer(HSIZE_T) :: out_size(2), maxsize(2), chunksize(2), totsize(2)
 integer(HSIZE_T) :: oset(2)
 character(len=*), parameter :: &
    this_sub_name = 'append_1d_real_hdf5'
 integer :: h5err
 logical :: l_exists
 integer :: ndims

  !check the presence of the requested dataset (actually checking the
  !link to it)
  out_size(1) = int(size(outdata,1), h5sz); out_size(2) = 1
  call h5Lexists_f(file_id, outname, l_exists, h5err)
  if (.not. l_exists) then
    !build the dataset
    maxsize = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    !create unlimited dataset
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err, &
                                                        maxdims=maxsize)
    !create chunking properties
    call h5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, h5err)
    !TODO:  set a sensible value for the chunk
    chunksize = int((/2,5/), h5sz)
    call h5Pset_chunk_f(prp_id, rank, chunksize, h5err)
    !create dataset with beginning dimensions
    call h5Dcreate_f(file_id, trim(outname), h5t_file_float, &
                     dspace_id, dset_id, h5err, dcpl_id=prp_id )
    !close the property
    call h5Pclose_f(prp_id, h5err)

    !call the write
    call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)
  else
    !expand dataset and append data
    !open dataset
    call h5Dopen_f(file_id, trim(outname), dset_id, h5err)
    !get dataspace from dataset
    call h5Dget_space_f(dset_id, dspace_id, h5err)
    !get dimensions of the dataset
    call h5Sget_simple_extent_ndims_f(dspace_id, ndims, h5err)
    call h5Sget_simple_extent_dims_f(dspace_id, totsize, maxsize, h5err)
    if ((ndims .ne. rank) .or. (totsize(1) .ne. out_size(1))) &
      call error(this_sub_name,this_mod_name, &
                'Wrong dimensions on the dataset about to be extended')
    call H5Sclose_f(dspace_id, h5err)
    !extend dataset
    oset(1) = 0; oset(2) = totsize(2)
    totsize(2) = totsize(2) +1
    call h5Dset_extent_f(dset_id, totsize, h5err)
    !select only the appropriate target and write
    call h5Dget_space_f(dset_id, dspace_id, h5err)
    call h5Screate_simple_f(rank, out_size, memspace_id, h5err)
    call h5Sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
                               oset, out_size, h5err)
    call h5Dwrite_f(dset_id, h5t_mem_float, outdata, out_size, h5err,  &
                    mem_space_id=memspace_id, file_space_id=dspace_id)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                             'Problems writing dataset '//outname)
    call h5Sclose_f(memspace_id, h5err)
  endif

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine append_1d_real_hdf5

!-----------------------------------------------------------------------
! ==== Read Reals ====
!-----------------------------------------------------------------------

!> Read a single real
subroutine read_real_hdf5(input, inputname, loc_id)
 real(wp), intent(out) :: input
 character(len=*), intent(in) :: inputname
 integer(HID_T), intent(in) :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_real_hdf5'
 integer(HSIZE_T) :: in_size(1)
 integer(HID_T) :: dset_id
 integer :: h5err

  in_size = int((/1/), h5sz)
  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)
  !release resources
  call h5Dclose_f(dset_id, h5err)

end subroutine read_real_hdf5

!-----------------------------------------------------------------------

!> Read a rank 1 array of reals
subroutine read_1d_real_hdf5(input, inputname, loc_id)
 real(wp), intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_real_hdf5'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if(in_size(1) .ne. int(size(input,1), h5sz)) &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_1d_real_hdf5

!-----------------------------------------------------------------------

!> Read a rank 2 array of reals
subroutine read_2d_real_hdf5(input, inputname, loc_id)
 real(wp), intent(out) :: input(:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_2d_real_hdf5'
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 2 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if((in_size(1) .ne. int(size(input,1), h5sz)) .or. &
     (in_size(2) .ne. int(size(input,2), h5sz))) &
       call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_2d_real_hdf5

!-----------------------------------------------------------------------

!> Read a rank 3 array of reals
subroutine read_3d_real_hdf5(input, inputname, loc_id)
 real(wp), intent(out) :: input(:,:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_3d_real_hdf5'
 integer(HSIZE_T) :: in_size(3), max_size(3)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 3) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 2 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if((in_size(1) .ne. int(size(input,1), h5sz)) .or. &
     (in_size(2) .ne. int(size(input,2), h5sz)) .or.&
     (in_size(3) .ne. int(size(input,3), h5sz))) &
       call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_3d_real_hdf5

!-----------------------------------------------------------------------
! ==== Read and allocate Reals ====
!-----------------------------------------------------------------------

!> Read and allocate a rank 1 array of reals
subroutine read_1d_real_hdf5_al(input, inputname, loc_id)
 real(wp), allocatable, intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_real_hdf5_al'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1)))
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_1d_real_hdf5_al

!-----------------------------------------------------------------------

!> Read and allocate a rank 2 array of reals
subroutine read_2d_real_hdf5_al(input, inputname, loc_id)
 real(wp), allocatable, intent(out) :: input(:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_2d_real_hdf5_al'
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1),in_size(2)))
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_2d_real_hdf5_al

!-----------------------------------------------------------------------

!> Read and allocate a rank 3 array of reals
subroutine read_3d_real_hdf5_al(input, inputname, loc_id)
 real(wp), allocatable, intent(out) :: input(:,:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_3d_real_hdf5_al'
 integer(HSIZE_T) :: in_size(3), max_size(3)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 3) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1),in_size(2),in_size(3)))
  !read
  call h5Dread_f(dset_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_3d_real_hdf5_al

!-----------------------------------------------------------------------
! ==== Write Integers ====
!-----------------------------------------------------------------------

!> Write a single integer value
!!
subroutine write_int_hdf5(outdata, outname, file_id)
 integer, intent(in)           :: outdata
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_int_hdf5'
 integer :: h5err
 logical :: l_exists

  !create a (single) scalar dataspace
  out_size(1) = 1
  call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Lexists_f(file_id, trim(outname), l_exists, h5err)

  if (.not. l_exists) then
    !create the dataset on the file
    call h5Dcreate_f(file_id, trim(outname), h5t_file_int, dspace_id, &
                     dset_id, h5err)
  else
    !open the dataset
    !TODO: perform a check on the dimensions
    call h5Dopen_f(file_id, trim(outname), dset_id, h5err)

  endif

  !write
  call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_int_hdf5

!-----------------------------------------------------------------------

!> Write a rank 1 array of integers
!!
subroutine write_1d_int_hdf5(outdata, outname, file_id)
 integer, intent(in)           :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_1d_int_hdf5'
 integer :: h5err

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_int, dspace_id, &
                   dset_id, h5err)
  !write
  call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_1d_int_hdf5

!-----------------------------------------------------------------------

!> Write a rank 2 array of integers
!!
subroutine write_2d_int_hdf5(outdata, outname, file_id)
 integer, intent(in)           :: outdata(:,:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=2
 integer(HSIZE_T) :: out_size(2)
 character(len=*), parameter :: &
    this_sub_name = 'write_2d_int_hdf5'
 integer :: h5err

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  out_size(2) = int(size(outdata,2), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_int, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_2d_int_hdf5

!-----------------------------------------------------------------------

!> Write a rank 3 array of integers
!!
subroutine write_3d_int_hdf5(outdata, outname, file_id)
 integer, intent(in)           :: outdata(:,:,:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=3
 integer(HSIZE_T) :: out_size(3)
 character(len=*), parameter :: &
    this_sub_name = 'write_3d_int_hdf5'
 integer :: h5err

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  out_size(2) = int(size(outdata,2), h5sz)
  out_size(3) = int(size(outdata,3), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !create the dataset on the file
  call h5Dcreate_f(file_id, trim(outname), h5t_file_int, dspace_id, &
                   dset_id, h5err)

  !write
  call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing dataset '//outname)

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_3d_int_hdf5

!-----------------------------------------------------------------------
! ==== Append Integers ====
!-----------------------------------------------------------------------

!> Subroutine to append a 1D integer array to an existing dataset, to
!! create a time series in a single file. If the dataset does not exist
!! (first time used) the dataset is created
subroutine append_1d_int_hdf5(outdata, outname, file_id)
 integer, intent(in)          :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: file_id

 integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: dset_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer(HID_T) :: prp_id
 integer, parameter :: rank=2
 integer(HSIZE_T) :: out_size(2), maxsize(2), chunksize(2), totsize(2)
 integer(HSIZE_T) :: oset(2)
 character(len=*), parameter :: &
    this_sub_name = 'append_1d_real_hdf5'
 integer :: h5err
 logical :: l_exists
 integer :: ndims

  !check the presence of the requested dataset (actually checking the
  !link to it)
  out_size(1) = int( size(outdata,1), h5sz)
  out_size(2) = int( 1, h5sz)
  call h5Lexists_f(file_id, outname, l_exists, h5err)
  if (.not. l_exists) then
    !build the dataset
    maxsize = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    !create unlimited dataset
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err, &
                                                        maxdims=maxsize)
    !create chunking properties
    call h5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, h5err)
    !TODO:  set a sensible value for the chunk
    chunksize = int((/2,5/), h5sz)
    call h5Pset_chunk_f(prp_id, rank, chunksize, h5err)
    !create dataset with beginning dimensions
    call h5Dcreate_f(file_id, trim(outname), h5t_file_int, &
                     dspace_id, dset_id, h5err, dcpl_id=prp_id )
    !close the property
    call h5Pclose_f(prp_id, h5err)
    !call the write
    call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                             'Problems writing dataset '//outname)
  else
    !expand dataset and append data
    !open dataset
    call h5Dopen_f(file_id, trim(outname), dset_id, h5err)
    !get dataspace from dataset
    call h5Dget_space_f(dset_id, dspace_id, h5err)
    !get dimensions of the dataset
    call h5Sget_simple_extent_ndims_f(dspace_id, ndims, h5err)
    call h5Sget_simple_extent_dims_f(dspace_id, totsize, maxsize, h5err)
    if ((ndims .ne. rank) .or. (totsize(1) .ne. out_size(1))) &
      call error(this_sub_name,this_mod_name, &
                'Wrong dimensions on the dataset about to be extended')
    call H5Sclose_f(dspace_id, h5err)
    !extend dataset
    oset(1) = 0; oset(2) = totsize(2)
    totsize(2) = totsize(2) +1
    call h5Dset_extent_f(dset_id, totsize, h5err)
    !select only the appropriate target and write
    call h5Dget_space_f(dset_id, dspace_id, h5err)
    call h5Screate_simple_f(rank, out_size, memspace_id, h5err)
    call h5Sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
                               oset, out_size, h5err)
    call h5Dwrite_f(dset_id, h5t_mem_int, outdata, out_size, h5err,  &
                    mem_space_id=memspace_id, file_space_id=dspace_id)
    if(h5err<0) call error(this_mod_name,this_sub_name, &
                             'Problems writing dataset '//outname)
    call h5Sclose_f(memspace_id, h5err)
  endif

  !release resources
  call h5Dclose_f(dset_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine append_1d_int_hdf5

!-----------------------------------------------------------------------
! ==== Read Integers ====
!-----------------------------------------------------------------------

!> Read a single integer
subroutine read_int_hdf5(input, inputname, loc_id)
 integer, intent(out) :: input
 character(len=*), intent(in) :: inputname
 integer(HID_T), intent(in) :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_int_hdf5'
 integer(HSIZE_T) :: in_size(1)
 integer(HID_T) :: dset_id
 integer :: h5err

  in_size = int((/1/), h5sz)
  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)
  !release resources
  call h5Dclose_f(dset_id, h5err)

end subroutine read_int_hdf5

!-----------------------------------------------------------------------

!> Read a rank 1 array of integers
subroutine read_1d_int_hdf5(input, inputname, loc_id)
 integer, intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_int_hdf5'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if(in_size(1) .ne. int(size(input,1), h5sz))  &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_1d_int_hdf5

!-----------------------------------------------------------------------

!> Read a rank 2 array of integers
subroutine read_2d_int_hdf5(input, inputname, loc_id)
 integer, intent(out) :: input(:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_2d_int_hdf5'
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 2 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if((in_size(1) .ne. int(size(input,1), h5sz)) .or. &
     (in_size(2) .ne. int(size(input,2), h5sz))) &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_2d_int_hdf5

!-----------------------------------------------------------------------

!> Read a rank 3 array of integers
subroutine read_3d_int_hdf5(input, inputname, loc_id)
 integer, intent(out) :: input(:,:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_3d_int_hdf5'
 integer(HSIZE_T) :: in_size(3), max_size(3)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 3) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 3 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if((in_size(1) .ne. int(size(input,1), h5sz)) .or. &
     (in_size(2) .ne. int(size(input,2), h5sz)) .or. &
     (in_size(3) .ne. int(size(input,3), h5sz))) &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_3d_int_hdf5

!-----------------------------------------------------------------------
! ==== Read and allocate Integers ====
!-----------------------------------------------------------------------

!> Read and allocate a rank 1 array of integers
subroutine read_1d_int_hdf5_al(input, inputname, loc_id)
 integer, allocatable, intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_int_hdf5_al'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1)))
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_1d_int_hdf5_al

!-----------------------------------------------------------------------

!> Read and allocate a rank 2 array of integers
subroutine read_2d_int_hdf5_al(input, inputname, loc_id)
 integer, allocatable, intent(out) :: input(:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_2d_int_hdf5_al'
 integer(HSIZE_T) :: in_size(2), max_size(2)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 2) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1),in_size(2)))
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_2d_int_hdf5_al

!-----------------------------------------------------------------------

!> Read and allocate a rank 3 array of integers
subroutine read_3d_int_hdf5_al(input, inputname, loc_id)
 integer, allocatable, intent(out) :: input(:,:,:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_3d_int_hdf5_al'
 integer(HSIZE_T) :: in_size(3), max_size(3)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: dset_id
 integer :: h5err

  !open the dataset
  call h5Dopen_f(loc_id, inputname, dset_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening dataset  '//inputname)
  !get its dataspace
  call h5Dget_space_f(dset_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 3) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  ! allocate memory data
  allocate(input(in_size(1),in_size(2),in_size(3)))
  !read
  call h5Dread_f(dset_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading dataset  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Dclose_f(dset_id, h5err)

end subroutine read_3d_int_hdf5_al

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ==== Write Attributes ====
!-----------------------------------------------------------------------

!> Write a single real value to an attribute
!!
subroutine write_real_hdf5_attr(outdata, outname, loc_id)
 real(wp), intent(in)           :: outdata
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: loc_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_real_hdf5_attr'
 integer :: h5err
 logical :: l_exists

  !create a scalar dataspace
  out_size(1) = 1
  call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(.not.l_exists) then
    !create the dataset on the file
    call h5Acreate_f(loc_id, trim(outname), h5t_file_float, dspace_id, &
                    attr_id, h5err)
  else
    call h5Aopen_f(loc_id, trim(outname), attr_id, h5err)
  endif

  !write
  call h5Awrite_f(attr_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_real_hdf5_attr

!-----------------------------------------------------------------------

!> Write a rank 1 array of reals to an attribute
!!
subroutine write_1d_real_hdf5_attr(outdata, outname, loc_id)
 real(wp), intent(in)           :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: loc_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_1d_real_hdf5_attr'
 integer :: h5err
 logical :: l_exists

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(l_exists) then
    call h5Adelete_f(loc_id, trim(outname), h5err)
  endif
  !create the dataset on the file
  call h5Acreate_f(loc_id, trim(outname), h5t_file_float, dspace_id, &
                   attr_id, h5err)
  !write
  call h5Awrite_f(attr_id, h5t_mem_float, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_1d_real_hdf5_attr

!-----------------------------------------------------------------------

!> Write an integer to an attribute
!!
subroutine write_int_hdf5_attr(outdata, outname, loc_id)
 integer, intent(in)           :: outdata
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: loc_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_1d_int_hdf5'
 integer :: h5err
 logical :: l_exists

  !create a scalar dataspace
  out_size(1) = 1
  call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(.not.l_exists) then
    !create the dataset on the file
    call h5Acreate_f(loc_id, trim(outname), h5t_file_int, dspace_id, &
                    attr_id, h5err)
  else
    call h5Aopen_f(loc_id, trim(outname), attr_id, h5err)
  endif

  !write
  call h5Awrite_f(attr_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_int_hdf5_attr

!-----------------------------------------------------------------------

!> Write a rank 1 array of integers to an attribute
!!
subroutine write_1d_int_hdf5_attr(outdata, outname, loc_id)
 integer, intent(in)           :: outdata(:)
 character(len=*), intent(in)  :: outname
 integer(HID_T), intent(in)    :: loc_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_1d_int_hdf5'
 integer :: h5err
 logical :: l_exists

  !create a 1D dataspace
  out_size(1) = int(size(outdata,1), h5sz)
  call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(l_exists) then
    call h5Adelete_f(loc_id, trim(outname), h5err)
  endif
  !create the dataset on the file
  call h5Acreate_f(loc_id, trim(outname), h5t_file_int, dspace_id, &
                   attr_id, h5err)
  !write
  call h5Awrite_f(attr_id, h5t_mem_int, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
end subroutine write_1d_int_hdf5_attr

!-----------------------------------------------------------------------

!> Write a string to an attribute
!!
subroutine write_string_hdf5_attr(outdata, outname, loc_id)
 character(len=*), intent(in), target :: outdata
 character(len=*), intent(in)         :: outname
 integer(HID_T), intent(in)           :: loc_id

 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_string_hdf5_attr'
 integer :: h5err
 integer ::strlen
 logical :: l_exists

  strlen = len(outdata)
  out_size(1) = 1

  !create the string datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tcopy_f(h5t_file_char, filetype_id, h5err)

  if(strlen .gt. 0) then !create the string type and a scalar space

    call h5Tset_size_f(memtype_id, int(strlen,HSIZE_T),h5err)
    call h5Tset_size_f(filetype_id, int(strlen,HSIZE_T),h5err)

    !create a (single) scalar dataspace
    call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  else !create a 1 rank space of lenth 0

    out_size = 0
    call h5Screate_simple_f(rank, out_size, dspace_id, h5err)

  endif

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(l_exists) then
    call h5Adelete_f(loc_id, trim(outname), h5err)
  endif
  !create the dataset on the file
  call h5Acreate_f(loc_id, trim(outname), filetype_id, dspace_id, &
                   attr_id, h5err)
  !write
  call h5Awrite_f(attr_id, memtype_id, outdata, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
end subroutine write_string_hdf5_attr

!-----------------------------------------------------------------------


!> Write a logical to an attribute
!!
!! WARNING: the hdf5 support for boolean data types is doubtful, so for the
!! moment a single character "T" or "F" is being printed
subroutine write_logical_hdf5_attr(outdata, outname, loc_id)
 logical, intent(in)                  :: outdata
 character(len=*), intent(in)         :: outname
 integer(HID_T), intent(in)           :: loc_id

 character(len=1) :: strout
 integer(HID_T) :: dspace_id!, memspace_id
 integer(HID_T) :: attr_id
 !integer(HID_T) :: filetype_id, memtype_id
 integer, parameter :: rank=1
 integer(HSIZE_T) :: out_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'write_logical_hdf5_attr'
 integer :: h5err
 logical :: l_exists

  out_size(1) = 1

  if (outdata) then
    strout = 'T'
  else
    strout = 'F'
  endif

  !create the string datatypes
  !call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  !call h5Tset_size_f(memtype_id, int(strlen,HSIZE_T),h5err)
  !call h5Tcopy_f(h5t_file_char, filetype_id, h5err)
  !call h5Tset_size_f(filetype_id, int(strlen,HSIZE_T),h5err)

  !create a (single) scalar dataspace
  call h5Screate_f(H5S_SCALAR_F, dspace_id, h5err)

  !check if the dataset on file exists
  call h5Aexists_f(loc_id, trim(outname), l_exists, h5err)

  if(l_exists) then
    call h5Adelete_f(loc_id, trim(outname), h5err)
  endif
  !create the dataset on the file
  call h5Acreate_f(loc_id, trim(outname), h5t_file_char, dspace_id, &
                   attr_id, h5err)
  !write
  call h5Awrite_f(attr_id, h5t_mem_char, strout, out_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems writing attribute '//outname)

  !release resources
  call h5Aclose_f(attr_id, h5err)
  call h5Sclose_f(dspace_id, h5err)
  !call h5Tclose_f(memtype_id, h5err)
  !call h5Tclose_f(filetype_id, h5err)
end subroutine write_logical_hdf5_attr

!-----------------------------------------------------------------------
! ==== Read Attributes ====
!-----------------------------------------------------------------------

!> Read a single real value from an attribute
!!
subroutine read_real_hdf5_attr(input, inputname, loc_id)
 real(wp), intent(out)           :: input
 character(len=*), intent(in)  :: inputname
 integer(HID_T), intent(in)    :: loc_id

 integer(HID_T) :: attr_id
 integer(HSIZE_T) :: in_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'read_real_hdf5_attr'
 integer :: h5err

  !create a scalar dataspace
  in_size = int((/1/), h5sz)

  !open the attribute
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening attribute  '//inputname)
  !read
  call h5Aread_f(attr_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading attribute  '//inputname)
  !release resources
  call h5Aclose_f(attr_id, h5err)

end subroutine read_real_hdf5_attr

!-----------------------------------------------------------------------

!> Read a rank 1 array of reals from an attribute
!!
subroutine read_1d_real_hdf5_attr(input, inputname, loc_id)
 real(wp), intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_real_hdf5_attr'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: attr_id
 integer :: h5err

  !open the dataset
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening attribute  '//inputname)
  !get its dataspace
  call h5Dget_space_f(attr_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Attribute rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if(in_size(1) .ne. int(size(input,1), h5sz)) &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Aread_f(attr_id, h5t_mem_float, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading attribute  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Aclose_f(attr_id, h5err)

end subroutine read_1d_real_hdf5_attr

!-----------------------------------------------------------------------

!> Read an integer from an attribute
!!
subroutine read_int_hdf5_attr(input, inputname, loc_id)
 integer, intent(out) :: input
 character(len=*), intent(in) :: inputname
 integer(HID_T), intent(in) :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_int_hdf5_attr'
 integer(HSIZE_T) :: in_size(1)
 integer(HID_T) :: attr_id
 integer :: h5err

  in_size = int((/1/), h5sz)
  !open the attribute
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening attribute  '//inputname)
  !read
  call h5Aread_f(attr_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading attribute  '//inputname)
  !release resources
  call h5Aclose_f(attr_id, h5err)

end subroutine read_int_hdf5_attr

!-----------------------------------------------------------------------

!> Read a rank 1 array of integers from an attribute
!!
subroutine read_1d_int_hdf5_attr(input, inputname, loc_id)
 integer, intent(out) :: input(:)
 character(len=*), intent(in)       :: inputname
 integer(HID_T), intent(in)         :: loc_id

 character(len=*), parameter :: &
    this_sub_name = 'read_1d_int_hdf5_attr'
 integer(HSIZE_T) :: in_size(1), max_size(1)
 integer :: in_rank
 integer(HID_T) :: filespace_id
 integer(HID_T) :: attr_id
 integer :: h5err

  !open the attribute
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening attribute  '//inputname)
  !get its dataspace
  call h5Aget_space_f(attr_id, filespace_id, h5err)
  !get dataspace rank and dimensions
  call h5Sget_simple_extent_ndims_f(filespace_id, in_rank, h5err)
  if (in_rank .ne. 1) call error(this_mod_name,this_sub_name, &
     'Dataset rank different from the requested 1 for data ' &
     //trim(inputname))
  call h5Sget_simple_extent_dims_f(filespace_id, in_size, &
                                   max_size, h5err)
  if(in_size(1) .ne. int(size(input,1), h5sz))  &
      call  error(this_mod_name, this_sub_name, &
      'Dimension of read file and target array differs, for array '//inputname)
  !read
  call h5Aread_f(attr_id, h5t_mem_int, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading attribute  '//inputname)

  !release resources
  call h5Sclose_f(filespace_id, h5err)
  call h5Aclose_f(attr_id, h5err)

end subroutine read_1d_int_hdf5_attr

!-----------------------------------------------------------------------

!> Read a string from an attribute
!!
subroutine read_string_hdf5_attr(input, inputname, loc_id)
 character(len=*), intent(out)        :: input
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 !integer(HID_T) :: dspace_id, memspace_id
 integer(HID_T) :: attr_id
 integer(HID_T) :: filetype_id, memtype_id
 integer(HSIZE_T) :: in_size(1)
 character(len=*), parameter :: &
    this_sub_name = 'read_string_hdf5_attr'
 integer :: h5err
 integer(HSIZE_T) ::strlen

  in_size(1) = 1
  !open the attribute
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  !get the type
  call h5Aget_type_f(attr_id, filetype_id, h5err)
  !get the size
  call h5Tget_size_f(filetype_id, strlen, h5err)

  !Check the size
  if (int(strlen) .gt. len(input)) call error(this_mod_name, this_sub_name, &
     'The string'//trim(inputname)//' about to be read is shorter than the&
     &data on the file')

  !create the memory datatypes
  call h5Tcopy_f(h5t_mem_char, memtype_id, h5err)
  call h5Tset_size_f(memtype_id, strlen,h5err)

  !Prepare the local string: it might be filled with random data, and
  !since we read only the first part is necessary to clear it to avoid
  !having garbage at the end of the string
  input = ''
  !read
  !WARNING: no dimensions check for the moment
  call h5Aread_f(attr_id, memtype_id, input, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                           'Problems reading attribute '//inputname)

  !release resources
  call h5Tclose_f(memtype_id, h5err)
  call h5Tclose_f(filetype_id, h5err)
  call h5Aclose_f(attr_id, h5err)
end subroutine read_string_hdf5_attr

!-----------------------------------------------------------------------


!> Read a logical from an attribute
!!
!! WARNING: the hdf5 support for boolean data types is doubtful, so for the
!! moment a single character "T" or "F" is being read
subroutine read_logical_hdf5_attr(input, inputname, loc_id)
 logical, intent(out)                 :: input
 character(len=*), intent(in)         :: inputname
 integer(HID_T), intent(in)           :: loc_id

 character(len=1) :: strin

 character(len=*), parameter :: &
    this_sub_name = 'read_logical_hdf5_attr'

 integer(HSIZE_T) :: in_size(1)
 integer(HID_T) :: attr_id
 integer :: h5err


  in_size = int((/1/), h5sz)
  !open the attribute
  call h5Aopen_f(loc_id, inputname, attr_id, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems opening attribute  '//inputname)
  !read
  call h5Aread_f(attr_id, h5t_mem_char, strin, in_size, h5err)
  if(h5err<0) call error(this_mod_name,this_sub_name, &
                               'Problems reading attribute  '//inputname)
  !release resources
  call h5Aclose_f(attr_id, h5err)

  select case(strin)

   case('t','T')
    input = .true.

   case('f','F')
    input = .false.

   case default
    call error(this_mod_name,this_sub_name, &
              'Invalid logical hdf5 read while readin atribute  '//inputname)

  end select

end subroutine read_logical_hdf5_attr

!-----------------------------------------------------------------------
end module mod_hdf5_io
