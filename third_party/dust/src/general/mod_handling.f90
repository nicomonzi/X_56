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

!> Module to handle basic execution of the code
module mod_handling

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

!-----------------------------------------------------------------------

use iso_fortran_env, only: &
  output_unit, error_unit

!$ use omp_lib, only: &
!$   omp_get_wtime

!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------

! Module interface

public :: &
  dust_abort, &
  pure_abort, &
  dust_time, &
  error,   &
  internal_error,   &
  warning, &
  info, &
  printout, &
  unit_stdout, &
  unit_errout, &
  t_realtime, &
  new_file_unit, &
  check_preproc, &
  check_basename, &
  check_file_exists

private

!-----------------------------------------------------------------------

! Module types and parameters

integer, parameter ::  erru = error_unit
integer, parameter ::  outu = output_unit
integer, parameter :: unit_stdout = outu
integer, parameter :: unit_errout = erru

!character(len=*), parameter :: nl = achar(10) ! new-line

integer, parameter :: t_realtime = &
  !$ kind(0.0d0) ! default double: omp_get_wtime()
  !$ integer, parameter :: omp_dummy1 = &
  kind(0.0) ! default real: for the standard cpu_time


! File handling
integer, parameter :: &
  low_limit = 10, &
  top_limit = 99
integer :: actual_unit_number = low_limit




character(len=*), parameter :: &
  this_mod_name = 'mod_handling'

interface error
  module procedure error, error_multiline, error_main
end interface
interface warning
  module procedure warning, warning_multiline
end interface
interface info
  module procedure info, info_multiline
end interface
interface printout
  module procedure printout
end interface

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

!> Stop the execution
subroutine dust_abort()

  stop
  
end subroutine dust_abort

!-----------------------------------------------------------------------

!> Stop the execution, can be called in a pure procedure
pure subroutine pure_abort()

  ! Something to (hopefully) trigger a core dump
  real :: x
  x = x/0.0

end subroutine pure_abort

!-----------------------------------------------------------------------


!> Get time (only time differences are meaningful)
!!
!! Should work both within OMP and in normal execution
function dust_time() result(t)
 real(t_realtime) :: t

  !$ t = omp_get_wtime()

!$ contains

 !$ subroutine omp_dummy2(t)
 !$ real(omp_dummy1), intent(out) :: t
 intrinsic :: cpu_time

  call cpu_time(t)
 !$ end subroutine omp_dummy2

end function dust_time

!-----------------------------------------------------------------------

!> Error: abort the execution, printing a message of the abort location
!! and reason
recursive subroutine error(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'error'


  gfortran_bug_text(1) = trim(text)
  call put_msg( 1 , format_out_string(                            &
      !'ERROR in "', caller , caller_mod , '"!' , (/trim(text)/) ) )
      'ERROR in "', caller , caller_mod , '"!' , gfortran_bug_text ) )

  call dust_abort()

end subroutine error

!-----------------------------------------------------------------------

!> Internal Error: abort the execution, printing a witty message of the
!! abort location and reason.
!! To be used for internal checks that should never fail
recursive subroutine internal_error(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: message_of_the_day = &
    nl//'What you experienced never happened. Please report this error in &
    &order to let us dispatch our highly trained space raccoons to overwrite &
    &this timeline'
 character(len=len_trim(text)+len(message_of_the_day)) :: another_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'internal_error'


  gfortran_bug_text(1) = trim(text)
  another_text(1) = gfortran_bug_text(1)//message_of_the_day
  !call put_msg( 1 , format_out_string(                            &
  !    'Interna ERROR in "', caller , caller_mod , '"!' , &
  !    gfortran_bug_text//nl//'This should have never happened, please &
  !    &report this error so that  a team of professionals could travel &
  !    &back in time and fix this' ) )

  call put_msg( 1 , format_out_string(                            &
      !'ERROR in "', caller , caller_mod , '"!' , (/trim(text)/) ) )
      'ERROR in "', caller , caller_mod , '"!' , another_text ) )
  call dust_abort()

end subroutine internal_error

!-----------------------------------------------------------------------

!> Error: abort the execution, printing a message of the abort location
!! and reason, without the caller module
recursive subroutine error_main(caller,text)
 character(len=*), intent(in) :: &
   caller, text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'error'


  gfortran_bug_text(1) = trim(text)
  call put_msg( 1 , format_out_string_main(                            &
      !'ERROR in "', caller , caller_mod , '"!' , (/trim(text)/) ) )
      'ERROR in "', caller , '"!' , gfortran_bug_text ) )

  call dust_abort()

end subroutine error_main


!-----------------------------------------------------------------------

!> Error: abort the execution, printing a message of the abort location
!! and reason, with multiline message
subroutine error_multiline(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text(:)

 character(len=*), parameter :: &
   this_sub_name = 'error'

  call put_msg( 1 , format_out_string(                  &
      'ERROR in "', caller , caller_mod , '"!' , text ) )

  call dust_abort()

end subroutine error_multiline

!-----------------------------------------------------------------------

!> Warning: print a scary warning, with location and reason, but do not
!! abort the execution
subroutine warning(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'warning'


  gfortran_bug_text(1) = trim(text)
  call put_msg( 1 , format_out_string(                              &
      !'WARNING in "', caller , caller_mod , '"!' , (/trim(text)/) ) )
      'WARNING in "', caller , caller_mod , '"!' , gfortran_bug_text ) )

end subroutine warning

!-----------------------------------------------------------------------

!> Warning: print a scary warning, with location and reason, but do not
!! abort the execution
subroutine warning_multiline(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text(:)

 character(len=*), parameter :: &
   this_sub_name = 'error'


  call put_msg( 1 , format_out_string(                    &
      'WARNING in "', caller , caller_mod , '"!' , text ) )

end subroutine warning_multiline

!-----------------------------------------------------------------------

!> Info: print a detailed info message, with information on location
subroutine info(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'info'

  gfortran_bug_text(1) = trim(text)
  call put_msg( 2 , format_out_string(                             &
      !'INFO from "', caller , caller_mod , '":' , (/trim(text)/) ) )
      'INFO from "', caller , caller_mod , '":' , gfortran_bug_text ) )

end subroutine info

!-----------------------------------------------------------------------

!> Info: print a detailed info message, with information on location
subroutine info_multiline(caller,caller_mod,text)
 character(len=*), intent(in) :: &
   caller, caller_mod, text(:)

 character(len=*), parameter :: &
   this_sub_name = 'error'

  call put_msg( 2 , format_out_string(                   &
      'INFO from "', caller , caller_mod , '":' , text ) )

end subroutine info_multiline

!-----------------------------------------------------------------------

!> Print a message to standard output without any additional information
subroutine printout(text)
 character(len=*), intent(in) :: &
    text

 ! compiler bug -- gfortran bug
 ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77649
 character(len=len_trim(text)) :: gfortran_bug_text(1)
 character(len=*), parameter :: &
   this_sub_name = 'info'

  gfortran_bug_text(1) = trim(text)
  call put_msg( 2 , gfortran_bug_text(1)  )

end subroutine printout

!-----------------------------------------------------------------------

!> Make a string and count its length
function format_out_string(msg,caller,caller_mod,end, text ) result(s)
 character(len=*), intent(in) :: msg, caller, caller_mod, end
 character(len=*), intent(in) :: text(:)
 integer :: i
 character( &
   len = 1+len(msg)+len(caller)+14+len(caller_mod)+len(end)    &
        +1+sum((/( 2+len_trim(text(i))+1 , i=1,size(text) )/)) &
          ) :: s

  !DEBUG:
  ! See ifort bug
  ! DPD200411297
  ! https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/629432
  ! and also
  ! https://software.intel.com/en-us/comment/1878150#comment-1878150
  write(s,'(a,*(a))') &
      nl//msg // caller // '", in module "' // caller_mod // end &
      // nl, ( "  "//trim(text(i))//nl ,i=1,size(text))

end function format_out_string

!-----------------------------------------------------------------------

!> Make a string and count its length, without module
function format_out_string_main(msg,caller,end, text ) result(s)
 character(len=*), intent(in) :: msg, caller, end
 character(len=*), intent(in) :: text(:)
 integer :: i
 character( &
   len = 1+len(msg)+len(caller)+len(end)    &
        +1+sum((/( 2+len_trim(text(i))+1 , i=1,size(text) )/)) &
          ) :: s

  ! See ifort bug
  ! DPD200411297
  ! https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/629432
  ! and also
  ! https://software.intel.com/en-us/comment/1878150#comment-1878150
  write(s,'(a,*(a))') &
      nl//msg // caller // end &
      // nl, ( "  "//trim(text(i))//nl ,i=1,size(text))

end function format_out_string_main

!-----------------------------------------------------------------------

!> Write a message to the selcted channel
subroutine put_msg( msg_type , msg )
 integer, intent(in) :: msg_type
 character(len=*), intent(in) :: msg

 character(len=*), parameter :: &
   this_sub_name = 'put_msg'

  msg_type_case: select case(msg_type)

   case(1) ! error

      write(erru,'(a)') msg
      !without flush call, sometimes, gfortran does not output to the
      !correct file if the output is redirected
      call flush(erru)

   case(2) ! info

      write(outu,'(a)') msg
      !without flush call, sometimes, gfortran does not output to the
      !correct file if the output is redirected
      call flush(outu)

   case default
    call error(this_sub_name,this_mod_name,         &
      'Internal module error: unknown message type.')
  end select msg_type_case

contains

 pure subroutine append_buffered_text(idx,buf,msg)
  character(len=*), intent(in) :: msg
  integer,          intent(inout) :: idx
  character(len=*), intent(inout) :: buf

   idx = len_trim(buf) + 1
   if( idx+len(msg)-1 .gt. len(buf) ) then ! reset the buffer
     buf = ''
     idx = 1
   endif
   buf(idx:) = msg ! msg is truncated if too long
 end subroutine append_buffered_text

end subroutine put_msg

!-----------------------------------------------------------------------

!> Return a new, not busy, file unit for input/output
subroutine new_file_unit(fu,ierr)
 integer, intent(out) :: fu, ierr

 character(len=*), parameter :: &
   this_sub_name = 'new_file_unit'
 logical :: busy
 integer :: previous_unit_number

  !$omp critical(omp_new_file_unit_critical)
  previous_unit_number = actual_unit_number
  busy = .true.
  do
   actual_unit_number = low_limit+mod(actual_unit_number+1-low_limit, &
                                      top_limit+1-low_limit)
   inquire(unit=actual_unit_number,opened=busy)
   if((busy.eqv..false.).or. & ! either a free unit or no available units
      (actual_unit_number.eq.previous_unit_number)) exit
  enddo

  if(busy.eqv..false.) then ! we have been able to find an available unit
    ierr = 0
    fu = actual_unit_number
  else
    ierr = 1
    fu = -1
    call warning(this_sub_name,this_mod_name, &
                 'No available file units.')
  endif

  !$omp end critical(omp_new_file_unit_critical)

end subroutine new_file_unit

!-----------------------------------------------------------------------

!> Check if the basename provided is legitimate
!!
!! Attempts to create a test file with such basename, if it fails,
!! aborts with error
subroutine check_basename(basename,sub_name,mod_name)
 character(len=*), intent(in) :: basename
 character(len=*), intent(in) :: sub_name
 character(len=*), intent(in), optional :: mod_name
 integer :: estat, cstat

  !try to print a test file in the provided basename
  call execute_command_line('touch '//trim(basename)//'.try > /dev/null 2>&1 ', &
                                          exitstat=estat,cmdstat=cstat)

  if( estat.ne.0 .or. cstat.ne.0) then
    !got a problem with the basename
    call error(trim(sub_name),'Problems creating output with basename: '&
         & //trim(basename)//', possibly invalid path; please check if folder exists &
         & and is writable.')

  else
    !the basename is good, remove the test file
    call execute_command_line('rm '//trim(basename)//'.try', &
                                          exitstat=estat,cmdstat=cstat)
  endif

end subroutine check_basename

!-----------------------------------------------------------------------

!> Check if one file exists (before opening it), if not returns an error
!!
subroutine check_file_exists(filename,sub_name,mod_name)
 character(len=*), intent(in) :: filename
 character(len=*), intent(in) :: sub_name
 character(len=*), intent(in), optional :: mod_name
 logical :: exists

  !Inquire file existence
  inquire(file=trim(filename), exist=exists)

  if( .not. exists) then
    !got a problem with the basename
    if(present(mod_name)) then
      call error(trim(sub_name),trim(mod_name),'The file "'//trim(filename) &
      //'" about to be read does not exist, check that input paths and &
      &commands are correct.')
    else
      call error(trim(sub_name),'The file "'//trim(filename) &
      //'" about to be read does not exist, check that input paths and &
      &commands are correct.')
    endif
  endif

end subroutine check_file_exists

!-----------------------------------------------------------------------

!> Check if the geometry files exists, if not, try to generate it
!! calling the preprocessor
subroutine check_preproc(filename)
 character(len=*), intent(in) :: filename

 character(len=*), parameter :: this_sub_name = 'check_preproc'
 character(len=max_char_len) :: msg
 logical :: exists
 integer :: estat, cstat

  !1) Check if the solver geometry input file exists
  inquire(file=trim(filename), exist=exists)

  if(exists) then
    !2) exists, return
    return
  else
    !3) does not exist, try to call the preprocessor
    write(msg,*) '!!!! The geometry input file ',trim(filename),' was not found,&
    & trying to call the preprocessor with default parameters to create it'
    call printout(msg)
    call execute_command_line('./dust_pre dust_pre.in '//trim(filename), &
                                                  exitstat=estat,cmdstat=cstat)

    !4) check the execution
    if(cstat .ne. 0) call error(this_sub_name, this_mod_name, 'Call to the &
                     &preprocessor not executed, system error')
    if(estat .ne. 0) call error(this_sub_name, this_mod_name, 'Call to the &
                     &preprocessor failed, stopping')

    !5) check that the file was actually created
    inquire(file=trim(filename), exist=exists)

    if(.not.exists) then
      write(msg,*) 'The execution of the preprocessor was succesful, however &
           & the input file was not created. Must stop the execution now.'
      call error(this_sub_name, this_mod_name, msg)
    endif

  endif

end subroutine check_preproc

!-----------------------------------------------------------------------

end module mod_handling

