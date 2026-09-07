!!=====================================================================
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
!! This file was originally part of Flexi, and has been (minimally) modified
!! to suit the needs of DUST
!!
!!=====================================================================
!
!
!Original copyright:
!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of flexi, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! flexi is free software: you can redistribute it and/or modify it under the terms of the gnu General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! flexi is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty
! of merchantability or fitness for A particular purpose. See the gnu General Public License v3.0 for more details.
!
! You should have received a copy of the gnu General Public License along with flexi. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!
! attention:
! The routines 'clear_formatting', 'set_formatting', 'get_escape_sequence' and 'split_string' are copied from the fortran output
! library (foul).
! Copyright and license see below.
! Full version of the foul-library can be found here:
!   http://foul.sourceforge.net
!
!------------------------------------------------------------
! foul - The Fortran Output Library
!------------------------------------------------------------
! Provides routines enabling Fortran programs to:
!
! - Write formatted console output using ansi escape codes
!   (see http://en.wikipedia.org/wiki/ANSI_escape_code
!    for more information)
! - Convert numbers to strings in a variety of
!   finely-controlled ways, including fully trimmed
!   of whitespace characters
! - Time program execution and other processes with the
!   highest accuracy provided by the system
!------------------------------------------------------------
! Copyright (C) 2010-2011 by Philipp Emanuel Weidmann
! E-Mail: philipp.weidmann@gmx.de
!------------------------------------------------------------
! This library is free software; you can redistribute it
! and/or modify it under the terms of the gnu General Public
! License as published by the Free Software Foundation;
! version 3 of the License.
!
! This library is distributed in the hope that it will be
! useful, but without any warranty; without even the implied
! warranty of merchantability or fitness for A particular
! purpose.
!
! See the gnu General Public License for more details.
!------------------------------------------------------------


!> Routines for performing operations on strings, which are not covered by
!! ISO_VARYING_STRING.
!!
!! The routines 'clear_formatting', 'set_formatting', 'get_escape_sequence'
!! and 'split_string' are copied from the fortran output library (foul).
!! Full version of the foul-library can be found here:
!!   http://foul.sourceforge.net

module MOD_StringTools

use mod_handling, only: &
  error, warning, info, unit_stdout, printout

use mod_param, only: &
  max_char_len, wp

use MOD_ISO_VARYING_STRING

!----------------------------------------------------------------------

implicit none
private

interface LowCase
  module procedure LowCase
  module procedure LowCase_overwrite
end interface

interface stricmp
  module procedure stricmp
end interface

interface StripSpaces
  module procedure StripSpaces
end interface

interface inttostr
  module procedure inttostr
end interface

interface isint
  module procedure isint
end interface

interface isreal
  module procedure isreal
end interface

interface isempty
  module procedure isempty
end interface

interface set_formatting
  module procedure set_formatting
end interface

interface clear_formatting
  module procedure clear_formatting
end interface

public :: LowCase
public :: stricmp
public :: StripSpaces
public :: inttostr
public :: isint
public :: isreal
public :: isempty
public :: set_formatting
public :: clear_formatting
public :: IsInList, strip_mult_appendix

logical :: use_escape_codes = .true.  !< If set to .false., output will consist only of standard text, allowing the
                                      !< escape characters to be switched off in environments which don't support them.
public :: use_escape_codes

character(len=*), parameter :: this_mod_name = 'mod_stringtools'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Transform upper case letters in "Str1" into lower case letters, result is "Str2" (in place version)
subroutine LowCase_overwrite(Str1)
! modules
implicit none
! input / output variables
character(len=*),intent(inout) :: Str1  !< Input/output string
! local variables
integer                    :: iLen,nLen,Upper
character(len=*),parameter :: lc='abcdefghijklmnopqrstuvwxyz'
character(len=*),parameter :: uc='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
logical                    :: HasEq
HasEq=.false.
nLen=LEN_TRIM(Str1)
do iLen=1,nLen
  ! Transformation stops at "="
  if(Str1(iLen:iLen).eq.'=') HasEq=.true.
  Upper=index(uc,Str1(iLen:iLen))
  if ((Upper > 0).and. .not. HasEq) then
    Str1(iLen:iLen) = lc(Upper:Upper)
  end if
end do
end subroutine LowCase_overwrite

!----------------------------------------------------------------------

!> Transform upper case letters in "Str1" into lower case letters, result is "Str2"
subroutine LowCase(Str1,Str2)
! modules
implicit none
! input / output variables
character(len=*),intent(in)  :: Str1  !< Input string
character(len=*),intent(out) :: Str2  !< Output string, lower case letters only
! local variables
integer                    :: iLen,nLen,Upper
character(len=*),parameter :: lc='abcdefghijklmnopqrstuvwxyz'
character(len=*),parameter :: uc='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
logical                    :: HasEq
HasEq=.false.
Str2=Str1
nLen=LEN_TRIM(Str1)
do iLen=1,nLen
  ! Transformation stops at "="
  if(Str1(iLen:iLen).eq.'=') HasEq=.true.
  Upper=index(uc,Str1(iLen:iLen))
  if ((Upper > 0).and. .not. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
end do
end subroutine LowCase

!----------------------------------------------------------------------

!> Case insensitive string comparison
function stricmp(a, b)
! modules
implicit none
! input/output variables
character(len=*),intent(in) :: a,b !< strings to compare with each other
! local variables
logical            :: stricmp
character(len=255) :: alow
character(len=255) :: blow
call LowCase(a, alow)
call LowCase(b, blow)
stricmp = (trim(alow).eq.trim(blow))
end function stricmp

!----------------------------------------------------------------------

function IsInList(str, str_list, pos) result(isin)
 character(len=*), intent(in) :: str
 character(len=*), allocatable, intent(in) :: str_list(:)
 integer, optional :: pos
 logical :: isin

 integer :: i

  if (present(pos)) pos = 0

  isin = .false.
  if(allocated(str_list)) then
   do i = lbound(str_list,1), ubound(str_list,1)
     if(stricmp(trim(str),trim(str_list(i)))) then
       if (present(pos))  pos = i
       isin = .true.
       return
     endif
   enddo
  endif

end function

!----------------------------------------------------------------------

!> Removes all whitespace from a string
subroutine StripSpaces(string)
! modules
implicit none
! input / output variables
character(len=*),intent(inout) :: string  !< input string
! local variables
integer :: stringLen
integer :: last, actual
stringLen = len(string)
last = 1
actual = 1
do while (actual < stringLen)
  if(string(last:last) == ' ') then
    actual = actual + 1
    string(last:last) = string(actual:actual)
    string(actual:actual) = ' '
  else
    last = last + 1
    if (actual < last) &
        actual = last
  endif
end do
end subroutine

!----------------------------------------------------------------------

!> Converts integer to string
pure function inttostr(value)
integer,intent(in)  :: value
character(len=255)  :: inttostr
write(inttostr,"(I20)") value
end function inttostr

!----------------------------------------------------------------------

!> Checks if a string is an integer
pure function isint(value)
character(len=*),intent(in)  :: value
logical                        :: isint
! local variables
integer                    :: i,stat
read(value, *, iostat=stat) i
isint=(stat.eq.0)
end function isint

!----------------------------------------------------------------------

!> Checks if a string is a real
pure function isreal(value)
character(len=*),intent(in)  :: value
logical                      :: isreal
! local variables
integer  :: stat
real(wp) :: i
read(value, *, iostat=stat) i
isreal=(stat.eq.0)
end function isreal

!----------------------------------------------------------------------

!> Checks if a string is a empty
pure function isempty(value)
character(len=*),intent(in)  :: value
logical                      :: isempty
isempty = ( LEN_TRIM(value) == 0 )
end function isempty

!----------------------------------------------------------------------

!> Splits the supplied string along a delimiter.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
subroutine split_string(string, delimiter, substrings, substring_count)
! modules
implicit none
! input / output variables
character (len = *), intent(in)  :: string          !< Variable-length character string that is to be split
character (len = *),           intent(in)  :: delimiter       !< Character along which to split
character (len = *), intent(out) :: substrings(*)   !< Array of substrings generated by split operation
integer,             intent(out) :: substring_count !< Number of substrings generated
! local variables
integer :: start_position, end_position, len_del
start_position  = 1
substring_count = 0
len_del = len(trim(delimiter))

do
  end_position = index(string(start_position:), delimiter)
  substring_count = substring_count + 1

  if (end_position == 0) then
    substrings(substring_count) = string(start_position:)
    exit
  else
    substrings(substring_count) = string(start_position : start_position + end_position - 2)
    start_position = start_position + end_position + len_del-1
  end if
end do
end subroutine split_string

!----------------------------------------------------------------------

subroutine strip_mult_appendix(string, rem, delimiter)
 character(len=*), intent(in) :: string
 character(len=*), intent(in) :: delimiter
 character(len=*), intent(out) :: rem

 character(len=max_char_len) :: substrings(5)
 integer :: nstr, striplen, nlen
 character(len=*), parameter :: this_sub_name = 'strip_mult_appendix'

 call split_string(trim(string), trim(delimiter), substrings, nstr)


 if (nstr .gt. 3) then
   call warning(this_sub_name, this_mod_name, 'More than one delimiter &
   &"'//trim(delimiter)//'" was found in string "'//trim(string)//'". Assuming&
   & that the delimiters previous to the last were part of an input string.&
   & This is a dengerous situation and could lead to unexpected behaviours, the&
   & user is encouraged to avoid the use of the aforementioned delimiter in&
   & input strings.')
 endif

 if (.not.isint(trim(substrings(nstr))) .and. nstr .gt. 1) then
   call warning(this_sub_name, this_mod_name, 'One or more instances of the &
   &delimiter "'//trim(delimiter)//'" were found in strin "'//trim(string)//'"&
   & without delimiting a multiple string numeral. &
   & This is a dengerous situation and could lead to unexpected behaviours, the&
   & user is encouraged to avoid the use of the aforementioned delimiter in&
   & input strings.')
 endif

 if (isint(trim(substrings(nstr)))) then
   striplen = len(trim(substrings(nstr)))
   striplen = striplen + len(trim(delimiter))
 else
   striplen = 0
 endif

 nlen = len(trim(string))-striplen
 rem = string(1:nlen)

end subroutine strip_mult_appendix

!----------------------------------------------------------------------

!> Generates an ansi escape sequence from the supplied style string.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
subroutine get_escape_sequence(style_string, escape_sequence)
! modules
implicit none
! input / output variables
character(len=*),intent(in)   :: style_string    !< String describing which styles to set (separated by space)
                                                 !< see source code for supported styles
character(len=16),intent(out) :: escape_sequence !< escape_sequence: ansi escape sequence generated from the specified styles
integer           :: i
character(len=32) :: style_substrings(16)
integer           :: style_substring_count
! Start sequence with command to clear any previous attributes
escape_sequence = char(27) // '[0'

call split_string(trim(style_string), ' ', style_substrings, style_substring_count)

do i = 1, style_substring_count
  call LowCase(style_substrings(i))

  select case (trim(style_substrings(i)))
  case ('bright')
    escape_sequence = trim(escape_sequence) // ';1'
  case ('faint')
    escape_sequence = trim(escape_sequence) // ';2'
  case ('italic')
    escape_sequence = trim(escape_sequence) // ';3'
  case ('underline')
    escape_sequence = trim(escape_sequence) // ';4'
  case ('blink_slow')
    escape_sequence = trim(escape_sequence) // ';5'
  case ('blink_fast')
    escape_sequence = trim(escape_sequence) // ';6'
  case ('black')
    escape_sequence = trim(escape_sequence) // ';30'
  case ('red')
    escape_sequence = trim(escape_sequence) // ';31'
  case ('green')
    escape_sequence = trim(escape_sequence) // ';32'
  case ('yellow')
    escape_sequence = trim(escape_sequence) // ';33'
  case ('blue')
    escape_sequence = trim(escape_sequence) // ';34'
  case ('magenta')
    escape_sequence = trim(escape_sequence) // ';35'
  case ('cyan')
    escape_sequence = trim(escape_sequence) // ';36'
  case ('white')
    escape_sequence = trim(escape_sequence) // ';37'
  case ('background_black')
    escape_sequence = trim(escape_sequence) // ';40'
  case ('background_red')
    escape_sequence = trim(escape_sequence) // ';41'
  case ('background_green')
    escape_sequence = trim(escape_sequence) // ';42'
  case ('background_yellow')
    escape_sequence = trim(escape_sequence) // ';43'
  case ('background_blue')
    escape_sequence = trim(escape_sequence) // ';44'
  case ('background_magenta')
    escape_sequence = trim(escape_sequence) // ';45'
  case ('background_cyan')
    escape_sequence = trim(escape_sequence) // ';46'
  case ('background_white')
    escape_sequence = trim(escape_sequence) // ';47'
  end select
end do

! Append end of sequence marker
escape_sequence = trim(escape_sequence) // 'm'
end subroutine get_escape_sequence

!----------------------------------------------------------------------

!> Sets output formatting to the supplied styles.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
subroutine set_formatting(style_string)
! modules
implicit none
! input / output variables
character (len = *), intent(in) :: style_string !< String describing which styles to set (separated by space).
                                                !< See get_escape_sequence for supported styles.
! local variables
character(len=16) :: escape_sequence
character         :: escape_sequence_array(16)
equivalence         (escape_sequence, escape_sequence_array)
character(len=16) :: format_string
character(len=16) :: istring
integer           :: i
if (use_escape_codes) then
  call get_escape_sequence(trim(style_string), escape_sequence)

  write(istring,*) LEN_TRIM(escape_sequence)
  format_string = '(' // trim(istring) // 'A1)'

  write(UNIT_stdOut, trim(format_string), advance="no") (escape_sequence_array(i), i = 1, LEN_TRIM(escape_sequence))
end if
end subroutine set_formatting

!----------------------------------------------------------------------

!> Resets output formatting to normal.
!> This function is copied from the fortran output library (foul). For the full version of this library see:
!>   http://foul.sourceforge.net
subroutine clear_formatting()
implicit none
if (use_escape_codes) then
  ! Clear all previously set styles
  write(UNIT_stdOut, '(3A1)', advance="no") (/ char(27), '[', 'm' /)
end if
end subroutine clear_formatting

!----------------------------------------------------------------------

end module MOD_StringTools
