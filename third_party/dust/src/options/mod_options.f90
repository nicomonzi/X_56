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
!! It was originally part of Flexi, and has been modified to suit the
!! needs of DUST
!!
!! Authors:
!!          Matteo Tugnoli
!!=====================================================================

! Original copyright:
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

!> Option-classes for values that are read from the parameter file
!! (integer,logical, real, string; each single or array).


!TODO: format all to the standards of Dust
module MOD_Options

use mod_handling , only: UNIT_StdOut, error

use mod_param, only: wp, max_char_len

!use MOD_StringTools

use MOD_StringTools, only: &
  !LowCase         ,&
  stricmp         !,&
  !StripSpaces     ,&
  !inttostr        ,&
  !isint           ,&
  !set_formatting   ,&
  !clear_formatting

!use MOD_ISO_VARYING_STRING
use MOD_ISO_VARYING_STRING, only: &
  varying_string, &
  assignment(=), &
  !operator(//) , &
  !operator(==) , &
  !operator(/=) , &
  !operator(<)  , &
  !operator(<=) , &
  !operator(>=) , &
  !operator(>)  , &
  !adjustl      , &
  !adjustr      , &
  char         , &
  !iachar       , &
  !ichar        , &
  !index        , &
  !len          , &
  len_trim     , &
  !lge          , &
  !lgt          , &
  !lle          , &
  !llt          , &
  !repeat       , &
  !scan         , &
  !trim         , &
  !verify       , &
  var_str      , &
  !get          , &
  !put          , &
  !put_line     , &
  !extract      , &
  !insert       , &
  !remove       , &
  replace      , &
  split
!----------------------------------------------------------------------

implicit none

public  :: option, IntOption, IntArrayOption, IntFromStringOption, &
           LogicalOption, LogicalArrayOption, RealOption, RealArrayOption, &
           StringOption, SubOption

private

character(len=*), parameter :: this_mod_name = 'mod_options'
!----------------------------------------------------------------------

!> Genereal, abstract option
type  :: option
  !< pointer to next option, used for a linked list of options MATTEO: not used
  !class(option),pointer :: next

  !< name of the option, case-insensitive (part before '=' in parameter file)
  character(len=max_char_len)    :: name

  !< comment in parameter file, after '!' character
  character(len=1000)   :: description

  !< section to which the option belongs. Not mandatory
  character(len=max_char_len)    :: section

  !< Is option set? default false. Becomes true, if found in parameter file
  logical               :: isSet

  !< Has option a default value? Default false. True if a default value is
  !! given in CreateXXXOption routine
  logical               :: hasDefault

  !< Is option multiple?  Default false.
  !! Indicates if an option can occur multiple times in parameter file
  logical               :: multiple

  !< Is it an additional option? Default false.
  !! Additional options are the multiple options after the first one
  logical               :: additional

  !< Has option been recovered in the code? Default false.
  !! Indicates if the option is already used (get... call) and therefore is
  !! no longer available in the list of parameters
  logical               :: isRemoved

contains

  !< function used to print option for a default parameter file
  procedure :: print

  !< function used to print the value
  procedure :: printValue

  !< function that parses a string from the parameter file to fill the value of
  !! the option
  procedure :: parse

  !< function that parses a string from the parameter file to fill the value of
  !! the option
  procedure :: parseReal

  !< function to compare case-insensitive a string with the name of this option
  procedure :: nameequals

  !< function that returns the string length of the name
  procedure :: getnamelen

  !< function that returns the string length required to print the value
  procedure :: getvaluelen

end type option

!----------------------------------------------------------------------

!> Integer Option
type, extends(option) :: IntOption
  integer :: value
  integer :: default_value
end type IntOption

!> Integer from String Option
type, extends(option) :: IntFromStringOption
  character(len=max_char_len)              :: value
  character(len=max_char_len)              :: default_value
  integer,allocatable             :: intList(:)
  character(len=max_char_len),allocatable  :: strList(:)
  integer                         :: listIndex
  logical                         :: foundInList = .false.
  integer                         :: maxLength=0
end type IntFromStringOption

!> Integer Array Option
type, extends(option) :: IntArrayOption
  integer,allocatable :: value(:)
  integer,allocatable :: default_value(:)
end type IntArrayOption

!> Logical Option
type, extends(option) :: LogicalOption
  logical :: value
  logical :: default_value
end type LogicalOption

!> Logical Array Option
type, extends(option) :: LogicalArrayOption
  logical,allocatable :: value(:)
  logical,allocatable :: default_value(:)
end type LogicalArrayOption

!> Real Option
type, extends(option) :: RealOption
  real(wp)    :: value
  real(wp)    :: default_value
  !< number of digits, the value has in parameter file
  !! - negative: -number of digits in exponential representation
  !! - 0 means not given
  integer :: digits = 0
end type RealOption

!> Real Array Option
type, extends(option) :: RealArrayOption
  real(wp),allocatable    :: value(:)
  real(wp),allocatable    :: default_value(:)
  !< number of digits, the value has in parameter file
  !! - negative: -number of digits in exponential representation
  !! - 0 means not given
  integer,allocatable :: digits(:)
end type RealArrayOption

!> String Option
type, extends(option) :: StringOption
  character(len=max_char_len) :: value
  character(len=max_char_len) :: default_value
end type StringOption

!> Sub-option: just an empty container of options
type, extends(option) :: SubOption
  character(len=max_char_len) :: value
  character(len=max_char_len) :: default_value
end type SubOption

!> String Array Option
!type,public,extends(option) :: StringArrayOption
  !character(len=max_char_len),allocatable  :: value(:)
!end type StringArrayOption

interface getstrlenreal
  module procedure getstrlenreal
end interface

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Compares name with the name of the option (case-insensitive)
function nameequals(this, name)
 class(option),intent(in)    :: this
 !< incoming name, which is compared with the name of this option
 character(len=*),intent(in) :: name
 logical            :: nameequals
 
  !nameequals = stricmp(this%name, name)
  
  ! compatibility with legacy naming convention 
  nameequals = stricmp(char(Replace(var_str(this%name),"_","",Every=.true.)), &
                      &char(Replace(var_str(name),"_","",Every=.true.)))

end function nameequals

!----------------------------------------------------------------------

!> return string-length of name
function getnamelen(this)
 class(option),intent(in) :: this
 !< length of option name
 integer                  :: getnamelen

  getnamelen = LEN_TRIM(this%name)

end function getnamelen

!----------------------------------------------------------------------

function getvaluelen(this)
 class(option),intent(in)    :: this
 !< string length
 integer                     :: getvaluelen

 integer           :: i
 character(len=50) :: tmp

  getvaluelen = 0  ! default return value

  ! only if option is set of has default value, we have to find the length of
  ! this value
  if ((this%isSet)) then
    ! each class has a different type, which requires different commands
    ! to get the string-length of the value
    select type (this)
    class is (IntOption)
      write(tmp,"(I0)") this%value
      getvaluelen = LEN_TRIM(tmp)
    class is (LogicalOption)
      getvaluelen = 1
    class is (RealOption)
      getvaluelen = getstrlenreal(this%value,this%digits)
    class is (StringOption)
      getvaluelen = LEN_TRIM(this%value)
    class is (IntFromStringOption)
      getvaluelen = this%maxLength
    class is (IntArrayOption)
      getvaluelen = 3 ! '(/ '
      do i=1,size(this%value)
        write(tmp,"(I0)") this%value(i)
        getvaluelen = getvaluelen + LEN_TRIM(tmp)
      end do
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%value)-1)
      ! ' /)'
      getvaluelen = getvaluelen + 3
    class is (LogicalArrayOption)
      getvaluelen = 3 ! '(/ '
      ! each value needs only one character
      getvaluelen = getvaluelen + size(this%value)
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%value)-1)
      getvaluelen = getvaluelen + 3 ! ' /)'
    class is (RealArrayOption)
      getvaluelen = 3 ! '(/ '
      do i=1,size(this%value)
        getvaluelen = getvaluelen + getstrlenreal(this%value(i), this%digits(i))
      end do
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%value)-1)
      getvaluelen = getvaluelen + 3 ! ' /)'
    !class is (StringArrayOption)
      !getvaluelen = 3 ! '(/ '
      !do i=1,size(this%value)
        !getvaluelen = getvaluelen + LEN_TRIM(this%value(i))
      !end do
      ! ', ' between array elements
      !getvaluelen = getvaluelen + 2*(size(this%value)-1)
      !getvaluelen = getvaluelen + 3 ! ' /)'
    class is (SubOption)
      getvaluelen = 2
      !TODO: fix this
    class default
      stop 'Unknown type'
    end select
  else if ((this%hasDefault)) then
    ! each class has a different type, which requires different commands
    ! to get the string-length of the value
    select type (this)
    class is (IntOption)
      write(tmp,"(I0)") this%default_value
      getvaluelen = LEN_TRIM(tmp)
    class is (LogicalOption)
      getvaluelen = 1
    class is (RealOption)
      getvaluelen = getstrlenreal(this%default_value,this%digits)
    class is (StringOption)
      getvaluelen = LEN_TRIM(this%default_value)
    class is (IntFromStringOption)
      getvaluelen = this%maxLength
    class is (IntArrayOption)
      getvaluelen = 3 ! '(/ '
      do i=1,size(this%default_value)
        write(tmp,"(I0)") this%default_value(i)
        getvaluelen = getvaluelen + LEN_TRIM(tmp)
      end do
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%default_value)-1)
      ! ' /)'
      getvaluelen = getvaluelen + 3
    class is (LogicalArrayOption)
      getvaluelen = 3 ! '(/ '
      ! each value needs only one character
      getvaluelen = getvaluelen + size(this%default_value)
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%default_value)-1)
      getvaluelen = getvaluelen + 3 ! ' /)'
    class is (RealArrayOption)
      getvaluelen = 3 ! '(/ '
      do i=1,size(this%default_value)
        getvaluelen = getvaluelen + getstrlenreal(this%default_value(i), this%digits(i))
      end do
      ! ', ' between array elements
      getvaluelen = getvaluelen + 2*(size(this%default_value)-1)
      getvaluelen = getvaluelen + 3 ! ' /)'
    !class is (StringArrayOption)
      !getvaluelen = 3 ! '(/ '
      !do i=1,size(this%default_value)
        !getvaluelen = getvaluelen + LEN_TRIM(this%default_value(i))
      !end do
      ! ', ' between array elements
      !getvaluelen = getvaluelen + 2*(size(this%default_value)-1)
      !getvaluelen = getvaluelen + 3 ! ' /)'
    class is (SubOption)
      getvaluelen = 2
      !TODO: fix this
    class default
      stop 'Unknown type'
    end select
  end if
end function getvaluelen

!----------------------------------------------------------------------

!> Returns length of a real represented as string with a given number of digits
function getstrlenreal(value,digits)
 !< real value to print
 real(wp),intent(in)         :: value
 !< number of digits (if < 1, then print as scientific)
 integer,intent(in)      :: digits
 !< length of real as string
 integer                 :: getstrlenreal
 ! local variables
 character(len=20)       :: fmtDigits
 character(len=50)       :: tmp

  if (digits.ge.1) then ! floating point representation
    write(fmtDigits,*) digits
    write(tmp,'(F0.'//fmtDigits//')') value
    if (index(tmp,'.').eq.1) tmp = '0'//tmp(1:49)
  else if (digits.le.-1) then ! scientific (exponential) representation
    write(fmtDigits,*) -digits
    write(tmp,'(E24.'//fmtDigits//')') value
  else ! digits not given
    write(tmp,'(E24.19)') value
  end if
  getstrlenreal = len(trim(adjustl(tmp)))
end function getstrlenreal

!----------------------------------------------------------------------

!> print option
subroutine print(this, maxNameLen, maxValueLen, mode)
!< option to print
 class(option),intent(in)    :: this
 !< max string length of name
 integer,intent(in)          :: maxNameLen
 !< max string length of name
 integer,intent(in)          :: maxValueLen
 !< 0: during readin, 1: default parameter file, 2: markdown
 integer,intent(in)          :: mode

 character(len=20)    :: fmtName
 character(len=20)    :: fmtValue
 type(VARYING_STRING) :: comment,headNewline,headSpace
 integer              :: length

  write(fmtName,*) maxNameLen
  write(fmtValue,*) maxValueLen
  ! print name
  if (mode.eq.0) then
    write(fmtName,*) maxNameLen
    write(UNIT_stdOut,'(a3)', advance='no')  " | "
    !call set_formatting("blue")
    write(UNIT_stdOut,"(a"//fmtName//")", advance='no') trim(this%name)
    !call clear_formatting()
  else
    write(UNIT_StdOut,"(A" // adjustl(fmtName) // ")",advance='no')  &
                                                         this%name(:maxNameLen)
  end if

  ! print delimiter between name and value
  select case(mode)
  case(0)
    write(UNIT_stdOut,'(a3)', advance='no')  " | "
  case(1)
    write(UNIT_StdOut,"(A3)",advance='no') " = "
  case(2)
    write(UNIT_StdOut,"(A3)",advance='no') "   "
  end select

  ! print value
  if ((mode.eq.0).or.(this%hasDefault)) then
    call this%printValue(maxValueLen)
  else
    write(UNIT_StdOut, "(A"//fmtValue//")", advance='no') ""
  end if


  ! print default/custom or print comment
  if (mode.eq.0) then
    ! print default/custom
    if (this%isSet) then
      write(UNIT_stdOut,"(a3)", advance='no') ' | '
      !call set_formatting("green")
      write(UNIT_stdOut,'(a7)', advance='no')  "*custom"
      !call clear_formatting()
      write(UNIT_stdOut,"(a3)") ' | '
    else
      write(UNIT_stdOut,"(a3)", advance='no') ' | '
      !call set_formatting("red")
      write(UNIT_stdOut,'(a7)', advance='no')  "default"
      !call clear_formatting()
      write(UNIT_stdOut,"(a3)") ' | '
    end if
  else
    ! print comment: this is complicated, since it includes line breaks for
    ! long comments line breaks are inserted at spaces and at '\n' characters
    comment = trim(this%description)
    write(fmtValue,*) maxNameLen + maxValueLen + 4
    length = 0
    write (UNIT_StdOut,'(A)',advance='no') " " ! insert space after value
    ! loop until comment is empty (split by newline)
    do while (LEN_TRIM(comment) .gt. 0)
      ! split comment at first newline. After the split:
      !   - comment contains remaining part after newline
      !   - headNewline contains part before newline or ==comment, !
      ! if there is no newline character in the comment anymore
      call split(comment, headNewline, "\n")
      ! loop until headNewline is empty (split by words)
      do while(LEN_TRIM(headNewline) .gt. 0)
        ! split comment at first space. After the split:
        !   - headNewline contains remaining part after the space
        !   - headSpace contains part before space or ==headNewline,
        ! if there is no space character in the headNewline anymore
        call split(headNewline, headSpace, " ")
        ! if word in headSpace does not fit into actual line -> insert newline
        if (length+LEN_TRIM(headSpace).gt.50) then
          write (UNIT_StdOut,*) ''
          write(UNIT_StdOut, "(A"//fmtValue//")", advance='no') ""
          length = 0 ! reset length of line
        end if
        ! insert word in headSpace and increase length of actual line
        if ((length.eq.0).and.(mode.eq.1)) then
          write(UNIT_StdOut,'(A2)',advance='no') "! "
        end if
        write (UNIT_StdOut,'(A)',advance='no') char(headSpace)//" "
        length = length + LEN_TRIM(headSpace)+1
      end do
      ! insert linebreak due to newline character in comment
      write(UNIT_StdOut,*) ''
      write(UNIT_StdOut, "(A"//fmtValue//")", advance='no') ""
      length = 0
    end do
    ! insert empty line after each option
    write(UNIT_StdOut,*) ''
  end if
end subroutine print

!----------------------------------------------------------------------

!> print value of an option
subroutine printValue(this,maxValueLen)
 !< option to print
 class(option),intent(in)    :: this
 !< max string length of name
 integer,intent(in)          :: maxValueLen

 character(len=20)            :: fmtValue
 character(len=20)            :: fmtDigits
 character(len=maxValueLen)   :: intFromStringOutput
 integer                      :: i,length

  write(fmtValue,*) maxValueLen
  select type (this)
  class is (IntOption)
    write(UNIT_StdOut,"(I"//fmtValue//")",advance='no') this%value
  class is (LogicalOption)
    write(UNIT_StdOut,"(L"//fmtValue//")",advance='no') this%value
  class is (RealOption)
    if (this%digits.ge.1) then ! floating point representation
      write(fmtDigits,*) this%digits
      write(UNIT_stdOut,'(F'//fmtValue//'.'//fmtDigits//')',advance='no') &
                                                                     this%value
    else if (this%digits.le.-1) then ! scientific (exponential) representation
      write(fmtDigits,*) -this%digits
      write(UNIT_stdOut,'(E'//fmtValue//'.'//fmtDigits//')',advance='no') &
                                                                     this%value
    else ! digits not given
      write(UNIT_stdOut,'(E'//fmtValue//'.19)',advance='no') this%value
    end if
  class is (StringOption)
    if (trim(this%value).eq."") then
      write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') '-'
    else
      write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') trim(this%value)
    end if
  class is (IntFromStringOption)
    if (trim(this%value).eq."") then
      write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') '-'
    else
      ! IntFromStringOption: Print in the format string (integer) if the
      ! given value is found in the mapping, otherwise print just the integer
      if (this%foundInList) then
        write(intFromStringOutput,"(A,A,I0,A)") &
                     trim(this%strList(this%listIndex)), &
                                         ' (', this%intList(this%listIndex), ')'
        write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') &
                                                      trim(intFromStringOutput)
      else
        write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') trim(this%value)
      end if
    end if
  class is (IntArrayOption)
    length=this%getvaluelen()
    if (maxValueLen - length.gt.0) then
      write(fmtValue,*) (maxValueLen - length)
      write(UNIT_stdOut,'('//fmtValue//'(" "))',advance='no')
    end if
    write(UNIT_StdOut,"(A3)",advance='no') "(/ "
    do i=1,size(this%value)
      write(fmtValue,'(I0)') this%value(i)
      write(fmtValue,*) LEN_TRIM(fmtValue)
      write(UNIT_StdOut,"(I"//fmtValue//")",advance='no') this%value(i)
      if (i.ne.size(this%value)) then
        write(UNIT_StdOut,"(A2)",advance='no') ", "
      end if
    end do
    write(UNIT_StdOut,"(A3)",advance='no') " /)"
  class is (LogicalArrayOption)
    length=this%getvaluelen()
    if (maxValueLen - length.gt.0) then
      write(fmtValue,*) (maxValueLen - length)
      write(UNIT_stdOut,'('//fmtValue//'(" "))',advance='no')
    end if
    write(UNIT_StdOut,"(A3)",advance='no') "(/ "
    do i=1,size(this%value)
      write(UNIT_StdOut,"(L1)",advance='no') this%value(i)
      if (i.ne.size(this%value)) then
        write(UNIT_StdOut,"(A2)",advance='no') ", "
      end if
    end do
    write(UNIT_StdOut,"(A3)",advance='no') " /)"
  class is (RealArrayOption)
    length=this%getvaluelen()
    if (maxValueLen - length.gt.0) then
      write(fmtValue,*) (maxValueLen - length)
      write(UNIT_stdOut,'('//fmtValue//'(" "))',advance='no')
    end if
    write(UNIT_stdOut,'(a3)',advance='no') '(/ '
    do i=1,size(this%value)
      write(fmtValue,*) getstrlenreal(this%value(i), this%digits(i))
        if (this%digits(i).ge.1) then ! floating point representation
        write(fmtDigits,*) this%digits(i)
        write(UNIT_stdOut,'(F'//fmtValue//'.'//fmtDigits//',A3)',advance='no')&
                                                                  this%value(i)
      ! scientific (exponential) representation
      else if (this%digits(i).le.-1) then
        write(fmtDigits,*) -this%digits(i)
        write(UNIT_stdOut,'(E'//fmtValue//'.'//fmtDigits//')',advance='no') &
                                                                  this%value(i)
      else ! digits not given
        write(UNIT_stdOut,'(E'//fmtValue//'.19,A3)',advance='no') this%value(i)
      end if
      if (i.ne.size(this%value)) then
        write(UNIT_stdOut,'(A2)',advance='no') ', '
      end if
    end do
    write(UNIT_stdOut,'(a3)',advance='no') ' /)'
  !###
  ! todo: Causes internal compiler error with gnu 6+ due to compiler bug
  ! (older gnu and Intel,Cray work). Uncomment as unused.
  !###
  !class is (StringArrayOption)
    !length=this%getvaluelen()
    !if (maxValueLen - length.gt.0) then
      !write(fmtValue,*) (maxValueLen - length)
      !write(UNIT_stdOut,'('//fmtValue//'(" "))',advance='no')
    !end if
    !write(UNIT_StdOut,"(A3)",advance='no') "(/ "
    !do i=1,size(this%value)
      !write(fmtValue,*) LEN_TRIM(this%value(i))
      !write(UNIT_StdOut,"(A"//fmtValue//")",advance='no') this%value(i)
      !if (i.ne.size(this%value)) then
        !write(UNIT_StdOut,"(A2)",advance='no') ", "
      !end if
    !end do
    !write(UNIT_StdOut,"(A3)",advance='no') " /)"
  class is(SubOption)
  !TODO: do something
  class default
    stop
  end select
end subroutine printValue

!----------------------------------------------------------------------

!> parse value from string 'rest_in'.
!!
!!This subroutine is used to readin values from the parameter file.
subroutine parse(this, rest_in, is_default)
 class(option)               :: this     !< class(option)
 character(len=*),intent(in) :: rest_in  !< string to parse
 logical, optional, intent(in) :: is_default
 character(len=max_char_len)  :: tmp,tmp2,rest
 integer             :: count,i,stat
 integer,allocatable :: inttmp(:)
 logical,allocatable :: logtmp(:)
 real(wp),allocatable    :: realtmp(:)
 logical :: def
 character(len=*), parameter :: this_sub_name = 'parse'
 !character(len=max_char_len),allocatable :: strtmp(:)

  def = .false.
  if(present(is_default)) then
    if(is_default) def = .true.
  endif


  stat=0
  ! Replace brackets
  rest=Replace(rest_in,"(/"," ",Every=.true.)
  rest=Replace(rest,   "/)"," ",Every=.true.)
  rest = trim(rest)

  select type (this)
  class is (IntOption)
    if(.not.def) then
      read(rest, *,iostat=stat) this%value
    else
      read(rest, *,iostat=stat) this%default_value
    endif
  class is (LogicalOption)
    if(.not.def) then
      read(rest, *,iostat=stat) this%value
    else
      read(rest, *,iostat=stat) this%default_value
    endif
  class is (RealOption)
    if(.not.def) then
      call this%parseReal(rest, this%value, this%digits)
    else
      call this%parseReal(rest, this%default_value, this%digits)
    endif
  class is (StringOption)
    if(.not.def) then
      read(rest, "(A)",iostat=stat) this%value
    else
      read(rest, "(A)",iostat=stat) this%default_value
    endif
  class is (IntFromStringOption)
    if(.not.def) then
      read(rest, "(A)",iostat=stat) this%value
    else
      read(rest, "(A)",iostat=stat) this%default_value
    endif
  class is (IntArrayOption)
    ! Array options are complicated, since we do not know a priori how long
    ! the array will be.
    ! Therefore we must read entry by entry and always increase the array
    ! size by reallocation.
    ! This is done using the temporary arrays inttmp/realtmp and with the
    ! MOVE_ALLOC routine.
    tmp2 = trim(adjustl(rest))
    if (allocated(this%value)) deallocate(this%value)
    count = 0
    allocate(this%value(count))
    do ! loop until no more entry
      i = index(trim(tmp2), ',')
      ! store text of tmp2 until next , in tmp
      if (i.gt.0) then
        tmp = tmp2(1:i-1)
      else
        tmp = tmp2
      end if
      ! remove entry and trim the result
      tmp2 = tmp2(i+1:)
      tmp2 = trim(adjustl(tmp2))
      ! increase the number of entries
      count = count + 1
      ! allocate temporary array and copy content of this%value to the first
      ! entries aftwards use MOVE_ALLOC to move it back to this%value
      ! (this deallocates the temporary array)
      allocate(inttmp(count))
      inttmp(:size(this%value)) = this%value
      call MOVE_ALLOC(inttmp, this%value) ! inttmp gets deallocated
      ! finally this%value has the correct size and we can readin the entry
      read(tmp, *,iostat=stat) this%value(count)
      if (i.eq.0) exit
    end do
    if(def) call move_alloc(this%value, this%default_value)
  class is (LogicalArrayOption)
    ! comments see IntArrayOption!
    tmp2 = trim(adjustl(rest))
    if (allocated(this%value)) deallocate(this%value)
    count = 0
    allocate(this%value(count))
    do ! loop until no more entry
      i = index(trim(tmp2), ',')
      ! store text of tmp2 until next , in tmp
      if (i.gt.0) then
        tmp = tmp2(1:i-1)
      else
        tmp = tmp2
      end if
      ! remove entry and trim the result
      tmp2 = tmp2(i+1:)
      tmp2 = trim(adjustl(tmp2))
      ! increase the number of entries
      count = count + 1
      ! allocate temporary array and copy content of this%value to the first
      ! entries aftwards use MOVE_ALLOC to move it back to this%value
      ! (this deallocates the temporary array)
      allocate(logtmp(count))
      logtmp(:size(this%value)) = this%value
      call MOVE_ALLOC(logtmp, this%value) ! logtmp gets deallocated
      ! finally this%value has the correct size and we can readin the entry
      read(tmp, *,iostat=stat) this%value(count)
      if (i.eq.0) exit
    end do
    if(def) call move_alloc(this%value, this%default_value)
  class is (RealArrayOption)
    ! comments see IntArrayOption!
    tmp2 = trim(adjustl(rest))
    if (allocated(this%value)) deallocate(this%value)
    if (allocated(this%digits)) deallocate(this%digits)
    count = 0
    count = 0
    allocate(this%value(count))
    allocate(this%digits(count))
    do
      i = index(trim(tmp2), ',')
      if (i.gt.0) then
        tmp = tmp2(1:i-1)
      else
        tmp = tmp2
      end if
      tmp2 = tmp2(i+1:)
      tmp2 = trim(adjustl(tmp2))
      count = count + 1
      allocate(realtmp(count))
      allocate(inttmp(count))
      realtmp(:size(this%value)) = this%value
      inttmp(:size(this%digits)) = this%digits
      call MOVE_ALLOC(realtmp, this%value)  ! realtmp gets deallocated
      call MOVE_ALLOC(inttmp,  this%digits) ! realtmp gets deallocated
      this%digits(count) = 0

      call this%parseReal(tmp, this%value(count), this%digits(count))

      if (i.eq.0) exit
    end do
    if(def) call move_alloc(this%value, this%default_value)
  !class is (StringArrayOption)
    !! comments see IntArrayOption!
    !tmp2 = trim(adjustl(rest))
    !if (allocated(this%value)) deallocate(this%value)
    !count = 0
    !allocate(this%value(count))
    !do ! loop until no more entry
      !i = index(trim(tmp2), ',')
      !! store text of tmp2 until next , in tmp
      !if (i.gt.0) then
        !tmp = tmp2(1:i-1)
      !else
        !tmp = tmp2
      !end if
      !! remove entry and trim the result
      !tmp2 = tmp2(i+1:)
      !tmp2 = trim(adjustl(tmp2))
      !! increase the number of entries
      !count = count + 1
      !! allocate temporary array and copy content of this%value to the first
      !! entries aftwards use MOVE_ALLOC to move it back to this%value
      !! (this deallocates the temporary array)
      !allocate(strtmp(count))
      !strtmp(:size(this%value)) = this%value
      !call MOVE_ALLOC(strtmp, this%value) ! strtmp gets deallocated
      !! finally this%value has the correct size and we can readin the entry
      !read(tmp, *,iostat=stat) this%value(count)
      !if (i.eq.0) exit
    !end do

  class is (SubOption)
    read(rest, "(A)",iostat=stat) this%value
    if (trim(this%value) .ne. '{') &
    call error(this_sub_name, this_mod_name, 'Error while parsing the option "'&
               //trim(rest_in)//' It was meant to be a su-option and sould have&
               & only "{" on its line to indicate the opening of a sub option &
               & block')
  class default
    call error(this_sub_name, this_mod_name, 'Unknown option class, this should&
    &have never happened, a team of professionals is on the way to eliminate &
    &the evidence')
  end select
  if(stat.gt.0)then
    call error(this_sub_name, this_mod_name, 'Error while parsing the option "'&
               //trim(this%name)//' " with value: '//trim(rest_in))
  end if

end subroutine parse

!----------------------------------------------------------------------

!> parse string to real and get the format of the number (floating,scientific)
subroutine parseReal(this,string_in, value, digits)
 class(option)                 :: this
 !< (in) string containing a real(wp) number
 character(len=max_char_len),intent(in) :: string_in
 !< (out) the converted real
 real(wp),intent(out)              :: value
 !< (out) the number of digits if floating representation, or -1 if scientific
 integer,intent(out)           :: digits

 integer            :: pos,posE,posMinus,stat
 character(len=max_char_len) :: string  ! left adjusted string
 character(len=*), parameter :: this_sub_name='parseReal'

  string = adjustl(string_in)
  posE=max(index(string, 'E'),index(string, 'e'))
  if (posE.gt.0) then
    ! exponential representation, will be printed as '0.xxEyy'.
    ! The number of digits (here: the length of 'xx') is needed.
    ! This number is the number of significant digits in the string,
    ! before and after the dot.
    pos = index(string, '.')
    if (pos.eq.0) then
      ! no dot in string, number of digits is length of string before
      ! the 'e' == posE-1
      digits = -(posE-1)
    else
      ! dot in string, number of digits is length of string before the 'e' -1
      ! for the dot == posE-2
      digits = -(posE-2)
    end if
    posMinus = index(string,'-')
    if (posMinus.eq.1) then
      ! string starts with a minus => decrease number of digits by one
      ! (here +1, since negative values for scientific format)
      digits = digits+1
      ! string is -0.xxx, meaning no significant stuff before the dot
      ! => remove the zero form the number of digits
      if (index(string,'0.').eq.2) digits = digits+1
    else
      ! string is 0.xxx, meaning no significant stuff before the dot
      ! => remove the zero form the number of digits
      if (index(string,'0.').eq.1) digits = digits+1
    end if

  else
    ! floating representation
    pos = index(string, '.')
    if (pos.eq.0) then
      digits = 1
    else
      digits = LEN_TRIM(string) - pos
      if (digits.eq.0) digits = 1
    end if
  end if
  read(string, *,iostat=stat) value
  if(stat.gt.0)then
    call error(this_sub_name, this_mod_name, 'Error while parsing the string "'&
               //string_in//'" to a real')
  end if

end subroutine parseReal

end module
