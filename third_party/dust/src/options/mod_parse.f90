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
!! It was originally part of Flexi, and has been extensively modified to suit
!! the needs of DUST
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

!> Module providing routines for reading Flexi parameter files.
!!
!! The whole structure to read options from the parameter file is as follows:
!!
!! All the options are stored in a linked list, which is defined as a class
!! which is exposed in this module and must be instantiated where employed.
!! In this way more intancies of parameters can co-exist
!!
!! The options are appended to this list via the type-bound procedures of the
!! t_parse class (CreateStringOption, CreateLogicalOption etc...),
!! which must be called before execution when needed. After calling all
!! the option creating procedures t_parse list contains all possible options
!! (with name, description, default value (optional)).
!!
!! After that the t_param\%read_options() routine is called,
!! which actually reads the options from the parameter file. Therefore the
!! parameter file is read line by line and each line is parsed for an option.
!! By this the values of the options, that are already in the linked list of
!! type 't_parse' are set.
!!
!! Now all the options are filled with the data from the parameter file and
!! can be accessed via the functions getint(array), getreal(array), ...
!! NOTE: this behaviour is likely to change in the future, in favour of
!! type-bound procedures also in this case
!!
!! A call of these functions then removes the specific option from the linked
!! list, such that every option can only be read once. This is necessary for
!! options with the same name, that occure multiple times in the parameter
!! file.
!!


!WARNING: while passing from module variable to external variable the class
!t_prms was added to dummy arguments in most of the functions, with intent(in),
!unless the compiler complaint. This anyway can cause unexpected behaviours...
module mod_parse

!TODO review all these

!use MOD_ISO_VARYING_STRING
use MOD_ISO_VARYING_STRING, only: &
  varying_string, &
  !assignment(=), &
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
  !len_trim     , &
  !lge          , &
  !lgt          , &
  !lle          , &
  !llt          , &
  !repeat       , &
  !scan         , &
  !trim         , &
  !verify       , &
  !var_str      , &
  get          , &
  !put          , &
  !put_line     , &
  !extract      , &
  !insert       , &
  !remove       , &
  replace      , &
  split

!use MOD_Options
use MOD_Options, only: &
    option, IntOption, IntArrayOption, IntFromStringOption, LogicalOption, &
    LogicalArrayOption, RealOption, RealArrayOption, &
    StringOption, SubOption

use MOD_StringTools, only: &
  LowCase, stricmp!, use_escape_codes!, set_formatting, clear_formatting

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, warning, info, unit_stdout, printout


!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------

public :: t_parse, ignoredParameters, finalizeparameters, &
          getstr, getlogical, getreal, getrealarray, getint, getintarray, &
          countoption, t_link, check_opt_consistency, getsuboption, &
          print_parse_debug

private

!-----------------------------------------------------------------------
!> Link for linked List
type :: t_link
  class(option), pointer :: opt => null()
  type(t_link), pointer   :: next => null()
  type(t_link), pointer   :: prev_read => null()
  type(t_link), pointer   :: next_read => null()
  type(t_parse), pointer   :: sub_list => null()
end type t_link

!> Class to store all options.
!! This is basically a linked list of options.
type :: t_parse
  type(t_link), pointer :: firstLink => null() !< first option in the list
  type(t_link), pointer :: lastLink  => null() !< last option in the list
  integer              :: maxNameLen          !< maximal string length of the name of an option in the list
  integer              :: maxValueLen         !< maximal string length of the value of an option in the list
  character(len=max_char_len)   :: actualSection = ""  !< actual section, to set section of an option, when inserted into list
  logical              :: printout_val !< print the value when collected in the code or not
  class(t_parse), pointer :: parent_parse => null()
contains
  procedure :: SetSection                 !< routine to set 'actualSection'
  procedure :: CreateOption               !< general routine to create a option and insert it into the linked list
  procedure :: CreateIntOption            !< routine to generate an integer option
  procedure :: CreateIntFromStringOption  !< routine to generate an integer option with a optional string representation
  procedure :: CreateLogicalOption        !< routine to generate an logical option
  procedure :: CreateRealOption           !< routine to generate an real option
  procedure :: CreateStringOption         !< routine to generate an string option
  procedure :: CreateIntArrayOption       !< routine to generate an integer array option
  procedure :: CreateLogicalArrayOption   !< routine to generate an logical array option
  procedure :: CreateRealArrayOption      !< routine to generate an real array option
  procedure :: CreateSubOption            !< routine to generate a sub-option block
  !procedure :: CreateStringArrayOption    !< routine to generate an string array option
  procedure :: CountOption_               !< function to count the number of options of a given name
  procedure :: read_options               !< routine that loops over the lines of a parameter files
                                          !< and calls read_option for every option. Outputs all unknow options
  procedure :: read_option                !< routine that parses a single line from the parameter file.
end type t_parse

interface IgnoredParameters
  module procedure IgnoredParameters
end interface

interface PrintDefaultParameterFile
  module procedure PrintDefaultParameterFile
end interface

interface CountOption
  module procedure CountOption
end interface

interface getint
  module procedure getint
end interface

interface getlogical
  module procedure getlogical
end interface

interface getreal
  module procedure getreal
end interface

interface getstr
  module procedure getstr
end interface

interface getintarray
  module procedure getintarray
end interface

interface getlogicalarray
  module procedure getlogicalarray
end interface

interface getrealarray
  module procedure getrealarray
end interface

interface getstrarray
  module procedure getstrarray
end interface

interface getintfromstr
  module procedure getintfromstr
end interface

interface getdescription
  module procedure getdescription
end interface

interface addStrListEntry
  module procedure addStrListEntry
end interface

interface FinalizeParameters
  module procedure FinalizeParameters
end interface

!public :: IgnoredParameters
!public :: PrintDefaultParameterFile
!public :: CountOption
!public :: getint
!public :: getlogical
!public :: getreal
!public :: getstr
!public :: getintarray
!public :: getlogicalarray
!public :: getrealarray
!public :: getstrarray
!public :: getdescription
!public :: getintfromstr
!public :: addStrListEntry
!public :: FinalizeParameters

!type(t_parse) :: prms
!public :: prms

  type :: t_str
     private
     character(len=max_char_len) :: chars
  end type t_str

!module variables

character(len=*), parameter :: this_mod_name = 'mod_parse'

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

!> Set actual section.
!!
!! All options created after calling this subroutine are
!! in this 'section'. The section property is only used to get nicer looking
!! parameter files when using --help or --markdown.
subroutine SetSection(this, section)
 class(t_parse),intent(inout) :: this    !< class(t_parse)
 character(len=*),intent(in)     :: section !< section to set

  this%actualSection = section
end subroutine SetSection

!-----------------------------------------------------------------------

!> General routine to create an option.
!!
!! Fills all fields of the option. Since the prms\%parse function is used to
!! set the value, this routine can be abstract for all  types of options.
subroutine CreateOption(this, opt, name, description, value, multiple)
 class(t_parse),intent(inout)       :: this          !< class(t_parse)
 class(option),intent(inout)           :: opt         !< option class
 character(len=*),intent(in)           :: name        !< option name
 character(len=*),intent(in)           :: description !< option description
 character(len=*),intent(in),optional  :: value       !< option value
 logical,intent(in),optional           :: multiple    !< marker if multiple option

 type(t_link), pointer :: newLink
 character(len=*), parameter :: this_sub_name = 'CreateOption'

  opt%hasDefault = present(value)
  if (opt%hasDefault) then
    call opt%parse(value,is_default=.true.)
  end if

  opt%multiple   = .false.
  opt%additional = .false.
  if (present(multiple)) opt%multiple = multiple
  if (opt%multiple.and.opt%hasDefault) then
    call error(this_sub_name, this_mod_name, &
        "A default value can not be given, when multiple=.true. in creation&
        &of option: '"//trim(name)//"'")
  end if

  opt%isSet = .false.
  opt%name = name
  opt%description = description
  opt%section = this%actualSection
  opt%isRemoved = .false.

  ! insert option into linked list
  if (.not. associated(this%firstLink)) then
    this%firstLink => constructor_Link(opt, this%firstLink)
    this%lastLink => this%firstLink
  else
    newLink => constructor_Link(opt, this%lastLink%next)
    this%lastLink%next => newLink
    this%lastLink => newLink
  end if
end subroutine CreateOption

!-----------------------------------------------------------------------

!> Create a sub-option, which is just a container for other options
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateSubOption(this, name, description, subparse, value, multiple)
 class(t_parse),intent(inout)      :: this           !< class(t_parse)
 character(len=*),intent(in)          :: name         !< option name
 character(len=*),intent(in)          :: description  !< option description
 type(t_parse), intent(out),pointer   :: subparse     !< Sublink for the suboptions
 character(len=*),intent(in),optional :: value        !< option value
 logical,intent(in),optional          :: multiple     !< marker if multiple option

 class(SubOption),allocatable,target :: subopt

  allocate(subopt)
  call this%CreateOption(subopt, name, description, value=value, &
                         multiple=multiple)
  !Create the sub-option parser and link it
  allocate(subparse)
  this%lastLink%sub_list => subparse

end subroutine CreateSubOption

!-----------------------------------------------------------------------

!> Create a new integer option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateIntOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this           !< class(t_parse)
 character(len=*),intent(in)          :: name         !< option name
 character(len=*),intent(in)          :: description  !< option description
 character(len=*),intent(in),optional :: value        !< option value
 logical,intent(in),optional          :: multiple     !< marker if multiple option

 class(IntOption),allocatable,target :: intopt

  allocate(intopt)
  call this%CreateOption(intopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateIntOption

!-----------------------------------------------------------------------

!> Create a new integer option with a optional string representation.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateIntFromStringOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(IntFromStringOption),allocatable,target :: intfromstropt

  allocate(intfromstropt)
  call this%CreateOption(intfromstropt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateIntFromStringOption

!-----------------------------------------------------------------------

!> Create a new logical option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateLogicalOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)         :: this        !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(LogicalOption),allocatable,target :: logicalopt

  allocate(logicalopt)
  call this%CreateOption(logicalopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateLogicalOption

!-----------------------------------------------------------------------

!> Create a new real option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateRealOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(RealOption),allocatable,target :: realopt

  allocate(realopt)
  call this%CreateOption(realopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateRealOption

!-----------------------------------------------------------------------

!> Create a new string option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateStringOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(StringOption),allocatable,target :: stringopt

  allocate(stringopt)
  call this%CreateOption(stringopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateStringOption

!-----------------------------------------------------------------------

!> Create a new integer array option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateIntArrayOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(IntArrayOption),allocatable,target :: intopt

  allocate(intopt)
  call this%CreateOption(intopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateIntArrayOption

!-----------------------------------------------------------------------

!> Create a new logical array option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateLogicalArrayOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(LogicalArrayOption),allocatable,target :: logicalopt

  allocate(logicalopt)
  call this%CreateOption(logicalopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateLogicalArrayOption

!-----------------------------------------------------------------------

!> Create a new real array option.
!!
!! Only calls the general prms\%createoption routine.
subroutine CreateRealArrayOption(this, name, description, value, multiple)
 class(t_parse),intent(inout)      :: this          !< class(t_parse)
 character(len=*),intent(in)          :: name        !< option name
 character(len=*),intent(in)          :: description !< option description
 character(len=*),intent(in),optional :: value       !< option value
 logical,intent(in),optional          :: multiple    !< marker if multiple option

 class(RealArrayOption),allocatable,target :: realopt

  allocate(realopt)
  call this%CreateOption(realopt, name, description, value=value, &
                         multiple=multiple)
end subroutine CreateRealArrayOption

!-----------------------------------------------------------------------

!> Create a new string array option. Only calls the general prms\%createoption routine.
!subroutine CreateStringArrayOption(this, name, description, value, multiple)
!class(t_parse),intent(inout)      :: this           !< class(t_parse)
!character(len=*),intent(in)          :: name           !< option name
!character(len=*),intent(in)          :: description    !< option description
!character(len=*),intent(in),optional :: value          !< option value
!logical,intent(in),optional          :: multiple       !< marker if multiple option
!class(StringArrayOption),allocatable,target :: stringopt
!allocate(stringopt)
!call this%CreateOption(stringopt, name, description, value=value, multiple=multiple)
!end subroutine CreateStringArrayOption

!-----------------------------------------------------------------------

!> Count number of occurrence of option with given name.
function CountOption_(this, name) result(count)
  class(t_parse),intent(inout)      :: this    !< class(t_parse)
  character(len=*),intent(in)       :: name    !< Search for this keyword in ini file
  integer                           :: count   !< number of found occurences of keyword

  type(t_link),pointer :: current
  count = 0
  ! iterate over all options and compare names
  current => this%firstLink
  
  do while (associated(current))
    if (current%opt%nameequals(name)) then
      if (current%opt%isSet) count = count + 1
    end if
    current => current%next
  end do
end function  CountOption_

!-----------------------------------------------------------------------

!> Insert a option in front of option with same name in the 'prms' linked list.
subroutine insertOption(first, opt)
 type(t_link),pointer,intent(inout) :: first !< first item in linked list
 class(option),intent(in)       :: opt   !< option to be inserted

 type(t_link),pointer :: current
 type(t_link),pointer :: newLink

  current =>  first
  do while (associated(current%next))
    if (.not.current%next%opt%nameequals(opt%name)) then
      exit
    end if
    current => current%next
  end do
  newLink => constructor_Link(opt, current%next)
  current%next => newLink
  !last thing: give back the last position of the list
  first => current
end subroutine insertOption

!-----------------------------------------------------------------------

recursive subroutine copy_parser(source_p, target_p)
 type(t_parse), pointer, intent(in)   :: source_p
 type(t_parse), pointer, intent(out)  :: target_p
 !type(t_parse), allocatable, intent(out)  :: target_p

 type(t_link), pointer :: newLink
 class(option), pointer  :: newopt
 type(t_link), pointer :: current_source, current_target

  !allocate(newParser)
  !target_p => newParser
  allocate(target_p, source=source_p)
  !allocate(target_p)
  target_p%firstLink => null(); target_p%lastLink => null();
  current_source => source_p%firstLink
  current_target => target_p%firstLink

  do while (associated(current_source))
    if (.not.current_source%opt%additional) then
      allocate(newopt, source=current_source%opt)

      newopt%isSet =.false.

      if (.not. associated(target_p%firstLink)) then
        target_p%firstLink => constructor_Link(newopt, target_p%firstLink)
        target_p%lastLink => target_p%firstLink
      else
        newLink => constructor_Link(newopt, target_p%lastLink%next)
        target_p%lastLink%next => newLink
        target_p%lastLink => newLink
      end if
      current_target => target_p%lastLink


      if(associated(current_source%sub_list)) then
        call copy_parser(current_source%sub_list, current_target%sub_list)
        !call copy_parser(current_source%sub_list, newParser)
        !current_target%sub_list => newParser
      endif
    endif

    newopt => null()
    current_source => current_source%next
  enddo

end subroutine copy_parser

!-----------------------------------------------------------------------

!> Read options from parameter file.
!!
!! Therefore the file is read line by line. After removing comments and all
!! white spaces each line is parsed in the prms\%read_option() routine.
!! Outputs all unknown options.
subroutine read_options(this, filename, printout_val)
 class(t_parse),intent(inout), target :: this     !< class(t_parse)
 character(len=*),intent(in)   :: filename !< name of file to be read
 logical, intent(in), optional :: printout_val !< set the print behaviour of the
                                           !! class, default true

 type(t_link), pointer  :: current
 type(t_link), pointer  :: prev_read
 class(t_parse), pointer :: actually_reading, root_reading
 class(t_parse), pointer :: sub_option
 integer               :: stat,iniUnit!,iExt
 type(Varying_String)  :: aStr,bStr
 character(len=max_char_len)    :: HelpStr
 logical               :: firstWarn=.true.,file_exists
 character(len=*), parameter :: this_sub_name = 'read_options'


  this%printout_val = .true.
  if(present(printout_val)) this%printout_val = printout_val


  inquire(file=trim(filename), exist=file_exists)
  if (.not.file_exists) then
    call error(this_sub_name, this_mod_name, &
         'The file '//trim(filename)//' does not exist')
  end if


  !call info(this_sub_name, this_mod_name, &
  !                  'Reading parameters from file "'//trim(filename)//'":')
  if (this%printout_val) then 
    call printout(nl//'Reading input parameters from file "'//&
                trim(filename)//'"')
  endif 
  
  ! Open parameter file for reading
  iniUnit= 100 !getfreeunit()
  open(unit=iniUnit,file=trim(filename),status='old',action='read', &
       access='sequential',iostat=stat)
  if(stat.ne.0) then
    call error(this_sub_name, this_mod_name, &
               'Not possible to read from file '//trim(filename))
  end if

  actually_reading => this
  root_reading => this
  this%parent_parse => this

  prev_read => null()
  ! infinte loop. Exit at eof
  do
    sub_option => null()
    ! read a line into 'aStr'
    call Get(iniUnit,aStr,iostat=stat)
    ! exit loop if eof
    if(IS_IOSTAT_END(stat)) exit
    if(.not.IS_IOSTAT_EOR(stat)) then
      call error(this_sub_name, this_mod_name, &
                 'Error during parameter file parsing')
    end if
    ! Remove comments with "!"
    call Split(aStr,bStr,"!")
    ! Remove comments with "#"
    call Split(bStr,aStr,"#")
    ! aStr may hold an option

    ! Remove blanks
    aStr=Replace(aStr," ","",Every=.true.)
    ! Replace brackets
    aStr=Replace(aStr,"(/"," ",Every=.true.)
    aStr=Replace(aStr,"/)"," ",Every=.true.)
    
    HelpStr = char(aStr)        

    ! If something remaind, this should be an option
    if (LEN_TRIM(HelpStr).gt.2) then
      ! read the option
      if (.not. actually_reading%read_option(HelpStr, sub_option, prev_read)) then
        if (firstWarn) then
          firstWarn=.false.
          !write(UNIT_StdOut, *) "WARNING: The following options in file "&
          !                              //trim(filename)//" are unknown!"
          !write(UNIT_StdOut,*) "   ", trim(HelpStr)                              
        end if
        ! write(UNIT_StdOut,*) "   ", trim(HelpStr)
      end if

      if (associated(sub_option)) then
        !we are getting inside a sub-option!
        root_reading => actually_reading
        actually_reading =>  sub_option
        sub_option%parent_parse => root_reading

      endif
    end if

    if(scan(HelpStr,'}').gt.0) then
      !we are getting out of a sub_option!
      actually_reading => root_reading
      root_reading => actually_reading%parent_parse

    endif

  end do
  if (.not.firstWarn) then
    !=!write(UNIT_StdOut,'(100("!"))')
  end if
  close(iniUnit)

  ! calculate the maximal string lenght of all option-names and option-values
  this%maxNameLen  = 0
  this%maxValueLen = 0
  current => this%firstLink
  do while (associated(current))
    this%maxNameLen = max(this%maxNameLen, current%opt%getnamelen())
    this%maxValueLen = max(this%maxValueLen, current%opt%getvaluelen())
    current => current%next
  end do

  ! check for colored output
  ! MATTEO: hardcode disabled for the moment
  ! use_escape_codes = getlogical(this, "ColoredOutput")
  !use_escape_codes = .false.
end subroutine read_options

!-----------------------------------------------------------------------

!> Parses one line of parameter file and sets the value of the specific option+
!! in the 'prms' linked list.
!!
!! Therefore it iterate over all entries of the linked list and compares the names.
function read_option(this, line, sub_option, prev_read) result(found)
 class(t_parse),intent(in)    :: this  !< class(t_parse)
 character(len=*),intent(in)  :: line  !< line to be parsed
 class(t_parse), pointer, intent(inout)  :: sub_option
 type(t_link), pointer, intent(inout) :: prev_read
 logical                      :: found !< marker if option found

 character(len=max_char_len)           :: name
 character(len=max_char_len)           :: rest
 type(t_link), pointer         :: current
 class(option),allocatable    :: newopt
 !type(t_parse),allocatable,target    :: newparse
 integer                      :: i
 character(len=*), parameter :: this_sub_name = 'read_option'

  found = .false.

  ! split at '='
  i = index(line, '=')
  if (i==0) return
  name = line(1:i-1)
  rest = line(i+1:)

  ! compatibility with legacy naming convention
  if(name .eq. "CompName" .or. name .eq. "GeometryFile" .or. name .eq. "MeshFileType" &
    & .or. name .eq. "Reference_Tag" .or. name .eq. "StartRes") then
    call printout(nl//'WARNING: the input file is using legacy syntax, which will be &
            & deprecated in the next release. Please update this file.'//nl)
  end if

  !default, we have not found a sub-option
  sub_option => null()

  current => this%firstLink
  ! iterate over all options and compare names
  do while (associated(current))
    ! compare name
    if (current%opt%nameequals(name)) then
      found = .true.
      if (current%opt%isSet) then
        if (.not.(current%opt%multiple)) then
          ! option already set, but is not a multiple option
          write(UNIT_StdOut,*) 'Option "', trim(name), &
                             '" is already set, but is not a multiple option!'
          stop
        else
          ! create new instance of multiple option
          allocate(newopt, source=current%opt)
          call newopt%parse(rest)
          newopt%isSet = .true.
          newopt%additional = .true.
          ! insert option
          call insertOption(current, newopt)
          !if it is a sub-option it is necessary to re-build all the chain
          select type(op=>current%opt)
          type is(SubOption)
            call copy_parser(current%sub_list, current%next%sub_list)
          end select
          !highly experimental
          current => current%next
          current%prev_read => prev_read
          if(associated(prev_read)) prev_read%next_read => current
          prev_read => current
          !treat the sub-option case
          select type (op => current%opt)
          class is(SubOption)
           sub_option => current%sub_list
          class default
          !do nothing
          end select
          return
        end if
      end if

      ! parse option
      call current%opt%parse(rest)
      current%opt%isSet = .true.
      current%prev_read => prev_read
      if(associated(prev_read)) prev_read%next_read => current
      prev_read => current

      !treat the sub-option case
      select type (op => current%opt)
       class is(SubOption)
        sub_option => current%sub_list
       class default
       !do nothing
      end select

      return
    end if
    current => current%next
  end do
end function read_option

!-----------------------------------------------------------------------

!subroutine check_arg_order(prms)
! class(t_parse),intent(in) :: prms  !< class(t_parse)
!
! type(t_link), pointer :: current
! type(t_link), pointer :: prev, next
!
! current => prms%firstLink
!
! do while(associated(current))
!
!  write(*,*) 'Parameter name: ',trim(current%opt%name)
!  prev => current%prev_read
!  if (associated(prev)) then
!    write(*,*) 'Previously read name: ',trim(prev%opt%name)
!  endif
!  next => current%next_read
!  if (associated(next)) then
!    write(*,*) 'Next read name: ',trim(next%opt%name)
!  endif
!
!
!  current => current%next
! enddo
!
!end subroutine

!> Check the consistency of a given option
!!
!! The given option is passed as a link, previously retrieved from a
!! get... subroutine. Then through previous or next optional inputs
!! is checked if the previous or next option actually introduced int the input
!! file is coherent with the required prev_opt or next_opt
subroutine check_opt_consistency(olink, prev, prev_opt, next, next_opt)
 type(t_link), intent(in) :: olink
 logical, intent(in), optional :: prev
 character(len=*), intent(in), optional :: prev_opt
 logical, intent(in), optional :: next
 character(len=*), intent(in), optional :: next_opt

 logical :: res
 character(len=*), parameter :: this_sub_name = 'check_opt_consistency'

 if(present(prev)) then
 if(prev .eqv. .true.) then

   if(.not.present(prev_opt)) call error(this_sub_name, this_mod_name, &
   'No prev_opt provided when asked to check previous argument')
   if (associated(olink%prev_read))  then
     res = olink%prev_read%opt%nameequals(trim(prev_opt))
     if (.not. res) call error(this_sub_name, this_mod_name, &
     'It is required that in the input file before the option '//&
     trim(olink%opt%name)//' the option '//trim(prev_opt)//' is introduced.&
     & However the option '//trim(olink%prev_read%opt%name)//' was found.')
   else
     call error(this_sub_name, this_mod_name, &
     'It is required that in the input file before the option '//&
     trim(olink%opt%name)//' the option '//trim(prev_opt)//' is introduced.&
     & However no option was found.')
   endif


 endif
 endif

 if(present(next)) then
 if(next .eqv. .true.) then

   if(.not.present(next_opt)) call error(this_sub_name, this_mod_name, &
   'No next_opt provided when asked to check next argument')

   res = olink%next_read%opt%nameequals(trim(next_opt))
   if (.not. res) call error(this_sub_name, this_mod_name, &
   'It is required that in the input file after the option '//&
   trim(olink%opt%name)//' the option '//trim(next_opt)//' is introduced.&
   & However the option '//trim(olink%next_read%opt%name)//' was found.')

 endif
 endif


end subroutine check_opt_consistency

subroutine print_parse_debug(prms)
 class(t_parse),intent(in) :: prms  !< class(t_parse)

 type(t_link),pointer :: current

  write(*,*) 'List of option: '
  current => prms%firstLink
  do while (associated(current))
    write(*,*) 'Option name: ',trim(current%opt%name)
    write(*,'(A)',advance='no') '       value: '
    call current%opt%printValue(20)
    write(*,*)
    write(*,*) '       is set? ',current%opt%isSet
    write(*,*) '       is removed? ',current%opt%isRemoved
    current => current%next
  end do


end subroutine print_parse_debug

!-----------------------------------------------------------------------

!> Output all parameters, which are defined but not set in the parameter file.
subroutine IgnoredParameters(prms)
 class(t_parse),intent(in) :: prms  !< class(t_parse)

 type(t_link), pointer :: current
 logical              :: prms_ignored

  !1) chech if some of the parameters were not collected
  prms_ignored = .false.
  current => prms%firstLink
  do while (associated(current))
    if (.not.current%opt%isRemoved) then
      prms_ignored = .true.
      exit
    end if
    current => current%next
  end do
  !2) Output the ignored parameters only if actually some were ignored
  if (prms_ignored) then
    current => prms%firstLink
    !call set_formatting("bright red")
    write(UNIT_StdOut,'(72("!"))')
    write(UNIT_StdOut,*) "warning: The following parameters were set in the &
    &input file, but are not currently employed :"
    do while (associated(current))
      if (.not.current%opt%isRemoved) then
        write(UNIT_StdOut,*) "   ", trim(current%opt%name)
      end if
      current => current%next
    end do
    write(UNIT_StdOut,'(72("!"))')
    !call clear_formatting()
  endif
end subroutine IgnoredParameters

!-----------------------------------------------------------------------

!> Print a default parameter file.
!!
!! The command line argument --help prints it in the format, that is used for
!! reading the parameter file. With --markdown one can print a default
!! parameter file in markdown format.
!! Also prints the descriptions of a single parameter or parameter sections,
!! if name corresponds to one of them.
subroutine PrintDefaultParameterFile(prms,markdown,name)
!TODO: remove module call from here
use MOD_StringTools ,only: stricmp
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 logical,intent(in)   :: markdown  !< marker whether markdown format is used for output
 character(len=max_char_len)   :: name      !< for this parameter help is printed. If empty print all.

 type(t_link), pointer   :: current
 class(option), pointer :: currentOpt
 integer                :: maxNameLen
 integer                :: maxValueLen
 integer                :: lineLen
 integer                :: spaceNameLen
 integer                :: spaceValueLen
 integer                :: mode
 character(len=max_char_len)     :: section = "-"
 character(len=max_char_len)     :: singlesection = ""
 character(len=max_char_len)     :: singleoption = ""
 character(len=20)      :: fmtLineLen
 character(len=20)      :: fmtName
 character(len=20)      :: fmtValue
 character(len=20)      :: fmtComment
 character(len=20)      :: fmtNamespace
 character(len=20)      :: fmtValuespace
 integer                :: i
 character(len=max_char_len)     :: intFromStringOutput
 character(len=max_char_len)     :: fmtIntFromStringLength
 character(len=max_char_len)     :: fmtStringIntFromString


  maxNameLen  = 0
  maxValueLen = 0
  current => prms%firstLink
  ! check if name is a section or a option
  do while (associated(current))
    if (stricmp(current%opt%section,name)) then
      singlesection = trim(name)
      exit
    end if
    if (current%opt%nameequals(name)) then
      singleoption = trim(name)
      singlesection = trim(current%opt%section)
      exit
    end if
    current => current%next
  end do

  ! if name is not specified, the complete parameter files needs to be printed
  if ((.not.markdown).and.(LEN_TRIM(name).eq.0)) then
    write(UNIT_StdOut,'(A80)')  "!===========================================&
      &===================================="
    write(UNIT_StdOut,'(A)')    "! Default Parameter File generated using  &
      &'dust --help' "
    write(UNIT_StdOut,'(A80)')  "!===========================================&
      &===================================="
  end if

  mode = 1
  if (markdown) then
    mode = 2
    write(UNIT_StdOut,'(A)') "## Parameterfile"
    write(UNIT_StdOut,'(A)') ""
  end if

  ! Find longest parameter name and length of the standard values
  current => prms%firstLink
  do while (associated(current))
    maxNameLen = max(maxNameLen, current%opt%getnamelen())
    maxValueLen = max(maxValueLen, current%opt%getvaluelen())
    current => current%next
  end do
  lineLen = maxNameLen + maxValueLen + 4 + 50
  spaceNameLen = maxNameLen - 9
  spaceValueLen = maxValueLen - 10
  write(fmtLineLen,*) lineLen
  write(fmtName,*)    maxNameLen
  write(fmtValue,*)   maxValueLen
  write(fmtComment,*) 50
  write(fmtNamespace,*) spaceNameLen
  write(fmtValuespace,*) spaceValueLen
  current => prms%firstLink
  do while (associated(current))
    if (stricmp(current%opt%section,'RecordPoints')) then
      current => current%next
      cycle
    end if
    if ((LEN_TRIM(singlesection).eq.0).or. &
        (stricmp(singlesection,current%opt%section))) then
      if (.not.stricmp(section,current%opt%section)) then
        section = current%opt%section
        if (markdown) then
          write(UNIT_StdOut,'('//fmtLineLen//'("-"))')
          write(UNIT_StdOut,'(A2,A,A2)')             "**",trim(section),"**"
          write(UNIT_StdOut,'('//fmtName//'("-")"--"A1)', advance='no')  " "
          write(UNIT_StdOut,'('//fmtValue//'("-")A1)', advance='no')     " "
          write(UNIT_StdOut,'('//fmtComment//'("-"))')
          write(UNIT_StdOut,'(A)', advance='no')              "**Variable**"
          write(UNIT_StdOut,'('//fmtNamespace//'(" "))', advance='no')
          write(UNIT_StdOut,'(A)', advance='no')               "**Default**"
          write(UNIT_StdOut,'('//fmtValuespace//'(" "))', advance='no')
          write(UNIT_StdOut,'(A)')                         "**Description**"
          write(UNIT_StdOut,'(A80)')                                     ""
        else
          write(UNIT_StdOut,'(A1,'//fmtLineLen//'("="))') "!"
          write(UNIT_StdOut,'(A2,A)') "! ", trim(section)
          write(UNIT_StdOut,'(A1,'//fmtLineLen//'("="))') "!"
        end if
      end if

      if ((LEN_TRIM(singleoption).eq.0).or. &
          (current%opt%nameequals(singleoption))) then
        call current%opt%print(maxNameLen, maxValueLen,mode)
      end if

      ! If help is called for a single IntFromStringOption,
      ! print the possible values of this parameter
      if (current%opt%nameequals(singleoption)) then
        currentOpt => current%opt
        select type(currentOpt)
        class is (IntFromStringOption)
          write(UNIT_StdOut,'(A)') 'Possible options for this parameter are:'
          write(fmtIntFromStringLength,*) currentOpt%maxLength   ! The biggest lenght of a named option
          write(fmtStringIntFromString,*) &
                                 "(A"//trim(fmtIntFromStringLength)//",A,I0,A)"
          do i=1,size(currentOpt%strList)
            ! Output is in the format string (integer)
            write(intFromStringOutput,trim(fmtStringIntFromString)) &
                  trim(currentOpt%strList(i)), ' (', currentOpt%intList(i), ')'
            write(UNIT_StdOut,'(A)') trim(intFromStringOutput)
          end do
        end select
      end if

      ! print ------ line at the end of a section in markdown mode
      if (associated(current%next).and.markdown) then
        if (.not.stricmp(section,current%next%opt%section)) then
          write(UNIT_StdOut,'('//fmtLineLen//'("-"))')
          write(UNIT_StdOut,*) ''
        end if
      end if
    end if
    current => current%next
  end do
end subroutine

!----------------------------------------------------------------------

!> Creates a new link to a option-object,
!! 'next' is the following link in the linked list
function constructor_Link(opt, next)
 type(t_link),pointer            :: constructor_Link  !< new link
 class(option),intent(in)       :: opt               !< option to be linked
 type(t_link),intent(in),pointer :: next              !< next link

  allocate(constructor_Link)
  constructor_Link%next => next
  allocate(constructor_Link%opt, source=opt)
end function constructor_Link

!----------------------------------------------------------------------

!> Count number of times a parameter is used within a file
!! in case of multiple parameters. This only calls the internal
!! function countoption_ of the parameters class.
function CountOption(prms,name) result(no)
  !TODO: remove this ugly call
  use MOD_Options
    class(t_parse),intent(inout) :: prms  !< class(t_parse)
    character(len=*),intent(in) :: name  !< parameter name
    integer                     :: no    !< number of parameters
      no = prms%CountOption_(name)
end function CountOption

!----------------------------------------------------------------------

!> General routine to get an option. This routine is called from getint,
!! getreal,getlogical,getstr to get the value a non-array  option.
subroutine GetGeneralOption(prms, value, name, proposal, olink)
!TODO: remove this horrible call
use MOD_Options
 class(t_parse),intent(in)           :: prms     !< class(t_parse)
 character(len=*),intent(in)          :: name     !< parameter name
 character(len=*),intent(in),optional :: proposal !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 class(*)                             :: value    !< parameter value

 type(t_link),pointer   :: current
 class(Option),pointer :: opt
 character(len=max_char_len)    :: proposal_loc
 character(len=*), parameter :: this_sub_name = 'GetGeneralOption'


  ! iterate over all options
  current => prms%firstLink
  do while (associated(current))
    ! if name matches option
    if (current%opt%nameequals(name).and.(.not.current%opt%isRemoved)) then
      opt => current%opt
      if(present(olink)) olink => current
      ! if proposal is present and the option is not set due to the parameter
      ! file, then return the proposal
      if ((present(proposal)).and.(.not.opt%isSet)) then
        proposal_loc = trim(proposal)
        call opt%parse(proposal_loc)
      else
        ! no proposal, no default and also not set in parameter file => abort
        if ((.not.opt%hasDefault).and.(.not.opt%isSet)) then
          call error(this_sub_name, this_mod_name, &
              "Required option '"//trim(name)//"' not set in parameter&
              &file and has no default value.")
          return
        end if
      end if
      if(opt%isSet) then
        ! copy value from option to result variable
        select type (opt)
        !class is (SubOption)
        !  select type(value)
        !  type is (t_parse)
        !    value => current%sub_list
        !  end select
        class is (IntOption)
          select type(value)
          type is (integer)
            value = opt%value
          end select
        class is (RealOption)
          select type(value)
          type is (real(wp))
            value = opt%value
          end select
        class is (LogicalOption)
          select type(value)
          type is (logical)
            value = opt%value
          end select
        class is (StringOption)
          select type(value)
          type is (t_str)
            value%chars = opt%value
          end select
        end select
      else
        ! copy value from default option to result variable
        select type (opt)
        !class is (SubOption)
        !  select type(value)
        !  type is (t_parse)
        !    value => current%sub_list
        !  end select
        class is (IntOption)
          select type(value)
          type is (integer)
            value = opt%default_value
          end select
        class is (RealOption)
          select type(value)
          type is (real(wp))
            value = opt%default_value
          end select
        class is (LogicalOption)
          select type(value)
          type is (logical)
            value = opt%default_value
          end select
        class is (StringOption)
          select type(value)
          type is (t_str)
            value%chars = opt%default_value
          end select
        end select
      endif
      !TODO: at the moment completely commented out, fix in another moment
      ! print option and value to stdout (if configured to do so)
      !if (prms%printout_val) &
      !  call opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
      ! remove the option from the linked list of all parameters
      current%opt%isRemoved = .true.
      return
    end if
    current => current%next
  end do
  write(*,*)'Option "'//trim(name)//'" is not defined in any DefineParameters... routine '//&
      'or already read (use get... routine only for multiple options more than once).'
  call abort
end subroutine GetGeneralOption

!----------------------------------------------------------------------

!> General routine to get an array option. This routine is called from
!! getintarray,getrealarray,getlogicalarray,getstrarray to get
!! the value an array option.
subroutine GetGeneralArrayOption(prms, value, name, no, proposal, olink)
!TODO: my god clear all this generic use
use MOD_Options
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name      !< parameter name
 integer,intent(in)                   :: no        !< size of array
 character(len=*),intent(in),optional :: proposal  !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 class(*)                             :: value(no) !< parameter value

 type(t_link),pointer   :: current
 class(Option),pointer :: opt
 character(len=max_char_len)    :: proposal_loc
 character(len=*), parameter :: this_sub_name = 'GetGeneralArrayOption'
 !integer               :: i

  ! iterate over all options
  current => prms%firstLink
  do while (associated(current))
    ! if name matches option
    if (current%opt%nameequals(name).and.(.not.current%opt%isRemoved)) then
      opt => current%opt
      if(present(olink)) olink => current
      ! if proposal is present and the option is not set due to the
      ! parameter file, then return the proposal
      if ((present(proposal)).and.(.not.opt%isSet)) then
        proposal_loc = trim(proposal)
        call opt%parse(proposal_loc)
      else
        ! no proposal, no default and also not set in parameter file => abort
        if ((.not.opt%hasDefault).and.(.not.opt%isSet)) then
          call error(this_sub_name, this_mod_name, "Required option '"// &
           trim(name)//"' not set in parameter file and has no default value.")
          return
        end if
      end if

      if(opt%isSet) then
        ! copy value from option to result variable
        select type (opt)
        class is (IntArrayOption)
          if (size(opt%value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (integer)
            value = opt%value
          end select
        class is (RealArrayOption)
          if (size(opt%value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (real(wp))
            value = opt%value
          end select
        class is (LogicalArrayOption)
          if (size(opt%value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (logical)
            value = opt%value
          end select
        !class is (StringArrayOption)
          !if (size(opt%value).ne.no) call Abort(__STAMP__,"Array size of
          !option '"//trim(name)//"' is not correct!")
          !select type(value)
          !type is (t_str)
            !do i=1,no
              !value(i)%chars = opt%value(i)
            !end do
          !end select
        end select
      else
        ! copy value from option to result variable
        select type (opt)
        class is (IntArrayOption)
          if (size(opt%default_value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (integer)
            value = opt%default_value
          end select
        class is (RealArrayOption)
          if (size(opt%default_value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (real(wp))
            value = opt%default_value
          end select
        class is (LogicalArrayOption)
          if (size(opt%default_value).ne.no) call error(this_sub_name, this_mod_name, &
                       "Array size of option '"//trim(name)//"' is not correct!")
          select type(value)
          type is (logical)
            value = opt%default_value
          end select
        !class is (StringArrayOption)
          !if (size(opt%default_value).ne.no) call Abort(__STAMP__,"Array size of
          !option '"//trim(name)//"' is not correct!")
          !select type(value)
          !type is (t_str)
            !do i=1,no
              !value(i)%chars = opt%default_value(i)
            !end do
          !end select
        end select
      endif
      ! print option and value to stdout
      !MATTEO: at the moment hardcoded NOT to print values
      !call opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
      ! remove the option from the linked list of all parameters
      current%opt%isRemoved = .true.
      return
    end if
    current => current%next
  end do
  call error(this_sub_name, this_mod_name, &
                'Option "'//trim(name)//'" is not defined in any &
                &DefineParameters... routine '//&
                 'or already read (use get... routine only for multiple &
                 &options more than once).')
end subroutine GetGeneralArrayOption

!----------------------------------------------------------------------

!> Get integer, where proposal is used as default value,
!! if the option was not set in parameter file
function getint(prms, name, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in) :: name              !< parameter name
 character(len=*),intent(in),optional :: proposal !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 integer                     :: value             !< parameter value

  value = -1
  call GetGeneralOption(prms, value, name, proposal, olink)
end function getint

!----------------------------------------------------------------------

!!> Get a sub option parse list
!function getsuboption(prms, name, proposal, olink) result(value)
! class(t_parse),intent(in) :: prms  !< class(t_parse)
! character(len=*),intent(in) :: name              !< parameter name
! character(len=*),intent(in),optional :: proposal !< reference value
! type(t_link), intent(out), pointer, optional :: olink   !< link to the option
! class(t_parse), pointer :: value             !< parameter value
!
!  value => null()
!
!  call GetGeneralOption(prms, value, name, proposal, olink)
!end function getsuboption

subroutine getsuboption(prms, name, value, olink) !result(value)
 class(t_parse),intent(in)           :: prms     !< class(t_parse)
 character(len=*),intent(in)          :: name     !< parameter name
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 type(t_parse), pointer, intent(out)         :: value    !< parameter value

 type(t_link),pointer   :: current
 class(Option),pointer :: opt
 !character(len=max_char_len)    :: proposal_loc
 character(len=*), parameter :: this_sub_name = 'GetSubOption'

  ! iterate over all options
  current => prms%firstLink
  do while (associated(current))
    ! if name matches option
    if (current%opt%nameequals(name).and.(.not.current%opt%isRemoved)) then
      opt => current%opt
      if(present(olink)) olink => current
      if ((.not.opt%hasDefault).and.(.not.opt%isSet)) then
          call error(this_sub_name, this_mod_name, &
              "Required option '"//trim(name)//"' not set in parameter&
              &file and has no default value.")
          return
      end if
      ! copy value from option to result variable
      select type (opt)
      class is (SubOption)
          value => current%sub_list
      class default
       write(*,*) 'No way, wrong kind'
      end select
      ! print option and value to stdout (if configured to do so)
      !if (prms%printout_val) &
      !  call opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
      ! remove the option from the linked list of all parameters
      current%opt%isRemoved = .true.
      return
    end if
    current => current%next
  end do
  write(*,*)'Option "'//trim(name)//'" is not defined in any DefineParameters... routine '//&
      'or already read (use get... routine only for multiple options more than once).'
  call abort
end subroutine getsuboption

!----------------------------------------------------------------------

!> Get logical, where proposal is used as default value,
!! if the option was not set in parameter file
function getlogical(prms, name, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in) :: name              !< parameter name
 character(len=*),intent(in),optional :: proposal !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 logical                     :: value             !< parameter value

  value = .false.
  call GetGeneralOption(prms, value, name, proposal, olink)
end function getlogical

!----------------------------------------------------------------------

!> Get real, where proposal is used as default value, if the option was not
!! set in parameter file
function getreal(prms, name, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name     !< parameter name
 character(len=*),intent(in),optional :: proposal !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 real(wp)                                 :: value    !< parameter value

  value = -1.0_wp
  call GetGeneralOption(prms, value, name, proposal, olink)
end function getreal

!----------------------------------------------------------------------

!> Get string, where proposal is used as default value, if the option was
!! not set in parameter file
function getstr(prms, name, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name     !< parameter name
 character(len=*),intent(in),optional :: proposal !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 character(len=max_char_len)                   :: value    !< parameter value

 type(t_str) :: tmp ! compiler bug workaround (gfortran 4.8.4)

  call GetGeneralOption(prms, tmp, name, proposal, olink)
  value = tmp%chars
end function getstr

!----------------------------------------------------------------------

!> Get integer array, where proposal is used as default value,
!! if the option was not set in parameter file
function getintarray(prms, name, no, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name      !< parameter name
 integer,intent(in)                   :: no        !< size of array
 character(len=*),intent(in),optional :: proposal  !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 integer                              :: value(no) !< array of integers

  value = -1
  call GetGeneralArrayOption(prms, value, name, no, proposal, olink)
end function getintarray

!----------------------------------------------------------------------

!> Get logical array, where proposal is used as default value,
!! if the option was not set in parameter file
function getlogicalarray(prms, name, no, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name      !< parameter name
 integer,intent(in)                   :: no        !< size of array
 character(len=*),intent(in),optional :: proposal  !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 logical                              :: value(no) !< array of logicals

  value = .false.
  call GetGeneralArrayOption(prms, value, name, no, proposal, olink)
end function getlogicalarray

!----------------------------------------------------------------------

!> Get real array, where proposal is used as default value,
!! if the option was not set in parameter file
function getrealarray(prms, name, no, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name      !< parameter name
 integer,intent(in)                   :: no        !< size of array
 character(len=*),intent(in),optional :: proposal  !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 real(wp)                                 :: value(no) !< array of reals

  value = -1.0_wp
  call GetGeneralArrayOption(prms, value, name, no, proposal, olink)
end function getrealarray

!----------------------------------------------------------------------

!> Get string array, where proposal is used as default value,
!! if the option was not set in parameter file
function getstrarray(prms, name, no, proposal, olink) result(value)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name      !< parameter name
 integer,intent(in)                   :: no        !< size of array
 character(len=*),intent(in),optional :: proposal  !< reference value
 type(t_link), intent(out), pointer, optional :: olink   !< link to the option
 character(len=max_char_len)                   :: value(no) !< array of strings

 type(t_str) :: tmp(no) ! compiler bug workaround (gfortran 4.8.4)
 integer      :: i

  call GetGeneralArrayOption(prms, tmp, name, no, proposal, olink)
  do i = 1, no
    value(i)=tmp(i)%chars
  end do ! i = 1, no
end function getstrarray

!----------------------------------------------------------------------

!> Get string array, where proposal is used as default value,
!! if the option was not set in parameter file
function getdescription(prms,name) result(description)
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)          :: name        !< parameter name
 character(len=1000)                  :: description !< description

 type(t_link),pointer :: current

  ! iterate over all options and compare names
  current => prms%firstLink
  do while (associated(current))
    if (current%opt%nameequals(name)) then
      description = current%opt%description
    end if
    current => current%next
  end do
end function getdescription

!----------------------------------------------------------------------

!> getint for options with string values.
!!
!! Requires a map that provides the link between the
!! possible integer values and the corresponding named values.
!! This map is set using the addStrListEntry routine during
!! parameter definition. If there is no named value to an option passed as
!! int a warning is returned.
function getintfromstr(prms, name) result(value)
!TODO: fix this use
use MOD_StringTools ,only: isint, stricmp
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)   :: name        !< parameter name
 integer                       :: value       !< return value

 type(t_link),pointer           :: current
 class(Option),pointer         :: opt
 integer                       :: i
 logical                       :: found
 integer                       :: listSize         ! current size of list
 character(len=*), parameter   :: this_sub_name = 'getintfromstr'

  ! iterate over all options and compare names
  current => prms%firstLink
  do while (associated(current))
    if (current%opt%nameequals(name)) then
      opt => current%opt
      select type (opt)
      class is (IntFromStringOption)
        ! Set flag indicating the given option has an entry in the mapping
        opt%foundInList = .true.
        ! Size of list with string-integer pairs
        listSize = size(opt%strList)
        ! Check if an integer has been specfied directly
        if (isint(opt%value)) then
          read(opt%value,*) value
          found=.false.
          ! Check if the integer is present in the list of possible integers
          do i=1,listSize
            if (opt%intList(i).eq.value)then
              found=.true.
              opt%listIndex = i ! Store index of the mapping
              exit
            end if
          end do
          ! If it is not found, print a warning and set the flag to later use
          ! the correct output format
          if(.not.found)then
            call warning(this_sub_name, this_mod_name, "No named option for &
            &parameter " //trim(name)// &
            " exists for this number, please ensure your input is correct.")
            opt%foundInList = .false.
          end if
          call opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
          return
        end if
        ! If a string has been supplied, check if this string exists in the
        ! list and set it's integer representation according to the
        ! mapping
        do i=1,listSize
          if (stricmp(opt%strList(i), opt%value)) then
            value = opt%intList(i)
            opt%listIndex = i ! Store index of the mapping
            call opt%print(prms%maxNameLen, prms%maxValueLen, mode=0)
            return
          end if
        end do
        call error(this_sub_name, this_mod_name, &
                     "Unknown value for option: "//trim(name))
      end select
    end if
    current => current%next
  end do
  call error(this_sub_name, this_mod_name, &
                     "Unknown value for option: "//trim(name))

end function getintfromstr

!----------------------------------------------------------------------

!> Add an entry to the mapping of string and integer values for the
!! StringToInt option.
subroutine addStrListEntry(prms, name,string_in,int_in)
!TODO: fix this use
!use MOD_Globals,     only: abort
 class(t_parse),intent(in) :: prms  !< class(t_parse)
 character(len=*),intent(in)    :: name      !< parameter name
 character(len=*),intent(in)    :: string_in !< (in) string used for the option value
 integer         ,intent(in)    :: int_in    !< (in) integer used internally for the option value

 type(t_link), pointer           :: current
 class(option),pointer          :: opt
 integer                        :: listSize         ! current size of list
 character(len=max_char_len),allocatable :: strListTmp(:)    ! temporary string list
 integer           ,allocatable :: intListTmp(:)    ! temporary integer list
 character(len=*), parameter    :: this_sub_name = 'addStrListEntry'

  ! iterate over all options and compare names
  current => prms%firstLink
  do while (associated(current))
    if (current%opt%nameequals(name)) then
      opt => current%opt
      select type (opt)
      class is (IntFromStringOption)
        ! Check if the arrays containing the string and integer values
        ! are already allocated
        if (.not.(allocated(opt%strList))) then
          ! This is the first call to addEntry, allocate the arrays
          ! with dimension one
          allocate(opt%strList(1))
          allocate(opt%intList(1))
          ! Store the values in the lists
          opt%strList(1) = trim(string_in)
          opt%intList(1) = int_in
          ! Save biggest length of string entry
          opt%maxLength = LEN_TRIM(string_in)
        else
          ! Subsequent call to addEntry, re-allocate the lists with
          ! one additional entry
          listSize = size(opt%strList)    ! opt size of the list
          ! store opt values in temporary arrays
          allocate(strListTmp(listSize))
          allocate(intListTmp(listSize))
          strListTmp = opt%strList
          intListTmp = opt%intList
          ! Deallocate and re-allocate the list arrays
          if(allocated(opt%strList)) deallocate(opt%strList)
          if(allocated(opt%intList)) deallocate(opt%intList)
          allocate(opt%strList(listSize+1))
          allocate(opt%intList(listSize+1))
          ! Re-write the old values
          opt%strList(1:listSize) = strListTmp
          opt%intList(1:listSize) = intListTmp
          ! Deallocate temp arrays
          if(allocated(strListTmp)) deallocate(strListTmp)
          if(allocated(intListTmp)) deallocate(intListTmp)
          ! Now save the actual new entry in the list
          opt%strList(listSize+1) = trim(string_in)
          opt%intList(listSize+1) = int_in
          ! Save biggest length of string entry
          opt%maxLength = max(opt%maxLength,LEN_TRIM(string_in))
        end if
        return
      class default
        call error(this_sub_name, this_mod_name, &
                     "Option is not of type IntFromString: "//trim(name))
      end select
    end if
    current => current%next
  end do
  call error(this_sub_name, this_mod_name, &
                     "Option not yet set: "//trim(name))

end subroutine addStrListEntry

!----------------------------------------------------------------------

!> Clear parameters list 'prms'.
recursive subroutine FinalizeParameters(prms)
 class(t_parse),intent(in) :: prms  !< class(t_parse)

 type(t_link), pointer         :: current, tmp

  current => prms%firstLink
  do while (associated(current%next))
    if(associated(current%sub_list)) call FinalizeParameters(current%sub_list)
    nullify(current%sub_list)
    deallocate(current%opt)
    nullify(current%opt)
    tmp => current%next
    deallocate(current)
    nullify(current)
    current => tmp
  end do
end subroutine FinalizeParameters

!----------------------------------------------------------------------

end module mod_parse
