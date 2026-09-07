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
!!=====================================================================

! Modified version of Rich Townsend's iso_varying_string.f90, minimal
! modifications have been made in order to make the internal string
! storage compatible with a null-terminated C string; the original api
! has not changed, it has been completed with a constructor from a
! pointer to a C null-terminated string and a function returning a C
! const char* pointer to an existing varying string. The iso_c_binding
! intrinsic module is now required.
!
!    Copyright 2003 Rich Townsend <rhdt@bartol.udel.edu>
!    Copyright 2011 Davide Cesari <dcesari69 at gmail dot com>
!
!    This file is part of FortranGIS.
!
! The original copyright notice follows:
! ******************************************************************************
! *                                                                            *
! * iso_varying_string.f90                                                     *
! *                                                                            *
! * Copyright (c) 2003, Rich Townsend <rhdt@bartol.udel.edu>                   *
! * All rights reserved.                                                       *
! *                                                                            *
! * Redistribution and use in source and binary forms, with or without         *
! * modification, are permitted provided that the following conditions are     *
! * met:                                                                       *
! *                                                                            *
! *  * Redistributions of source code must retain the above copyright notice,  *
! *    this list of conditions and the following disclaimer.                   *
! *  * Redistributions in binary form must reproduce the above copyright       *
! *    notice, this list of conditions and the following disclaimer in the     *
! *    documentation and/or other materials provided with the distribution.    *
! *                                                                            *
! * this software is provided by the copyright holders and contributors "as    *
! * is" and any express or implied warranties, including, but not limited to,  *
! * the implied warranties of merchantability and fitness for A particular     *
! * purpose are disclaimed. in no event shall the copyright owner or           *
! * contributors be liable for any direct, indirect, incidental, special,      *
! * exemplary, or consequential damages (including, but not limited to,        *
! * procurement of substitute goods or services; loss of use, data, or         *
! * profits; or business interruption) however caused and on any theory of     *
! * liability, whether in contract, strict liability, or tort (including       *
! * negligence or otherwise) arising in any way out of the use of this         *
! * software, even if advised of the possibility of such damage.               *
! *                                                                            *
! ******************************************************************************
!
! Author    : Rich Townsend <rhdt@bartol.udel.edu>
! Synopsis  : Definition of iso_varying_string module, conformant to the api
!             specified in iso/iec 1539-2:2000 (varying-length strings for
!             Fortran 95).
! Version   : 1.3-F
! Thanks    : Lawrie Schonfelder (bugfixes and design pointers), Walt Brainerd
!             (conversion to F).

module MOD_ISO_VARYING_STRING


use,intrinsic :: iso_c_binding

! No implicit typing

  implicit none

! Parameter definitions

  integer, parameter, private :: GET_BUFFER_LEN = 256

! Type definitions

  type, public :: varying_string
     private
     character(len=1), dimension(:), allocatable :: chars
  end type varying_string

! Interface blocks

  interface assignment(=)
     module procedure op_assign_CH_VS
     module procedure op_assign_VS_CH
  end interface assignment(=)

  interface operator(//)
     module procedure op_concat_VS_VS
     module procedure op_concat_CH_VS
     module procedure op_concat_VS_CH
  end interface operator(//)

  interface operator(==)
     module procedure op_eq_VS_VS
     module procedure op_eq_CH_VS
     module procedure op_eq_VS_CH
  end interface operator(==)

  interface operator(/=)
     module procedure op_ne_VS_VS
     module procedure op_ne_CH_VS
     module procedure op_ne_VS_CH
  end interface operator (/=)

  interface operator(<)
     module procedure op_lt_VS_VS
     module procedure op_lt_CH_VS
     module procedure op_lt_VS_CH
  end interface operator (<)

  interface operator(<=)
     module procedure op_le_VS_VS
     module procedure op_le_CH_VS
     module procedure op_le_VS_CH
  end interface operator (<=)

  interface operator(>=)
     module procedure op_ge_VS_VS
     module procedure op_ge_CH_VS
     module procedure op_ge_VS_CH
  end interface operator (>=)

  interface operator(>)
     module procedure op_gt_VS_VS
     module procedure op_gt_CH_VS
     module procedure op_gt_VS_CH
  end interface operator (>)

  interface adjustl
     module procedure adjustl_
  end interface adjustl

  interface adjustr
     module procedure adjustr_
  end interface adjustr

  interface char
     module procedure char_auto
     module procedure char_fixed
  end interface char

  interface iachar
     module procedure iachar_
  end interface iachar

  interface ichar
     module procedure ichar_
  end interface ichar

  interface index
     module procedure index_VS_VS
     module procedure index_CH_VS
     module procedure index_VS_CH
  end interface index

  interface len
     module procedure len_
  end interface len

  interface len_trim
     module procedure len_trim_
  end interface len_trim

  interface lge
     module procedure lge_VS_VS
     module procedure lge_CH_VS
     module procedure lge_VS_CH
  end interface lge

  interface lgt
     module procedure lgt_VS_VS
     module procedure lgt_CH_VS
     module procedure lgt_VS_CH
  end interface lgt

  interface lle
     module procedure lle_VS_VS
     module procedure lle_CH_VS
     module procedure lle_VS_CH
  end interface lle

  interface llt
     module procedure llt_VS_VS
     module procedure llt_CH_VS
     module procedure llt_VS_CH
  end interface llt

  interface repeat
     module procedure repeat_
  end interface repeat

  interface scan
     module procedure scan_VS_VS
     module procedure scan_CH_VS
     module procedure scan_VS_CH
  end interface scan

  interface trim
     module procedure trim_
  end interface trim

  interface verify
     module procedure verify_VS_VS
     module procedure verify_CH_VS
     module procedure verify_VS_CH
  end interface verify

  interface var_str
     module procedure var_str_
     module procedure var_str_c_ptr
  end interface var_str

  interface get
     module procedure get_
     module procedure get_unit
     module procedure get_set_VS
     module procedure get_set_CH
     module procedure get_unit_set_VS
     module procedure get_unit_set_CH
  end interface get

  interface put
     module procedure put_VS
     module procedure put_CH
     module procedure put_unit_VS
     module procedure put_unit_CH
  end interface put

  interface put_line
     module procedure put_line_VS
     module procedure put_line_CH
     module procedure put_line_unit_VS
     module procedure put_line_unit_CH
  end interface put_line

  interface extract
     module procedure extract_VS
     module procedure extract_CH
  end interface extract

  interface insert
     module procedure insert_VS_VS
     module procedure insert_CH_VS
     module procedure insert_VS_CH
     module procedure insert_CH_CH
  end interface insert

  interface remove
     module procedure remove_VS
     module procedure remove_CH
  end interface remove

  interface replace
     module procedure replace_VS_VS_auto
     module procedure replace_CH_VS_auto
     module procedure replace_VS_CH_auto
     module procedure replace_CH_CH_auto
     module procedure replace_VS_VS_fixed
     module procedure replace_CH_VS_fixed
     module procedure replace_VS_CH_fixed
     module procedure replace_CH_CH_fixed
     module procedure replace_VS_VS_VS_target
     module procedure replace_CH_VS_VS_target
     module procedure replace_VS_CH_VS_target
     module procedure replace_CH_CH_VS_target
     module procedure replace_VS_VS_CH_target
     module procedure replace_CH_VS_CH_target
     module procedure replace_VS_CH_CH_target
     module procedure replace_CH_CH_CH_target
  end interface

  interface split
     module procedure split_VS
     module procedure split_CH
  end interface split

  interface c_ptr_new
     module procedure c_ptr_new_VS
  end interface c_ptr_new


! Access specifiers

  public :: assignment(=)
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>=)
  public :: operator(>)
  public :: adjustl
  public :: adjustr
  public :: char
  public :: iachar
  public :: ichar
  public :: index
  public :: len
  public :: len_trim
  public :: lge
  public :: lgt
  public :: lle
  public :: llt
  public :: repeat
  public :: scan
  public :: trim
  public :: verify
  public :: var_str
  public :: get
  public :: put
  public :: put_line
  public :: extract
  public :: insert
  public :: remove
  public :: replace
  public :: split

  private :: op_assign_CH_VS
  private :: op_assign_VS_CH
  private :: op_concat_VS_VS
  private :: op_concat_CH_VS
  private :: op_concat_VS_CH
  private :: op_eq_VS_VS
  private :: op_eq_CH_VS
  private :: op_eq_VS_CH
  private :: op_ne_VS_VS
  private :: op_ne_CH_VS
  private :: op_ne_VS_CH
  private :: op_lt_VS_VS
  private :: op_lt_CH_VS
  private :: op_lt_VS_CH
  private :: op_le_VS_VS
  private :: op_le_CH_VS
  private :: op_le_VS_CH
  private :: op_ge_VS_VS
  private :: op_ge_CH_VS
  private :: op_ge_VS_CH
  private :: op_gt_VS_VS
  private :: op_gt_CH_VS
  private :: op_gt_VS_CH
  private :: adjustl_
  private :: adjustr_
  private :: char_auto
  private :: char_fixed
  private :: iachar_
  private :: ichar_
  private :: index_VS_VS
  private :: index_CH_VS
  private :: index_VS_CH
  private :: len_
  private :: len_trim_
  private :: lge_VS_VS
  private :: lge_CH_VS
  private :: lge_VS_CH
  private :: lgt_VS_VS
  private :: lgt_CH_VS
  private :: lgt_VS_CH
  private :: lle_VS_VS
  private :: lle_CH_VS
  private :: lle_VS_CH
  private :: llt_VS_VS
  private :: llt_CH_VS
  private :: llt_VS_CH
  private :: repeat_
  private :: scan_VS_VS
  private :: scan_CH_VS
  private :: scan_VS_CH
  private :: trim_
  private :: verify_VS_VS
  private :: verify_CH_VS
  private :: verify_VS_CH
  private :: var_str_
  private :: var_str_c_ptr
  private :: get_
  private :: get_unit
  private :: get_set_VS
  private :: get_set_CH
  private :: get_unit_set_VS
  private :: get_unit_set_CH
  private :: put_VS
  private :: put_CH
  private :: put_unit_VS
  private :: put_unit_CH
  private :: put_line_VS
  private :: put_line_CH
  private :: put_line_unit_VS
  private :: put_line_unit_CH
  private :: extract_VS
  private :: extract_CH
  private :: insert_VS_VS
  private :: insert_CH_VS
  private :: insert_VS_CH
  private :: insert_CH_CH
  private :: remove_VS
  private :: remove_CH
  private :: replace_VS_VS_auto
  private :: replace_CH_VS_auto
  private :: replace_VS_CH_auto
  private :: replace_CH_CH_auto
  private :: replace_VS_VS_fixed
  private :: replace_CH_VS_fixed
  private :: replace_VS_CH_fixed
  private :: replace_CH_CH_fixed
  private :: replace_VS_VS_VS_target
  private :: replace_CH_VS_VS_target
  private :: replace_VS_CH_VS_target
  private :: replace_CH_CH_VS_target
  private :: replace_VS_VS_CH_target
  private :: replace_CH_VS_CH_target
  private :: replace_VS_CH_CH_target
  private :: replace_CH_CH_CH_target
  private :: split_VS
  private :: split_CH

! Procedures

contains

!****

  elemental subroutine op_assign_CH_VS (var, exp)

    character(len=*), intent(out)    :: var
    type(varying_string), intent(in) :: exp

! Assign a varying string to a character string

    var = char(exp)

! Finish

    return

  end subroutine op_assign_CH_VS

!****

  elemental subroutine op_assign_VS_CH (var, exp)

    type(varying_string), intent(out) :: var
    character(len=*), intent(in)      :: exp

! Assign a character string to a varying string

    var = var_str(exp)

! Finish

    return

  end subroutine op_assign_VS_CH

!****

  elemental function op_concat_VS_VS (string_a, string_b) result (concat_string)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    type(varying_string)             :: concat_string

    integer                          :: len_string_a

! Concatenate two varying strings

    len_string_a = len(string_a)

    allocate(concat_string%chars(len_string_a+len(string_b)+1))
    concat_string%chars(:len_string_a) = string_a%chars(:len_string_a)
    concat_string%chars(len_string_a+1:) = string_b%chars(:)


! Finish

    return

  end function op_concat_VS_VS

!****

  elemental function op_concat_CH_VS (string_a, string_b) result (concat_string)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    type(varying_string)             :: concat_string

! Concatenate a character string and a varying
! string

    concat_string = op_concat_VS_VS(var_str(string_a), string_b)

! Finish

    return

  end function op_concat_CH_VS

!****

  elemental function op_concat_VS_CH (string_a, string_b) result (concat_string)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    type(varying_string)             :: concat_string

! Concatenate a varying string and a character
! string

    concat_string = op_concat_VS_VS(string_a, var_str(string_b))

! Finish

    return

  end function op_concat_VS_CH

!****

  elemental function op_eq_VS_VS (string_a, string_b) result (op_eq)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_eq

! Compare (==) two varying strings

    op_eq = char(string_a) == char(string_b)

! Finish

    return

  end function op_eq_VS_VS

!****

  elemental function op_eq_CH_VS (string_a, string_b) result (op_eq)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_eq

! Compare (==) a character string and a varying
! string

    op_eq = string_a == char(string_b)

! Finish

    return

  end function op_eq_CH_VS

!****

  elemental function op_eq_VS_CH (string_a, string_b) result (op_eq)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_eq

! Compare (==) a varying string and a character
! string

    op_eq = char(string_a) == string_b

! Finish

    return

  end function op_eq_VS_CH

!****

  elemental function op_ne_VS_VS (string_a, string_b) result (op_ne)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ne

! Compare (/=) two varying strings

    op_ne = char(string_a) /= char(string_b)

! Finish

    return

  end function op_ne_VS_VS

!****

  elemental function op_ne_CH_VS (string_a, string_b) result (op_ne)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ne

! Compare (/=) a character string and a varying
! string

    op_ne = string_a /= char(string_b)

! Finish

    return

  end function op_ne_CH_VS

!****

  elemental function op_ne_VS_CH (string_a, string_b) result (op_ne)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_ne

! Compare (/=) a varying string and a character
! string

    op_ne = char(string_a) /= string_b

! Finish

    return

  end function op_ne_VS_CH

!****

  elemental function op_lt_VS_VS (string_a, string_b) result (op_lt)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_lt

! Compare (<) two varying strings

    op_lt = char(string_a) < char(string_b)

! Finish

    return

  end function op_lt_VS_VS

!****

  elemental function op_lt_CH_VS (string_a, string_b) result (op_lt)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_lt

! Compare (<) a character string and a varying
! string

    op_lt = string_a < char(string_b)

! Finish

    return

  end function op_lt_CH_VS

!****

  elemental function op_lt_VS_CH (string_a, string_b) result (op_lt)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_lt

! Compare (<) a varying string and a character
! string

    op_lt = char(string_a) < string_b

! Finish

    return

  end function op_lt_VS_CH

!****

  elemental function op_le_VS_VS (string_a, string_b) result (op_le)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_le

! Compare (<=) two varying strings

    op_le = char(string_a) <= char(string_b)

! Finish

    return

  end function op_le_VS_VS

!****

  elemental function op_le_CH_VS (string_a, string_b) result (op_le)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_le

! Compare (<=) a character string and a varying
! string

    op_le = string_a <= char(string_b)

! Finish

    return

  end function op_le_CH_VS

!****

  elemental function op_le_VS_CH (string_a, string_b) result (op_le)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_le

! Compare (<=) a varying string and a character
! string

    op_le = char(string_a) <= string_b

! Finish

    return

  end function op_le_VS_CH

!****

  elemental function op_ge_VS_VS (string_a, string_b) result (op_ge)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ge

! Compare (>=) two varying strings

    op_ge = char(string_a) >= char(string_b)

! Finish

    return

  end function op_ge_VS_VS

!****

  elemental function op_ge_CH_VS (string_a, string_b) result (op_ge)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_ge

! Compare (>=) a character string and a varying
! string

    op_ge = string_a >= char(string_b)

! Finish

    return

  end function op_ge_CH_VS

!****

  elemental function op_ge_VS_CH (string_a, string_b) result (op_ge)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_ge

! Compare (>=) a varying string and a character
! string

    op_ge = char(string_a) >= string_b

! Finish

    return

  end function op_ge_VS_CH

!****

  elemental function op_gt_VS_VS (string_a, string_b) result (op_gt)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_gt

! Compare (>) two varying strings

    op_gt = char(string_a) > char(string_b)

! Finish

    return

  end function op_gt_VS_VS

!****

  elemental function op_gt_CH_VS (string_a, string_b) result (op_gt)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: op_gt

! Compare (>) a character string and a varying
! string

    op_gt = string_a > char(string_b)

! Finish

    return

  end function op_gt_CH_VS

!****

  elemental function op_gt_VS_CH (string_a, string_b) result (op_gt)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: op_gt

! Compare (>) a varying string and a character
! string

    op_gt = char(string_a) > string_b

! Finish

    return

  end function op_gt_VS_CH

!****

  elemental function adjustl_ (string) result (adjustl_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: adjustl_string

! Adjust the varying string to the left

    adjustl_string = adjustl(char(string))

! Finish

    return

  end function adjustl_

!****

  elemental function adjustr_ (string) result (adjustr_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: adjustr_string

! Adjust the varying string to the right

    adjustr_string = adjustr(char(string))

! Finish

    return

  end function adjustr_

!****

  elemental function len_ (string) result (length)

    type(varying_string), intent(in) :: string
    integer                          :: length

! Get the length of a varying string

    if(allocated(string%chars)) then
       length = size(string%chars)-1
    else
       length = 0
    endif

! Finish

    return

  end function len_

!****

  elemental function len_trim_ (string) result (length)

    type(varying_string), intent(in) :: string
    integer                          :: length

! Get the trimmed length of a varying string

    if(allocated(string%chars)) then
       length = LEN_TRIM(char(string))
    else
       length = 0
    endif

! Finish

    return

  end function len_trim_

!****

  pure function char_auto (string) result (char_string)

    type(varying_string), intent(in) :: string
    character(len=len(string))       :: char_string

    integer                          :: i_char

! Convert a varying string into a character string
! (automatic length)

    forall(i_char = 1:len(string))
       char_string(i_char:i_char) = string%chars(i_char)
    end forall

! Finish

    return

  end function char_auto

!****

  pure function char_fixed (string, length) result (char_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: length
    character(len=length)            :: char_string

! Convert a varying string into a character string
! (fixed length)

    char_string = char(string)

! Finish

    return

  end function char_fixed

!****

  elemental function iachar_ (c) result (i)

    type(varying_string), intent(in) :: c
    integer                          :: i

! Get the position in the iso 646 collating sequence
! of a varying string character

    i = ichar(char(c))

! Finish

    return

  end function iachar_

!****

  elemental function ichar_ (c) result (i)

    type(varying_string), intent(in) :: c
    integer                          :: i

! Get the position in the processor collating
! sequence of a varying string character

    i = ichar(char(c))

! Finish

    return

  end function ichar_

!****

  elemental function index_VS_VS (string, substring, back) result (i_substring)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a varying substring within a
! varying string

    i_substring = index(char(string), char(substring), back)

! Finish

    return

  end function index_VS_VS

!****

  elemental function index_CH_VS (string, substring, back) result (i_substring)

    character(len=*), intent(in)     :: string
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a varying substring within a
! character string

    i_substring = index(string, char(substring), back)

! Finish

    return

  end function index_CH_VS

!****

  elemental function index_VS_CH (string, substring, back) result (i_substring)

    type(varying_string), intent(in) :: string
    character(len=*), intent(in)     :: substring
    logical, intent(in), optional    :: back
    integer                          :: i_substring

! Get the index of a character substring within a
! varying string

    i_substring = index(char(string), substring, back)

! Finish

    return

  end function index_VS_CH

!****

  elemental function lge_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lge) two varying strings

    comp = (char(string_a) >= char(string_b))

! Finish

    return

  end function lge_VS_VS

!****

  elemental function lge_CH_VS (string_a, string_b) result (comp)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lge) a character string and a varying
! string

    comp = (string_a >= char(string_b))

! Finish

    return

  end function lge_CH_VS

!****

  elemental function lge_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (lge) a varying string and a character
! string

    comp = (char(string_a) >= string_b)

! Finish

    return

  end function lge_VS_CH

!****

  elemental function lgt_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lgt) two varying strings

    comp = (char(string_a) > char(string_b))

! Finish

    return

  end function lgt_VS_VS

!****

  elemental function lgt_CH_VS (string_a, string_b) result (comp)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lgt) a character string and a varying
! string

    comp = (string_a > char(string_b))

! Finish

    return

  end function lgt_CH_VS

!****

  elemental function lgt_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (lgt) a varying string and a character
! string

    comp = (char(string_a) > string_b)

! Finish

    return

  end function lgt_VS_CH

!****

  elemental function lle_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lle) two varying strings

    comp = (char(string_a) <= char(string_b))

! Finish

    return

  end function lle_VS_VS

!****

  elemental function lle_CH_VS (string_a, string_b) result (comp)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (lle) a character string and a varying
! string

    comp = (string_a <= char(string_b))

! Finish

    return

  end function lle_CH_VS

!****

  elemental function lle_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (lle) a varying string and a character
! string

    comp = (char(string_a) <= string_b)

! Finish

    return

  end function lle_VS_CH

!****

  elemental function llt_VS_VS (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (llt) two varying strings

    comp = (char(string_a) < char(string_b))

! Finish

    return

  end function llt_VS_VS

!****

  elemental function llt_CH_VS (string_a, string_b) result (comp)

    character(len=*), intent(in)     :: string_a
    type(varying_string), intent(in) :: string_b
    logical                          :: comp

! Compare (llt) a character string and a varying
! string

    comp = (string_a < char(string_b))

! Finish

    return

  end function llt_CH_VS

!****

  elemental function llt_VS_CH (string_a, string_b) result (comp)

    type(varying_string), intent(in) :: string_a
    character(len=*), intent(in)     :: string_b
    logical                          :: comp

! Compare (llt) a varying string and a character
! string

    comp = (char(string_a) < string_b)

! Finish

    return

  end function llt_VS_CH

!****

  elemental function repeat_ (string, ncopies) result (repeat_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: ncopies
    type(varying_string)             :: repeat_string

! Concatenate several copies of a varying string

    repeat_string = var_str(repeat(char(string), ncopies))

! Finish

    return

  end function repeat_

!****

  elemental function scan_VS_VS (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a varying string for occurrences of
! characters in a varying-string set

    i = scan(char(string), char(set), back)

! Finish

    return

  end function scan_VS_VS

!****

  elemental function scan_CH_VS (string, set, back) result (i)

    character(len=*), intent(in)     :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a character string for occurrences of
! characters in a varying-string set

    i = scan(string, char(set), back)

! Finish

    return

  end function scan_CH_VS

!****

  elemental function scan_VS_CH (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    character(len=*), intent(in)     :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Scan a varying string for occurrences of
! characters in a character-string set

    i = scan(char(string), set, back)

! Finish

    return

  end function scan_VS_CH

!****

  elemental function trim_ (string) result (trim_string)

    type(varying_string), intent(in) :: string
    type(varying_string)             :: trim_string

! Remove trailing blanks from a varying string

    trim_string = trim(char(string))

! Finish

    return

  end function trim_

!****

  elemental function verify_VS_VS (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a varying string for occurrences of
! characters in a varying-string set

    i = verify(char(string), char(set), back)

! Finish

    return

  end function verify_VS_VS

!****

  elemental function verify_CH_VS (string, set, back) result (i)

    character(len=*), intent(in)     :: string
    type(varying_string), intent(in) :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a character string for occurrences of
! characters in a varying-string set

    i = verify(string, char(set), back)

! Finish

    return

  end function verify_CH_VS

!****

  elemental function verify_VS_CH (string, set, back) result (i)

    type(varying_string), intent(in) :: string
    character(len=*), intent(in)     :: set
    logical, intent(in), optional    :: back
    integer                          :: i

! Verify a varying string for occurrences of
! characters in a character-string set

    i = verify(char(string), set, back)

! Finish

    return

  end function verify_VS_CH

!****

  elemental function var_str_ (char_) result (string)

    character(len=*), intent(in) :: char_
    type(varying_string)         :: string

    integer                      :: length
    integer                      :: i_char

! Convert a character string to a varying string

    length = len(char_)

    allocate(string%chars(length+1))

    forall(i_char = 1:length)
       string%chars(i_char) = char_(i_char:i_char)
    end forall
    string%chars(length+1) = char(0)

! Finish

    return

  end function var_str_

!****

  function var_str_c_ptr (char_c_ptr) result (string)

    type(c_ptr), intent(in) :: char_c_ptr
    type(varying_string)         :: string

    character(len=1),pointer :: char_(:)
    integer :: length

! Convert a character string to a varying string

    if (c_ASSOCIATED(char_c_ptr)) then

      call C_F_POINTER(char_c_ptr, char_, (/huge(1)-1/))

      do length = 1, size(char_)
        if (char_(length) == char(0)) exit
      enddo

      allocate(string%chars(length))
      string%chars(:) = char_(1:length)
      string%chars(length) = char(0) ! handle absurdus huge() case

    else

      string = var_str('')

    endif

! Finish

    return

  end function var_str_c_ptr

!****

  subroutine get_ (string, maxlen, iostat)

    type(varying_string), intent(out) :: string
    integer, intent(in), optional     :: maxlen
    integer, intent(out), optional    :: iostat

    integer                           :: n_chars_remain
    integer                           :: n_chars_read
    character(len=GET_BUFFER_LEN)     :: buffer
    integer                           :: local_iostat

! Read from the default unit into a varying string

    string = ""

    if(present(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = huge(1)
    endif

    read_loop : do

       if(n_chars_remain <= 0) return

       n_chars_read = min(n_chars_remain, GET_BUFFER_LEN)

       if(present(iostat)) then
          read(unit=*, fmt="(A)", advance="no", &
               iostat=iostat, size=n_chars_read) buffer(:n_chars_read)
          if(iostat < 0) exit read_loop
          if(iostat > 0) return
       else
          read(unit=*, fmt="(A)", advance="no", &
               iostat=local_iostat, size=n_chars_read) buffer(:n_chars_read)
          if(local_iostat < 0) exit read_loop
       endif

       string = string//buffer(:n_chars_read)
       n_chars_remain = n_chars_remain - n_chars_read

    end do read_loop

    string = string//buffer(:n_chars_read)

! Finish (end-of-record)

    return

  end subroutine get_

!****

  subroutine get_unit (unit, string, maxlen, iostat)

    integer, intent(in)               :: unit
    type(varying_string), intent(out) :: string
    integer, intent(in), optional     :: maxlen
    integer, intent(out), optional    :: iostat

    integer                           :: n_chars_remain
    integer                           :: n_chars_read
    character(len=GET_BUFFER_LEN)     :: buffer
    integer                           :: local_iostat

! Read from the specified unit into a varying string

    string = ""

    if(present(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = huge(1)
    endif

    read_loop : do

       if(n_chars_remain <= 0) return

       n_chars_read = min(n_chars_remain, GET_BUFFER_LEN)

       if(present(iostat)) then
          read(unit=unit, fmt="(A)", advance="no", &
               iostat=iostat, size=n_chars_read) buffer(:n_chars_read)
          if(iostat < 0) exit read_loop
          if(iostat > 0) return
       else
          read(unit=unit, fmt="(A)", advance="no", &
               iostat=local_iostat, size=n_chars_read) buffer(:n_chars_read)
          if(local_iostat < 0) exit read_loop
       endif

       string = string//buffer(:n_chars_read)
       n_chars_remain = n_chars_remain - n_chars_read

    end do read_loop

    string = string//buffer(:n_chars_read)

! Finish (end-of-record)

    return

  end subroutine get_unit

!****

  subroutine get_set_VS (string, set, separator, maxlen, iostat)

    type(varying_string), intent(out)           :: string
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

! Read from the default unit into a varying string,
! with a custom varying-string separator

    call get(string, char(set), separator, maxlen, iostat)

! Finish

    return

  end subroutine get_set_VS

!****

  subroutine get_set_CH (string, set, separator, maxlen, iostat)

    type(varying_string), intent(out)           :: string
    character(len=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

    integer                                     :: n_chars_remain
    character(len=1)                            :: buffer
    integer                                     :: i_set
    integer                                     :: local_iostat

! Read from the default unit into a varying string,
! with a custom character-string separator

    string = ""

    if(present(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = huge(1)
    endif

    if(present(separator)) separator = ""

    read_loop : do

       if(n_chars_remain <= 0) return

       if(present(iostat)) then
          read(unit=*, fmt="(A1)", advance="no", iostat=iostat) buffer
          if(iostat /= 0) exit read_loop
       else
          read(unit=*, fmt="(A1)", advance="no", iostat=local_iostat) buffer
          if(local_iostat /= 0) exit read_loop
       endif

       i_set = scan(buffer, set)

       if(i_set == 1) then
          if(present(separator)) separator = buffer
          exit read_loop
       endif

       string = string//buffer
       n_chars_remain = n_chars_remain - 1

    end do read_loop

! Finish

    return

  end subroutine get_set_CH

!****

  subroutine get_unit_set_VS (unit, string, set, separator, maxlen, iostat)

    integer, intent(in)                         :: unit
    type(varying_string), intent(out)           :: string
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

! Read from the specified unit into a varying string,
! with a custom varying-string separator

    call get(unit, string, char(set), separator, maxlen, iostat)

! Finish

    return

  end subroutine get_unit_set_VS

!****

  subroutine get_unit_set_CH (unit, string, set, separator, maxlen, iostat)

    integer, intent(in)                         :: unit
    type(varying_string), intent(out)           :: string
    character(len=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    integer, intent(in), optional               :: maxlen
    integer, intent(out), optional              :: iostat

    integer                                     :: n_chars_remain
    character(len=1)                            :: buffer
    integer                                     :: i_set
    integer                                     :: local_iostat

! Read from the default unit into a varying string,
! with a custom character-string separator

    string = ""

    if(present(maxlen)) then
       n_chars_remain = maxlen
    else
       n_chars_remain = huge(1)
    endif

    if(present(separator)) separator = ""

    read_loop : do

       if(n_chars_remain <= 0) return

       if(present(iostat)) then
          read(unit=unit, fmt="(A1)", advance="no", iostat=iostat) buffer
          if(iostat /= 0) exit read_loop
       else
          read(unit=unit, fmt="(A1)", advance="no", iostat=local_iostat) buffer
          if(local_iostat /= 0) exit read_loop
       endif

       i_set = scan(buffer, set)

       if(i_set == 1) then
          if(present(separator)) separator = buffer
          exit read_loop
       endif

       string = string//buffer
       n_chars_remain = n_chars_remain - 1

    end do read_loop

! Finish

    return

  end subroutine get_unit_set_CH

!****

  subroutine put_VS (string, iostat)

    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of
! the default unit

    call put(char(string), iostat)

! Finish

  end subroutine put_VS

!****

  subroutine put_CH (string, iostat)

    character(len=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a character string to the current record of
! the default unit

    if(present(iostat)) then
       write(unit=*, fmt="(A)", advance="no", iostat=iostat) string
    else
       write(unit=*, fmt="(A)", advance="no") string
    endif

! Finish

  end subroutine put_CH

!****

  subroutine put_unit_VS (unit, string, iostat)

    integer, intent(in)              :: unit
    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of
! the specified unit

    call put(unit, char(string), iostat)

! Finish

    return

  end subroutine put_unit_VS

!****

  subroutine put_unit_CH (unit, string, iostat)

    integer, intent(in)            :: unit
    character(len=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a character string to the current record of
! the specified unit

    if(present(iostat)) then
       write(unit=unit, fmt="(A)", advance="no", iostat=iostat) string
    else
       write(unit=unit, fmt="(A)", advance="no") string
    endif

! Finish

    return

  end subroutine put_unit_CH

!****

  subroutine put_line_VS (string, iostat)

    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of
! the default unit, terminating the record

    call put_line(char(string), iostat)

! Finish

    return

  end subroutine put_line_VS

!****

  subroutine put_line_CH (string, iostat)

    character(len=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a varying string to the current record of
! the default unit, terminating the record

    if(present(iostat)) then
       write(unit=*, fmt="(A,/)", advance="no", iostat=iostat) string
    else
       write(unit=*, fmt="(A,/)", advance="no") string
    endif

! Finish

    return

  end subroutine put_line_CH

!****

  subroutine put_line_unit_VS (unit, string, iostat)

    integer, intent(in)              :: unit
    type(varying_string), intent(in) :: string
    integer, intent(out), optional   :: iostat

! Append a varying string to the current record of
! the specified unit, terminating the record

    call put_line(unit, char(string), iostat)

! Finish

    return

  end subroutine put_line_unit_VS

!****

  subroutine put_line_unit_CH (unit, string, iostat)

    integer, intent(in)            :: unit
    character(len=*), intent(in)   :: string
    integer, intent(out), optional :: iostat

! Append a varying string to the current record of
! the specified unit, terminating the record

    if(present(iostat)) then
       write(unit=unit, fmt="(A,/)", advance="no", iostat=iostat) string
    else
       write(unit=unit, fmt="(A,/)", advance="no") string
    endif

! Finish

    return

  end subroutine put_line_unit_CH

!****

  elemental function extract_VS (string, start, finish) result (ext_string)

    type(varying_string), intent(in) :: string
    integer, intent(in), optional    :: start
    integer, intent(in), optional    :: finish
    type(varying_string)             :: ext_string

! Extract a varying substring from a varying string

    ext_string = extract(char(string), start, finish)

! Finish

    return

  end function extract_VS

!****

  elemental function extract_CH (string, start, finish) result (ext_string)

    character(len=*), intent(in)  :: string
    integer, intent(in), optional :: start
    integer, intent(in), optional :: finish
    type(varying_string)          :: ext_string

    integer                       :: start_
    integer                       :: finish_

! Extract a varying substring from a character string

    if(present(start)) then
       start_ = max(1, start)
    else
       start_ = 1
    endif

    if(present(finish)) then
       finish_ = min(len(string), finish)
    else
       finish_ = len(string)
    endif

    ext_string = var_str(string(start_:finish_))

! Finish

    return

  end function extract_CH

!****

  elemental function insert_VS_VS (string, start, substring) result (ins_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: ins_string

! Insert a varying substring into a varying string

    ins_string = insert(char(string), start, char(substring))

! Finish

    return

  end function insert_VS_VS

!****

  elemental function insert_CH_VS (string, start, substring) result (ins_string)

    character(len=*), intent(in)     :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: ins_string

! Insert a varying substring into a character string

    ins_string = insert(string, start, char(substring))

! Finish

    return

  end function insert_CH_VS

!****

  elemental function insert_VS_CH (string, start, substring) result (ins_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    character(len=*), intent(in)     :: substring
    type(varying_string)             :: ins_string

! Insert a character substring into a varying string

    ins_string = insert(char(string), start, substring)

! Finish

    return

  end function insert_VS_CH

!****

  elemental function insert_CH_CH (string, start, substring) result (ins_string)

    character(len=*), intent(in) :: string
    integer, intent(in)          :: start
    character(len=*), intent(in) :: substring
    type(varying_string)         :: ins_string

    integer                      :: start_

! Insert a character substring into a character
! string

    start_ = max(1, min(start, len(string)+1))

    ins_string = var_str(string(:start_-1)//substring//string(start_:))

! Finish

    return

  end function insert_CH_CH

!****

  elemental function remove_VS (string, start, finish) result (rem_string)

    type(varying_string), intent(in) :: string
    integer, intent(in), optional    :: start
    integer, intent(in), optional    :: finish
    type(varying_string)             :: rem_string

! Remove a substring from a varying string

    rem_string = remove(char(string), start, finish)

! Finish

    return

  end function remove_VS

!****

  elemental function remove_CH (string, start, finish) result (rem_string)

    character(len=*), intent(in)  :: string
    integer, intent(in), optional :: start
    integer, intent(in), optional :: finish
    type(varying_string)          :: rem_string

    integer                       :: start_
    integer                       :: finish_

! Remove a substring from a character string

    if(present(start)) then
       start_ = max(1, start)
    else
       start_ = 1
    endif

    if(present(finish)) then
       finish_ = min(len(string), finish)
    else
       finish_ = len(string)
    endif

    if(finish_ >= start_) then
       rem_string = var_str(string(:start_-1)//string(finish_+1:))
    else
       rem_string = string
    endif

! Finish

    return

  end function remove_CH

!****

  elemental function replace_VS_VS_auto (string, start, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring

    rep_string = replace(char(string), start, max(start, 1)+len(substring)-1, char(substring))

! Finish

    return

  end function replace_VS_VS_auto

!****

  elemental function replace_CH_VS_auto (string, start, substring) result (rep_string)

    character(len=*), intent(in)     :: string
    integer, intent(in)              :: start
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring

    rep_string = replace(string, start, max(start, 1)+len(substring)-1, char(substring))

! Finish

    return

  end function replace_CH_VS_auto

!****

  elemental function replace_VS_CH_auto (string, start, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    character(len=*), intent(in)     :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring

    rep_string = replace(char(string), start, max(start, 1)+len(substring)-1, substring)

! Finish

    return

  end function replace_VS_CH_auto

!****

  elemental function replace_CH_CH_auto (string, start, substring) result (rep_string)

    character(len=*), intent(in) :: string
    integer, intent(in)          :: start
    character(len=*), intent(in) :: substring
    type(varying_string)         :: rep_string

! Replace part of a character string with a character
! substring

    rep_string = replace(string, start, max(start, 1)+len(substring)-1, substring)

! Finish

    return

  end function replace_CH_CH_auto

!****

  elemental function replace_VS_VS_fixed (string, start, finish, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring

    rep_string = replace(char(string), start, finish, char(substring))

! Finish

    return

  end function replace_VS_VS_fixed

!****

!****

  elemental function replace_CH_VS_fixed (string, start, finish, substring) result (rep_string)

    character(len=*), intent(in)     :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    type(varying_string), intent(in) :: substring
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring

    rep_string = replace(string, start, finish, char(substring))

! Finish

    return

  end function replace_CH_VS_fixed

!****

  elemental function replace_VS_CH_fixed (string, start, finish, substring) result (rep_string)

    type(varying_string), intent(in) :: string
    integer, intent(in)              :: start
    integer, intent(in)              :: finish
    character(len=*), intent(in)     :: substring
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring

    rep_string = replace(char(string), start, finish, substring)

! Finish

    return

  end function replace_VS_CH_fixed

!****

  elemental function replace_CH_CH_fixed (string, start, finish, substring) result (rep_string)

    character(len=*), intent(in) :: string
    integer, intent(in)          :: start
    integer, intent(in)          :: finish
    character(len=*), intent(in) :: substring
    type(varying_string)         :: rep_string

    integer                      :: start_
    integer                      :: finish_

! Replace part of a character string with a character
! substring

    start_ = max(1, start)
    finish_ = min(len(string), finish)

    if(finish_ < start_) then
       rep_string = insert(string, start_, substring)
    else
       rep_string = var_str(string(:start_-1)//substring//string(finish_+1:))
    endif

! Finish

    return

  end function replace_CH_CH_fixed

!****

  elemental function replace_VS_VS_VS_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a varying
! substring, at a location matching a varying-
! string target

    rep_string = replace(char(string), char(target), char(substring), every, back)

! Finish

    return

  end function replace_VS_VS_VS_target

!****

  elemental function replace_CH_VS_VS_target (string, target, substring, every, back) result (rep_string)

    character(len=*), intent(in)     :: string
    type(varying_string), intent(in) :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a varying-
! string target

    rep_string = replace(string, char(target), char(substring), every, back)

! Finish

    return

  end function replace_CH_VS_VS_target

!****

  elemental function replace_VS_CH_VS_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    character(len=*), intent(in)     :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a character-
! string target

    rep_string = replace(char(string), target, char(substring), every, back)

! Finish

    return

  end function replace_VS_CH_VS_target

!****

  elemental function replace_CH_CH_VS_target (string, target, substring, every, back) result (rep_string)

    character(len=*), intent(in)     :: string
    character(len=*), intent(in)     :: target
    type(varying_string), intent(in) :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a varying
! substring, at a location matching a character-
! string target

    rep_string = replace(string, target, char(substring), every, back)

! Finish

    return

  end function replace_CH_CH_VS_target

!****

  elemental function replace_VS_VS_CH_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    type(varying_string), intent(in) :: target
    character(len=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring, at a location matching a varying-
! string target

    rep_string = replace(char(string), char(target), substring, every, back)

! Finish

    return

  end function replace_VS_VS_CH_target

!****

  elemental function replace_CH_VS_CH_target (string, target, substring, every, back) result (rep_string)

    character(len=*), intent(in)     :: string
    type(varying_string), intent(in) :: target
    character(len=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a character string with a character
! substring, at a location matching a varying-
! string target

    rep_string = replace(string, char(target), substring, every, back)

! Finish

    return

  end function replace_CH_VS_CH_target

!****

  elemental function replace_VS_CH_CH_target (string, target, substring, every, back) result (rep_string)

    type(varying_string), intent(in) :: string
    character(len=*), intent(in)     :: target
    character(len=*), intent(in)     :: substring
    logical, intent(in), optional    :: every
    logical, intent(in), optional    :: back
    type(varying_string)             :: rep_string

! Replace part of a varying string with a character
! substring, at a location matching a character-
! string target

    rep_string = replace(char(string), target, substring, every, back)

! Finish

    return

  end function replace_VS_CH_CH_target

!****

  elemental function replace_CH_CH_CH_target (string, target, substring, every, back) result (rep_string)

    character(len=*), intent(in)  :: string
    character(len=*), intent(in)  :: target
    character(len=*), intent(in)  :: substring
    logical, intent(in), optional :: every
    logical, intent(in), optional :: back
    type(varying_string)          :: rep_string

    logical                       :: every_
    logical                       :: back_
    type(varying_string)          :: work_string
    integer                       :: length_target
    integer                       :: i_target

! Handle special cases when len(target) == 0. Such
! instances are prohibited by the standard, but
! since this function is elemental, no error can be
! thrown. Therefore, it makes sense to handle them
! in a sensible manner

    if(len(target) == 0) then
       if(len(string) /= 0) then
          rep_string = string
       else
          rep_string = substring
       endif
       return
    end if

! Replace part of a character string with a character
! substring, at a location matching a character-
! string target

    if(present(every)) then
       every_ = every
    else
       every_ = .false.
    endif

    if(present(back)) then
       back_ = back
    else
       back_ = .false.
    endif

    rep_string = ""

    work_string = string

    length_target = len(target)

    replace_loop : do

       i_target = index(work_string, target, back_)

       if(i_target == 0) exit replace_loop

       if(back_) then
          rep_string = substring//extract(work_string, start=i_target+length_target)//rep_string
          work_string = extract(work_string, finish=i_target-1)
       else
          rep_string = rep_string//extract(work_string, finish=i_target-1)//substring
          work_string = extract(work_string, start=i_target+length_target)
       endif

       if(.not. every_) exit replace_loop

    end do replace_loop

    if(back_) then
       rep_string = work_string//rep_string
    else
       rep_string = rep_string//work_string
    endif

! Finish

    return

  end function replace_CH_CH_CH_target

!****

  elemental subroutine split_VS (string, word, set, separator, back)

    type(varying_string), intent(inout)         :: string
    type(varying_string), intent(out)           :: word
    type(varying_string), intent(in)            :: set
    type(varying_string), intent(out), optional :: separator
    logical, intent(in), optional               :: back

! Split a varying string into two verying strings

    call split_CH(string, word, char(set), separator, back)

! Finish

    return

  end subroutine split_VS

!****

  elemental subroutine split_CH (string, word, set, separator, back)

    type(varying_string), intent(inout)         :: string
    type(varying_string), intent(out)           :: word
    character(len=*), intent(in)                :: set
    type(varying_string), intent(out), optional :: separator
    logical, intent(in), optional               :: back

    logical                                     :: back_
    integer                                     :: i_separator

! Split a varying string into two verying strings

    if(present(back)) then
       back_ = back
    else
       back_ = .false.
    endif

    i_separator = scan(string, set, back_)

    if(i_separator /= 0) then

       if(back_) then
          word = extract(string, start=i_separator+1)
          if(present(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, finish=i_separator-1)
       else
          word = extract(string, finish=i_separator-1)
          if(present(separator)) separator = extract(string, start=i_separator, finish=i_separator)
          string = extract(string, start=i_separator+1)
       endif

    else

       word = string
       if(present(separator)) separator = ""
       string = ""

    endif

! Finish

    return

  end subroutine split_CH


  function c_ptr_new_VS(string) result(c_ptr_new)
  type(varying_string),intent(in),target :: string
  type(c_ptr) :: c_ptr_new

  c_ptr_new = C_LOC(string%chars(1))

  end function c_ptr_new_VS


end module MOD_ISO_VARYING_STRING

