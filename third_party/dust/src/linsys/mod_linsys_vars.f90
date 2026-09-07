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


!> Module to expose only the linsys type
!!
!! needed to break circular dependencies

module mod_linsys_vars


use mod_param, only: &
  wp

!----------------------------------------------------------------------

implicit none

public :: t_linsys

private

!----------------------------------------------------------------------
!> Type containing all the data necessary to solve the linear system
!! for the aerodynamic elements doublets intensities
type :: t_linsys

  !> Rank of the linear system
  integer :: rank

  !> Linear system matrix, will include Morino-kutta effect
  real(wp), allocatable :: A(:,:)

  !> Linear system matrix, will not include Morino-kutta effect
  real(wp), allocatable :: A_wake_free(:,:)

  !> Linear system right hand side
  real(wp), allocatable :: b(:)

  !> Linear system matrix (for pressure integral equation)
  real(wp), allocatable :: A_pres(:,:)

  !> Linear system right hand side (for pressure integral equation)
  real(wp), allocatable :: b_pres(:)

  !> LU solvers pivot permutation matrix (in vector form)
  integer, allocatable :: P(:)
  integer, allocatable :: P_pres(:)
  integer, allocatable :: P_wake_free(:) 
  !> Static part of the right hand side
  !!
  !! The right hand side contains contributions from all the surface panels.
  !! The contributions from surface panels that do not actually move can be
  !! calculated just once, than multiplied for the freestream velocity and
  !! finally the contribution of the moving panels is added
  real(wp), allocatable :: b_static(:,:)
  !> static part for integral equation (for pressure integral equation)
  real(wp), allocatable :: b_static_pres(:,:) ! not strictly required
  ! OSS -> b_static_pres = b_static(1:idSurfPan_static , 1:idSurfPan_static )
  ! OSS    with idSurfPan_static = idSurfPan(1:nstatic_SurfPan)

  !> Static part of the lifting line contribution
  !!
  !! The contribution of the lifting line is similar to wake contribution,
  !! however a part could be static on static and can be pre-computed
  real(wp), allocatable :: L_static(:,:)

  !> Static part of the actuator disk contribution
  !!
  !! The contribution of the actuator disk is similar to wake contribution,
  !! however a part could be static on static and can be pre-computed
  !real(wp), allocatable :: D_static(:,:)

  !> Result of the linear system solution (intensity of the panels doublets)
  real(wp), allocatable :: res(:)

  !> Result of the pressure linear system solution
  real(wp), allocatable :: res_pres(:)

  !> Results of the explicit part of the elements, has more than one column
  !! since holds previous records for extrapolation
  real(wp), allocatable :: res_expl(:,:)

  !> Number of static and moving panels
  integer :: nstatic, nmoving

  !> Number of static and moving surface panels
  integer :: nstatic_sp, nmoving_sp, n_sp

  !> Number of static and moving actuator disks
  !integer :: nstatic_ad, nmoving_ad, n_ad

  !> Number of static and moving explicit elements
  integer :: nstatic_expl, nmoving_expl

  !> Number of explicit elements
  integer :: n_expl

  !> Global id of surface panel elements ( copy of the structure in geo )
  integer, allocatable :: idSurfPan(:)

  !> Global id of surface panel elements ( copy of the structure in geo )
  integer, allocatable :: idSurfPanG2L(:)

  !> Skip dynamic update
  logical :: skip

  !> TL matrix: left hand side of linear doublet 
  real(wp), allocatable :: TL(:,:)

  !> TR matrix: right hand side of linear doublet
  real(wp), allocatable :: TR(:,:)

end type t_linsys

!----------------------------------------------------------------------

end module mod_linsys_vars
