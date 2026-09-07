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
!!          Davide Montagnani
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================



module mod_precice_rbf

use mod_param, only: &
  wp

use mod_math, only: &
  sort_vector_real

implicit none

private

public :: t_precice_rbf

!> Connectivity, indices and weights of the structural nodes driving
! the surface points
type :: t_rbf_conn
  !> Indices (surface to structure)
  integer,  allocatable :: ind(:,:)
  !> Weights (surface to structure)
  real(wp), allocatable :: wei(:,:)
end type t_rbf_conn

!> RBF coupling structures
type :: t_precice_rbf

  !> Coupling nodes
  real(wp), allocatable :: nodes(:,:)

  !> Global position of the coupling nodes
  real(wp), allocatable :: rrb(:,:)

  !> Orientation of the coupling nodes
  real(wp), allocatable :: rrb_rot(:,:)

  !> Grid nodes connectivity
  type(t_rbf_conn) :: nod
  !> Grid virtual nodes connectivity
  type(t_rbf_conn) :: nod_virtual
  !> Elem centers connectivity
  type(t_rbf_conn) :: cen
  !> Elem virtual centers connectivity
  type(t_rbf_conn) :: cen_virtual
  !> Ac stripe connectivity (for vl corrected) 
  type(t_rbf_conn) :: ctr_pt 
  !> Generic point 
  type(t_rbf_conn) :: point

  !> --- Parameters of rbf interpolation ---
  !> Number of points for transferring the motion from the structure to the
  ! surface, through a weighted average
  integer :: n_wei = 2

  !> Order of the norm used for computing distance-based weights
  real(wp) :: w_order = 1.0_wp

  contains

  procedure, pass(this) :: build_connectivity

end type t_precice_rbf

!> --------------------------------------------------------------
contains
!----------------------------------------------------------------

!----------------------------------------------------------------
!> Read local coordinates of surface nodes, rr, and build data for
! structure to surface interpolation of the motion
subroutine build_connectivity(this, aero_coord, coupling_node_rot)
  class(t_precice_rbf), intent(inout) :: this
  real(wp),             intent(in)    :: aero_coord(:,:)
  real(wp),             intent(in)    :: coupling_node_rot(3,3)
  
  real(wp), allocatable               :: diff_all(:,:), diff_all_transpose(:,:) 
  real(wp), allocatable               :: dist_all(:), mat_dist_all(:,:), wei_v(:), Wnorm(:,:)
  integer , allocatable               :: ind_v(:)
  integer                             :: np, ns
  integer                             :: ip, is

  !> Number of surface points
  np = size(aero_coord,2)

  !> Number of coupling nodes of the structure
  ns = size(this%nodes,2)

  !> === Surface nodes ===
  allocate(this%point%ind(this%n_wei, np)); this%point%ind = 0
  allocate(this%point%wei(this%n_wei, np)); this%point%wei = 0.0_wp

  allocate(diff_all(3,ns)); diff_all = 0.0_wp
  allocate(diff_all_transpose(ns,3)); diff_all_transpose = 0.0_wp
  allocate(dist_all(ns)); dist_all = 0.0_wp
  allocate(mat_dist_all(ns,ns)); mat_dist_all = 0.0_wp
  allocate(Wnorm(3,3)); Wnorm = 0.0_wp

  !> anisotropy matrix: section is rigid chordwise
  Wnorm(1,1) = 1e-3_wp
  Wnorm(2,2) = 1e+0_wp
  Wnorm(3,3) = 1e-3_wp

  !> From beam ref. sys to Dust ref. sys
  Wnorm = matmul(transpose(coupling_node_rot),(matmul(Wnorm,coupling_node_rot)))
  
  do ip = 1, np

    !> Distance of the surface nodes from the structural nodes
    do is = 1, ns
      diff_all(:,is)  = aero_coord(:,ip) - this%nodes(:,is)
    end do
    !> interpolation matrix
    ! [ns x ns] =                       [ns x 3]        *     [3 x 3] *  [3 x ns]
    mat_dist_all(:,:)   =    matmul(transpose(diff_all), matmul(Wnorm , diff_all)) 
    
    do is = 1, ns
      dist_all(is) = sqrt(mat_dist_all(is,is))
    end do

    call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

    wei_v = 1.0_wp / max( wei_v, 1e-9_wp )**this%w_order
    wei_v = wei_v / sum(wei_v)

    this%point%wei(:,ip) = wei_v
    this%point%ind(:,ip) = ind_v

  enddo

  !> Deallocate and cleaning
  if ( allocated(dist_all) )  deallocate(dist_all)
  if ( allocated(wei_v   ) )  deallocate(wei_v   )
  if ( allocated(ind_v   ) )  deallocate(ind_v   ) 

end subroutine build_connectivity


end module mod_precice_rbf
