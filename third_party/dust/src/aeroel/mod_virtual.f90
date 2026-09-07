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
!!          Alessandro Cocco
!!          Federico Gentile
!!          Matteo Dall'Ora
!!=========================================================================


!> Module to treat surface doublet + source panels

module mod_virtual

    use mod_param, only: &
      wp , nl, &
      prev_tri , next_tri , &
      prev_qua , next_qua , &
      pi, max_char_len
    
    use mod_handling, only: &
      error, warning, printout
    
    use mod_sim_param, only: &
      t_sim_param, sim_param
    
    use mod_math, only: &
      cross , compute_qr

!----------------------------------------------------------------------
    
implicit none

public :: c_elem_virtual, t_elem_virtual_p

private
    
!----------------------------------------------------------------------

type :: t_elem_virtual_p
    class(c_elem_virtual), pointer :: p
end type

!----------------------------------------------------------------------

type :: c_elem_virtual

  !> Center of the element
  real(wp) :: cen(3)

  !> id of the component to which it belongs
  integer :: comp_id

  !> Id of the element, TODO consider removing this
  integer :: id

  !> Number of vertexes
  integer :: n_ver

  !> Vertexes coordinates
  real(wp), allocatable :: ver(:,:)

  !> Id of the vertexes
  integer, allocatable :: i_ver(:)

  !> Element area
  real(wp) :: area

  !> Element normal vector
  real(wp) :: nor(3)

  !> Element normal vector (at previous dt, to compute dn_dt)
  real(wp) :: nor_old(3)

  !> Element normal vector (at previous dt, to compute dn_dt)
  real(wp) :: dn_dt(3)

  !> Element tangential vectors
  real(wp) :: tang(3,2) ! tangent unit vectors as in PANAIR

  !> Vector of each edge
  real(wp), allocatable :: edge_vec(:,:)

  !> Length of each edge
  real(wp), allocatable :: edge_len(:)

  !> Unit vector of each edge
  real(wp), allocatable :: edge_uni(:,:)

  !> Body velocity at the centre
  real(wp) :: ub(3)

  !> Orientation vector of the center 
  real(wp) :: ori(3)

  !> Orientation matrix of the center 
  real(wp) :: R_cen(3,3)
  
  !> Is the element moving during simulation?
  logical  :: moving
  
  !> Hinge motion
  ! Initialization to zero *** to do: restart??? ***
  !> Delta position, due to hinge motion of the element center
  real(wp) :: dcen_h(3)     ! = 0.0_wp
  real(wp) :: dcen_h_old(3) ! = 0.0_wp
  !> Delta velocity, due to hinge motion of the element center, evaluated
  ! with finite difference method: dvel_h = ( dcen_h - dcen_h_old ) / dt
  real(wp) :: dvel_h(3) ! = 0.0_wp
  real(wp) :: loc_ctr_pt(3)

  contains
  procedure, pass(this) :: calc_geo_data_virtual
  
end type c_elem_virtual

!----------------------------------------------------------------------

contains 

!> Calculate the geometrical quantities of a virtual panel
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! virtual panel
subroutine calc_geo_data_virtual(this,vert)
  class(c_elem_virtual), intent(inout) :: this
  real(wp), intent(in) :: vert(:,:)

  integer :: nsides, is
  real(wp):: nor(3), tanl(3)

  this%ver = vert
  nsides = this%n_ver

  ! center
  this%cen =  sum ( this%ver,2 ) / real(nsides,wp)

  ! unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross(this%ver(:,3) - this%ver(:,1) , &
                this%ver(:,4) - this%ver(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross(this%ver(:,3) - this%ver(:,2) , &
                this%ver(:,1) - this%ver(:,2)     )
  end if

  this%area = 0.5_wp * norm2(nor)
  this%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

  this%tang(:,1) = tanl / norm2(tanl)
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )

  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_tri(is)) - this%ver(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_qua(is)) - this%ver(:,is)
    end do
  end if

  ! edge: edge_len(:)
  do is = 1 , nsides
    this%edge_len(is) = norm2(this%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    this%edge_uni(:,is) = this%edge_vec(:,is) / this%edge_len(is)
  end do

end subroutine calc_geo_data_virtual

!----------------------------------------------------------------------

end module mod_virtual