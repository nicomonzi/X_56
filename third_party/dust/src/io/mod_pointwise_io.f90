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

!> Module to read meshes and geometrical entities expressed as a series of
!! points
!! Disclamer: this has no association with the homonymous grid-generating
!! software

! list of subroutines and functions
! - read_mesh_pointwise
! - read_points
! - read_lines
! - set_parser_pointwise

module mod_pointwise_io

use mod_param, only: &
  wp, max_char_len, nl, pi

use mod_math, only: &
  geoseries, geoseries_both, geoseries_hinge

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, &
  countoption, getsuboption, getintarray , finalizeParameters

use mod_spline, only: &
  t_spline, hermite_spline , deallocate_spline

use mod_parametric_io, only: &
  define_section , define_division

!----------------------------------------------------------------------

implicit none

public :: read_mesh_pointwise , read_mesh_pointwise_ll

private

character(len=*), parameter :: this_mod_name = 'mod_pointwise_io'

!----------------------------------------------------------------------
!> module variables


!----------------------------------------------------------------------
!> point type
type :: t_point
  integer                 :: id
  real(wp)                :: coord(3)
  character(max_char_len) :: airfoil
  character(max_char_len) :: airfoil_table
  real(wp)                :: chord       !
  real(wp)                :: theta       ! deg
  character(max_char_len) :: sec_nor_str
  logical                 :: flip_sec
  real(wp)                :: sec_nor(3)
  real(wp) , allocatable  :: xy(:,:)
end type t_point

!> line type
type :: t_line
  character(max_char_len) :: l_type
  integer                 :: end_points(2)
  integer                 :: nelems
  real(wp)                :: tension
  real(wp)                :: bias
  character(max_char_len) :: type_span   ! discretisation in span
  real(wp)                :: r_ib        ! geometric series parameters
  real(wp)                :: r_ob
  real(wp)                :: y_ref 
  real(wp)                :: leng
  integer                 :: neigh_line(2) = 0
  real(wp), allocatable   :: t_vec1(:) , t_vec2(:)
end type t_line

! !> refline_pt type: type containing info of the points on the reference line
! type :: t_refline_pt
!   real(wp)   :: r(3)     ! coordinate
!   real(wp)   :: t(3)     ! tangent (unit) vector to ref line
!   !> interpolation from inputs
!   integer    :: id_pt(2) ! neighboring points
!   real(wp)   :: s        ! curvilinear coord s\in(0,1) between the nodes
!   real(wp)   :: n(3)     ! unit vector normal to the plane where the
! end type t_refline_pt

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
!> read pointwise input file and build ee, rr arrays of connectivity
!  and point coordinates
subroutine read_mesh_pointwise (mesh_file, ee, rr, y_fountain, nelem_chord, ref_chord_fraction, ElType, &
                                npoints_chord_tot , nelem_span_tot, &
                                airfoil_list_actual, i_airfoil_e, normalised_coord_e, & 
                                aero_table_out, thickness) 


  character(len=*), intent(in)                                      :: mesh_file
  character, intent(in)                                             :: ElType
  real(wp), intent(in)                                              :: ref_chord_fraction
  integer, intent(in)                                               :: nelem_chord 
  integer  , allocatable, intent(out)                               :: ee(:,:)
  real(wp) , allocatable, intent(out)                               :: rr(:,:)
  real(wp) , intent(out)                                            :: y_fountain
  integer,  allocatable, intent(out), optional                      :: i_airfoil_e(:,:)
  real(wp), allocatable, intent(out), optional                      :: normalised_coord_e(:,:)
  real(wp), allocatable, intent(out), optional                      :: thickness(:,:) 
  character(len=max_char_len), allocatable , intent(out), optional  :: airfoil_list_actual(:) 
  integer  , intent(out), optional                                  :: npoints_chord_tot, nelem_span_tot
  logical, intent(out), optional                                    :: aero_table_out

  !> Parser and sub-parsers
  type(t_parse)                                                     :: pmesh_prs
  type(t_parse) , pointer                                           :: point_prs , line_prs
  logical                                                           :: aero_table
  character(max_char_len)                                           :: type_chord
  
  real(wp), allocatable                                             :: chord_fraction(:)
  real(wp)                                                          :: thickness_section                             
  !> Point and line structures
  type(t_point) , allocatable                                       :: points(:)
  type(t_line ) , allocatable                                       :: lines(:)
  !> geo series mesh
  real(wp)                                                          :: r , r_le , r_te , r_le_fix 
  real(wp)                                                          :: r_le_moving , r_te_fix , r_te_moving 
  real(wp)                                                          :: x_refinement
  !> ee, rr size
  integer                                                           :: nelem_chord_tot, npoint_chord_tot, npoint_span_tot
  integer                                                           :: ee_size , rr_size
  integer                                                           :: i_ch , i_sp , i_point , i_elem

  real(wp) , allocatable                                            :: ref_line_points(:,:)
  real(wp) , allocatable                                            :: ref_line_normal(:,:)
  integer  , allocatable                                            :: ref_line_interp_p(:,:)
  real(wp) , allocatable                                            :: ref_line_interp_s(:)
  real(wp) , allocatable                                            :: ref_line_interp_s_all(:)
  real(wp) , allocatable                                            :: s_in(:) , nor_in(:,:)

  real(wp) , allocatable                                            :: xy1(:,:) , xy2(:,:) , xy(:,:)
  real(wp) , allocatable                                            :: rr_s(:,:)
  real(wp)                                                          :: twist_rad , theta

  real(wp)                                                          :: w1 , w2
  integer                                                           :: i , i1 , i2, j, k, kk
  integer                                                           :: iAirfoil, nAirfoils 
  real(wp) , allocatable                                            :: airfoil_list_actual_s(:)

  real(wp)                                                          :: s_cen_e
  character(len=*), parameter                                       :: this_sub_name = 'read_mesh_pointwise'
  
  aero_table = .false.

  !> Prepare parser and sub-parsers
  ! pass "p" for <panel> or <vortlat> elems that contain "airfoil" field,
  ! while <liftlin>s contain "airfoil_table" field
  call set_parser_pointwise( 'p' , pmesh_prs , point_prs , line_prs )

  !> Enable option reading
  call pmesh_prs%read_options(mesh_file,printout_val=.true.)

  type_chord  = getstr(pmesh_prs,'type_chord')

  !> geo series mesh 
  r = getreal(pmesh_prs,'r')
  r_le = getreal(pmesh_prs,'r_le')
  r_te = getreal(pmesh_prs,'r_te')
  r_le_fix = getreal(pmesh_prs,'r_le_fix')
  r_le_moving = getreal(pmesh_prs,'r_le_moving')
  r_te_fix = getreal(pmesh_prs,'r_te_fix')
  r_te_moving = getreal(pmesh_prs,'r_te_moving')
  x_refinement = getreal(pmesh_prs,'x_refinement')
  !> end geo series mesh

  !> y fountain 
  y_fountain = getreal(pmesh_prs,'y_fountain') 

  aero_table = getlogical(pmesh_prs,     'airfoil_table_correction')
  
  !> Read points and lines
  call read_points ( ElType, pmesh_prs , point_prs , points , aero_table)
  call read_lines  (       pmesh_prs ,  line_prs , lines  , nelem_span_tot )
  
  !> Set dimensions of ee, rr mat
  if (ElType.eq.'p') then ;  nelem_chord_tot = 2 * nelem_chord
  else                    ;  nelem_chord_tot =     nelem_chord
  endif

  npoint_chord_tot = nelem_chord_tot + 1
  npoint_span_tot  = nelem_span_tot  + 1

  ee_size =  nelem_chord_tot *  nelem_span_tot
  rr_size = npoint_chord_tot * npoint_span_tot

  ! === connectivity matrix, ee ===
  allocate( ee( 4 , ee_size ) ) ; ee = 0
  do i_ch = 1 , nelem_chord_tot
    do i_sp = 1 , nelem_span_tot
      i_elem  =  nelem_chord_tot*( i_sp - 1 ) + i_ch
      i_point = npoint_chord_tot*( i_sp - 1 ) + i_ch
      ee(1,i_elem) = i_point + npoint_chord_tot         ! iPoint
      ee(2,i_elem) = i_point                            ! iPoint + 1
      ee(3,i_elem) = i_point + 1                        ! iPoint + npoint_chord_tot + 1
      ee(4,i_elem) = i_point + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
    end do
  end do

  ! === define the reference line ===
  !> Check line connectivity, re-order lines and define tangent vectors
  call sort_lines( lines )

  call sort_points( points )

  !> check input consistency (todo: improve this subroutine)
  call check_point_line_inputs ( points , lines )

  call fill_line_tan_vec( points , lines )

  !>
  call build_reference_line(mesh_file, npoint_span_tot   , &
                            points , lines    , &
                            ref_line_points   , &
                            ref_line_normal   , &
                            ref_line_interp_p , &
                            ref_line_interp_s , &
                            ref_line_interp_s_all , &
                            s_in , nor_in )

  !> update ref_line_normal
  call update_ref_line_normal( points , ref_line_normal   , &
                                        ref_line_interp_p , &
                                        ref_line_interp_s , &
                                        ref_line_interp_s_all , &
                                        s_in , nor_in    )


  ! === define the coordinates of the sections at all the input points ===
  !> check that the first and last node has no interp attribute
  if ( trim(points( lines(1)%end_points(1) )%airfoil) .eq. 'interp' ) then
    write(*,*) ' error in read_mesh_pointwise: "airfoil" input &
                &of the first point cannot be "interp". Stop ' ; stop
  end if
  if ( trim(points( lines(size(lines))%end_points(2) )%airfoil) .eq. 'interp' ) then
    write(*,*) ' error in read_mesh_pointwise: "airfoil" input &
              &of the last point cannot be "interp". Stop ' ; stop
  end if 

  allocate(chord_fraction(nelem_chord+1))
  call define_division(type_chord, nelem_chord, chord_fraction, &
                      r, r_le, r_te, r_le_fix, r_le_moving, r_te_fix, r_te_moving, x_refinement)

  if (present(thickness)) then 
    allocate(thickness(2,size(points)))
    thickness = 0.0_wp
  endif
  
  !> first point
  call define_section(points(1)%chord , trim(adjustl(points(1)%airfoil)) , &
                      points(1)%theta , ElType , nelem_chord             , &
                      type_chord , chord_fraction , ref_chord_fraction   , &
                      (/ 0.0_wp , 0.0_wp , 0.0_wp /) , points(1)%xy, thickness_section )
  if (present(thickness)) thickness(1,1) = thickness_section 
  if ( points(1)%flip_sec ) call flip_section( ElType , points(1)%xy )

  !> last point
  i = size(points)
  call define_section( points(i)%chord , trim(adjustl(points(i)%airfoil)) , &
                      points(i)%theta , ElType , nelem_chord             , &
                      type_chord , chord_fraction , ref_chord_fraction   , &
                      (/ 0.0_wp , 0.0_wp , 0.0_wp /) , points(i)%xy, thickness_section )
  if (present(thickness)) thickness(1,i) = thickness_section 
  if ( points(i)%flip_sec ) call flip_section( ElType , points(i)%xy ) 

  do i = 2 , size(points)-1

    if ( trim(points(i)%airfoil) .ne. 'interp' ) then
      call define_section( points(i)%chord , trim(adjustl(points(i)%airfoil)) , &
                          points(i)%theta , ElType , nelem_chord             , &
                          type_chord , chord_fraction , ref_chord_fraction   , &
                          (/ 0.0_wp , 0.0_wp , 0.0_wp /) , points(i)%xy, thickness_section )
      if (present(thickness)) thickness(1,i) = thickness_section
      if ( points(i)%flip_sec ) call flip_section( ElType , points(i)%xy   )

    else ! points of the section must be interpolated

    !> find the input points and the weights to be used
    i1 = i-1 ; i2 = i+1
    do while ( trim(points(i1)%airfoil) .eq. 'interp' )
      i1 = i1-1
    end do
    do while ( trim(points(i2)%airfoil) .eq. 'interp' )
      i2 = i2+1
    end do

    w1 = ( s_in(i2) - s_in(i ) ) / ( s_in(i2) - s_in(i1) )
    w2 = ( s_in(i ) - s_in(i1) ) / ( s_in(i2) - s_in(i1) )

    if ( allocated(xy1) ) deallocate(xy1)
    if ( allocated(xy2) ) deallocate(xy2)

    !> compute the xy coordinates of the neighboring points w/o interp attribute
    !> interpolate the shape only:
    !> -----> unitary chord, theta = 0.0, ref_chord fraction = 0.0_wp
    call define_section( 1.0_wp , trim(adjustl(points(i1)%airfoil))  , &
                          0.0_wp , ElType , nelem_chord              , & 
                          type_chord , chord_fraction , 0.0_wp       , &
                          (/ 0.0_wp , 0.0_wp , 0.0_wp /) , xy1, thickness_section )
    call define_section( 1.0_wp , trim(adjustl(points(i2)%airfoil))  , &
                          0.0_wp , ElType , nelem_chord              , &
                          type_chord , chord_fraction , 0.0_wp       , &
                          (/ 0.0_wp , 0.0_wp , 0.0_wp /) , xy2, thickness_section )
    if ( points(i1)%flip_sec ) call flip_section( ElType , xy1 )
    if ( points(i2)%flip_sec ) call flip_section( ElType , xy2 )

    ! linear interpolation (weighted sum)
    if ( allocated(xy) ) deallocate(xy)
    allocate( xy( size(xy1,1) , size(xy1,2) ) )

    xy = xy1 * w1 + xy2 * w2

    ! transformations: 1. translation, 2. scaling, 3. rotation
    twist_rad = points(i)%theta * pi / 180.0_wp
    !> 1. translation
    xy(1,:) = xy(1,:) - ref_chord_fraction
    !> 2. scaling
    xy      = xy * points(i)%chord
    !> 3. rotation
    xy = matmul( reshape( (/ cos(twist_rad),-sin(twist_rad) , &
                              sin(twist_rad), cos(twist_rad) /) , (/2,2/) ) , xy )
    
    allocate(points(i)%xy(size(xy,1),size(xy,2))) ; points(i)%xy = xy

    end if

  end do

  ! === define the desired (output) points of the geometry ===
  allocate( rr( 3 , rr_size ) ) ; rr = 0.0_wp
  allocate( rr_s ( 2 , size(points(1)%xy,2) ) )  

  do i = 1 , npoint_span_tot

    rr_s = points( ref_line_interp_p(i,1) )%xy * &
                                          ( 1.0_wp - ref_line_interp_s(i) ) + &
            points( ref_line_interp_p(i,2) )%xy  *  ref_line_interp_s(i)

    i1 = 1 + ( i-1 ) * npoint_chord_tot
    i2 =       i     * npoint_chord_tot

    ! rotation
    theta = atan2( ref_line_normal(i,3) , ref_line_normal(i,2) )

    rr(1,i1:i2) = rr_s(1,:)
    rr(2,i1:i2) =-rr_s(2,:) * sin(theta)
    rr(3,i1:i2) = rr_s(2,:) * cos(theta) 

    ! traslation
    rr(1,i1:i2) = rr(1,i1:i2) + ref_line_points(i,1)
    rr(2,i1:i2) = rr(2,i1:i2) + ref_line_points(i,2)
    rr(3,i1:i2) = rr(3,i1:i2) + ref_line_points(i,3)

  end do


  deallocate( rr_s )
  if (ElType .eq. 'v' .and. aero_table) then
    !> airfoil_list_actual: save only the actual airfoil, not "interp" inputs
    nAirfoils = 0
    do i = 1 , size(points)
      if ( trim(points(i)%airfoil) .ne. 'interp' ) nAirfoils = nAirfoils + 1
    end do

    allocate(airfoil_list_actual(  nAirfoils))
    allocate(airfoil_list_actual_s(nAirfoils))

    iAirfoil = 0

    do i = 1 , size(points)
      if ( trim(points(i)%airfoil) .ne. 'interp' ) then

        iAirfoil  = iAirfoil + 1
        airfoil_list_actual(iAirfoil) = trim( points(i)%airfoil_table )
        airfoil_list_actual_s(iAirfoil) = s_in( i ) ! curvilinear coord. of a sec.
                                                    ! where an airfoil is defined

        !> check  if the file containing the .c81 table exists
        call check_file_exists(airfoil_list_actual(iAirfoil), this_sub_name, &
                              this_mod_name)
      end if
    end do
  

    !> i_airfoil_e, normalised_coord_e
    allocate(i_airfoil_e       (2,nelem_span_tot)) ; i_airfoil_e        = 0
    allocate(normalised_coord_e(2,nelem_span_tot)) ; normalised_coord_e = 0.0_wp
    ! Find normalised_coord_e for ll interpolation ---
    j = 1
    do i = 1 , nelem_span_tot ! loop over the elements in the spanwise direction

      s_cen_e = 0.5_wp * ( ref_line_interp_s_all(i) + ref_line_interp_s_all(i+1) )
    
      ! scan airfoil_list_actual_s and update j, index of the lower bound for interpolation
      if  ( s_cen_e .gt. airfoil_list_actual_s(j+1) ) then
        j = j + 1
        if ( j .eq. size(airfoil_list_actual_s) ) then
          write(*,*) ' error in read_mesh_pointwise: &
                        &out of bounds while scanning line sectons; stop ' ; stop
        endif
      end if

      ! check if the elem belongs to more than two sections of the reference line: if .t. -> error
      if ( ( ref_line_interp_s_all(i  ) .lt. airfoil_list_actual_s(j  ) ) .and. &
            ( ref_line_interp_s_all(i+1) .gt. airfoil_list_actual_s(j+1) ) ) then
            
        write(*,*) '  warning in read_mesh_pointwise: '
        write(*,*) '  element i =', i , ' belongs to more than 2 sections: '
        write(*,*) '  centre of the el.            : ' , s_cen_e
        write(*,*) '  ref_line_interp_s_all(',i,':',i+1,') : ' , ref_line_interp_s_all(i:i+1)
        write(*,*) '  airfoil_list_actual_s(',j,':',j+1,') : ' , airfoil_list_actual_s(j:j+1)
      end if

      ! i_airfoil_e and normalised_coord_e (of the edges of the element) for ll interpolation
      i_airfoil_e(1,i) = j
      i_airfoil_e(2,i) = j + 1
      
      k = 1
      kk = size(airfoil_list_actual_s)
      if ( kk .le. k ) then
        write(*,*) ' error in read_mesh_pointwise: invalid airfoil table bounds. stop ' ; stop
      end if
      normalised_coord_e(1,i) = ( ref_line_interp_s_all(i)   - airfoil_list_actual_s(k) ) / &
                                ( airfoil_list_actual_s(kk) - airfoil_list_actual_s(k) )
      normalised_coord_e(2,i) = ( ref_line_interp_s_all(i+1) - airfoil_list_actual_s(k) ) / &
                                ( airfoil_list_actual_s(kk) - airfoil_list_actual_s(k) ) 
    end do
  end if
  ! optional output ----
  npoints_chord_tot = npoint_chord_tot
  ! optional output ----
  if (present(aero_table_out)) then
    aero_table_out = aero_table
  end if

  call finalizeParameters(pmesh_prs)


end subroutine read_mesh_pointwise

!----------------------------------------------------------------------
!> read pointwise mesh of ll elements.
! The implementation of this routine may appear not so clean, as a result
! of paste-and-copy of some sections of the routines:
! - read_mesh_pointwise -> build_reference_line
! - read_mesh_ll (io/mod_ll_io) -> airfoil_table as Input,
!                                  i_airtoil_e ,  normalised_coord_e ,
!                                  chord_p , theta_p
subroutine read_mesh_pointwise_ll(mesh_file, ee, rr, y_fountain, ElType, &
                        airfoil_list_actual, nelem_span_list, &
                        i_airfoil_e , normalised_coord_e    , &
                        npoints_chord_tot , nelem_span_tot  , &
                        chord_p,theta_p,theta_e)

  character(len=*), intent(in)                              :: mesh_file
  character, intent(in)                                     :: ElType
  integer  , allocatable, intent(out)                       :: ee(:,:)
  real(wp) , allocatable, intent(out)                       :: rr(:,:)
  character(len=max_char_len), allocatable , intent(out)    :: airfoil_list_actual(:)
  integer  , allocatable, intent(out)                       :: nelem_span_list(:)
  integer  , allocatable, intent(out)                       :: i_airfoil_e(:,:)
  real(wp) , allocatable, intent(out)                       :: normalised_coord_e(:,:)
  integer  ,              intent(out)                       :: npoints_chord_tot, nelem_span_tot
  real(wp) , allocatable, intent(out)                       :: chord_p(:),theta_p(:),theta_e(:)
  real(wp) ,              intent(out)                       :: y_fountain

  real(wp)                                                  :: s_cen_e 
  logical                                                   :: aero_table 
  ! parser and sub-parsers
  type(t_parse)                                             :: pmesh_prs
  type(t_parse), pointer                                    :: point_prs, line_prs
  !>
  integer                                                   :: nelem_chord
  logical                                                   :: mesh_flat 

  !> point and line structures
  type(t_point) , allocatable                               :: points(:)
  type(t_line ) , allocatable                               :: lines(:)

  integer                                                   :: nelem_chord_tot, npoint_chord_tot, npoint_span_tot
  integer                                                   :: ee_size, rr_size
  integer                                                   :: i_ch, i_sp, i_point, i_elem

  !>
  real(wp), allocatable                                     :: ref_line_points(:,:)
  real(wp), allocatable                                     :: ref_line_normal(:,:)
  integer , allocatable                                     :: ref_line_interp_p(:,:)
  real(wp), allocatable                                     :: ref_line_interp_s(:)
  real(wp), allocatable                                     :: ref_line_interp_s_all(:)
  real(wp), allocatable                                     :: s_in(:), nor_in(:,:)
  real(wp), allocatable                                     :: airfoil_list_actual_s(:)

  real(wp), allocatable                                     :: rr_s(:,:)
  real(wp)                                                  :: twist_rad, theta

  !>
  integer                                                   :: nAirfoils, iAirfoil
  integer                                                   :: i, i1, i2, j, k, kk

  character(len=*), parameter                               :: this_sub_name = 'read_mesh_pointwise_ll'

  aero_table = .false. 
  ! === Reference line, by points and lines as in read_mesh_pointwise ===
  !> Prepare parser and sub-parsers
  ! pass "p" for <panel> or <vortlat> elems that contain "airfoil" field,
  ! while <liftlin>s contain "airfoil_table" field
  call set_parser_pointwise( 'l' , pmesh_prs , point_prs , line_prs )

  !> Enable option reading
  call pmesh_prs%read_options(mesh_file,printout_val=.true.)

  nelem_chord     = 1 ! nelem_chord = 1 and no type_chord for liftlin
  nelem_chord_tot = 1
  if ( eltype .ne. 'l' ) then
    write(*,*) ' Error in mod_pointwise_io: read_mesh_pointwise_ll. '
    write(*,*) ' expected ElType = "l", but got: ', Eltype
    write(*,*) ' Stop. ' ; stop
  end if

  y_fountain = getreal(pmesh_prs, 'y_fountain') 

  !> Option for flat meshes: panels are not rotated, but normals and
  !! tangent vectors are rotated with the twist
  mesh_flat = getlogical(pmesh_prs,'mesh_flat')

  if ( mesh_flat .and. trim(ElType) .ne. 'l' ) then
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
        &flat mesh option is available only for lifting line elements.')
  end if 

  !> Read points and lines
  call read_points ( 'l' , pmesh_prs , point_prs , points , aero_table)
  call read_lines  (       pmesh_prs ,  line_prs , lines  , nelem_span_tot )

  npoint_chord_tot = nelem_chord_tot + 1
  npoint_span_tot  = nelem_span_tot  + 1

  !> === build ee, rr arrays (and all the remaining variables) ===
  do i = 1 , size(points)

  end do
  ee_size = nelem_span_tot
  rr_size = npoint_span_tot * npoint_chord_tot

  ! Directly build connectivity matrix
  allocate(ee(4,ee_size)) ; ee = 0
  do i_ch = 1,nelem_chord_tot
    do i_sp = 1,nelem_span_tot
      i_elem  =  nelem_chord_tot*(i_sp-1) + i_ch
      i_point = npoint_chord_tot*(i_sp-1) + i_ch
      ee(1,i_elem) = i_point + npoint_chord_tot         ! iPoint
      ee(2,i_elem) = i_point                            ! iPoint + 1
      ee(3,i_elem) = i_point + 1                        ! iPoint + npoint_chord_tot + 1
      ee(4,i_elem) = i_point + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
    enddo
  enddo

  ! === define the reference line ===
  !> Check line connectivity, re-order lines and define tangent vectors
  call sort_lines( lines )

  call sort_points( points )

  !> check input consistency (todo: improve this subroutine)
  call check_point_line_inputs ( points , lines )

  !>
  call fill_line_tan_vec( points , lines )

  !>
  call build_reference_line(mesh_file, npoint_span_tot   , &
                            points , lines    , &
                            ref_line_points   , &
                            ref_line_normal   , &
                            ref_line_interp_p , &
                            ref_line_interp_s , &
                            ref_line_interp_s_all , &
                            s_in , nor_in )

  !> update ref_line_normal
  call update_ref_line_normal( points , ref_line_normal   , &
                              ref_line_interp_p , &
                              ref_line_interp_s , &
                              ref_line_interp_s_all , &
                              s_in , nor_in    )

  !> output ~ read_mesh_ll(): nelem_span_list
  allocate(nelem_span_list(size(lines)))
  do i = 1 , size(lines)
    nelem_span_list(i) = lines(i)%nelems
  end do

  allocate(rr   (3,rr_size)) ; rr      = 0.0_wp
  allocate(chord_p(npoint_span_tot)) ; chord_p = 0.0_wp
  allocate(theta_p(npoint_span_tot)) ; theta_p = 0.0_wp
  allocate(theta_e(nelem_span_tot )) ; theta_e = 0.0_wp  


  ! === define the coordinates of the sections at all the input points ===
  !> check that the first and last node has no interp attribute
  if ( trim(points( lines(1)%end_points(1) )%airfoil) .eq. 'interp' ) then
    write(*,*) ' error in read_mesh_pointwise: "airfoil" input &
                &of the first point cannot be "interp". Stop ' ; stop
  end if
  if ( trim(points( lines(size(lines))%end_points(2) )%airfoil) .eq. 'interp' ) then
    write(*,*) ' error in read_mesh_pointwise: "airfoil" input &
                &of the last point cannot be "interp". Stop ' ; stop
  end if


  ! === fields on the input points ===
  do i = 1 , size(points)

    twist_rad = points(i) % theta * pi / 180.0_wp
    allocate( points(i)%xy ( 2 , 2 ) )
    ! Rotate the section around the reference line with the twist angle
    !
    ! For flat meshes the section is not rotated, but the
    ! normal/tangent vector are

    points(i)%xy(1,1) = 0.0_wp
    points(i)%xy(2,1) = 0.0_wp
    if ( mesh_flat ) then
      points(i)%xy(1,2) = points(i)%chord
      points(i)%xy(2,2) = 0.0_wp
    else
      points(i)%xy(1,2) = cos(twist_rad)*points(i)%chord
      points(i)%xy(2,2) = -sin(twist_rad)*points(i)%chord
    endif

    if ( points(i)%flip_sec ) points(i)%xy(2,:) = -points(i)%xy(2,:)

  end do


  ! === fields on the output points ===
  allocate(rr_s(2,2)) ; rr_s = 0.0_wp
  do i = 1 , npoint_span_tot

    !> chord
    chord_p(i) = points( ref_line_interp_p(i,1) )%chord * &
                                              ( 1.0_wp - ref_line_interp_s(i) ) + &
                 points( ref_line_interp_p(i,2) )%chord * ref_line_interp_s(i)
    !> theta [rad]
    theta_p(i) = points( ref_line_interp_p(i,1) )%theta * &
                                              ( 1.0_wp - ref_line_interp_s(i) ) + &
                 points( ref_line_interp_p(i,2) )%theta * ref_line_interp_s(i)
    theta_p(i) = theta_p(i) * pi/180.0_wp

    !> rr
    rr_s = points( ref_line_interp_p(i,1) )%xy * &
                                        ( 1.0_wp - ref_line_interp_s(i) ) + &
           points( ref_line_interp_p(i,2) )%xy  *  ref_line_interp_s(i)

    i1 = 1 + ( i-1 ) * npoint_chord_tot
    i2 =       i     * npoint_chord_tot

    ! rotation
    theta = atan2( ref_line_normal(i,3) , ref_line_normal(i,2) )

    rr(1,i1:i2) = rr_s(1,:)
    rr(2,i1:i2) =-rr_s(2,:) * sin(theta)
    rr(3,i1:i2) = rr_s(2,:) * cos(theta)

    ! traslation
    rr(1,i1:i2) = rr(1,i1:i2) + ref_line_points(i,1)
    rr(2,i1:i2) = rr(2,i1:i2) + ref_line_points(i,2)
    rr(3,i1:i2) = rr(3,i1:i2) + ref_line_points(i,3)

  end do

  ! For flat meshes the section is not rotated, but the
  ! normal/tangent vector are

  if ( mesh_flat ) then
    do i = 1,nelem_span_tot
      theta_e(i) = 0.5_wp*(theta_p(i)+theta_p(i+1))
    enddo
  endif

  ! === from the output of the read_mesh_pointwise-like section of the  ===
  ! === routine to the desired format (~ read_mesh_ll)                  ===
  ! ref_line_interp_p, red_line_interp_s -> i_airfoil_e , normalised_coord_e
  ! ...

  !> airfoil_list_actual: save only the actual airfoil, not "interp" inputs
  nAirfoils = 0
  do i = 1 , size(points)
    if ( trim(points(i)%airfoil) .ne. 'interp' ) nAirfoils = nAirfoils + 1
  end do

  allocate(airfoil_list_actual(  nAirfoils)) 
  allocate(airfoil_list_actual_s(nAirfoils))

  iAirfoil = 0

  do i = 1 , size(points)
    if ( trim(points(i)%airfoil) .ne. 'interp' ) then

      iAirfoil  = iAirfoil + 1
      airfoil_list_actual(iAirfoil) = trim( points(i)%airfoil )
      airfoil_list_actual_s(iAirfoil) = s_in( i ) ! curvilinear coord. of a sec.
                                                  ! where an airfoil is defined
      !> check  if the file containing the .c81 table exists
      call check_file_exists(airfoil_list_actual(iAirfoil), this_sub_name, &
                            this_mod_name)
    end if
  end do


  !> i_airfoil_e, normalised_coord_e
  allocate(i_airfoil_e       (2,nelem_span_tot)) ; i_airfoil_e        = 0
  allocate(normalised_coord_e(2,nelem_span_tot)) ; normalised_coord_e = 0.0_wp

  ! Find normalised_coord_e for ll interpolation ---
  j = 1
  do i = 1 , nelem_span_tot ! loop over the elements in the spanwise direction

    s_cen_e = 0.5_wp * ( ref_line_interp_s_all(i) + ref_line_interp_s_all(i+1) )

    ! scan airfoil_list_actual_s and update j, index of the lower bound for interpolation
    if  ( s_cen_e .gt. airfoil_list_actual_s(j+1) ) then
      j = j + 1
      if ( j .eq. size(airfoil_list_actual_s) ) then
        write(*,*) ' error in read_mesh_pointwise_ll: &
                      &out of bounds while scanning line sectons; stop ' ; stop
      endif
    end if

    ! check if the elem belongs to more than two sections of the reference line: if .t. -> error
    if ( ( ref_line_interp_s_all(i  ) .lt. airfoil_list_actual_s(j  ) ) .and. &
          ( ref_line_interp_s_all(i+1) .gt. airfoil_list_actual_s(j+1) ) ) then

      write(*,*) ' error in read_mesh_pointwise_ll: '
      write(*,*) '  element i =', i , ' belongs to more than 2 sections: '
      write(*,*) '  centre of the el.            : ' , s_cen_e
      write(*,*) '  ref_line_interp_s_all(',i,':',i+1,') : ' , ref_line_interp_s_all(i:i+1)
      write(*,*) '  airfoil_list_actual_s(',j,':',j+1,') : ' , airfoil_list_actual_s(j:j+1)
    end if

    ! i_airfoil_e and normalised_coord_e (of the edges of the element) for ll interpolation
    i_airfoil_e(1,i) = j
    i_airfoil_e(2,i) = j + 1 
    k = 1
    kk = size(airfoil_list_actual_s)
    if ( kk .le. k ) then
      write(*,*) ' error in read_mesh_pointwise_ll: invalid airfoil table bounds. stop ' ; stop
    end if
    normalised_coord_e(1,i) = ( ref_line_interp_s_all(i)   - airfoil_list_actual_s(k) ) / &
                                ( airfoil_list_actual_s(kk) - airfoil_list_actual_s(k) )
    normalised_coord_e(2,i) = ( ref_line_interp_s_all(i+1) - airfoil_list_actual_s(k) ) / &
                                ( airfoil_list_actual_s(kk) - airfoil_list_actual_s(k) ) 

  end do

  ! optional output ----
  npoints_chord_tot = npoint_chord_tot

  call finalizeParameters(pmesh_prs)

end subroutine read_mesh_pointwise_ll

!----------------------------------------------------------------------
!> build reference line
subroutine build_reference_line(mesh_file, npoint_span_tot   , points, lines     , &
                                ref_line_points   , ref_line_normal   , &
                                ref_line_interp_p , ref_line_interp_s , &
                                ref_line_interp_s_all ,                 &
                                s_in , nor_in )

  integer                     , intent(in)    :: npoint_span_tot
  character(len=*), intent(in)                :: mesh_file
  type(t_point)               , intent(inout) :: points(:)
  type(t_line )               , intent(inout) :: lines( :)
  real(wp)      , allocatable , intent(out)   :: ref_line_points(:,:)
  real(wp)      , allocatable , intent(out)   :: ref_line_normal(:,:)
  integer       , allocatable , intent(out)   :: ref_line_interp_p(:,:)
  real(wp)      , allocatable , intent(out)   :: ref_line_interp_s(:)
  real(wp)      , allocatable , intent(out)   :: ref_line_interp_s_all(:)
  real(wp)      , allocatable , intent(out)   :: s_in(:)
  real(wp)      , allocatable , intent(out)   :: nor_in(:,:)

  real(wp)      , allocatable                 :: ref_line_spline_s(:)

  integer       , allocatable                 :: ip(:,:)
  real(wp)      , allocatable                 :: s_in_1(:) , nor_in_1(:,:)

  type(t_spline) :: spl

  integer :: i , i1 , i2 , j , n

  allocate(ref_line_points(  npoint_span_tot,3))
  allocate(ref_line_normal(  npoint_span_tot,3))
  allocate(ref_line_interp_p(npoint_span_tot,2))
  allocate(ref_line_spline_s(npoint_span_tot  ))     ; ref_line_spline_s     = 0.0_wp
  allocate(ref_line_interp_s(npoint_span_tot  ))     ; ref_line_interp_s     = 0.0_wp
  allocate(ref_line_interp_s_all(npoint_span_tot  )) ; ref_line_interp_s_all = 0.0_wp

  allocate(  s_in(size(points)  )) ;   s_in = 0.0_wp
  allocate(nor_in(size(points),3)) ; nor_in = 0.0_wp


  !> Starting point
  ref_line_points(1,:) = points( lines(1)%end_points(1) ) % coord
  i2 = 1
  do i = 1 , size(lines)
    ! end points
    i1 = i2
    i2 = i1 + lines(i)%nelems
    ref_line_points(i2,:) = points( lines(i)%end_points(2) ) % coord

    if (trim(lines(i)%l_type) .eq. 'straight') then

      !> points for interpolation
      do j = i1 , i2
        ref_line_interp_p(j,:) = lines(i)%end_points(:)
      end do

      call straight_line( mesh_file                              , &
                          points( lines(i)%end_points(1) )%coord , &
                          points( lines(i)%end_points(2) )%coord , &
                          lines(i)%nelems                        , &
                          lines(i)%type_span                     , &
                          lines(i)%r_ib                          , &
                          lines(i)%r_ob                          , &
                          lines(i)%y_ref                         , &
                          ref_line_points(i1:i2,:)               , &
                          ref_line_normal(i1:i2,:)               , &
                          ref_line_interp_s(i1:i2)               , &
                          lines(i)%leng  )

      !> s_all
      if ( i .eq. 1 ) then
        ref_line_interp_s_all(i1:i2) = ref_line_interp_s(i1:i2) * lines(i)%leng
      else
        ref_line_interp_s_all(i1:i2) = ref_line_interp_s(i1:i2) * lines(i)%leng + &
            ref_line_interp_s_all(i1)
      end if

      !> s
      if ( allocated(s_in_1) ) deallocate(s_in_1) ; allocate(s_in_1(2))
      s_in_1 = (/ 0.0_wp , 1.0_wp /)
      s_in(lines(i)%end_points(2)) = &
              s_in(lines(i)%end_points(1)) + s_in_1(2) * lines(i)%leng 
      
      !> nor
      if ( allocated(nor_in_1) ) deallocate(nor_in_1) ; allocate(nor_in_1(2,3))
      nor_in_1(1,:) = ref_line_normal(i1,:)
      nor_in_1(2,:) = ref_line_normal(i1,:)
      nor_in_1(1,:) = nor_in_1(1,:) / norm2(nor_in_1(1,:))
      nor_in_1(2,:) = nor_in_1(2,:) / norm2(nor_in_1(2,:))

      nor_in(lines(i)%end_points(1),:) = nor_in_1(1,:)
      nor_in(lines(i)%end_points(2),:) = nor_in_1(2,:)

      deallocate(s_in_1 , nor_in_1)


    else if ( trim(lines(i)%l_type) .eq. 'spline'   ) then

      !> build t_spline structure
      n = lines(i)%end_points(2) - lines(i)%end_points(1) + 1
      allocate(spl%rr( n , 3 ))
      allocate(spl%d0( 3 )) ; spl%d0 = lines(i)%t_vec1
      allocate(spl%d1( 3 )) ; spl%d1 = lines(i)%t_vec2
      spl%e_bc = (/ 'derivative' , 'derivative' /)

      do j = 1 , n
        spl%rr( j , : ) = points( lines(i)%end_points(1)-1+j )%coord
      end do

      if ( allocated(  ip)     ) deallocate(  ip    ) ; allocate(  ip(i2-i1+1,2))
      if ( allocated(  s_in_1) ) deallocate(  s_in_1) ; allocate(  s_in_1(n))
      if ( allocated(nor_in_1) ) deallocate(nor_in_1) ; allocate(nor_in_1(n,3))

      !> compute ref_line_points on the spline
      call hermite_spline( mesh_file, spl, lines(i)%nelems           , &
                                lines(i)%tension          , &
                                lines(i)%bias             , &
                                lines(i)%type_span        , &
                                lines(i)%r_ib             , &
                                lines(i)%r_ob             , &
                                lines(i)%y_ref            , &
                                ref_line_points(i1:i2,:)  , &
                                ref_line_normal(i1:i2,:)  , &
                                ip                        , &
                                ref_line_interp_s(i1:i2)  , &
                                ref_line_spline_s(i1:i2)  , &
                                lines(i)%leng             , &
                                s_in_1 , nor_in_1 )
      !> s_all
      if ( i .eq. 1 ) then
        ref_line_interp_s_all(i1:i2) = ref_line_spline_s(i1:i2) * lines(i)%leng
      else
        ref_line_interp_s_all(i1:i2) = ref_line_spline_s(i1:i2) * lines(i)%leng + &
            ref_line_interp_s_all(i1)
      end if

      !> s
      s_in( lines(i)%end_points(1): lines(i)%end_points(2) ) = &
                      s_in_1 + s_in(lines(i)%end_points(1))

      !> nor
      nor_in( lines(i)%end_points(1):lines(i)%end_points(2),:) = nor_in_1

      !> from ip to ref_line_interp_p
      do j = 1 , i2-i1+1
          ref_line_interp_p(i1+j-1,:) = lines(i)%end_points(1) - 1 + ip(j,:)
      end do

      if ( allocated(       ip ) ) deallocate(       ip )
      if ( allocated(   s_in_1 ) ) deallocate(   s_in_1 )
      if ( allocated( nor_in_1 ) ) deallocate( nor_in_1 )

      call deallocate_spline( spl )

    else
      write(*,*) ' error in build_reference_line(): '
      write(*,*) ' lines(',i,')%l_type = ' , lines(i)%l_type
      write(*,*) ' but l_type must be either spline or straight.'
      write(*,*) ' Stop. '; stop
    end if

  end do


end subroutine build_reference_line

!----------------------------------------------------------------------
!> straight_line subdivision
subroutine straight_line( mesh_file, r1, r2, nelems, type_span, r_ib, r_ob, y_ref, &  
                          rr, nor, s, leng )
  character(len=*)        , intent(in)    :: mesh_file
  real(wp)                , intent(in)    :: r1(3) , r2(3)
  integer                 , intent(in)    :: nelems
  character(max_char_len) , intent(in)    :: type_span
  real(wp)                , intent(in)    :: r_ib, r_ob, y_ref
  real(wp)                , intent(inout) :: rr(:,:)
  real(wp)                , intent(inout) ::nor(:,:)
  real(wp)                , intent(inout) ::  s(:)
  real(wp)                , intent(out)   :: leng
  real(wp), allocatable                   :: division(:), divisionIB(:), divisionOB(:) 
  character(len=*), parameter             :: this_sub_name = 'straight_line'
  real(wp) :: nor_v(3)

  integer :: i


  nor_v = r2 - r1
  leng = norm2(nor_v)
  nor_v = nor_v / leng
  
  do i = 1, nelems + 1

    !> rr
    if ( trim(type_span) .eq. 'uniform' ) then ! uniform spacing
      rr(i,:) = r1 + real(i - 1,wp)/real(nelems,wp) * ( r2 - r1 )
    elseif ( trim(type_span) .eq. 'cosine' ) then ! cosine  spacing in span
      rr(i,:) = 0.5_wp * ( r1 + r2 ) - &
                0.5_wp * ( r2 - r1 ) * &
                cos( real(i-1,wp)*pi/ real(nelems,wp) )
    elseif ( trim(type_span) .eq. 'cosine_ob' ) then ! cosine  spacing in span: outboard refinement
      rr(i,:) = r1 + &
                ( r2 - r1 ) * &
                sin( 0.5_wp*real(i-1,wp)*pi/ real(nelems,wp) )
    elseif ( trim(type_span) .eq. 'cosine_ib' ) then ! cosine  spacing in span: inboard refinement
      rr(i,:) = r2 - &
                ( r2 - r1 ) * &
                cos( 0.5_wp*real(i-1,wp)*pi/ real(nelems,wp) )
    elseif ( trim(type_span) .eq. 'equalarea' ) then ! equalarea spacing in span 
      rr(i,2) = sqrt(real(i-1,wp)/real(nelems,wp))*(r2(2) - r1(2)) + r1(2) 
      rr(i,1) = r1(1) + (rr(i,2) - r1(2))*(r2(1) - r1(1))/(r2(2) - r1(2))
      rr(i,3) = r1(3) + (rr(i,2) - r1(2))*(r2(3) - r1(3))/(r2(2) - r1(2))
    elseif ( trim(type_span) .eq. 'geoseries') then ! geometric series spacing in span
      call geoseries_both(r1(2), r2(2), nelems, y_ref, r_ib, r_ob, division) 
      rr(i,2) = division(i) 
      rr(i,1) = r1(1) + (rr(i,2) - r1(2))*(r2(1) - r1(1))/(r2(2) - r1(2))
      rr(i,3) = r1(3) + (rr(i,2) - r1(2))*(r2(3) - r1(3))/(r2(2) - r1(2))
    elseif ( trim(type_span) .eq. 'geoseries_ob') then ! geometric series spacing in span: outboard refinement
      call geoseries(r1(2), r2(2), nelems, r_ob, divisionIB, divisionOB) 
      rr(i,2) = divisionOB(i)  
      rr(i,1) = r1(1) + (rr(i,2) - r1(2))*(r2(1) - r1(1))/(r2(2) - r1(2))
      rr(i,3) = r1(3) + (rr(i,2) - r1(2))*(r2(3) - r1(3))/(r2(2) - r1(2))  
    elseif ( trim(type_span) .eq. 'geoseries_ib') then ! geometric series spacing in span: inboard refinement
      call geoseries(r1(2), r2(2), nelems, r_ib, divisionIB, divisionOB)  
      rr(i,2) = divisionIB(i) 
      rr(i,1) = r1(1) + (rr(i,2) - r1(2))*(r2(1) - r1(1))/(r2(2) - r1(2))
      rr(i,3) = r1(3) + (rr(i,2) - r1(2))*(r2(3) - r1(3))/(r2(2) - r1(2))
    else
      write(*,*) ' Mesh file   : ' , trim(mesh_file)
      write(*,*) ' type_span   : ' , trim(type_span)
      call error(this_sub_name, this_mod_name, 'Incorrect input: &
            & type_span must be equal to uniform, cosine, cosine_ib, cosine_ob, equalarea, geoseries, geoseries_ib, geoseries_ob.')
    end if

    !> nor
    nor(i,:) = nor_v

    !> curvilinear coordinate, s
    s(i) = norm2(rr(i,:)-r1) / norm2(r2-r1)

  end do


end subroutine straight_line

!----------------------------------------------------------------------
!> sort points:
subroutine sort_points( points )
  type(t_point) , intent(inout) :: points(:)
  type(t_point) , allocatable   :: points_tmp(:)

  integer :: n_points
  integer :: i , j

  n_points = size(points)

  allocate(points_tmp(n_points)) ; points_tmp = points

  do i = 1 , n_points
    do j = 1 , n_points
      if ( points_tmp(j)%id .eq. i ) then
        points( i ) = points_tmp( j )
      end if
    end do
  end do

  deallocate(points_tmp)

end subroutine sort_points

!----------------------------------------------------------------------
!> sort lines: first line has the first point with id = 1
!  and then chain last-to-first points until the last line
subroutine sort_lines( lines )
  type(t_line ) , intent(inout) :: lines( :)
  type(t_line ) , allocatable :: lines_tmp(:)
  integer :: n_lines
  integer :: i , i1 , i2 , i2_old , il , n_1

  n_lines = size(lines)

  allocate(lines_tmp(n_lines))
  lines_tmp = lines

  !> some preliminary checks:
  ! -> only one line has end_points(1) = 1 -> first line
  n_1 = 0
  do i = 1 , n_lines
    if ( lines_tmp(i)%end_points(1) .eq. 1 ) then
      n_1 = n_1 + 1 ; i1 = i
    end if
  end do
  if ( n_1 .ne. 1 ) then
    write(*,*) ' error in sort_lines: only 1 line must have'
    write(*,*) ' end_points(1) = 1. Stop ' ; stop
  end if
  lines(1) = lines_tmp(i1)
  ! -> no multiple start/end points
  ! ...

  !> find the chain of lines
  i1 = lines(1)%end_points(1) ; i2 = lines(1)%end_points(2)
  do i = 2 , n_lines

    il = 0 ; i2_old = i2 ! re-initialisation

    do while ( ( i2 .eq. i2_old ) .and. ( il .lt. n_lines ) )
      il = il + 1
      if (      lines_tmp(il)%end_points(1) .eq. i2_old ) then
        i1 =    lines_tmp(il)%end_points(1)
        i2 =    lines_tmp(il)%end_points(2)
      end if

    end do
    !>
    if ( i2 .eq. i2_old ) then
      write(*,*) ' error in sort_lines: broken chain. stop ' ; stop
    end if

    lines(i) = lines_tmp(il)

  end do


end subroutine sort_lines

!----------------------------------------------------------------------
!> flip section
subroutine flip_section ( eltype , xy )
  character, intent(in)    :: eltype
  real(wp) , intent(inout) :: xy(:,:)
  real(wp) , allocatable   :: xy_tmp(:,:)

  integer :: i

  allocate( xy_tmp(size(xy,1),size(xy,2)) )

  if ( eltype .eq. 'p' ) then

    do i = 1 , size(xy,2)
      xy_tmp(1,i) =  xy(1,size(xy,2)-i+1)
      xy_tmp(2,i) = -xy(2,size(xy,2)-i+1)
    enddo

  else if ( eltype .eq. 'v' ) then

    do i = 1 , size(xy,2)
      xy_tmp(1,i) =  xy(1,i)
      xy_tmp(2,i) = -xy(2,i)
    enddo

  else
    write(*,*) ' Error in flip_section. Eltype must be either &
                &"p" or "v". Stop.' ; stop
  end if

  xy = xy_tmp

  deallocate( xy_tmp )

end subroutine flip_section

!----------------------------------------------------------------------
!> check input consistency
subroutine check_point_line_inputs( points , lines )
  type(t_point) , intent(inout) :: points(:)
  type(t_line ) , intent(inout) :: lines( :)

  integer :: i , n_lines

  n_lines = size(lines)

  !> all the points are present
  if ( size(points) .lt. lines(size(lines))%end_points(2) ) then
    write(*,*) ' n. points .lt. last line end point id. stop ' ; stop
  end if

  !> if first or last lines are Splines, the user must provide the
  !  initial or final direction
  if ( ( trim(lines(1)%l_type) .eq. 'spline' ) .and. &
      ( .not. allocated(lines(1)%t_vec1) ) ) then
    write(*,*) ' error in check_point_lines_input: the first line is a Spline and &
                &has no TangentVec1 provided as an input. Stop ' ; stop
  end if

  if ( ( trim(lines(n_lines)%l_type) .eq. 'spline' ) .and. &
      ( .not. allocated(lines(n_lines)%t_vec2) ) ) then
    write(*,*) ' error in check_point_lines_input: the last line is a spline and &
                &has no TangentVec2 provided as an input. Stop. ' ; stop
  end if

  !> no consecutive splines are allowed so far
  do i = 1 , n_lines-1
    if ( ( trim(lines(i  )%l_type) .eq. 'spline' ) .and. &
        ( trim(lines(i+1)%l_type) .eq. 'spline' ) ) then
      write(*,*) ' error in check_point_lines_input: two consecutive splines are not &
                  &allowed so far. '
      write(*,*) ' In future release, this could not be true anymore: '
      write(*,*) ' two consecutive splines could be joined together if no'
      write(*,*) ' TangentVec is provided, or two consecutive spline with'
      write(*,*) ' TangentVec provided could be treated as well. So far, Stop.'
      stop
    end if
  end do


end subroutine check_point_line_inputs

!----------------------------------------------------------------------
!> update_ref_line_nor
subroutine update_ref_line_normal( points , ref_line_normal   , &
                                            ref_line_interp_p , &
                                            ref_line_interp_s , &
                                            ref_line_interp_s_all , &
                                            s_in , nor_in    )
  type(t_point) , intent(inout) :: points(:)
  integer       , intent(in)    :: ref_line_interp_p(:,:)
  real(wp)      , intent(in)    :: ref_line_interp_s(:)
  real(wp)      , intent(in)    :: ref_line_interp_s_all(:)
  real(wp)      , intent(inout) :: ref_line_normal(:,:)
  real(wp)      , intent(in)    :: s_in(:)
  real(wp)      , intent(in)    :: nor_in(:,:)

  real(wp) :: nor1(3) , nor2(3) , nor0(3)
  character(max_char_len) :: str1 , str2
  real(wp) , parameter :: tol = 1e-6_wp

  real(wp) :: ds

  real(wp) :: w1 , w2 , w0
  integer :: i


  !> set points(...)%sec_nor
  do i = 1 , size(points)
    if ( trim(points(i)%sec_nor_str) .eq. 'reference_line' ) then
      points(i)%sec_nor = nor_in(i,:)
    elseif ( trim(points(i)%sec_nor_str) .eq. 'y_axis'    ) then
      points(i)%sec_nor = (/ 0.0_wp , 1.0_wp , 0.0_wp /)
    elseif ( trim(points(i)%sec_nor_str) .eq. 'y_axis_neg' ) then
      points(i)%sec_nor = (/ 0.0_wp ,-1.0_wp , 0.0_wp /)
    elseif ( trim(points(i)%sec_nor_str) .eq. 'vector' ) then
      ! % sec_nor assigned during reading in read_points()
    else
      write(*,*) ' Error in mod_pointwise_io.f90: '
      write(*,*) ' points(',i,')%sec_nor_str : ' , trim(points(i)%sec_nor_str)
      write(*,*) ' while the possible inputs are: "reference_line" (default), '
      write(*,*) ' "y_axis", "y_axis_neg" , "vector". Stop. ' ; stop
    endif
  end do


  !> set ref_line_normal vector
  do i = 1 , size(ref_line_normal,1)

    str1 = trim( points( ref_line_interp_p(i,1) ) %sec_nor_str )
    str2 = trim( points( ref_line_interp_p(i,2) ) %sec_nor_str )

    if (     trim(str1) .eq. 'y_axis'         ) then
      nor1 = (/ 0.0_wp , 1.0_wp , 0.0_wp /)
    elseif ( trim(str1) .eq. 'y_axis_neg'      ) then
      nor1 = (/ 0.0_wp ,-1.0_wp , 0.0_wp /)
    elseif ( trim(str1) .eq. 'vector'        ) then
      nor1 = points( ref_line_interp_p(i,1) )%sec_nor
    elseif ( trim(str1) .eq. 'reference_line' ) then
      nor1 = ref_line_normal(i,:)
    end if
    if (     trim(str2) .eq. 'y_axis'         ) then
      nor2 = (/ 0.0_wp , 1.0_wp , 0.0_wp /)
    elseif ( trim(str2) .eq. 'y_axis_neg'      ) then
      nor2 = (/ 0.0_wp ,-1.0_wp , 0.0_wp /)
    elseif ( trim(str2) .eq. 'vector'        ) then
      nor2 = points( ref_line_interp_p(i,2) )%sec_nor
    elseif ( trim(str2) .eq. 'reference_line' ) then
      nor2 = ref_line_normal(i,:)
    end if

    nor1 = nor1 / norm2(nor1)
    nor2 = nor2 / norm2(nor2)

    if ( norm2( nor1+nor2 ) .gt. tol ) then

      w1 = ( s_in( ref_line_interp_p(i,2) ) - ref_line_interp_s_all(i)       ) / &
          ( s_in( ref_line_interp_p(i,2) ) - s_in( ref_line_interp_p(i,1) ) )
      w2 = ( ref_line_interp_s_all(i)       - s_in( ref_line_interp_p(i,1) ) ) / &
          ( s_in( ref_line_interp_p(i,2) ) - s_in( ref_line_interp_p(i,1) ) )

      ref_line_normal(i,:) = nor1 * w1 +  nor2 * w2

    else

      nor0 = points( ref_line_interp_p(i,2) ) % coord - &
              points( ref_line_interp_p(i,1) ) % coord

      !> non-dimensional coord \in (-1,1)
      ds = (ref_line_interp_s_all(i) - &
            0.5_wp * ( s_in( ref_line_interp_p(i,2) ) + s_in( ref_line_interp_p(i,1) ) ) ) * &
            2.0_wp / ( s_in( ref_line_interp_p(i,2) ) - s_in( ref_line_interp_p(i,1) ) )

      w1 = 0.5_wp * ds * ( ds-1.0_wp)
      w2 = 0.5_wp * ds * ( ds+1.0_wp)
      w0 = (1.0_wp - ds)*(1.0_wp + ds)

      ref_line_normal(i,:) = nor1 * w1 +  nor2 * w2 + w0 * nor0

    end if

    !> normalisation
    ref_line_normal(i,:) = ref_line_normal(i,:) / norm2(ref_line_normal(i,:))

  end do

end subroutine update_ref_line_normal

!----------------------------------------------------------------------
!>
subroutine fill_line_tan_vec( points , lines )
  type(t_point) , intent(inout) :: points(:)
  type(t_line ) , intent(inout) :: lines( :)

  integer :: n_lines , i

  n_lines = size(lines)  

  ! first straight lines
  do i = 1 , n_lines 

    if ( trim(lines(i)%l_type) .eq. 'straight' ) then

      allocate( lines(i)%t_vec1(3) , lines(i)%t_vec2(3) ) ;

      lines(i)%t_vec1 = points( lines(i)%end_points(2) ) %coord - &
                        points( lines(i)%end_points(1) ) %coord

      if ( norm2(lines(i)%t_vec1) .lt. 1.0e-9_wp ) then
        write(*,*) ' error in fill_line_tan_vec ... stop ' ; stop
      end if
    
      lines(i)%t_vec1 = lines(i)%t_vec1 / norm2(lines(i)%t_vec1)
      lines(i)%t_vec2 = lines(i)%t_vec1

    end if

  end do

  !> then splines
  do i = 1 , n_lines

    if ( trim(lines(i)%l_type) .eq. 'spline' ) then

      if ( .not. allocated(lines(i)%t_vec1) ) then

        if  ( i .eq. 1 ) then  ! check
          write(*,*) ' error in fill_line_tan_vec: '
          write(*,*) ' first line is a spline w/o tangentVec1 input. '
          write(*,*) ' Some minor modifications in mod_pointwise_io.f90 '
          write(*,*) ' should be enough to have "free" ends.'
          write(*,*) ' stop. ' ; stop
        end if

        allocate(lines(i)%t_vec1(3))
        lines(i)%t_vec1 = lines(i-1)%t_vec2
      
      end if

      if ( .not. allocated(lines(i)%t_vec2) ) then

        if  ( i .eq. n_lines ) then  ! check
          write(*,*) ' error in fill_line_tan_vec: '
          write(*,*) ' last line is a spline w/o tangentVec2 input '
          write(*,*) ' Some minor modifications in mod_pointwise_io.f90 '
          write(*,*) ' should be enough to have "free" ends.'
          write(*,*) ' stop. ' ; stop

        end if

        allocate(lines(i)%t_vec2(3))
        lines(i)%t_vec2 = lines(i+1)%t_vec1

      end if

    end if

  end do


end subroutine fill_line_tan_vec

!----------------------------------------------------------------------
!> fill point structure
subroutine read_points ( eltype , pmesh_prs , point_prs , points, aero_table )
  character                   , intent(in)    :: eltype
  logical                     , intent(in)    :: aero_table
  type(t_parse)               , intent(inout) :: pmesh_prs
  type(t_parse) , pointer     , intent(inout) :: point_prs
  type(t_point) , allocatable , intent(out)   :: points(:)

  integer :: nPoints , i 
  
  ! === Read point groups ===
  nPoints = countoption(pmesh_prs,'point') ; allocate( points(nPoints) ) 
  ! loop over Point groups
  do i = 1 , nPoints

    call getsuboption(pmesh_prs, 'Point', point_prs)
    points(i)%id = getint(point_prs, 'Id')
    points(i)%coord = getrealarray(point_prs, 'Coordinates',3)
    if ((eltype .eq. 'p' )) then
      points(i)%airfoil     = getstr(point_prs, 'airfoil')
    elseif ((eltype .eq. 'v' )) then
      if (aero_table) then
        points(i)%airfoil_table = getstr( point_prs, 'airfoil_table')
        points(i)%airfoil = getstr(point_prs, 'airfoil')
      else
        points(i)%airfoil = getstr(point_prs, 'airfoil')
      endif 
    else if (eltype .eq. 'l') then
      points(i)%airfoil = getstr(point_prs, 'airfoil_table')
    else
      write(*,*) ' Error in read_points (Internal error: some of the programmers &
                  &did it bad): eltype must be either "p", "v", "l". Stop' ; stop
    end if
    points(i)%chord = getreal(point_prs, 'chord')
    points(i)%theta = getreal(point_prs, 'twist')
    points(i)%sec_nor_str = getstr(point_prs, 'section_normal')
    if ( trim(points(i)%sec_nor_str) .eq. 'vector' ) then
      points(i)%sec_nor = getrealarray(  point_prs, 'section_normal_vector',3)
    end if
    !> flipSection for 'p' or 'v'
    if ((eltype .eq. 'p') .or. (eltype .eq. 'v') .or. (eltype .eq. 'l')) then
      points(i)%flip_sec = getlogical(point_prs, 'flip_section', 'F')
    else
      points(i)%flip_sec = .false.
    end if
  
    point_prs => null()

  end do

end subroutine read_points

!----------------------------------------------------------------------
!> fill line structure
subroutine read_lines ( pmesh_prs , line_prs , lines  , nelems_span_tot)
  type(t_parse) ,               intent(inout) :: pmesh_prs
  type(t_parse) , pointer     , intent(inout) ::  line_prs
  type(t_line ) , allocatable , intent(out)   :: lines(:)
  integer                     , intent(out)   :: nelems_span_tot

  integer :: nLines , i

  ! === Read line groups ===
  nLines  = countoption(pmesh_prs,'Line') ; allocate( lines(nLines) )

  nelems_span_tot = 0
  do i = 1 , nLines

    call getsuboption( pmesh_prs , 'Line' , line_prs )

    lines(i)%l_type     = getstr(      line_prs , 'type'   )
    lines(i)%end_points = getintarray( line_prs , 'end_points' , 2 )
    lines(i)%nelems     = getint(      line_prs , 'nelem_line' )
    lines(i)%type_span  = getstr(      line_prs , 'type_span')
    !> geometric series parameters
    if ( trim(lines(i)%type_span) .eq. 'geoseries' ) then
      lines(i)%r_ib = getreal( line_prs , 'r_ib' )
      lines(i)%r_ob = getreal( line_prs , 'r_ob' )
      lines(i)%y_ref = getreal( line_prs , 'y_refinement' )
    elseif ( trim(lines(i)%type_span) .eq. 'geoseries_ob' ) then
      lines(i)%r_ob = getreal( line_prs , 'r_ob' )
    elseif ( trim(lines(i)%type_span) .eq. 'geoseries_ib' ) then
      lines(i)%r_ib = getreal( line_prs , 'r_ib' )
    else 
      lines(i)%r_ib = 0.0_wp
      lines(i)%r_ob = 0.0_wp
      lines(i)%y_ref = 0.0_wp
    end if

    !> tension , bias parameters
    if ( trim(lines(i)%l_type) .eq. 'spline' ) then
      lines(i)%tension = getreal(line_prs , 'tension')
      lines(i)%bias    = getreal(line_prs , 'bias'   )
    else
      lines(i)%tension = 0.0_wp
      lines(i)%bias    = 0.0_wp
    end if

    !> allocate tvec for spline if they are provided as a input
    if ( trim(lines(i)%l_type) .eq.'spline' ) then

      !> allocate tangent vec1 if specified as an input
      if ( countoption( line_prs , 'tangent_vec_1' ) .eq. 0 ) then
        ! do nothing: tan vec must be inherited
      elseif ( countoption( line_prs , 'tangent_vec_1' ) .eq. 1 ) then
        allocate(lines(i)%t_vec1(3))
        lines(i)%t_vec1 = getrealarray(line_prs , 'tangent_vec_1' , 3 )
      else
        write(*,*) ' Error in read_lines(). Provided more than one &
                    &tangent_vec_1 as an input. Stop' ; stop
      end if

      !> allocate tangent vec2 if specified as an input
      if ( countoption( line_prs , 'tangent_vec_2' ) .eq. 0 ) then
        ! do nothing: tan vec must be inherited
      elseif ( countoption( line_prs , 'tangent_vec_2' ) .eq. 1 ) then
        allocate(lines(i)%t_vec2(3))
        lines(i)%t_vec2 = getrealarray(line_prs , 'tangent_vec_2' , 3 )
      else
        write(*,*) ' Error in read_lines(). Provided more than one &
                    &tangent_vec_2 as an input. Stop' ; stop
      end if

    end if

    nelems_span_tot = nelems_span_tot + lines(i)%nelems

    line_prs => null()

  end do

end subroutine read_lines

!----------------------------------------------------------------------
!> prepare all the parser (and sub-parsers) options
subroutine set_parser_pointwise( eltype , pmesh_prs , point_prs , line_prs )
  character               , intent(in)  :: eltype
  type(t_parse)           , intent(out) :: pmesh_prs
  type(t_parse) , pointer , intent(out) :: point_prs , line_prs

  if ( ( eltype .eq. 'p' ) .or. ( eltype .eq. 'v' ) ) then

    call pmesh_prs%CreateStringOption('type_chord', &
                    'type of chord-wise division: uniform, cosine, &
                    &cosine_le, cosine_te, geoseries, geoseries_le, &
                    &geoseries_te, geoseries_hi', 'uniform', & 
                    multiple=.false.)
      
      !> geometric series parameters 
    call pmesh_prs%CreateRealOption( 'r', 'growth ratio of the elements at edge', '0.125', &
                    multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_le', 'growth ratio of the elements at leading edge', &
                    '0.142', multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_te', 'growth ratio of the elements at trailing edge', &
                    '0.067', multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_le_fix', 'growth ratio of the elements at leading edge &
                    & fixed part', '0.125', multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_te_fix', 'growth ratio of the elements at trailing edge &
                    & fixed part', '0.142', multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_le_moving', 'growth ratio of the elements at leading edge &
                    & moving part', '0.142', multiple=.false.)

    call pmesh_prs%CreateRealOption( 'r_te_moving', 'growth ratio of the elements at trailing edge & 
                    & moving part', '0.142', multiple=.false.) 

    call pmesh_prs%CreateRealOption( 'r_te_moving', 'growth ratio of the elements at trailing edge & 
                    & moving part', '0.1', multiple=.false.) 

    call pmesh_prs%CreateRealOption( 'x_refinement', 'chordwise station to which the refinement start', &
                    '0.5', multiple=.false.)

    call pmesh_prs%CreateLogicalOption('airfoil_table_correction', &
                  'include presence of aerodynamic .c81 for corrections', &
                  'F', &
                  multiple=.false.)
  end if

  ! === Point subparser ===
  call pmesh_prs%CreateSubOption('point','Point group',point_prs, &
                  multiple=.true.)

  call point_prs%CreateIntOption('Id', 'Point Id.' )
  call point_prs%CreateRealArrayOption('Coordinates', &
                'coordinates of the points used to define the comp' ) 
                
  if ( ( eltype .eq. 'p' ) .or. ( eltype .eq. 'v' ) ) then
    call point_prs%CreateStringOption( 'airfoil', 'section airfoil' ) 
    call point_prs%CreateStringOption(  'airfoil_table', 'section airfoil' )
  elseif ( eltype .eq. 'l') then
    call point_prs%CreateStringOption(  'airfoil_table', 'section airfoil' )
  endif

  if ( eltype .eq. 'l' ) then
    call pmesh_prs%CreateLogicalOption('mesh_flat', 'flat mesh yes/no','F' )
  endif
  call point_prs%CreateRealOption(      'chord', 'section chord' )
  call point_prs%CreateRealOption(      'twist', 'section twist angle' )
  call point_prs%CreateStringOption('section_normal', &
                'normal vector (str) of the plane section containing the airfoil &
                &points', 'reference_line' ) ! default y-axis
  call point_prs%CreateRealArrayOption('section_normal_vector', &
                'normal vector of the plane section containing the airfoil' )
  if ( ( eltype .eq. 'p' ) .or. ( eltype .eq. 'v' ) .or. ( eltype .eq. 'l' ) ) then
    call point_prs%CreateLogicalOption('flip_section', &
                  'flip section definition, e.g. for box wing configurations' , &
                  'F')
  end if
  ! == Remove fountain (for hover condition) == 
  call pmesh_prs%CreateRealOption('y_fountain', &
                'Upper y-coordinates beyond which the wake particles are not released', '0.0' ) 
  ! === Line sub-parser ===
  call pmesh_prs%CreateSubOption('Line','Line group',line_prs, &
                multiple=.true.)
  ! --- line group ---
  call line_prs%CreateStringOption(       'type', &
                'type of the line connecting points: Straight or Spline' )

  call line_prs%CreateIntArrayOption('end_points', &
                'list of point id.s belonging to the line' )

  call line_prs%CreateIntOption(        'nelem_line', &
                'n. spanwise elems of the line section' )

  call line_prs%CreateRealOption('tension', &
                'tension factor of the spline', '0.0' )

  call line_prs%CreateRealOption('bias', &
                'bias factor for the spline', '0.0' )

  call line_prs%CreateStringOption('type_span', 'type of span-wise division: &
                &uniform, cosine, cosine_ib, cosine_ob, equalarea, geoseries, geoseries_ib, geoseries_ob', &
                'uniform' ) ! default

  !> geo series parameters
  call line_prs%CreateRealOption( 'r_ib', 'growth ratio of the elements inboard')

  call line_prs%CreateRealOption( 'r_ob', 'growth ratio of the elements at outboard')

  call line_prs%CreateRealOption( 'y_refinement', 'spanwise station to which the refinement start') 

  !> TangentVec1,2: in the code:
  ! straight lines: useles input -> tan vec computed
  ! spline lines  : either assigned or inherited from neighbouring lines
  !
  call line_prs%CreateRealArrayOption('tangent_vec_1' , &
                'assign tangent vector to the line at its first point' )
  call line_prs%CreateRealArrayOption('tangent_vec_2' , &
                'assign tangent vector to the line at its first point' )

end subroutine set_parser_pointwise

end module mod_pointwise_io
