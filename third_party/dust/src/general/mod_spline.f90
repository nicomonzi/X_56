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
!!=========================================================================

module mod_spline

use mod_param, only: &
  wp, max_char_len, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_math, only: &
  linspace, geoseries, geoseries_both

use, intrinsic :: ieee_arithmetic

implicit none

public :: t_spline , hermite_spline , hermite_spline_profile, catmull_rom_chain, & 
          deallocate_spline

private

character(len=*), parameter :: this_mod_name = 'mod_spline' 
! ----------------------------------------------------------------------
!> spline type
type :: t_spline
  real(wp) , allocatable  :: rr(:,:)  ! points to be interpolate ( n_d , n_points)
  real(wp) , allocatable  :: dp(:,:)  ! derivative at int. points( n_d , n_points)
  real(wp) , allocatable  :: t(:)     ! array of the param describing the spline
  real(wp) , allocatable  :: dt(:)    ! dt between two points of the spline
  character(max_char_len) :: e_bc(2)  ! bc: 'der', 'free'
  real(wp) , allocatable  :: d0(:)    ! derivative at the first point
  real(wp) , allocatable  :: d1(:)    ! derivative at the end point
end type t_spline

! ----------------------------------------------------------------------
contains
! ----------------------------------------------------------------------
!> build hermite_spline
subroutine hermite_spline ( mesh_file, spl, nelems , par_tension , par_bias      , &
                            type_span, r_ib, r_ob, y_refinement,  & 
                            rr_spl , nn_spl , &
                            ip_spl, ss_spl , &
                            spl_spl, &
                            leng, s_in , nor_in )
  type(t_spline)               , intent(inout) :: spl
  integer                      , intent(in)    :: nelems
  real(wp)                     , intent(in)    :: par_tension
  real(wp)                     , intent(in)    :: par_bias
  character(len=*)             , intent(in)    :: type_span
  real(wp)                     , intent(in)    :: r_ib, r_ob, y_refinement
  character(len=*)             , intent(in)    :: mesh_file
  real(wp)                     , intent(out)   :: rr_spl(:,:)
  real(wp)                     , intent(out)   :: nn_spl(:,:)
  integer                      , intent(out)   :: ip_spl(:,:)
  real(wp)                     , intent(out)   :: ss_spl(:)
  real(wp)                     , intent(out)   :: spl_spl(:)
  real(wp)                     , intent(out)   :: leng
  real(wp)                     , intent(out)   ::   s_in(:)
  real(wp)                     , intent(out)   :: nor_in(:,:)
  character(len=*), parameter                  :: this_sub_name = 'hermite_spline'
  integer :: n_d , n_points , n_splines
  real(wp) , allocatable :: ll(:)  ! useless with this def of %t

  real(wp) , allocatable :: spl_s(:), divisionIB(:), divisionOB(:)

  ! spline reconstruction
  integer , parameter :: n_points1 = 100  ! now hardcoded (pass as a default input)
  real(wp) , allocatable :: tv1(:) , s(:) , s_spl(:)
  real(wp) , allocatable :: rr(:,:)

  integer :: i_r , i , j

  !> scale derivatives
  spl%d0 = spl%d0
  spl%d1 = spl%d1


  !> check input consistency
  if ( nelems .ne. size(rr_spl,1)-1 ) then
    write(*,*) ' error in hermite_spline : '
    write(*,*) '  nelems .ne. size(rr_spl,1)-1 . stop ' ; stop
  end if
  if ( nelems .ne. size(nn_spl,1)-1 ) then
    write(*,*) ' error in hermite_spline : '
    write(*,*) '  nelems .ne. size(nn_spl,1)-1 . stop ' ; stop
  end if

  !> ...
  if ( allocated(spl% t) ) deallocate(spl% t)
  if ( allocated(spl%dt) ) deallocate(spl%dt)
  if ( allocated(spl%dp) ) deallocate(spl%dp)


  !> dimensions
  n_d      = size(spl%rr,2)
  n_points = size(spl%rr,1) ; n_splines = n_points - 1

  ! === parameter vectors %t, %dt ===
  ! todo: check the sensitivity of the spline w.r.t. the param array t
  allocate( spl%t ( n_points  ) ) ; spl%t  = 0.0_wp
  allocate( spl%dt( n_splines ) ) ; spl%dt = 0.0_wp
  allocate( ll    ( n_splines ) ) ; ll     = 0.0_wp
  spl%t( 1) = 0.0_wp
  ll(    1) = norm2( spl%rr(2,:) - spl%rr(1,:) )
  spl%t( 2) = ll(1)
  spl%dt(1) = spl%t(2) - spl%t(1)
  do i = 2 , n_splines
    ll(    i  ) = norm2( spl%rr(i+1,:) - spl%rr(i,:) ) + ll(i-1)
    spl%t( i+1) = ll(i)
    spl%dt(i  ) = spl%t(i+1) - spl%t(i)
  end do
  deallocate(ll)

  allocate(spl%dp(n_points,n_d)) ; spl%dp = 0.0_wp

  !> following Paul Bourke description of splines
  !> http://paulbourke.net/miscellaneous/interpolation/

  !> === inner points ===
  do i = 2 , n_points - 1

    spl%dp(i,:) = 0.5_wp * ( 1.0_wp - par_tension ) * &
               ( ( spl%rr(i  ,:) -spl%rr(i-1,:) ) * ( 1.0_wp + par_bias ) + &
                 ( spl%rr(i+1,:) -spl%rr(i  ,:) ) * ( 1.0_wp - par_bias ) )

  end do

  !> === end points ===
  !-> first point
  if (     trim( spl%e_bc(1) ) .eq. 'derivative' ) then

    spl%dp(1,:) = spl%d0 / norm2(spl%d0) * norm2( spl%dp(2,:) )

  elseif ( trim( spl%e_bc(1) ) .eq. 'free'       ) then
    spl%dp(1,:) = 0.5_wp * ( - spl%dp(2,:) + 3.0_wp / spl%dt(1) * &
                            ( spl%rr(2,:) - spl%rr(1,:) ) ) ;
  else
    write(*,*) ' hermite_spline. Wrong b.c. '
    write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
  end if

  !-> last point
  if (     trim( spl%e_bc(2) ) .eq. 'derivative' ) then
    spl%dp(n_points,:) = spl%d1 / norm2(spl%d1) * norm2(spl%dp(n_points-1,:))
  elseif ( trim( spl%e_bc(2) ) .eq. 'free'       ) then
    spl%dp(n_points,:) = 0.5_wp * ( - spl%dp(n_points-1,:) + 3.0_wp / spl%dt(n_splines) * &
                                    ( spl%rr(n_points  ,:) - spl%rr(n_points-1,:) ) )
  else
    write(*,*) ' hermite_spline. Wrong b.c. '
    write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
  end if


  ! === tangent vector to the reference line ===
  nor_in = 0.0_wp
  do i = 1 , n_points
    nor_in(i,:) = spl%dp(i,:)
    nor_in(i,:) = nor_in(i,:) / norm2(nor_in(i,:))
  end do


  ! === spline (non-uniform param) ===
  allocate(rr( n_points1*n_splines + 1,n_d))
  allocate(tv1(n_points1))
  allocate( s( n_points1*n_splines + 1 ) )   ! overall curvilinear coord.
  allocate(spl_s(n_points)) ; spl_s = 0.0_wp ! overall curv. coord of the input points
  i_r = 0
  do i = 1 , n_splines

    if ( i .ne. n_splines ) then
      tv1 = (/ ( real(j,wp) , j = 0,n_points1-1 ) /) / real(n_points1,wp)
    else
      if ( allocated(tv1) ) deallocate(tv1)
      allocate(tv1(n_points1+1))
      tv1 = (/ ( real(j,wp) , j = 0,n_points1   ) /) / real(n_points1,wp)
    end if

    do j = 1 , size(tv1)
      i_r = i_r + 1  ! update index
      !> point coords
      rr(i_r,:) = hermite_p1(tv1(j)) * spl%rr(i  ,:) + &
                  hermite_p2(tv1(j)) * spl%rr(i+1,:) + &
                  hermite_d1(tv1(j)) * spl%dp(i  ,:) + &
                  hermite_d2(tv1(j)) * spl%dp(i+1,:)
      !> curvilinear coord.
      if ( i_r .ne. 1 ) then
        s(i_r) = s(i_r-1) + norm2(rr(i_r,:) - rr(i_r-1,:) )
      else
        s(i_r) = 0.0_wp  ! first point: s = 0
      end if
    end do

    spl_s(i+1) = s(i_r)

  end do

  !> s \in [ 0 , 1 ]
  leng = s(size(s))
  s     = s / s(size(s))

  s_in = spl_s           ! curvilinear coord of the input points
  spl_s = spl_s / spl_s(size(spl_s)) ! and its non-dimensionalisation


  deallocate(tv1)

  ! === spline ===
  !  allocate(rr_spl(nelems+1,3)) <- passed as an input
  allocate(s_spl(nelems+1)) ! curvilinear coord of the output points
  if ( trim(type_span) .eq. 'uniform' ) then
    s_spl = (/ ( real(j,wp) , j = 0, nelems ) /) / real(nelems,wp)
  elseif ( trim(type_span) .eq. 'cosine' ) then
    s_spl = 0.5_wp * ( 1.0_wp - cos( (/ (real(j,wp) , j = 0, nelems) /)*pi/ real(nelems,wp) ) )
  elseif ( trim(type_span) .eq. 'cosine_ob' ) then
    s_spl = sin( 0.5_wp*(/(real(j,wp) , j = 0, nelems)/)*pi/ real(nelems,wp) )
  elseif ( trim(type_span) .eq. 'cosine_ib' ) then
    s_spl = 1.0_wp - cos(0.5_wp*(/(real(j,wp) , j = 0, nelems)/)*pi/ real(nelems,wp) )
  elseif ( trim(type_span) .eq. 'equalarea' ) then
    s_spl = sqrt((/(real(j,wp) , j = 0, nelems)/)/real(nelems,wp)) 
  elseif ( trim(type_span) .eq. 'geoseries' ) then
    call geoseries_both(0.0_wp, 1.0_wp, nelems, y_refinement, r_ib, r_ob, s_spl)  
  elseif ( trim(type_span) .eq. 'geoseries_ob' ) then
    call geoseries(0.0_wp, 1.0_wp, nelems, r_ob, divisionIB, s_spl)
  elseif ( trim(type_span) .eq. 'geoseries_ib' ) then 
    call geoseries(0.0_wp, 1.0_wp, nelems, r_ib, s_spl, divisionOB)
  else 
    write(*,*) ' Mesh file   : ' , trim(mesh_file)
    write(*,*) ' type_span   : ' , trim(type_span)
    call error(this_sub_name, this_mod_name, 'Incorrect input: &
          & type_span must be equal to uniform, cosine, cosineIB, cosineOB.')
  end if

  !> first and last points
  rr_spl(       1,:) = spl%rr(       1,:)
  rr_spl(nelems+1,:) = spl%rr(n_points,:)

  nn_spl(       1,:) = rr(2,:) - rr(1,:)
  nn_spl(       1,:) = nn_spl(1,:) / norm2(nn_spl(1,:))
  nn_spl(nelems+1,:) = rr(size(rr,1),:) - rr(size(rr,1)-1,:)
  nn_spl(nelems+1,:) = nn_spl(nelems+1,:) / norm2(nn_spl(nelems+1,:))

  ip_spl(       1,:) = (/ 1 , 2 /)
  ip_spl(nelems+1,:) = (/ n_points-1, n_points /)
  ss_spl(       1)   = 0.0_wp
  ss_spl(nelems+1)   = 1.0_wp

  do i = 2 , nelems

    ! interp from the finest discretisation
    do j = 1 , size(s)
      if ( s(j) .gt. s_spl(i)  ) then

        if ( j .eq. 1 ) then
          write(*,*) ' error in hermite_spline, during the definition &
                      &of the points on the spline. Stop ' ; stop
        end if

        rr_spl(i,:) = ( s(j) - s_spl(i)   )/( s(j)-s(j-1) ) * rr(j-1,:) + &
                      ( s_spl(i) - s(j-1) )/( s(j)-s(j-1) ) * rr(j  ,:)

        nn_spl(i,:) = rr(j,:) - rr(j-1,:)
        nn_spl(i,:) = nn_spl(i,:) / norm2(nn_spl(i,:))

        exit

      end if
    end do

    ! find neighbouring input points (for input interpolation)
    do j = 1 , n_points
      if ( spl_s(j) .gt. s_spl(i) ) then

        ip_spl(i,:) = (/ j-1 , j /)
        ss_spl(i)   = ( s_spl(i) - spl_s(j-1) )/( spl_s(j)-spl_s(j-1) )

        exit

      end if
    end do

    ! spl_spl
    spl_spl(i) = ( s_spl(i) - spl_s(1) ) / ( spl_s(size(spl_s)) - spl_s(1) )

  end do


end subroutine hermite_spline

subroutine hermite_spline_profile(x, y, xq, tang_start, tang_end, yq) 
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  ! Inputs/Outputs declaration
  real(wp), intent(in)    :: x(:)          ! Input x coordinates database
  real(wp), intent(in)    :: y(:)          ! Input y coordinates database
  real(wp), intent(in)    :: xq(:)         ! Query x coordinates for interpolation
  real(wp), intent(in)    :: tang_start    ! Prescribed tangent slope at the start point x(1)
  real(wp), intent(in)    :: tang_end      ! Prescribed tangent slope at the end point x(end)
  real(wp), intent(inout) :: yq(:)         ! Output interpolated y coordinates 
  
  ! Local variables
  real(wp), allocatable   :: m(:), h(:), delta(:)
  real(wp)                :: t, h00, h01, h10, h11 
  integer                 :: i, j, n, nq
  integer                 :: k_start

  n  = size(x)
  nq = size(xq)

  ! Memory allocation for slopes and step sizes
  allocate(m(n));        m = 0.0_wp
  allocate(h(n-1));      h = 0.0_wp
  allocate(delta(n-1));  delta = 0.0_wp

  ! 1. Calculate interval step sizes (h) and secant slopes (delta)
  do i = 1, n - 1
     h(i) = max(x(i+1) - x(i), 1e-12_wp)
     delta(i) = (y(i+1) - y(i)) / h(i)
  end do

  ! 2. Assign boundary slopes at endpoints
  m(1) = tang_start
  m(n) = tang_end

  ! 3. Compute internal slopes using PCHIP (Fritsch-Carlson) monotonicity limiter
  !    This prevents undershoot/overshoot oscillations on non-uniform point distributions
  do i = 2, n - 1
     if (delta(i-1) * delta(i) > 0.0_wp) then
        ! Weighted harmonic mean for non-uniform grid spacing
        m(i) = (h(i-1) + h(i)) / ((2.0_wp * h(i) + h(i-1)) / delta(i-1) + &
                                  (h(i) + 2.0_wp * h(i-1)) / delta(i))
     else
        ! Set tangent to zero at local extrema to preserve monotonicity (avoid bulges)
        m(i) = 0.0_wp
     end if
  end do

  ! 4. Perform Hermite cubic interpolation for each query point (xq)
  k_start = 1
  do i = 1, nq
     ! Out-of-bounds protection (clamping at the edges)
     if (xq(i) <= x(1)) then
        yq(i) = y(1)
        cycle
     else if (xq(i) >= x(n)) then
        yq(i) = y(n)
        cycle
     end if

     ! Search for the interval [x(j), x(j+1)] containing xq(i)
     do j = k_start, n - 1
        if ((xq(i) >= x(j)) .and. (xq(i) <= x(j+1))) then
           k_start = j  ! Optimization for sorted xq arrays
           
           ! Normalized interval coordinate t in [0, 1]
           t = (xq(i) - x(j)) / h(j) 

           ! Cubic Hermite basis functions
           h00 =  2.0_wp*t**3 - 3.0_wp*t**2 + 1.0_wp
           h10 =  t**3 - 2.0_wp*t**2 + t
           h01 = -2.0_wp*t**3 + 3.0_wp*t**2
           h11 =  t**3 - t**2

           ! Interpolated ordinate calculation
           yq(i) = h00 * y(j) + &
                   h10 * h(j) * m(j) + &
                   h01 * y(j+1) + &
                   h11 * h(j) * m(j+1)

           ! Safety check for NaN values
           if (ieee_is_nan(yq(i))) yq(i) = 0.0_wp
           exit
        end if
     end do 
  end do

  ! Cleanup local allocations
  deallocate(m, h, delta) 
end subroutine hermite_spline_profile

! Calculate Catmull-Rom for a sequence of initial points and return the combined curve.
! param points: Base points from which the quadruples for the algorithm are taken
! param num_points: The number of points to include in each curve segment
! return: The chain of all points (points of all segments)
subroutine catmull_rom_chain(points, num_points, chain_points)
  real(wp),              intent(in)  :: points(:,:)
  integer,               intent(in)  :: num_points
  real(wp), allocatable, intent(out) :: chain_points(:,:) 

  integer, parameter    :: quadruple_size = 4
  real(wp), allocatable :: points_quadruples(:,:,:) 
  integer               :: idx_segment_start, j, i, num_segments 
  
  num_segments = size(points, 2) - 3 
  allocate(points_quadruples(2, quadruple_size, num_segments)); points_quadruples = 0.0_wp 
  do idx_segment_start = 1, num_segments
    points_quadruples(:, :, idx_segment_start) = points(:, idx_segment_start:idx_segment_start + quadruple_size-1)
  end do

  allocate(chain_points(2,num_points*num_segments)); chain_points = 0.0_wp
  j = 1 
  do i = 1, num_segments
    chain_points(:, j:j + num_points-1) = catmull_rom_spline( points_quadruples(:, 1, i), &
                                                              points_quadruples(:, 2, i), &
                                                              points_quadruples(:, 3, i), &
                                                              points_quadruples(:, 4, i), &
                                                              num_points)
    j = j + num_points 
  end do 
  !> cleanup
  deallocate(points_quadruples)

end subroutine catmull_rom_chain 

! Compute the points in the spline segment
! param P0, P1, P2, and P3: The (x,y) point pairs that define the Catmull-Rom spline
! param num_points: The number of points to include in the resulting curve segment
! param alpha: 0.5 for the centripetal spline (default), 0.0 for the uniform spline, 1.0 for the chordal spline.
! return: The points
function catmull_rom_spline(P0, P1, P2, P3, num_points) result(spline_points) 
  real(wp), intent(in) :: P0(:), P1(:), P2(:), P3(:)
  integer, intent(in)  :: num_points

  real(wp)             :: t0, t1, t2, t3
  real(wp), allocatable :: A1(:,:), A2(:,:), A3(:,:), B1(:,:), B2(:,:) 
  real(wp), allocatable :: spline_points(:,:), t(:)
  integer               :: i, j
  

  t0 = 0.0_wp
  t1 = tj(t0, P0, P1)
  t2 = tj(t1, P1, P2)
  t3 = tj(t2, P2, P3)
  
  allocate(t(num_points)); t = 0.0_wp 
  allocate(A1(2,num_points)); A1 = 0.0_wp
  allocate(A2(2,num_points)); A2 = 0.0_wp
  allocate(A3(2,num_points)); A3 = 0.0_wp
  allocate(B1(2,num_points)); B1 = 0.0_wp
  allocate(B2(2,num_points)); B2 = 0.0_wp
  allocate(spline_points(2,num_points)); spline_points = 0.0_wp
  
  call linspace(t1, t2, t) 
  
  do i = 1, 2
    do j = 1, num_points
    A1(i,j) = (t1 - t(j))/(t1 - t0)*P0(i)   + (t(j) - t0)/(t1 - t0)*P1(i) 
    A2(i,j) = (t2 - t(j))/(t2 - t1)*P1(i)   + (t(j) - t1)/(t2 - t1)*P2(i)
    A3(i,j) = (t3 - t(j))/(t3 - t2)*P2(i)   + (t(j) - t2)/(t3 - t2)*P3(i)
    B1(i,j) = (t2 - t(j))/(t2 - t0)*A1(i,j) + (t(j) - t0)/(t2 - t0)*A2(i,j) 
    B2(i,j) = (t3 - t(j))/(t3 - t1)*A2(i,j) + (t(j) - t1)/(t3 - t1)*A3(i,j) 
    spline_points(i,j) = (t2 - t(j))/(t2 - t1)*B1(i,j) + (t(j) - t1)/(t2 - t1)*B2(i,j) 
    enddo 
  enddo 
  !> cleanup
  deallocate(A1, A2, A3, B1, B2, t)
  
end function catmull_rom_spline

! Calculate t0 to t4. Then only calculate points between P1 and P2.
! Reshape linspace so that we can multiply by the points P0 to P3
! and get a point for each value of t.
function tj(ti, pii, pjj) result(t) 
  real(wp), intent(in) :: ti
  real(wp), intent(in) :: pii(:)
  real(wp), intent(in) :: pjj(:) 
  
  real(wp) :: t, dx, dy, l, xi, yi, xj, yj 
  real(wp), parameter :: alpha = 0.5_wp !> centripetal spline (default), 
                                        !  0.0 for the uniform spline, 
                                        !  1.0 for the chordal spline.

  xi = pii(1)
  yi = pii(2)
  xj = pjj(1)
  yj = pjj(2)

  dx = xj - xi
  dy = yj - yi
  l = sqrt(dx**2.0_wp + dy**2.0_wp)
  t = ti + l**alpha

end function

! ----------------------------------------------------------------------
!> hermite functions
!> h1: h1(0)=1 , h1(1)=0 , h1'(0)=0 , h1'(1)=0
function hermite_p1(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h =   2.0_wp * t**3.0_wp - 3.0_wp * t**2.0_wp + 1.0_wp

end function hermite_p1

!> h2: h2(0)=0 , h2(1)=1 , h2'(0)=0 , h2'(1)=0
function hermite_p2(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = - 2.0_wp * t**3.0_wp + 3.0_wp * t**2.0_wp

end function hermite_p2

!> h3: h3(0)=0 , h3(1)=0 , h3'(0)=1 , h3'(1)=0
function hermite_d1(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = t**3.0_wp - 2.0_wp * t**2.0_wp + t

end function hermite_d1

!> h4: h4(0)=0 , h4(1)=0 , h4'(0)=0 , h4'(1)=1
function hermite_d2(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = t**3.0_wp - t**2.0_wp

end function hermite_d2

! ----------------------------------------------------------------------
!> deallocate t_spline structure
subroutine deallocate_spline( spl )
  type(t_spline) , intent(inout) :: spl

  if ( allocated(spl%rr  ) ) deallocate(spl%rr  )
  if ( allocated(spl%dp  ) ) deallocate(spl%dp  )
  if ( allocated(spl%t   ) ) deallocate(spl%t   )
  if ( allocated(spl%dt  ) ) deallocate(spl%dt  )
  if ( allocated(spl%d0  ) ) deallocate(spl%d0  )
  if ( allocated(spl%d1  ) ) deallocate(spl%d1  )

end subroutine deallocate_spline

! ----------------------------------------------------------------------


end module mod_spline


