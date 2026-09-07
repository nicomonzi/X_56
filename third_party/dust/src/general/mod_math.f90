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

!> Module containing some mathematical utilities
module mod_math

use mod_param, only: &
  wp, eps, max_char_len

use mod_handling, only: &
  error, warning

implicit none

public :: dot, cross , linear_interp , compute_qr, &
          rotation_vector_combination, sort_vector_real, & 
          unique, unique_row, infinite_plate_spline, tessellate, & 
          vec2mat, mat2vec, invmat, invmat_banded, linspace, &
          geoseries, geoseries_both, geoseries_hinge

private

interface linear_interp
  module procedure  linear_interp_vector, &
                    linear_interp_array
end interface linear_interp

character(len=*), parameter :: this_mod_name='mod_math'

contains

pure function dot(a, b)
  real(wp)                           :: dot
  real(wp), dimension(:), intent(in) :: a, b

  dot = sum(a * b)

end function dot


! ----------------------------------------------------------------------

pure function cross(a, b)
  real(wp), dimension(3) :: cross
  real(wp), dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

! ----------------------------------------------------------------------

subroutine linear_interp_vector(val_vec , t_vec , t , val)
  real(wp), intent(in)        :: val_vec(:)
  real(wp), intent(in)        :: t_vec(:)
  real(wp), intent(in)        :: t
  real(wp), intent(out)       :: val
  integer                     :: it , nt
  character(len=*), parameter :: this_sub_name='linear_interp_vector'

  nt = size(t_vec)
  ! Check dimensions -----
  if ( size(val_vec) .ne. nt ) then
    call error(this_sub_name, this_mod_name, 'Different sizes for x and y &
                                  &data vector provided for interpolation')
  end if
  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .le. minval(t_vec)-eps ) then
    write(*,*) t, minval(t_vec)
    call warning(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is lower than the minimum of the interpolation data')
  end if
  if ( t .ge. maxval(t_vec)+eps ) then
    write(*,*) t, maxval(t_vec)
    call warning(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is higher than the maximum of the interpolation data')
  end if

  do it = 1 , nt-1

    if ( ( t .ge. t_vec( it ) ) .and. ( t .le. t_vec( it+1 ) ) ) then

      val = val_vec(it) + (t-t_vec(it))/(t_vec(it+1)-t_vec(it)) * &
                                  ( val_vec(it+1)-val_vec(it) )

    end if

  end do

end subroutine linear_interp_vector

! ----------------------------------------------------------------------

subroutine linear_interp_array( val_arr , t_vec , t , val )
  real(wp) , intent(in) :: val_arr(:,:)
  real(wp) , intent(in) :: t_vec(:)
  real(wp) , intent(in) :: t
  real(wp) , allocatable , intent(out) :: val(:)

  integer :: it , nt
  character(len=*), parameter :: this_sub_name='linear_interp_array'

  nt = size(t_vec)
  allocate(val(size(val_arr,1)))

  ! Check dimensions -----
  if ( size(val_arr,2) .ne. nt ) then
    call error(this_sub_name, this_mod_name, 'Different sizes for x and y &
                                  &data vector provided for interpolation')
  end if

  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .le. minval(t_vec)-eps ) then
    call warning(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is lower than the minimum of the interpolation data')
  end if
  if ( t .ge. maxval(t_vec)+eps ) then
    call warning(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is higher than the maximum of the interpolation data')
  end if

  do it = 1 , nt-1

    if ( ( t .ge. t_vec( it ) ) .and. ( t .le. t_vec( it+1 ) ) ) then

      val = val_arr(:,it) + (t-t_vec(it))/(t_vec(it+1)-t_vec(it)) * &
                                  ( val_arr(:,it+1)-val_arr(:,it) )

    end if
  end do

end subroutine linear_interp_array

! ----------------------------------------------------------------------
subroutine invmat(A) 
  real(wp), intent(inout) :: A(:,:)
  real(wp), allocatable   :: work(:)
  integer                 :: info
  integer, allocatable    :: ipiv(:)
  character(len=*), parameter :: this_sub_name= 'invmat'

  !> invert the A matrix
  allocate(ipiv(size(A,1)))
  allocate(work(size(A,1)))
#if (DUST_PRECISION == 1)
  call sgetrf(size(A,1), size(A,1), A, size(A,1), ipiv, info)  
  call sgetri(size(A,1), A, size(A,1), ipiv, work, size(A,1), info)
#elif(DUST_PRECISION == 2)
  call dgetrf(size(A,1), size(A,1), A, size(A,1), ipiv, info)  
  call dgetri(size(A,1), A, size(A,1), ipiv, work, size(A,1), info)
#endif
  deallocate(ipiv, work) 
end subroutine invmat
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
subroutine invmat_banded(A, n)
  integer, intent(in)                     :: n   ! n is the A size 
  real(wp), intent(inout)                 :: A(n,n)
  
  integer                                 :: info, nb
  real(wp)                                :: work(n)
  integer                                 :: i1, i2, j1, j2
  real(wp), allocatable                   :: Anz(:,:)
  integer, allocatable                    :: ipiv(:)
  
  character(len=max_char_len)             :: msg
  character(len=*), parameter             :: this_sub_name = 'invmat'
  !> Compute the inverse of the non-zero blocks of the matrix
  !> Has to explore whole matrix because it could be banded, ex:
  !> |0101|
  !> |0000|
  !> |0101|
  !> 
  
  ! (i1,j1) and (i2,j2) are the extremes of the non-zero block
  i2 = 0
  j2 = 0

!do i1=1,n
!write(*,*) A(i1,:)
!enddo
!write(*,*) "++++"
  do i1 = 1,n
      do j1 = 1,n
        ! check if we are still inside a previous block
        if(j1 .ge. j2+1 .or. i1 .ge. i2+1) then
          if (abs(A(i1,j1)) .ge. 1e-16_wp) then
            i2 = i1
            
            ! find the end of the block
            do while(abs(A(i2,j1)) .ge. 1e-16_wp .and. i2 .lt. n)
              i2 = i2 + 1
            enddo
            j2 = j1
            do while(abs(A(i2-1,j2)) .ge. 1e-16_wp .and. j2 .lt. n)
              j2 = j2 + 1
            enddo
        
            if (i2 .ne. n) i2 = i2-1
            if (j2 .ne. n) j2 = j2-1
            nb = i2-i1+1
            
            allocate(Anz(nb,nb)); Anz = 0.0_wp
            ! extract the block
            Anz = A(i1:i2,j1:j2)
            allocate(ipiv(nb)); ipiv = 0
! write(*,*) Anz
            !> Factorize the block
#if (DUST_PRECISION==1)
            call sgetrf(nb, nb, Anz, nb, ipiv, info)
            if (info /= 0) then
              write(msg,*) 'sgetrf failed with info = ', info
              call error(this_sub_name, this_mod_name, trim(msg))
            end if
#elif(DUST_PRECISION==2)
            call dgetrf(nb, nb, Anz, nb, ipiv, info)
            if (info /= 0) then
              write(msg,*) 'dgetrf failed with info = ', info
              call error(this_sub_name, this_mod_name, trim(msg))
            end if
#endif /*DUST_PRECISION*/

            !> Compute the inverse of the non-zero block
#if (DUST_PRECISION==1)
            call sgetri(nb, Anz, nb, ipiv, work, nb, info)
            if (info /= 0) then
              write(msg,*) 'sgetri failed with info = ', info
              call error(this_sub_name, this_mod_name, trim(msg))
            end if
#elif(DUST_PRECISION==2)
            call dgetri(nb, Anz, nb, ipiv, work, nb, info)
            if (info /= 0) then
              write(msg,*) 'dgetri failed with info = ', info
              call error(this_sub_name, this_mod_name, trim(msg)) 
            end if
#endif /*DUST_PRECISION*/
        
            !> Update the matrix with the inverse of the non-zero block
            A(i1:i2, j1:j2) = Anz
            
            deallocate(Anz)  
            deallocate(ipiv)      
          endif
        
        ! skip because still inside a previous block
        else
          cycle
        endif
          
      enddo !j1
  enddo !i1
end subroutine invmat_banded


! ----------------------------------------------------------------------
subroutine compute_qr ( A , Q , R )
  real(wp) , intent(inout) ::  A(:,:)
  real(wp) , allocatable , intent(inout) :: Q(:,:) , R(:,:)

  integer :: m , n , i_m

  ! lapack dgeqrf routine ----
  real(wp) , allocatable :: tau(:) , work(:)
  integer :: lwork , info

  ! tmp matrices to get Q matrix -----
  real(wp) , allocatable :: H(:,:) , v(:,:) , eye(:,:)
  character(len=*), parameter :: this_sub_name='compute_qr'


  ! input check and warnings
  if ( allocated(Q) ) then
    call warning(this_sub_name, this_mod_name, ' Q was already allocated.&
                                          & Deallocated and re-allocated')
    deallocate(Q)
  end if
  if ( allocated(R) ) then
    call warning(this_sub_name, this_mod_name, ' R was already allocated.&
                                          & Deallocated and re-allocated')
    deallocate(R)
  end if

  m = size(A,1)
  n = size(A,2)

  ! qr factorisation of matrix B
  allocate(   tau(min(m,n)) ) ;    tau = 0.0_wp
  allocate(  work(    n   ) ) ;   work = 0.0_wp
  lwork = n       ! <-- its size should be .ge. n*nb
                  ! with nb = optimal blocksize (???)

#if (DUST_PRECISION == 1)
  call sgeqrf( m , n , A , m , tau , work , lwork , info )
#elif(DUST_PRECISION == 2)
  call dgeqrf( m , n , A , m , tau , work , lwork , info )
#endif /*DUST_PRECISION*/

  ! build Q , R matrices
  ! allocataion and initialisation to 0.0 (useless for Q)
  allocate( Q(m,m) , R(m,n) ) ; R = 0.0_wp

  ! R upper triangular matrix (initialised to zero!)
  do i_m = 1 , m
    R( i_m , i_m : size(A,2) ) = A( i_m , i_m : size(A,2) )
  end do

  ! Q square unitary matrix
  allocate( H(m,m) , v(m,1) ) ! ; H = 0.0_wp  ; v(:,1) = 0.0_wp
  allocate( eye(m,m) ) ; eye = 0.0_wp ; 
  
  do i_m = 1, m 
    eye(i_m,i_m) = 1.0_wp 
  end do

  Q = eye ! Initialisation
  do i_m = 1 , min(m,n)
    v(:,1) = 0.0_wp ; v(i_m,1) = 1.0_wp ; v(i_m+1:m,1) = A(i_m+1:m,i_m) ;
    H = eye - tau(i_m) * matmul( v , transpose(v) )
    Q = matmul( Q , H )
  end do

  deallocate(H,v,eye)
  deallocate(tau,work)

end subroutine compute_qr


!> convert an orientation matrix into an orientation vector 
subroutine mat2vec(R, theta)
  real(wp), intent(in)  :: R(3,3)
  real(wp), intent(out) :: theta(3) 
  
  real(wp)              :: n_sintheta(3), ax_R(3)
  real(wp)              :: trace_R, sintheta, costheta, ntheta 
  real(wp)              :: RR(3, 3) 

  !> take the axial of the rotation matrix 
  RR = (R - transpose(R))/2.0_wp 
  ax_R = (/RR(3, 2), RR(1, 3), RR(2, 1)/) 

  !> take the trace of the rotation matrix 
  trace_R = R(1, 1) + R(2, 2) + R(3, 3) 

  n_sintheta = ax_R
  sintheta = norm2(n_sintheta)
  costheta = (trace_R - 1.0_wp)/2.0_wp
  ntheta = atan2(sintheta, costheta) 

  if (sintheta .lt.  1e-6_wp) then
    theta = n_sintheta
  else 
    theta = (ntheta/sintheta)*n_sintheta
  endif 

end subroutine

!> convert an orientation vector into an orientation matrix
subroutine vec2mat(theta, R)
  real(wp), intent(in)  :: theta(3)
  real(wp), intent(out) :: R(3,3) 
  
  real(wp)              :: t, a, b
  real(wp)              :: thcross(3,3), eye(3,3)
  integer               :: i 
  
  eye = 0.0_wp
  do i = 1, 3
    eye(i,i)  =  1.0_wp
  end do

  t = norm2(theta)
  if (t .ge. 1e-6_wp) then 
    a = sin(t)/t 
    b = (1 - cos(t)) / t**2 
    thcross(1,:) = (/0.0_wp,  -theta(3), theta(2)/)
    thcross(2,:) = (/theta(3), 0.0_wp,  -theta(1)/)
    thcross(3,:) = (/-theta(2), theta(1), 0.0_wp /)

    R = eye + a*thcross + b*matmul(thcross,thcross)
    
  else 
    R = eye  
  endif 


end subroutine
!----------------------------------------------------------------
!> Combination of two rotations, provided as rotation vectors r1, r2
! x: original vector
! y: rotated vector
! R1: rotation matrix from rotation vector r1, r1 -> R1
! R2: rotation matrix from rotation vector r2, r2 -> R2
!   y = R1 * R2 * x = R * x
! Inputs: r1, r2
! Outputs: r, th, n
! R <- r = th * n, with |n| = 1
!
subroutine rotation_vector_combination( r1, r2, r, th, n )
  real(wp),              intent(in)  :: r1(3), r2(3)
  real(wp),              intent(out) :: r(3)
  real(wp),              intent(out) :: n(3)
  real(wp),              intent(out) :: th

  real(wp) :: n1(3), n2(3), q1(4), q2(4), qs(4)
  real(wp) :: th1, th2, sin_t_2, cos_t_2

  if ( ( size(r1) .ne. 3 ) .or. ( size(r2) .ne. 3 ) ) then
    write(*,*) ' Wrong input dimensions in rotation_vector_combinations. Stop '
  end if

  th1 = norm2(r1)
  if ( th1 .ne. 0.0_wp ) then
    n1 = r1 / th1
    q1 = (/ cos(th1/2), n1(1)*sin(th1/2), &
                        n1(2)*sin(th1/2), &
                        n1(3)*sin(th1/2) /)
  else
    q1 = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  end if
  th2 = norm2(r2)
  if ( th2 .ne. 0.0_wp ) then
    n2 = r2 / th2
    q2 = (/ cos(th2/2), n2(1)*sin(th2/2), &
                        n2(2)*sin(th2/2), &
                        n2(3)*sin(th2/2) /)
  else
    q2 = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  end if

  !> Quaternion qs = q1 * q2 = cos(th/2) + sin(th/2) n · (i,j,k)
  qs = (/ &
        q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4) , &
        q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3) , &
        q1(1)*q2(3) + q1(3)*q2(1) + q1(4)*q2(2) - q1(2)*q2(4) , &
        q1(1)*q2(4) + q1(4)*q2(1) + q1(2)*q2(3) - q1(3)*q2(2)   &
        /)

  !> Rotation vector
  sin_t_2 = norm2(qs(2:4))
  if ( sin_t_2 .eq. 0.0_wp ) then
    n = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    cos_t_2 = 1.0_wp
    th = 0.0_wp
  else
    n = qs(2:4) / sin_t_2
    cos_t_2 = qs(1)
    th = 2.0_wp * atan2( sin_t_2, cos_t_2 )
  end if

  r = th * n


end subroutine rotation_vector_combination


! ----------------------------------------------------------------------

subroutine sort_vector_real( vec, nel, sort, ind )
  real(wp), intent(inout)               :: vec(:)
  integer , intent(in)                  :: nel
  real(wp), allocatable, intent(out)    :: sort(:)
  integer , allocatable, intent(out)    :: ind(:)

  real(wp)                              :: maxv
  integer                               :: i

  allocate(sort(nel)); sort = 0.0_wp
  allocate(ind(nel)); ind = 0

  maxv = maxval(vec)
  do i = 1, nel
    sort(i) = minval(vec, 1)
    ind(i) = minloc(vec, 1)
    vec(ind(i)) = maxv + 0.1_wp 
  end do

end subroutine sort_vector_real

subroutine unique(vec, vec_unique, tol)  
  real(wp), intent(inout)               :: vec(:) 
  real(wp), allocatable, intent(out)    :: vec_unique(:)
  real(wp), intent(in)                  :: tol

  real(wp), allocatable                 :: vec_sort(:) 
  real(wp), allocatable                 :: vec_tmp(:)
  integer, allocatable                  :: ind(:)
  integer                               :: nel, i, j, u 
  
  nel = size(vec)

  call sort_vector_real(vec, nel, vec_sort, ind)

  i = 1;   j = 2;   u = 1

  allocate(vec_tmp(nel)); vec_tmp = 0.0_wp 

  vec_tmp(1) = vec_sort(1)
  do while (j .le. nel)
    if (abs(vec_sort(j) - vec_sort(i)) .le. tol) then  
      j = j + 1
    else
      u = u + 1
      vec_tmp(u) = vec_sort(j) 
      i = j 
      j = j + 1   
    endif 
  enddo
  deallocate(vec_sort)
  allocate(vec_unique(u)); vec_unique = 0.0_wp 
  vec_unique = vec_tmp(1:u) 

end subroutine unique  

subroutine unique_row(vec, vec_unique, tol)  
  real(wp), intent(inout)                  :: vec(:,:) 
  real(wp), intent(in)                  :: tol
  real(wp), allocatable, intent(out)    :: vec_unique(:,:)
  
  real(wp), allocatable                 :: vec_sort(:) 
  real(wp), allocatable                 :: vec_tmp(:,:)
  real(wp), allocatable                 :: vec_sort_row(:,:)
  integer,  allocatable                 :: ind(:)
  integer                               :: nel, i, j, u 
  !> column vector 
  nel = size(vec, 1) 
  allocate(vec_sort_row(nel, size(vec,2))); vec_sort_row = 0.0_wp

  call sort_vector_real(vec(:,1), nel, vec_sort, ind)

  vec_sort_row = vec(ind, :) 

  i = 1;   j = 2;   u = 1

  allocate(vec_tmp(nel, size(vec,2))); vec_tmp = 0.0_wp 

  vec_tmp(1, :) = vec_sort_row(1,:) 
  do while (j .le. nel)
    if (abs(vec_sort(j) - vec_sort(i)) .le. tol) then  
      j = j + 1
    else
      u = u + 1
      vec_tmp(u,:) = vec_sort_row(j, :) 
      i = j 
      j = j + 1   
    endif 
  enddo
  
  deallocate(vec_sort, vec_sort_row)

  allocate(vec_unique(u, size(vec,2))); vec_unique = 0.0_wp 
  vec_unique = vec_tmp(1:u, :)  

end subroutine unique_row 


! ----------------------------------------------------------------------
! RBF interpolation weight matrix for 2D structures (plates)
subroutine infinite_plate_spline(pos_interp, pos_ref, W)
  real(wp), intent(in)                  :: pos_interp(:,:) ! positions to be interpolated into
  real(wp), intent(in)                  :: pos_ref(:,:) ! positions of the reference points
  real(wp), allocatable, intent(inout)  :: W(:,:) ! weight matrix
  
  integer                               :: n_r, n_i, i, j
  real(wp), allocatable                 :: R_r(:,:), R_i(:,:), Z_r(:,:), Z_ir(:,:), Y_r(:,:)
  real(wp), allocatable                 :: eye(:,:)                    
  real(wp)                              :: nrm 
  
  n_r = size(pos_ref,2)
  n_i = size(pos_interp,2)
  
  allocate(R_r(n_r,4))
  allocate(R_i(n_i,4))
  ! matrixes R [1|x|y|z]
  R_r(:,1) = 1.0_wp
  R_r(:,2:4) = transpose(pos_ref)
  R_i(:,1) = 1.0_wp
  R_i(:,2:4) = transpose(pos_interp) 

  allocate(Z_r(n_r,n_r))
  do i=1,n_r
    do j=1,n_r
      if (i .eq. j) then
        Z_r(i,j) = 0.0_wp
      else
        nrm = norm2(pos_ref(:,i)-pos_ref(:,j))
        !dist = pos_ref(:,i)-pos_ref(:,j)
        !nrm = norm2(dist*matmul(Wnorm,dist))
        Z_r(i,j) = nrm**2*log(nrm)
      endif
    enddo
  enddo
  allocate(Z_ir(n_i,n_r))
  do i=1,n_i
    do j=1,n_r
      nrm = norm2(pos_interp(:,i)-pos_ref(:,j))
      !dist = pos_interp(:,i)-pos_ref(:,j)
      !nrm = norm2(dist*matmul(Wnorm,dist))
      if (nrm .le. 1e-10_wp) then
        Z_ir(i,j) = 0.0_wp
      else
        Z_ir(i,j) = nrm**2*log(nrm)
      endif
    enddo
  enddo  

  ! inverse matrix (NB the inverse is overwritten into Z_r)
  call invmat(Z_r)

  allocate(Y_r(4,4))  
  
  ! inverse matrix (NB the inverse is overwritten into Y_r)
  Y_r = matmul(transpose(R_r),matmul(Z_r,R_r))
  Y_r = Y_r + 1e-6_wp
  call invmat(Y_r)

  allocate(eye(n_r,n_r))
  eye = 0.0_wp
  do i = 1, n_r
    eye(i,i)  =  1.0_wp
  end do

  allocate(W(n_i,n_r))
  W = matmul((matmul(R_i,matmul(Y_r,transpose(R_r))) + &
      matmul(Z_ir,(eye-matmul(Z_r,matmul(R_r,matmul(Y_r,transpose(R_r))))))),Z_r)

  deallocate(R_i)
  deallocate(R_r)
  deallocate(Z_r, Z_ir)
  deallocate(Y_r, eye)
  
end subroutine infinite_plate_spline

! ----------------------------------------------------------------------
! Tessellate a given quadrilateral into triangles
! Firstly divides the sides in multiples of the shortest side and creates 
! triangles with the resulting points
! The following n_steps steps are midpoint tessellations
! Returns centres of triangles and radius of each circumscribed circle
!
! Midpoint tessellation based on: 
!   https://lindenreidblog.com/2017/12/03/simple-mesh-tessellation-triangulation-tutorial/
!   and adapted to C by Stefano Colli
subroutine tessellate(vertices, n_steps, tol, cen_sbprt, rad_sbprt, n_sbprt, dry_run_in)
  real(wp), intent(in)               :: vertices(:,:) ! coordinates of quadri vertices
  integer, intent(in)                :: n_steps ! number of midpoint triangulation steps
  real(wp), intent(in)               :: tol ! tolerance for quadri initial division
  real(wp), allocatable, intent(out) :: cen_sbprt(:,:), rad_sbprt(:) ! centre coords and radius of subparticles
  integer, intent(out)               :: n_sbprt ! number of generated subparts
  logical, intent(in), optional      :: dry_run_in
  
  real(wp), allocatable :: tex(:,:,:) ! stores the tessellation 
  integer               :: i, j, n_tri_tot
  logical               :: dry_run
  
  ! variables for initial division
  real(wp)              :: lengths(4), min_side_len, succ(3,4), mid(3), abc(3,3)
  real(wp), allocatable :: division_points(:,:), init_tri(:,:,:)
  integer               :: n_div(4), n_init_tri, i_min_side, ip, first_v(4), a, b, c, abcmax
  
  ! variables for triforce
  real(wp) :: mid_ab(3), mid_bc(3), mid_ca(3)
  integer  ::  n_tri_used, idx_glob, idx_loc, offset
  
  if (present(dry_run_in)) then
    dry_run = dry_run_in
  else
    dry_run = .false.
  endif
  
  ! Initial division
  succ = cshift(vertices,1,2) ! array of successive points to vertices
  lengths = norm2(succ-vertices,1) ! length of each side of quadri
  min_side_len = minval(lengths) ! minimum side length
  i_min_side = minloc(lengths,1) ! which side is minimum
  ! compute how many times the min side is into each side
  n_div = max(floor(lengths/min_side_len+tol),1) ! number of subdivision for each side
  
  n_sbprt = (sum(n_div*2)-4)*4**n_steps
  
  if (dry_run) then ! return only n_sbprt *CAUTION* other arrays are not allocated
    return
  endif  
    
  ! arrays to hold points and triangles for initial subdivision
  allocate(division_points(3,sum(n_div)), init_tri(3,3,sum(n_div*2)-4))
  
  ! start subdivision: cycle over all sides and if it is longer than min side
  ! divide it in parts equal to the min side
  ip = 1 ! points index
  do i = 1,4 ! for all 4 sides
    first_v(i) = ip ! index of first vertex of side
    division_points(:,ip) = vertices(:,i) ! first point is the vertex
    ip = ip + 1
    if (n_div(i) .gt. 1) then ! if side is too long
      do j = 1,n_div(i)-1 ! divide it 
        ! linear combination of side vertices
        division_points(:,ip) = (real(j,wp)/real(n_div(i),wp))*succ(:,i) + &
                          (1.0_wp-real(j,wp)/real(n_div(i),wp))*vertices(:,i)
        ip = ip +1
      enddo
    endif
  enddo
  ip = ip-1 ! ip now is the total number of points

  ! Starting from min side, connect the points alternately creating triangles
  ! and actually divide this triangles in half (to get a clean behaviour with squares)
  n_init_tri = 0 ! counter for tri
  ! first point of min side
  b = first_v(i_min_side)
  ! second point of min side
  a = modd(b+1,ip)

  do i = 1,ip-2
    ! third point, alternate between opposite sides as we move along
    c = mod(i,2)*modd(a+1,ip) + mod(i-1,2)*modd(a-1,ip)
    
    ! coordinates of the points of the triangle
    abc(:,1) = division_points(:,a)
    abc(:,2) = division_points(:,b)
    abc(:,3) = division_points(:,c)
    
    ! which side of the tri is max      
    abcmax = maxloc(norm2(cshift(abc,1,2)-abc,1),1)
    
    ! divide in half max side
    mid = (abc(:,abcmax)+abc(:,modd(abcmax+1,3)))/2
    
    ! split the tri in two tris with the midpoint
    n_init_tri = n_init_tri+1
    init_tri(:,1,n_init_tri) = mid
    init_tri(:,2,n_init_tri) = abc(:,modd(abcmax+2,3))
    init_tri(:,3,n_init_tri) = abc(:,modd(abcmax,3))
    
    n_init_tri = n_init_tri+1
    init_tri(:,1,n_init_tri) = mid
    init_tri(:,2,n_init_tri) = abc(:,modd(abcmax+2,3))
    init_tri(:,3,n_init_tri) = abc(:,modd(abcmax+1,3))
    
    a = b
    b = c
  enddo  
  
  ! compute the total number of tri in the tessellation
  n_tri_tot = n_init_tri
  ! inital tris + n_steps of triforce
  do i=1,n_steps
    n_tri_tot = n_tri_tot + 4**i*n_init_tri
  end do
  
  ! allocate array of tris for the tessellation
  ! new tris will be appended to it
  allocate(tex(3,3,n_tri_tot))
  tex = 0.0_wp
  ! fill it with the tris we already have
  do i = 1,n_init_tri
    tex(:,:,i) = init_tri(:,:,i)
  enddo
  
  deallocate(division_points, init_tri)
  
  ! Midpoint division
  ! each tri is now divided in 4 by taking its sides' midpoints
  ! process is repeated n_steps times
  n_tri_used = 0  
  do i=1,n_steps
    offset = 1
    do j=1,n_init_tri*4**(i-1)
      ! idx_glob is the index for where the new tris are going in the tex array
      ! idx_loc is the index for which tri to divide in 4
      idx_glob = n_tri_used + n_init_tri*4**(i-1) + offset
      idx_loc = n_tri_used + j
            
      mid_ab = (tex(:,1,idx_loc)+tex(:,2,idx_loc))/2
      mid_bc = (tex(:,2,idx_loc)+tex(:,3,idx_loc))/2
      mid_ca = (tex(:,3,idx_loc)+tex(:,1,idx_loc))/2

      tex(:,1,idx_glob) = tex(:,1,idx_loc) ! tri.a
      tex(:,2,idx_glob) = mid_ab
      tex(:,3,idx_glob) = mid_ca
      
      tex(:,1,idx_glob+1) = tex(:,2,idx_loc) ! tri.b
      tex(:,2,idx_glob+1) = mid_ab
      tex(:,3,idx_glob+1) = mid_bc
      
      tex(:,1,idx_glob+2) = tex(:,3,idx_loc) ! tri.c
      tex(:,2,idx_glob+2) = mid_bc
      tex(:,3,idx_glob+2) = mid_ca
      
      tex(:,1,idx_glob+3) = mid_ab
      tex(:,2,idx_glob+3) = mid_bc
      tex(:,3,idx_glob+3) = mid_ca
              
      offset = offset+4
    enddo
    n_tri_used = n_tri_used+n_init_tri*4**(i-1)
  enddo
              
  ! compute cen and area
  ! only the tris generated in the last step are used
  allocate(cen_sbprt(3,n_init_tri*4**n_steps), rad_sbprt(n_init_tri*4**n_steps))

  if(size(rad_sbprt,1) .ne. n_sbprt) then
    write(*,*) "ERROR n_sbprt", size(rad_sbprt,1), n_sbprt
    stop
  endif
  
  do i = 1, n_init_tri*4**n_steps
    cen_sbprt(:,i) = (tex(:,1,n_tri_tot-i+1)+tex(:,2,n_tri_tot-i+1)+tex(:,3,n_tri_tot-i+1))/3.0_wp
    rad_sbprt(i) = circumradius(tex(:,:,n_tri_tot-i+1))
  enddo

  deallocate(tex)
  
end subroutine


! ----------------------------------------------------------------------
! Computes radius of circumscribed circle of a triangle 
real(wp) function circumradius(vertices)
  
  real(wp), intent(in) :: vertices(3,3) ! triangle vertices, xyz for three points
  real(wp) :: a, b, c, p, area
  
  a = norm2(vertices(:,1)-vertices(:,2))
  b = norm2(vertices(:,2)-vertices(:,3))
  c = norm2(vertices(:,3)-vertices(:,1))
  
  p = (a+b+c)/2
  
  area = sqrt(p*(p-a)*(p-b)*(p-c)) ! Heron's formula for area
  
  if (area .ge. 1e-9_wp) then
    circumradius = a*b*c/(4*area)    
  else
    circumradius = 1e-9_wp
  endif
    
  return 
  
end function circumradius

! ----------------------------------------------------------------------
subroutine linspace(from, to, array)
  real(wp), intent(in)    :: from, to
  real(wp), intent(out)   :: array(:)

  real(wp)                :: range
  integer                 :: n, i
  
  n = size(array)
  range = to - from

  if (n .eq. 0) return

  if (n .eq. 1) then
    array(1) = from
    return
  end if

  do i=1, n
    array(i) = from + range * real(i - 1, wp)/real(n - 1,wp)
  end do

end subroutine
! ----------------------------------------------------------------------

subroutine geoseries(start_x, end_x, n_elem, r, dcsi_le, dcsi_te) 
  real(wp), intent(in)               :: start_x, end_x
  integer, intent(in)                :: n_elem
  real(wp), intent(in)               :: r
  real(wp), allocatable, intent(out) :: dcsi_le(:), dcsi_te(:) 

  real(wp), allocatable              :: dl(:)
  real(wp)                           :: r_d
  integer :: i
  character(len=*), parameter :: this_sub_name = 'geoseries' 

  allocate(dl(n_elem)); dl = 0.0_wp
  if (.not. allocated(dcsi_te)) allocate(dcsi_te(n_elem + 1)); dcsi_te = 0.0_wp 
  if (.not. allocated(dcsi_le)) allocate(dcsi_le(n_elem + 1)); dcsi_le = 0.0_wp
  
  !> check series convergence 
  if ((r .ge. 1.0_wp) .or. (r .le. 0.0_wp)) then 
    call error(this_mod_name, this_sub_name, 'insert a growth value between 0. and 1.') 
  end if
  
  r_d = 1-r ! growth ratio (intended as decreasing ratio) 
  ! get size of the first element knowing the sum of the geometric series
  dl(1) = (1.0_wp-r_d)/(1.0_wp - r_d**real(n_elem,wp));
  do i = 2, n_elem
    dl(i) = r_d*dl(i-1) ! geometric series-> get length of element  
  end do 
  do i = 1, n_elem
    dcsi_te(i + 1) = dcsi_te(i) + dl(i); !-> position  
  end do

  dcsi_le = (1-dcsi_te)
  !> flip dcsi_le
  dcsi_le = dcsi_le(size(dcsi_le):1:-1)  
  !> rescale in the interval 
  dcsi_te = dcsi_te*(end_x - start_x) + start_x 
  dcsi_le = dcsi_le*(end_x - start_x) + start_x
  !> cleanup 
  if (allocated(dl)) deallocate(dl)

end subroutine geoseries 

subroutine geoseries_both(start_x, end_x, n_elem, xh, r_le, r_te, dcsi)
  real(wp),               intent(in)    :: start_x, end_x
  integer,                intent(in)    :: n_elem
  real(wp),               intent(in)    :: xh
  real(wp),               intent(in)    :: r_le, r_te
  real(wp), allocatable,  intent(out)   :: dcsi(:) 

  real(wp), allocatable                 :: dcsi_le(:), dcsi_te(:) 
  real(wp), allocatable                 :: dcsi_le_(:), dcsi_te_(:) !> dummy
  real(wp)                              :: mid_point, start_x_le, end_x_le, start_x_te, end_x_te 
  integer                               :: n_elem_le, n_elem_te
  character(len=*), parameter           :: this_sub_name = 'geoseries_both' 
  
  allocate(dcsi(n_elem+1)); dcsi = 0.0_wp 

  n_elem_le = ceiling(real(n_elem,wp)*(1.0_wp - xh)) 
  n_elem_te = n_elem - n_elem_le
  mid_point = (start_x + end_x)*xh
  start_x_le = start_x
  end_x_le = mid_point
  start_x_te = mid_point
  end_x_te = end_x
  call geoseries(start_x_le, end_x_le, n_elem_le, r_le, dcsi_le, dcsi_te_)
  call geoseries(start_x_te, end_x_te, n_elem_te, r_te, dcsi_le_, dcsi_te)
  dcsi = (/dcsi_le(1:(size(dcsi_le)-1)), dcsi_te/) 
  !> cleanup
  if (allocated(dcsi_le))  deallocate(dcsi_le)
  if (allocated(dcsi_te))  deallocate(dcsi_te)
  if (allocated(dcsi_le_)) deallocate(dcsi_le_)
  if (allocated(dcsi_te_)) deallocate(dcsi_te_)

end subroutine geoseries_both

subroutine geoseries_hinge(start_x, end_x, n_elem, xh, &
                          r_le_aft, r_te_aft, r_le_fore, r_te_fore, dcsi)
  real(wp),               intent(in)     :: start_x, end_x
  integer,                intent(in)     :: n_elem
  real(wp),               intent(in)     :: xh
  real(wp),               intent(in)     :: r_le_aft, r_te_aft, r_le_fore, r_te_fore
  real(wp), allocatable,  intent(out)    :: dcsi(:)

  real(wp), parameter                 :: x_mid = 0.5_wp !> put as input? 
  real(wp)                            :: start_x_aft, end_x_aft, start_x_fore, end_x_fore
  real(wp), allocatable               :: dcsi_aft(:), dcsi_fore(:) 
  integer                             :: n_elem_aft, n_elem_fore 
  character(len=*), parameter         :: this_sub_name = 'geoseries_hinge' 
  
  allocate(dcsi(n_elem+1)); dcsi = 0.0_wp 
  n_elem_aft = ceiling(real(n_elem,wp)*xh) 
  n_elem_fore = n_elem - n_elem_aft
  start_x_aft = start_x
  end_x_aft = xh
  start_x_fore = xh
  end_x_fore = end_x
  call geoseries_both(start_x_aft, end_x_aft, n_elem_aft, x_mid, r_le_aft, r_te_aft, dcsi_aft);
  call geoseries_both(start_x_fore, end_x_fore, n_elem_fore, x_mid, r_le_fore, r_te_fore, dcsi_fore);

  dcsi = (/dcsi_aft(1:(size(dcsi_aft)-1)), dcsi_fore/)

  !> cleanup
  if (allocated(dcsi_aft))  deallocate(dcsi_aft)
  if (allocated(dcsi_fore)) deallocate(dcsi_fore)  

end subroutine geoseries_hinge


! Modified modulo function to cycle arrays 
! Substitutes 0 with last element
! Eg (1,2,3) +1 mod 3  -> (2,0,1)
!    (1,2,3) +1 modd 3 -> (2,3,1)
integer function modd(a,b)

    integer, intent(in) :: a,b
    
    modd = mod(a,b)
    if (modd .eq. 0) modd = b
    return 
  
end function





end module mod_math
