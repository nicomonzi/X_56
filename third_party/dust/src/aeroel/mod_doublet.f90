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
!!          Alessandro Cocco 
!!=========================================================================

!> Module containing subroutines about doublets distributions.
!!
!! Note that this module does not expose an aerodynamic element of its own,
!! but collects utilities for other aerodynamic elements that contain
!! a doublet distribution
module mod_doublet

use mod_param, only: &
  wp, nl, &
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi, eps, max_char_len

use mod_handling, only: &
  error, warning, printout

use mod_aeroel, only: &
  c_pot_elem

use mod_math, only: &
  cross

use mod_sim_param, only: &
  sim_param

implicit none

public ::    initialize_doublet,  &
          potential_calc_doublet, &
          velocity_calc_doublet,  &
          gradient_calc_doublet,  & 
          linear_potential_calc_doublet

private

real(wp) :: ff_ratio
real(wp) :: eps_dou
real(wp) :: r_Rankine
real(wp) :: r_cutoff

character(len=*), parameter :: this_mod_name = 'mod_doublet' 
contains

!----------------------------------------------------------------------

!> Subroutine to populate the module variables from input
!!
subroutine initialize_doublet()

  ff_ratio  = sim_param%FarFieldRatioDoublet
  eps_dou   = sim_param%DoubletThreshold
  r_Rankine = sim_param%RankineRad
  r_cutoff  = sim_param%CutoffRad

end subroutine initialize_doublet

!----------------------------------------------------------------------

!> compute AIC of panel <this>, on the control point <pos>
!! Compute Omega_ik: the solid angle seen from the point <pos>_i (passive)
!!   looking to the panel <this>_k (active). This can be interpreted as the
!!   velocity potential in <pos>_i induced by a (-4*pi)-intensity surface
!!   doublet on the panel <this>_k.               ^---------
!!
!! Omega_ik = int_{S_k} { n \cdot (r_i - r) / |r_i - r|^3 }
!!
!! The relation with unitary surface doublet D_ik is:
!!   D_ik = -(1/4pi) * Omega_ik
!!
subroutine potential_calc_doublet(this, dou, pos)
  class(c_pot_elem), intent(in) :: this
  real(wp), intent(out)         :: dou
  real(wp), intent(in)          :: pos(:)
  real(wp), dimension(3)        :: e3
  real(wp)                      :: zQ, ap1nxam1n
  real(wp), dimension(3)        :: ei , ap1 , am1
  real(wp)                      :: ap1n , am1n , sinB , cosB , beta
  integer                       :: indm1
  ! QUAD elem as 2 TRIA elems
  real(wp) :: radius
  integer :: i1

  radius = norm2(pos - this%cen)
  ! unit normal
  e3 = this%nor
  ! Control point (Q): distance (normal proj) of the point <pos> from the panel <this>
  zQ = dot_product((pos - this%cen),e3)

  if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1)
    dou = zQ * this%area / radius**3.0_wp
  else
    !> initialisation
    dou = 0.0_wp
    !> old implmentation that did not distinguish TRIA from QUAD
    do i1 = 1 , this%n_ver

      if ( this%n_ver .eq. 3 ) then
        indm1 = prev_tri(i1)
      else if ( this%n_ver .eq. 4 ) then
        indm1 = prev_qua(i1)
      end if
      
      ! doublet  -----
      ! it is possible to use ver, instead of ver_p for the doublet
      ei = - ( pos - this%ver (:,i1) ) ; ei = ei / norm2(ei)
      ap1 =   this%edge_vec(:,i1   ) - ei * dot_product( ei, this%edge_vec(:,i1   ) )
      am1 = - this%edge_vec(:,indm1) + ei * dot_product( ei, this%edge_vec(:,indm1) )
      ap1n= norm2(ap1)
      am1n= norm2(am1)
      ap1nxam1n = ap1n*am1n
      if (ap1nxam1n .lt. eps) then !> avoid singularities 
        sinB = 0.0_wp
        cosB = 1.0_wp
      else
        sinB = dot_product(ei, cross(am1, ap1)) / ap1nxam1n
        cosB = dot_product(am1, ap1) / ap1nxam1n
      end if
      beta = atan2( sinB , cosB )
      dou = dou + beta
    end do

    ! Correct the result to obtain the solid angle (from Gauss-Bonnet theorem)
    !TODO: use "sign" here and check the results
    if     ( dou .lt. -real(this%n_ver-2,wp)*pi + eps ) then
      dou = dou + real(this%n_ver-2,wp) * pi
    elseif ( dou .gt. +real(this%n_ver-2,wp)*pi - eps ) then
      dou = dou - real(this%n_ver-2,wp) * pi
    end if

  end if

end subroutine potential_calc_doublet

!----------------------------------------------------------------------
!> compute velocity AIC of panel <this>, on the control point <pos>
!! Biot-Savart law. Regular kernel: linear core (Rankine vortex)
!! Compute the velocity induced by a vortex ring (equivalent to a constant
!!  intentsity surface doublet) with intensity 4*pi. <------
!!
!! V^{vr}_ik = grad_{r_i} { - int_{S_k} { n \cdot (r_i - r) / |r_i - r|^3 } }
!!
subroutine velocity_calc_doublet(this, v_dou, pos)
  class(c_pot_elem), intent(in) :: this
  real(wp), intent(out)         :: v_dou(3)
  real(wp), intent(in)          :: pos(:)
  real(wp)                      :: phix , phiy , pdou
  integer                       :: indp1 
  real(wp)                      :: av(3) , hv(3)
  real(wp)                      :: ai    , hi
  real(wp)                      :: R1 , R2
  real(wp)                      :: radius_v(3)
  real(wp)                      :: radius , rati
  real(wp)                      :: r_Ran
  integer                       :: i1
  ! WARNING:
  ! to introduce the FLAT-PANEL APPROXIMATION, the projection of the nodes
  ! on the mean palen should be used ...
  ! ti = this%edge_uni(i1)
  ! si = this%edge_len(i1)
  v_dou = 0.0_wp
  radius_v = pos - this%cen
  radius   = norm2(radius_v)
  
  if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1)

    rati = 3.0_wp * this%area * dot_product(radius_v, this%nor) / radius**5.0_wp

    phix = rati * dot_product(radius_v, this%tang(:,1))
    phiy = rati * dot_product(radius_v, this%tang(:,2))
    pdou = this%area * ( - radius**2.0_wp + 3.0_wp*dot_product(radius_v, this%nor)**2.0_wp ) / radius**5.0_wp

    ! vdou = matmul((/tang(:,:), nor/), (/phix, phiy, pdou/)) 
    v_dou(1) = this%tang(1,1)*phix + this%tang(1,2)*phiy + this%nor(1)* pdou
    v_dou(2) = this%tang(2,1)*phix + this%tang(2,2)*phiy + this%nor(2)* pdou
    v_dou(3) = this%tang(3,1)*phix + this%tang(3,2)*phiy + this%nor(3)* pdou

  else

    ! === Original DUST regularisation ===
    do i1 = 1 , this%n_ver

      !> index of next vertex (/2, 3, 4, 1/)
      if ( this%n_ver .eq. 3 ) then
        indp1 = next_tri(i1)
      else if ( this%n_ver .eq. 4 ) then
        indp1 = next_qua(i1)
      end if
      
      !> use this%ver instead of its projection this%verp
      av = pos-this%ver(:,i1)
      ai = dot_product(av, this%edge_uni(:,i1))
      R1 = norm2(av)
      R2 = norm2(pos - this%ver(:,indp1))
      hv = av - ai*this%edge_uni(:,i1)
      hi = norm2(hv)
      
      ! === TEST THRESHOLDS ===
      !> Relative threshold ---
      if ( hi .gt. this%edge_len(i1)*r_Rankine ) then
        v_dou = v_dou + ( (this%edge_len(i1)-ai)/r2 + ai/r1 )/(hi**2.0_wp) * &
                        cross(this%edge_uni(:,i1),hv)
      else
        ! === TEST THRESHOLDS ===
        !> Relative threshold ---
        if ( ( R1 .gt. this%edge_len(i1)*r_cutoff ) .and. &! avoid singularity
            ( R2 .gt. this%edge_len(i1)*r_cutoff )   ) then
          r_Ran = r_Rankine * this%edge_len(i1)
          v_dou = v_dou + ((this%edge_len(i1)-ai)/R2 + ai/R1)/(r_Ran**2.0_wp)* &
                          cross(this%edge_uni(:,i1),hv)
        end if
      end if
    end do
  end if

end subroutine velocity_calc_doublet

!----------------------------------------------------------------------

!> Compute gradient of the velcoity field using the Rosenhead kernel
!>
subroutine gradient_calc_doublet(this, grad_dou, pos)
  class(c_pot_elem), intent(in) :: this
  real(wp), intent(out) :: grad_dou(3,3)
  real(wp), intent(in) :: pos(:)

  integer :: indp1 , i1 , i2 , i
  real(wp) :: R1(3) , R2(3) , a1(3) , a2(3) , l(3) , a
  real(wp) :: R1v(3,1) , R2v(3,1) , a1v(3,1) , a2v(3,1) , lv(3,1)
  real(wp) :: lx(3,3) , aa1(3,3) , aa2(3,3) , ax1(3,3) , ax2(3,3) , al1(3,3) , al2(3,3)
  real(wp) :: del , a2del2

  del = sim_param % RankineRad
  !del = sim_param % VortexRad

  grad_dou = 0.0_wp

  do i = 1 , this%n_ver

    if ( this%n_ver .eq. 3 ) then
      indp1 = next_tri(i)
    else if ( this%n_ver .eq. 4 ) then
      indp1 = next_qua(i)
    end if

    i1 = i
    i2 = indp1

    l = this%edge_uni(:,i)
    lv(:,1) = l

    R1 = pos-this%ver(:,i1) ;   a1 = cross( l , R1 )
    R2 = pos-this%ver(:,i2) ;   a2 = cross( l , R2 )
    a = norm2(a1)  
    a2del2 = a**2.0_wp + del**2.0_wp

    R1v(:,1) = R1 ;  a1v(:,1) = a1
    R2v(:,1) = R2 ;  a2v(:,1) = a2

    lx(:,1) = (/  0.0_wp ,  l(3)   , -l(2)   /)
    lx(:,2) = (/ -l(3)   ,  0.0_wp ,  l(1)   /)
    lx(:,3) = (/  l(2)   , -l(1)   ,  0.0_wp /)

    aa1 = matmul( a1v , transpose(a1v) )
    ax1 = matmul( a1v , transpose(R1v) )
    al1 = matmul( a1v , transpose( lv) )
    aa2 = matmul( a2v , transpose(a2v) )
    ax2 = matmul( a2v , transpose(R2v) )
    al2 = matmul( a2v , transpose( lv) )

    grad_dou = grad_dou &
      + 1.0_wp / ( a2del2**1.5_wp * norm2(R1) ) * &
      ( (   a * lx &
          + ( 1.0_wp/a - 3.0_wp*a/a2del2 ) * matmul(aa1,lx) ) * sum(l*R1) &
        + a * ( al1 - ax1 * sum( l * R1 ) / norm2(R1)**2.0_wp ) )  &
        - 1.0_wp / ( a2del2**1.5_wp * norm2(R2) ) * &
      ( (   a * lx &
          + ( 1.0_wp/a - 3.0_wp*a/a2del2 ) * matmul(aa2,lx) ) * sum(l*R2) &
      + a * ( al2 - ax2 * sum( l * R2 ) / norm2(R2)**2.0_wp ) )

  end do

end subroutine gradient_calc_doublet

!----------------------------------------------------------------------
subroutine linear_potential_calc_doublet(this, TL, TR, pos) 
  class(c_pot_elem), intent(in) :: this
  real(wp), intent(in)          :: pos(:)
  real(wp), intent(out)         :: TL, TR  

  real(wp)                      :: xQ, yQ, zQ
  real(wp)                      :: R1 , R2, Q, v
  real(wp), dimension(3)        :: Qp
  real(wp)                      :: den, sum_xy, phi 
  real(wp), allocatable         :: verp(:,:), cosTi(:), sinTi(:)
  integer                       :: indp1, i1
  character(len=max_char_len)   :: msg(3)  
  character(len=*), parameter   :: this_sub_name='linear_potential_calc_doublet'
  
  ! From Newman 'Distributions of sources and normal dipoles over a quadrilateral panel' Sec. 4
  
  ! Control point (Q): distance (normal proj) of the point <pos> from the panel <this>
  xQ = dot_product(pos-this%cen, this%tang(:,1))
  yQ = dot_product(pos-this%cen, this%tang(:,2))
  ! Control point (Q)
  zQ = dot_product((pos-this%cen) , this%nor)

  Qp = pos - zQ * this%nor 

  ! projection of the vertices on the mean plane 
  allocate(verp(3,this%n_ver)); verp = 0.0_wp 
  allocate(cosTi(this%n_ver), sinTi(this%n_ver)); cosTi = 0.0_wp; sinTi = 0.0_wp 

  do i1 = 1 , this%n_ver
    cosTi(i1) = dot_product( this%edge_uni(:,i1), this%tang(:,1))
    sinTi(i1) = dot_product( this%edge_uni(:,i1), this%tang(:,2))
    verp(:,i1) = this%ver(:,i1) - this%nor * &
                dot_product( (this%ver(:,i1) - this%cen ), this%nor )
  enddo 

  call potential_calc_doublet(this, phi, pos)  
  sum_xy = 0.0_wp;  

  do i1 = 1 , this%n_ver

    if ( this%n_ver .eq. 3 ) then
      indp1 = next_tri(i1)
    else if ( this%n_ver .eq. 4 ) then
      indp1 = next_qua(i1)
    end if

    R1 = norm2(pos - verp(:,i1))
    R2 = norm2(pos - verp(:,indp1))
    den = R1+R2-this%edge_len(i1)

    if(den .lt. 1e-6_wp .and. sim_param%debug_level .ge. 5) then
      write(msg(1),'(A)') 'Too small denominator in &
      &source computation with point projection, using actual &
      &points instead.'
      write(msg(2),'(A,F12.6,F12.6,F12.6)') 'Computing sources on point: ',&
      pos(1),pos(2),pos(3)
      write(msg(3),'(A)')'This is most likely due to severely warped &
      &quadrilateral elements adjacent to small elements.'//nl//&
      &'      === CHECK MESH QUALITY! ==='
      call warning(this_sub_name, this_mod_name, msg)
    endif

    Q = log((R1+R2+this%edge_len(i1))/(den))
    v = dot_product(cross(Qp - verp(:,i1), this%edge_vec(:,i1)), this%nor) / this%edge_len(i1)
    sum_xy = sum_xy + cosTi(i1)*(v*Q*sinTi(i1) - (R2 - R1)*cosTi(i1)) 

  end do
  
  TR = -(phi*xq*yq + zq*sum_xy)
  TL = -phi + TR 

  deallocate(verp, cosTi, sinTi)

end subroutine linear_potential_calc_doublet


end module mod_doublet
