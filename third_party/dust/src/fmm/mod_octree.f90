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
!!          Matteo Tugnoli
!!=========================================================================


!> Module to handle the octree grid
module mod_octree

use mod_param, only: &
  wp, nl, pi, one_4pi, max_char_len

use mod_math, only: &
  cross

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, internal_error

use mod_aeroel, only: &
  t_elem_p, t_pot_elem_p, c_pot_elem, c_elem

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_vortline, only: &
  t_vortline

use mod_multipole, only: &
  t_multipole, t_polyexp, t_ker_der, set_multipole_param

use mod_wind, only: &
  variable_wind
!----------------------------------------------------------------------

implicit none

public :: initialize_octree, destroy_octree, sort_particles, t_octree, &
          calculate_multipole, apply_multipole, apply_multipole_panels, & 
          apply_multipole_optimized 
private

!> Pointer to a cell
type :: t_cell_p
  type(t_cell), pointer :: p => null()
end type

!----------------------------------------------------------------------

!> Type containing all the data relative to a single octree cell
type :: t_cell

  !> Is the cell a leaf?
  logical :: leaf

  !> Is the cell a branch (active and not a leaf)
  logical :: branch

  !> Is the cell active?
  logical :: active

  !> Pointer to the parent cell
  type(t_cell_p) :: parent

  !> Pointers to the children cells
  type(t_cell_p) :: children(8)

  !> Pointers to the neighbouring cells
  type(t_cell_p) :: neighbours(-1:1,-1:1,-1:1)

  !> Pointers to far leaf neighbours
  type(t_cell_p), allocatable :: leaf_neigh(:)

  !> Pointers to the cells in the interaction list
  !! (i.e. child of the parent siblings which are not neighbouring)
  type(t_cell_p), allocatable :: interaction_list(:)

  !> Indices in the implicit cartesian frame of the cell inside the layer
  integer :: cart_index(3)

  !> Level of the cell in the layers
  !integer :: level

  !> Number of particles contained in each cell
  integer :: npart

  !> Pointers to all the particles contained in the cell
  type(t_vortpart_p), allocatable :: cell_parts(:)

  !> Number of panels (centres) contained in the cell
  integer :: npan

  !> Poiners to all the panels contained in the cell
  type(t_pot_elem_p), allocatable :: cell_pans(:)

  !> Is it near to a solid wall?
  logical :: near_wall

  !> Center of the cell
  real(wp) :: cen(3)

  !> Multipole data relative to the cell
  type(t_multipole) :: mp

  !> Average vorticity magnitude
  real(wp) :: ave_vortmag

end type

!----------------------------------------------------------------------

!> Type containing the pointers to a layer in the octree
type :: t_cell_layer

  !> Size of all the cubic cells in the layer
  real(wp) :: cell_size

  !> All the cells contained in the layer, the main allocation happens here
  type(t_cell), allocatable :: lcells(:,:,:)

  !> All the possible derivatives of the kernel among all the possible
  !! different relative locations for interaction at the current level
  type(t_ker_der), allocatable :: ker_der(:,:,:)

  !length (positive and negative) of kernel derivatives in each direction
  integer :: ker_der_len(3)

end type


!----------------------------------------------------------------------

!> Type containing the whole octree structure
type :: t_octree

  !> maximum and minimum of the coordinates of the whole octree
  real(wp) :: xmin(3), xmax(3)

  !> number of level 1 cubic boxes in each direction
  integer :: nbox(3)

  !> number of levels of octree division
  integer :: nlevels

  !> Level at which the solid walls are considered (for particles redistribution)
  integer :: lvl_solid

  !> total number of cells (sum of cells at each level)
  integer :: ncells_tot

  !> number of cells in each direction, at each level
  integer, allocatable :: ncl(:,:)

  !> all the actual cells, in a 1d array
  !type(t_cell), allocatable :: cells(:)

  !> cell pointers organized in layers, each layer corresponding to
  !! a level of the octree
  type(t_cell_layer), allocatable :: layers(:)

  !> pointer array to all the leaves
  type(t_cell_p), allocatable :: leaves(:)

  !> number of leaves
  integer :: nleaves

  !> Degree of multipole expansion
  integer :: degree

  !> Polynomial expansion tools
  type(t_polyexp) :: pexp
  type(t_polyexp) :: pexp_der

  !> Regularized kernel radius
  real(wp) :: delta

  !> Time spent in multipole calculations
  real(t_realtime) :: t_mp

  !> Time spent in leaves calculations
  real(t_realtime) :: t_lv

end type

!----------------------------------------------------------------------
! Interfaces
interface push_ptr
  module procedure push_cell_ptr_scalar, push_cell_ptr_vec
  module procedure push_part_ptr_scalar, push_part_ptr_vec
  module procedure push_pan_ptr_scalar,  push_pan_ptr_vec
end interface

interface pull_ptr
  module procedure pull_part_ptr_scalar
end interface

!----------------------------------------------------------------------
character(len=*), parameter :: this_mod_name='mod_octree'
character(len=max_char_len) :: msg
real(t_realtime) :: t1 , t0
real(t_realtime) :: t11, t00

integer :: min_part_4_cell

integer, parameter :: allocation_block = 32
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the octree
subroutine initialize_octree(box_length, nbox, origin, nlevels, min_part, &
                              degree, delta, octree)
  real(wp), intent(in) :: box_length
  integer,  intent(in) :: nbox(3)
  real(wp), intent(in) :: origin(3)
  integer,  intent(in) :: nlevels
  integer,  intent(in) :: min_part
  integer,  intent(in) :: degree
  real(wp), intent(in) :: delta
  type(t_octree), intent(out), target :: octree

  type(t_cell_p) :: inter_list(208) !all the possible interactions
  integer :: l, i,j,k, ic,jc,kc, ip_list
  integer :: imax, jmax, kmax
  integer :: p, child
  integer :: indx(3)
  character(len=*), parameter :: this_sub_name = 'initialize_octree'

  octree%xmin = origin
  octree%xmax = origin + real(nbox,wp)*box_length
  octree%nbox = nbox
  octree%nlevels = nlevels
  octree%delta = delta

  !Check the particles bounding
  if (any(octree%xmax-sim_param%particles_box_max .lt. 0.0_wp) .or. &
      any(sim_param%particles_box_min-octree%xmin .lt. 0.0_wp)) then
    where(sim_param%particles_box_max .gt. octree%xmax) &
          sim_param%particles_box_max  =   octree%xmax
    where(sim_param%particles_box_min .lt. octree%xmin) &
          sim_param%particles_box_min  =   octree%xmin
    call warning(this_sub_name, this_mod_name, 'Particles box bigger than octree box, particles box resized to the octree box')
  endif

  octree%lvl_solid = sim_param%lvl_solid

  call set_multipole_param(sim_param%use_vs)

  min_part_4_cell = min_part

  octree%degree = degree
  !polynomial tools at the multipole degree
  call octree%pexp%set_degree(degree)
  !for the derivatives is needed 2*degree+1 degree
  call octree%pexp_der%set_degree(2*degree+2)

  !Count all the cells contained in the octree
  octree%ncells_tot = product(nbox)
  do l = 2,nlevels
    octree%ncells_tot = octree%ncells_tot + product(nbox)*8**(l-1)
  enddo

  !Allocate all the possible leaves (to avoid pushing the pointer vector many
  ! times), then only the first nleaves will be employed
  allocate(octree%leaves(product(nbox)*8**(nlevels-1)))

  if(sim_param%use_dyn_layers) then
    allocate(octree%layers(sim_param%NMaxOctreeLevels))
  else
    allocate(octree%layers(nlevels))
  endif

  !compute the number of cells in each direction, for each layer
  allocate(octree%ncl(3,sim_param%NMaxOctreeLevels))
  do l=1,nlevels
    octree%ncl(1,l) = nbox(1)*2**(l-1);
    octree%ncl(2,l) = nbox(2)*2**(l-1);
    octree%ncl(3,l) = nbox(3)*2**(l-1);
  enddo

  ! == Prepare the hierarchical tree

  !first level
  octree%layers(1)%cell_size = box_length
  allocate(octree%layers(1)%lcells(nbox(1),nbox(2),nbox(3)))

  !cycle on all the cells on the first level
  do k = 1,nbox(3); do j = 1,nbox(2); do i = 1,nbox(1)
    !Set the quantities for each cell
    octree%layers(1)%lcells(i,j,k)%cart_index = (/i,j,k/)
    !octree%layers(1)%lcells(i,j,k)%level = 1
    octree%layers(1)%lcells(i,j,k)%cen = octree%xmin + &
                (real((/i,j,k/),wp)-0.5_wp)*(octree%xmax-octree%xmin) &
                /real(nbox,wp)
  enddo; enddo; enddo

  !dive down  all the other levels
  do l = 2,nlevels
    !Set the quantities for the layer
    octree%layers(l)%cell_size = box_length/real((2**(l-1)),wp)
    allocate(octree%layers(l)%&
                        lcells(octree%ncl(1,l),octree%ncl(2,l),octree%ncl(3,l)))
    !cycle on the parents (i.e. all the cells at the upper level)
    do k = 1,octree%ncl(3,l-1)
    do j = 1,octree%ncl(2,l-1)
    do i = 1,octree%ncl(1,l-1)

      p = 1 !reset the counter for the parent to set the children
      do kc = 1,2; do jc = 1,2; do ic = 1,2
        !point the cell in the ordered layer
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cart_index &
                              = (/2*(i-1)+ic, 2*(j-1)+jc, 2*(k-1)+kc/)
    !    octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%level = l

        !set the parent
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%parent%p &
                                      => octree%layers(l-1)%lcells(i,j,k)
    !    octree%layers(l)%lcells(i,j,k)%level = l
        !set the present cell as children of the parent
        octree%layers(l-1)%lcells(i,j,k)%children(p)%p => &
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)

        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cen = &
                octree%xmin + &
                (real((/2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc/),wp)-0.5_wp)* &
                (octree%xmax-octree%xmin)/real(nbox*2**(l-1),wp)

        !update the counters
        p = p+1

      enddo; enddo; enddo !children ic,jc,kc
    enddo; enddo; enddo !parents i,j,k
  enddo

  !PROFILE
  t0 = dust_time()
  ! == Set the neighbours and perform initializations
  do l = 1,nlevels

    imax=octree%ncl(1,l); jmax=octree%ncl(2,l); kmax=octree%ncl(3,l);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      !cycle on all the possible neighbours
      do kc = -1,1; do jc = -1,1; do ic = -1,1

        if( ((i+ic).ge.1 .and. (i+ic).le.imax) .and. &
            ((j+jc).ge.1 .and. (j+jc).le.jmax) .and. &
            ((k+kc).ge.1 .and. (k+kc).le.kmax) .and. &
            .not.(ic.eq.0 .and. jc.eq.0 .and. kc.eq.0) ) then
          !if it is inside the domain, set the pointer to the neighbour
          octree%layers(l)%lcells(i,j,k)%neighbours(ic,jc,kc)%p =>  &
              octree%layers(l)%lcells(i+ic,j+jc,k+kc)
        endif

      enddo; enddo; enddo !neighbours ic,jc,kc

      !initialize few things for each cell
      !(most of these actions are performed at the beginning of each timestep
      ! however rely on de-allocating and re-allocating without check for
      ! speed, so things are allocated to zero here)
      allocate(octree%layers(l)%lcells(i,j,k)%cell_parts(0))
      allocate(octree%layers(l)%lcells(i,j,k)%leaf_neigh(0))
      octree%layers(l)%lcells(i,j,k)%npart = 0
      allocate(octree%layers(l)%lcells(i,j,k)%cell_pans(0))
      octree%layers(l)%lcells(i,j,k)%npan = 0
      octree%layers(l)%lcells(i,j,k)%active = .false.
      octree%layers(l)%lcells(i,j,k)%leaf = .false.
      allocate(octree%layers(l)%lcells(i,j,k)%interaction_list(0) )
      !call octree%layers(l)%lcells(i,j,k)%mp%init(octree%pexp)

    enddo; enddo; enddo !layer cells i,j,k
  enddo
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Initialized neighbours in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)



  ! == Set the interaction list
  !first level: interaction list is all the elements which are not neighbours
  imax=nbox(1); jmax=nbox(2); kmax=nbox(3);
  !cycle on the elements on the level
  do k = 1,kmax; do j = 1,jmax; do i = 1,imax
    !cycle again on all the cells on first level
    do kc = 1,kmax; do jc = 1,jmax; do ic = 1,imax
          indx = (/ic,jc,kc/)
          if(any(abs(indx - octree%layers(1)%lcells(i,j,k)%cart_index)&
                                                                    .gt.1)) then
            !If it is not a neighbour, push it in the interaction list
            call push_ptr(octree%layers(1)%lcells(i,j,k)%interaction_list, &
            octree%layers(1)%lcells(ic,jc,kc))
          endif
    enddo; enddo; enddo !layer cells i,j,k
  enddo; enddo; enddo !layer cells i,j,k

  !PROFILE
  t0 = dust_time()
  !second level and downwards: the interaction list are all the children of
  !parent which are not a neighbour of the cell
  do l = 2,nlevels
    !cycle on the elements on the level
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      !cycle on all the neighbours of the parent
      ip_list = 0
!DIR$ IVDEP
      do kc = -1,1; do jc = -1,1; do ic = -1,1
        if(associated(octree%layers(l)%lcells(i,j,k)%parent%p%&
                                                  neighbours(ic,jc,kc)%p)) then
          !if the neighbour of the parent exists, cycle on all the childs
          do child = 1,8
            !Index of the child of the parent neighbour
            indx = octree%layers(l)%lcells(i,j,k)%parent%p%&
                            neighbours(ic,jc,kc)%p%children(child)%p%cart_index
            if(any(abs(indx - octree%layers(l)%lcells(i,j,k)%cart_index)&
                                                                    .gt.1)) then
              !If the parent neighbour child is not a neighbour =)
              !push it in the interaction list
              !call push_ptr(octree%layers(l)%lcells(i,j,k)%interaction_list, &
              !       octree%layers(l)%lcells(i,j,k)%parent%p%&
              !       neighbours(ic,jc,kc)%p%children(child)%p)
              ip_list = ip_list + 1
              inter_list(ip_list)%p => octree%layers(l)%lcells(i,j,k)%parent%p%&
                      neighbours(ic,jc,kc)%p%children(child)%p

            endif
          enddo !child
        endif
      enddo; enddo; enddo !parents neighbours ic,jc,kc


        call push_ptr(octree%layers(l)%lcells(i,j,k)%interaction_list, &
                inter_list(1:ip_list))

    enddo; enddo; enddo !layer cells i,j,k
  enddo
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Initialized interaction lists in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)


  !pre-build all the kernel derivatives for cell-cell interactions
  !first level can have many cells, must be treated separately

  octree%layers(1)%ker_der_len = max(octree%nbox,(/3,3,3/))
  allocate(octree%layers(1)%ker_der&
       (-octree%layers(1)%ker_der_len(1):octree%layers(1)%ker_der_len(1),&
        -octree%layers(1)%ker_der_len(2):octree%layers(1)%ker_der_len(2),&
        -octree%layers(1)%ker_der_len(3):octree%layers(1)%ker_der_len(3)))

  do k=-octree%layers(1)%ker_der_len(3),octree%layers(1)%ker_der_len(3)
  do j=-octree%layers(1)%ker_der_len(2),octree%layers(1)%ker_der_len(2)
  do i=-octree%layers(1)%ker_der_len(1),octree%layers(1)%ker_der_len(1)
    if(any(abs((/i,j,k/)) .ge. 2)) then
      call octree%layers(1)%ker_der(i,j,k)%&
        compute_der(octree%layers(1)%cell_size*real((/-i,-j,-k/),wp), &
        octree%delta, octree%pexp_der)
    endif
  enddo; enddo; enddo


  ! In all the other levels there can't be more than +-3 cells of distance
  ! in the interaction list
  do l=2,nlevels
    allocate(octree%layers(l)%ker_der(-3:3,-3:3,-3:3))
    !cycle on all the possible interactions (interaction list cannot go farther
    !than +-3 cells in the same layer)
    do k=-3,3; do j=-3,3; do i=-3,3
      if(any(abs((/i,j,k/)) .ge. 2)) then
        call octree%layers(l)%ker_der(i,j,k)%&
        compute_der(octree%layers(l)%cell_size*real((/-i,-j,-k/),wp), &
        octree%delta, octree%pexp_der)
      endif
    enddo; enddo; enddo
  enddo

end subroutine initialize_octree

!----------------------------------------------------------------------
!> Initialize the octree
subroutine add_layer(octree)
 type(t_octree), intent(inout), target :: octree

 type(t_cell_p) :: inter_list(208) !all the possible interactions
 integer :: ll, i,j,k, ic,jc,kc, ip_list!, il
 integer :: imax, jmax, kmax
 integer :: p, child
 integer :: indx(3)
 !type(t_cell_p), allocatable :: leaves_temp(:)
 !type(t_cell_layer), allocatable :: layers_temp(:)

  octree%nlevels = octree%nlevels+1
  ll = octree%nlevels
  write(msg,'(A,I0,A)') 'Adding an octree level!! Now there are ',ll,&
                          ' levels.'
  if(sim_param%debug_level.ge.5) call printout(msg)


  !add the last level to the number of cells
  octree%ncells_tot = octree%ncells_tot + product(octree%nbox)*8**(ll-1)

  !extend leaves vector
  !TODO: check a way to do it less ugly with move allocs
  !allocate(leaves_temp(product(octree%nbox)*8**(octree%nlevels-1)))
  deallocate(octree%leaves)
  allocate(octree%leaves(product(octree%nbox)*8**(octree%nlevels-1)))
  !do il = 1,size(octree%leaves)
  !  octree%leaves(il) = leaves_temp(il)
  !enddo
  !deallocate(leaves_temp)

  !extend the layers vector
  !TODO: check a way to do it less ugly with move allocs
  !allocate(layers_temp(octree%nlevels))
  !do il = 1,size(octree%layers)
  !  layers_temp(il) = octree%layers(il)
  !enddo
  !deallocate(octree%layers); allocate(octree%layers(octree%nlevels))
  !do il = 1,ll-1
  !  octree%layers(il) = layers_temp(il)
  !enddo
  !deallocate(layers_temp)

  !compute the number of cells in each direction, for last layer
  octree%ncl(1,ll) = octree%nbox(1)*2**(ll-1);
  octree%ncl(2,ll) = octree%nbox(2)*2**(ll-1);
  octree%ncl(3,ll) = octree%nbox(3)*2**(ll-1);

  ! == Prepare the hierarchical tree

  !Set the quantities for the  last layer
  octree%layers(ll)%cell_size = octree%layers(ll-1)%cell_size/2.0_wp
  allocate(octree%layers(ll)%&
                  lcells(octree%ncl(1,ll),octree%ncl(2,ll),octree%ncl(3,ll)))
  !cycle on the parents (i.e. all the cells at the upper level)
  do k = 1,octree%ncl(3,ll-1)
  do j = 1,octree%ncl(2,ll-1)
  do i = 1,octree%ncl(1,ll-1)

    p = 1 !reset the counter for the parent to set the children
    do kc = 1,2; do jc = 1,2; do ic = 1,2
      !point the cell in the ordered layer
      octree%layers(ll)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cart_index &
                            = (/2*(i-1)+ic, 2*(j-1)+jc, 2*(k-1)+kc/)
    !  octree%layers(ll)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%level = ll

      !set the parent
      octree%layers(ll)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%parent%p &
                                    => octree%layers(ll-1)%lcells(i,j,k)
      ! octree%layers(ll)%lcells(i,j,k)%level = ll
      !set the present cell as children of the parent
      octree%layers(ll-1)%lcells(i,j,k)%children(p)%p => &
      octree%layers(ll)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)

      octree%layers(ll)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cen = &
              octree%xmin + &
              (real((/2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc/),wp)-0.5_wp)* &
              (octree%xmax-octree%xmin)/real(octree%nbox*2**(ll-1),wp)

      !update the counters
      p = p+1

    enddo; enddo; enddo !children ic,jc,kc
  enddo; enddo; enddo !parents i,j,k

  ! == Set the neighbours and perform initializations
  imax=octree%ncl(1,ll); jmax=octree%ncl(2,ll); kmax=octree%ncl(3,ll);
  !cycle on the elements on the level
  do k = 1,kmax; do j = 1,jmax; do i = 1,imax
    !cycle on all the possible neighbours
    do kc = -1,1; do jc = -1,1; do ic = -1,1

      if( ((i+ic).ge.1 .and. (i+ic).le.imax) .and. &
          ((j+jc).ge.1 .and. (j+jc).le.jmax) .and. &
          ((k+kc).ge.1 .and. (k+kc).le.kmax) .and. &
          .not.(ic.eq.0 .and. jc.eq.0 .and. kc.eq.0) ) then
        !if it is inside the domain, set the pointer to the neighbour
        octree%layers(ll)%lcells(i,j,k)%neighbours(ic,jc,kc)%p =>  &
            octree%layers(ll)%lcells(i+ic,j+jc,k+kc)
      endif

    enddo; enddo; enddo !neighbours ic,jc,kc

    !initialize few things for each cell
    !(most of these actions are performed at the beginning of each timestep
    ! however rely on de-allocating and re-allocating without check for
    ! speed, so things are allocated to zero here)
    allocate(octree%layers(ll)%lcells(i,j,k)%cell_parts(0))
    octree%layers(ll)%lcells(i,j,k)%npart = 0
    allocate(octree%layers(ll)%lcells(i,j,k)%cell_pans(0))
    octree%layers(ll)%lcells(i,j,k)%npan = 0
    octree%layers(ll)%lcells(i,j,k)%active = .true.
    octree%layers(ll)%lcells(i,j,k)%leaf = .false.
    allocate(octree%layers(ll)%lcells(i,j,k)%interaction_list(0) )
    call octree%layers(ll)%lcells(i,j,k)%mp%init(octree%pexp)

  enddo; enddo; enddo !layer cells i,j,k



  ! == Set the interaction list
  !the interaction list are all the children of
  !parent which are not a neighbour of the cell
  !cycle on the elements on the level
  do k=1,octree%ncl(3,ll); do j=1,octree%ncl(2,ll); do i = 1,octree%ncl(1,ll)
    !cycle on all the neighbours of the parent
    ip_list = 0
    do kc = -1,1; do jc = -1,1; do ic = -1,1
      if(associated(octree%layers(ll)%lcells(i,j,k)%parent%p%&
                                                neighbours(ic,jc,kc)%p)) then
        !if the neighbour of the parent exists, cycle on all the childs
        do child = 1,8
          !Index of the child of the parent neighbour
          indx = octree%layers(ll)%lcells(i,j,k)%parent%p%&
                          neighbours(ic,jc,kc)%p%children(child)%p%cart_index
          if(any(abs(indx - octree%layers(ll)%lcells(i,j,k)%cart_index)&
                                                                  .gt.1)) then
            !If the parent neighbour child is not a neighbour =)
            !push it in the interaction list
            ip_list = ip_list + 1
            inter_list(ip_list)%p => octree%layers(ll)%lcells(i,j,k)%parent%p%&
                    neighbours(ic,jc,kc)%p%children(child)%p

          endif
        enddo !child
      endif
    enddo; enddo; enddo !parents neighbours ic,jc,kc

      call push_ptr(octree%layers(ll)%lcells(i,j,k)%interaction_list, &
              inter_list(1:ip_list))

  enddo; enddo; enddo !layer cells i,j,k


  !pre-build all the kernel derivatives for cell-cell interactions
  allocate(octree%layers(ll)%ker_der(-3:3,-3:3,-3:3))
  !cycle on all the possible interactions (interaction list cannot go farther
  !than +-3 cells in the same layer)
  do k=-3,3; do j=-3,3; do i=-3,3
    if(any(abs((/i,j,k/)) .ge. 2)) then
      call octree%layers(ll)%ker_der(i,j,k)%&
      compute_der(octree%layers(ll)%cell_size*real((/-i,-j,-k/),wp), &
      octree%delta, octree%pexp_der)
    endif
  enddo; enddo; enddo

end subroutine add_layer

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_octree(octree)
  type(t_octree), intent(out) :: octree


end subroutine

!----------------------------------------------------------------------

!> Sort particles inside the octree grid
subroutine sort_particles(wparts, n_prt, elem, octree)
  type(t_vortpart), intent(inout), target :: wparts(:)
  integer, intent(inout) :: n_prt
  type(t_pot_elem_p), intent(in) :: elem(:)
  type(t_octree), intent(inout), target :: octree

  integer :: ip, ipp
  integer :: idx(3)
  real(wp) :: csize
  integer :: l, i,j,k, child, ic,jc,kc
  integer :: ll, nl
  logical :: got_leaves, got_branches
  integer :: nprt, iq
  real(wp) :: redistr(3), vort(3)
  character(len=*), parameter :: this_sub_name = 'sort_particles'

  !PROFILE
  t0 = dust_time()

  ll = octree%nlevels
  csize = octree%layers(ll)%cell_size


  !set all the cells to empty
  do l = 1,ll
    !cycle on the elements on the level
!$omp parallel do collapse(3) private(i,j,k) schedule(dynamic)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      call reset_cell(octree%layers(l)%lcells(i,j,k))
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  enddo

  !to stay on the safe side, nullify all the leaves pointers
!$omp parallel do private(l)
  do l=1,size(octree%leaves)
    nullify(octree%leaves(l)%p)
  enddo
!$omp end parallel do

  !cycle on all the particles
  !TODO: make this parallel too, but pay attention to the race conditions
  ipp = 0
  do ip=1,n_prt
    ipp = ipp + 1
    do while(wparts(ipp)%free)
        ipp = ipp +1
        if(ipp .gt. size(wparts)) &
        call internal_error(this_sub_name, this_mod_name, &
        'not enough non-free particles while sorting')
    enddo
    !check in which cell at the lowest level it is located
    idx = ceiling((wparts(ipp)%cen-octree%xmin)/csize)

    !check that we are not sorting things outside the octree
    if(any(idx.le.0).or.any(idx.gt.shape(octree%layers(ll)%lcells))) then
      write(msg,*) 'Sorted particle resulted outside octree box at: ',&
                    wparts(ipp)%cen
      call error(this_sub_name, this_mod_name, msg)
    endif

    !add the particle to the lowest level
    octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npart = &
                    octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npart + 1
    call push_ptr(octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%cell_parts, &
                  wparts(ipp))
    octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%ave_vortmag =  &
      octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%ave_vortmag + &
      wparts(ipp)%mag

  enddo


  !Sort also all the elements
  if(sim_param%use_pr .or. sim_param%use_fmm_pan) then
    do ip=1,size(elem)
      !check in which cell at the lowest level it is located
      idx = ceiling((elem(ip)%p%cen-octree%xmin)/csize)

      !check that we are not sorting things outside the octree
      if(any(idx.le.0).or.any(idx.gt.shape(octree%layers(ll)%lcells))) then
        write(msg,*) 'Sorted element resulted outside octree box at: ',&
                      elem(ip)%p%cen
        if(sim_param%use_fmm_pan) then
          call error(this_sub_name, this_mod_name, msg)
        else
          if(sim_param%debug_level.ge.5) &
            call warning(this_sub_name, this_mod_name, msg)
          cycle
        endif
      endif

      !add the particle to the lowest level
      octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npan = &
                      octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npan + 1
      call push_ptr(octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%cell_pans&
                    , elem(ip)%p)

    enddo
    !From the second to last level upwards, add the panels
    do l = ll-1,1,-1
      !cycle on the elements on the level
!!!!!$omp parallel do collapse(3) private(i,j,k,child)
      do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)

        !cycle on the children, gather the number of particles and if are
        !leaves
        do child = 1,8
          octree%layers(l)%lcells(i,j,k)%npan = &
                octree%layers(l)%lcells(i,j,k)%npan + &
                octree%layers(l)%lcells(i,j,k)%children(child)%p%npan
          !push all the panels pointers
          call push_ptr(octree%layers(l)%lcells(i,j,k)%cell_pans, &
          octree%layers(l)%lcells(i,j,k)%children(child)%p%cell_pans)
        enddo !child
      enddo; enddo; enddo !layer cells i,j,k
!!!!!!$omp end parallel do
    enddo
  endif




  ! --- Particles Redistribution ---
  if(sim_param%use_pr) then

    !Check if the cells contain a solid boundary at level lvl_solid
    l = octree%lvl_solid
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      !Check if the cell has solid boundaries in it
      if (octree%layers(l)%lcells(i,j,k)%npan .gt. 0) &
        call set_near_wall(octree%layers(l)%lcells(i,j,k))
      !Check also that no neighbour has solid boundaries
      do kc = -1,1; do jc = -1,1; do ic = -1,1
        if(associated(octree%layers(l)%lcells(i,j,k)%&
                                                neighbours(ic,jc,kc)%p)) then
          if(octree%layers(l)%lcells(i,j,k)%neighbours(ic,jc,kc)%p%npan.gt.0) &
            call set_near_wall(octree%layers(l)%lcells(i,j,k))
        endif
      enddo; enddo; enddo !layer cells in,jn,kn
    enddo; enddo; enddo !layer cells i,j,k

    !On the bottom level remove the too small particles
    do k=1,octree%ncl(3,ll); do j=1,octree%ncl(2,ll); do i = 1,octree%ncl(1,ll)
      if(.not. octree%layers(ll)%lcells(i,j,k)%near_wall) then
        nprt = octree%layers(ll)%lcells(i,j,k)%npart
        if (nprt .gt. 0) then
          octree%layers(ll)%lcells(i,j,k)%ave_vortmag =  &
              octree%layers(ll)%lcells(i,j,k)%ave_vortmag/real(nprt,wp)
        endif

        ip = 1
        !do ip=1,nprt
        do while (ip.lt.nprt+1)
          if(octree%layers(ll)%lcells(i,j,k)%cell_parts(ip)%p%mag .lt. &
             1.0_wp/sim_param%part_redist_ratio *&
              octree%layers(ll)%lcells(i,j,k)%ave_vortmag) then
            !the particle is significantly smaller than the average in the cell
            !redistribute it
            redistr = (octree%layers(ll)%lcells(i,j,k)%cell_parts(ip)%p%mag*&
                        octree%layers(ll)%lcells(i,j,k)%cell_parts(ip)%p%dir)&
                        / real(nprt,wp)
            !free
            octree%layers(ll)%lcells(i,j,k)%cell_parts(ip)%p%free = .true.
            call pull_ptr(octree%layers(ll)%lcells(i,j,k)%cell_parts, ip)
            octree%layers(ll)%lcells(i,j,k)%npart = &
                                      octree%layers(ll)%lcells(i,j,k)%npart - 1
            nprt = nprt -1
            !lower also the global counter
            n_prt = n_prt -1

            do iq=1,nprt
                if(.not. octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%free) then
                  vort =  octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%mag*&
                          octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%dir
                  vort = vort + redistr
                  octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%mag = norm2(vort)
                  if(norm2(vort) .gt. 0.0_wp) then
                    octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%dir = &
                                                     vort/norm2(vort)
                  else
                    octree%layers(ll)%lcells(i,j,k)%cell_parts(iq)%p%dir = 0.0_wp
                  endif
                endif !not free particle
            enddo !iq, rest of the particles
          else
            ip = ip + 1
          endif !too small particle
        enddo !ip, particles in cell
      endif !near wall
    enddo; enddo; enddo !layer cells i,j,k
  endif !particles redistribution

  nl = 0
  !Bottom level: just check if are leaves
!$omp parallel do collapse(3) private(i,j,k)
  do k=1,octree%ncl(3,ll); do j=1,octree%ncl(2,ll); do i = 1,octree%ncl(1,ll)
    if( (octree%layers(ll)%lcells(i,j,k)%npart .ge. min_part_4_cell) .or. &
         ll .eq. 1) then !if last level is first, all are leaves
      !mark as leaf
      call set_leaf(octree%layers(ll)%lcells(i,j,k), octree%leaves, nl, &
                                                             octree%pexp)
    endif
  enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do


  !From the second to last level upwards
  do l = ll-1,1,-1
    !cycle on the elements on the level
!$omp parallel do collapse(3) private(i,j,k,got_leaves, got_branches, child)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)

      !cycle on the children, gather the number of particles and if are
      !leaves
      got_leaves = .false.
      got_branches = .false.
      do child = 1,8
        if(octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf) &
                                                          got_leaves = .true.
        if(octree%layers(l)%lcells(i,j,k)%children(child)%p%branch) &
                                                        got_branches = .true.
        octree%layers(l)%lcells(i,j,k)%npart = &
              octree%layers(l)%lcells(i,j,k)%npart + &
              octree%layers(l)%lcells(i,j,k)%children(child)%p%npart
        !push all the particles pointers
        call push_ptr(octree%layers(l)%lcells(i,j,k)%cell_parts, &
        octree%layers(l)%lcells(i,j,k)%children(child)%p%cell_parts)
      enddo !child

      !check how we need to flag the cell
      if(got_leaves .or. got_branches) then
        !Some of the children are leaves or branches, all which is not
        !set becomes leaf, the cell becomes branch
        do child = 1,8
          if(.not. octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf &
          .and. .not. octree%layers(l)%lcells(i,j,k)%children(child)&
                                                          %p%branch) then
            call set_leaf(octree%layers(l)%lcells(i,j,k)%children(child)%p, &
                      octree%leaves, nl, octree%pexp)
          endif
        enddo
        call set_branch(octree%layers(l)%lcells(i,j,k), octree%pexp)
      else
        !None of the children are leaves nor branches
        !Particles have alerady been inherited, just flag the rest inactive
        do child = 1,8
          call set_inactive(octree%layers(l)%lcells(i,j,k)%children(child)%p)
        enddo
        !then check if itself it is a leaf
        if( (octree%layers(l)%lcells(i,j,k)%npart .ge. min_part_4_cell) &
        .or. l .eq. 1) then !if last level is first, all are leaves
          call set_leaf(octree%layers(l)%lcells(i,j,k), &
                      octree%leaves, nl, octree%pexp)
        endif
      endif
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  enddo

  !Last dive down the levels to set the leaf neighbours
  do l = 2,ll
    !cycle on the elements on the level
!$omp parallel do collapse(3) private(i,j,k,ic,kc,jc) schedule(dynamic)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
    if(octree%layers(l)%lcells(i,j,k)%active) then
      !1) get the leaf neighbours from the parent if it exists
      !if(associated(octree%layers(l)%lcells(i,j,k)%parent%p)) then
        call push_ptr(octree%layers(l)%lcells(i,j,k)%leaf_neigh, &
                      octree%layers(l)%lcells(i,j,k)%parent%p%leaf_neigh)
      !endif
      !2) If one of my parent neighbours is a leaf, push it into my leaf neigh
      do kc = -1,1; do jc = -1,1; do ic = -1,1
        if(associated(octree%layers(l)%lcells(i,j,k)%parent%p%&
                                              neighbours(ic,jc,kc)%p)) then
        if(octree%layers(l)%lcells(i,j,k)%parent%p%&
                                          neighbours(ic,jc,kc)%p%leaf) then
          call push_ptr(octree%layers(l)%lcells(i,j,k)%leaf_neigh, &
            octree%layers(l)%lcells(i,j,k)%parent%p%neighbours(ic,jc,kc)%p)
        endif
        endif
      enddo; enddo; enddo !neighbour cells ic,jc,kc
    endif !active
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  enddo

  octree%nleaves = nl

  !PROFILE
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Sorted particles in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)
end subroutine sort_particles


!----------------------------------------------------------------------

!> Calculate the multipole and local expansion in all the cells
subroutine calculate_multipole(part,octree)
 type(t_vortpart_p), intent(in), target :: part(:)
 type(t_octree), intent(inout) :: octree

 integer :: i, j, k, lv, l, child, il
 integer :: ll
 integer :: idx_diff(3)
 real(t_realtime) :: tsta , tend

  tsta = dust_time()
  ll = octree%nlevels
  !reset multipole data (might be done elsewhere)
  ! amended: now the multipole data is reset when assigned leafs/branches
  !do l = ll,1,-1
  !  !cycle on the elements on the level
  !  do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
  !       call octree%layers(l)%lcells(i,j,k)%mp%reset
  !  enddo; enddo; enddo !layer cells i,j,k
  !enddo


  !For all the leaves, calculate the multipole coefficients
  !PROFILE
  t0 = dust_time()
!$omp parallel do private(lv) schedule(guided)
  do lv = 1, octree%nleaves
    call octree%leaves(lv)%p%mp%leaf_M(&
              octree%leaves(lv)%p%cen, &
              octree%leaves(lv)%p%cell_parts, &
              octree%leaves(lv)%p%npart, &
              octree%pexp  )
  enddo
!$omp end parallel do
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated multipoles in leaves in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)


  !layer by layer (except the last, in which are only leaves or inactive)
  !perform multipole to multipole
  !PROFILE
  t00 = dust_time()
  do l = ll-1,1,-1
    !cycle on the elements on the level
    !PROFILE
    t0 = dust_time()
!$omp parallel do collapse(3) private(i,j,k,child)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      if(octree%layers(l)%lcells(i,j,k)%active .and. &
         octree%layers(l)%lcells(i,j,k)%branch) then
         !if it is active and is a branch perform the M2M from all the
         !children
        do child = 1,8
          call octree%layers(l)%lcells(i,j,k)%mp%M2M(&
               octree%layers(l)%lcells(i,j,k)%cen, &
               octree%layers(l)%lcells(i,j,k)%children(child)%p%mp, &
               octree%layers(l)%lcells(i,j,k)%children(child)%p%cen, &
               octree%pexp )
        enddo !child
      endif
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  enddo !l
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated M2M in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)


  !from the top, interact in the interaction list, then pass to the children
  !(if not already at the lowest level)
  do l = 1,ll

  !PROFILE
  t0 = dust_time()
    !cycle on the elements on the level
!$omp parallel do collapse(3) private(i,j,k,il,idx_diff) schedule(dynamic,16)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      if(octree%layers(l)%lcells(i,j,k)%active) then
         !if it is active perform M2L with all the interaction list
        do il = 1,size(octree%layers(l)%lcells(i,j,k)%interaction_list)
          idx_diff = octree%layers(l)%lcells(i,j,k)%interaction_list(il)%p%cart_index -&
                  (/i,j,k/)
          if(octree%layers(l)%lcells(i,j,k)%interaction_list(il)%p%active) then
          call octree%layers(l)%lcells(i,j,k)%mp%M2L(&
               octree%layers(l)%ker_der(idx_diff(1),idx_diff(2),idx_diff(3)), &
               octree%pexp, &
               octree%pexp_der, &
               octree%layers(l)%lcells(i,j,k)%interaction_list(il)%p%mp)
          endif
        enddo !il
      endif
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,I0,A,F9.3,A)') 'Calculated M2L layer ',l,' in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)

  !PROFILE
  t0 = dust_time()
!$omp parallel do collapse(3) private(i,j,k,child) schedule(dynamic)
    do k=1,octree%ncl(3,l); do j=1,octree%ncl(2,l); do i = 1,octree%ncl(1,l)
      if(octree%layers(l)%lcells(i,j,k)%active .and. l.lt.ll) then
        !if active and not in the bottom layer, bring down the local expansion
        ! to the children
        do child = 1,8
          if(octree%layers(l)%lcells(i,j,k)%children(child)%p%active) then
            call octree%layers(l)%lcells(i,j,k)%children(child)%p%mp%L2L(&
                 octree%layers(l)%lcells(i,j,k)%children(child)%p%cen, &
                 octree%layers(l)%lcells(i,j,k)%mp, &
                 octree%layers(l)%lcells(i,j,k)%cen, &
                 octree%pexp )
          endif
        enddo
      endif
    enddo; enddo; enddo !layer cells i,j,k
!$omp end parallel do
  !PROFILE
  t1 = dust_time()
  write(msg,'(A,I0,A,F9.3,A)') 'Calculated L2L layer ',l,' in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)
  enddo !l
  tend = dust_time()
  octree%t_mp = tend-tsta
  !PROFILE
  t11 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated Multipole in all layers in: ' , t11 - t00,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)

end subroutine calculate_multipole

!----------------------------------------------------------------------

!> Apply the local expansion and calculate interaction, and also calculate
! the interaction from other elements (not particles)
subroutine apply_multipole(part, octree, elem, wpan, wrin, wvort)
  type(t_vortpart_p), intent(in), target :: part(:)
  type(t_octree), intent(inout)          :: octree
  type(t_pot_elem_p), intent(in)         :: elem(:)
  type(t_pot_elem_p), intent(in)         :: wpan(:)
  type(t_pot_elem_p), intent(in)         :: wrin(:)
  class(t_vortline), intent(in)          :: wvort(:)

  integer                     :: i, j, k, lv, ip, ipp, m, ie, iln
  real(wp)                    :: Rnorm2, vel(3), pos(3), v(3), stretch(3), stretch_alone(3), str(3), alpha(3), dir(3)
  real(wp), allocatable       :: velv(:,:), stretchv(:,:), rotuv(:,:)
  real(wp)                    :: grad(3,3)
  real(t_realtime)            :: tsta , tend
  real(wp)                    :: turbvisc, ave_ros
  real(wp)                    :: rotu(3), ru(3)

  real(wp) :: grad_elem(3,3) , stretch_elem(3)

  tsta = dust_time()
  !for all the leaves apply the local expansion and then local interactions
  t0 = dust_time()
  
!$omp parallel do private(lv, ip, ie, vel, pos, m, i, j, k, ipp, Rnorm2, &
!$omp& v, stretch, stretch_alone, str, grad, alpha, dir, velv, stretchv, ave_ros, &
!$omp& grad_elem, stretch_elem, iln, rotuv, rotu, ru) schedule(dynamic)
  do lv = 1, octree%nleaves
    !I am on a leaf, cycle on all the particles inside the leaf
    allocate(velv(3,octree%leaves(lv)%p%npart),&
              stretchv(3,octree%leaves(lv)%p%npart), &
              rotuv(3,octree%leaves(lv)%p%npart))

    ave_ros = 0.0_wp

!$omp simd private(vel, stretch, stretch_alone, str, grad, pos, alpha, dir, m) reduction(+:ave_ros)
    do ip = 1,octree%leaves(lv)%p%npart
    if(.not. octree%leaves(lv)%p%cell_parts(ip)%p%free) then
      vel = 0.0_wp
      stretch = 0.0_wp
      stretch_alone = 0.0_wp
      grad = 0.0_wp
      rotu = 0.0_wp
      pos = octree%leaves(lv)%p%cell_parts(ip)%p%cen
      alpha = octree%leaves(lv)%p%cell_parts(ip)%p%mag * &
              octree%leaves(lv)%p%cell_parts(ip)%p%dir
      dir = octree%leaves(lv)%p%cell_parts(ip)%p%dir

      !first apply the local multipole expansion
      do m = 1,size(octree%leaves(lv)%p%mp%b,2)
        !velocity
        vel = vel + &
              octree%leaves(lv)%p%mp%b(:,m)* &
              product((pos-octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
        if(sim_param%use_vs) then
          !stretching
          grad = grad + octree%leaves(lv)%p%mp%c(:,:,m)* &
                product((pos-octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
          !NOTE: grad_ij = d u_j / d x_i
        endif
      enddo
      velv(:,ip) = vel
      if(sim_param%use_vs)  then
        !str = matmul(alpha, grad)
        str = matmul(alpha, transpose(grad))
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
        stretch = stretch + str
        !remove the parallel component
!       stretch = stretch + str - sum(str*dir)*dir
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
        stretchv(:,ip) = stretch
        if(sim_param%use_divfilt) then
          rotu(1) = grad(2,3) - grad(3,2)
          rotu(2) = grad(3,1) - grad(1,3)
          rotu(3) = grad(1,2) - grad(2,1)
          rotuv(:,ip) = rotu
        endif
      endif
      if(sim_param%use_vd)  then
        ave_ros = ave_ros + sum((0.5_wp*(grad+transpose(grad)))**2)
      endif
    endif
    enddo !ip in the leave
!$omp end simd

    if(sim_param%use_vd)  then
      if(octree%leaves(lv)%p%npart .gt. 0) then
        ave_ros = ave_ros/real(octree%leaves(lv)%p%npart,wp)
      endif
      !if(sim_param%use_tv) then
      !  turbvisc = (0.11_wp*sim_param%VortexRad)**2 * &
      !             sqrt(2.0_wp*ave_ros)
      !else
      !  turbvisc = 0.0_wp
      !endif

    endif

    do ip = 1,octree%leaves(lv)%p%npart
    if(.not. octree%leaves(lv)%p%cell_parts(ip)%p%free) then

      vel = velv(:,ip)
      stretch = stretchv(:,ip)
      stretch_alone = stretchv(:,ip)
      grad = 0.0_wp
      rotu  = rotuv(:,ip)
      pos = octree%leaves(lv)%p%cell_parts(ip)%p%cen
      alpha = octree%leaves(lv)%p%cell_parts(ip)%p%mag * &
              octree%leaves(lv)%p%cell_parts(ip)%p%dir
      dir = octree%leaves(lv)%p%cell_parts(ip)%p%dir
      if(sim_param%use_tv) then
        !turbvisc = (0.11_wp*sim_param%VortexRad)**2 * &
        !           sqrt(2.0_wp*ave_ros)
        octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = &
                      (0.11_wp*octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex)**2 * &
                      sqrt(2.0_wp*ave_ros)
      else
        octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = 0.0_wp
      endif
      !octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = ! turbvisc


      !then interact with all the neighbouring cell particles
      do k=-1,1; do j=-1,1; do i=-1,1
        if(associated(octree%leaves(lv)%p%neighbours(i,j,k)%p)) then
        if(octree%leaves(lv)%p%neighbours(i,j,k)%p%active) then
          do ipp = 1,octree%leaves(lv)%p%neighbours(i,j,k)%p%npart

            call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                  %compute_vel(pos, v) 
            vel = vel +v*one_4pi
            if(sim_param%use_vs) then
              call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                  %compute_stretch(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, str)
              stretch = stretch + str*one_4pi
              stretch_alone = stretch_alone + str*one_4pi
              if(sim_param%use_divfilt) then
                call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                    %compute_rotu(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex,  ru)
                rotu = rotu + ru*one_4pi
              endif
            endif
            if(sim_param%use_vd) then
              call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                  %compute_diffusion(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                  octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
              stretch = stretch +str*(2.0_wp*sim_param%nu_inf + octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)!turbvisc)
                ! 21/12/2023 Added factor 2 (see Winckelmans)
            endif
          enddo
        endif
        endif
      enddo; enddo; enddo

      !interact with the leaf neighbours
      do iln = 1,size(octree%leaves(lv)%p%leaf_neigh)
        do ipp = 1,octree%leaves(lv)%p%leaf_neigh(iln)%p%npart

          call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
              %compute_vel(pos, v)
          vel = vel +v*one_4pi
          if(sim_param%use_vs) then
            call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
              %compute_stretch(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex,  str)
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            stretch = stretch + str*one_4pi
            stretch_alone = stretch_alone + str*one_4pi
!           !removed the parallel component
!           stretch = stretch +(str - sum(str*dir)*dir)*one_4pi
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            if(sim_param%use_divfilt) then
              call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
                %compute_rotu(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, ru)
              rotu = rotu + ru*one_4pi
           endif
         endif
         if(sim_param%use_vd) then
           call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
              %compute_diffusion(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
           stretch = stretch + str*(2.0_wp*sim_param%nu_inf+octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)!turbvisc)
            ! 21/12/2023 Added factor 2 (see Winckelmans)
         endif
       enddo
      enddo

      !finally interact with the particles inside the cell
      do ipp = 1,octree%leaves(lv)%p%npart
        if (ipp .ne. ip) then
          call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_vel(pos, &
                                                         v)
          vel = vel +v*one_4pi
          if(sim_param%use_vs) then
            call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_stretch(pos, &
                  alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex,  str)
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            stretch = stretch + str*one_4pi
            stretch_alone = stretch_alone + str*one_4pi
!           !> removed the parallel component ( proj to avoid numerical instability ? )
!           stretch = stretch +(str - sum(str*dir)*dir)*one_4pi
! === VORTEX STRETCHING: AVOID NUMERICAL INSTABILITIES ? ===
            if(sim_param%use_divfilt) then
              call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_rotu(pos, &
                    alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, ru)
              rotu = rotu + ru*one_4pi
            endif
          endif

          if(sim_param%use_vd) then
            call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_diffusion(pos, &
                     !                             alpha, str)
                     alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex,&
                     octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
            stretch = stretch + str*(2.0_wp*sim_param%nu_inf+octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)!turbvisc)
             ! 21/12/2023 Added factor 2 (see Winckelmans)
          endif

        endif
      enddo

      !Calculate the interactions with all the other elements
      !calculate the influence of the solid bodies
      do ie=1,size(elem)
        call elem(ie)%p%compute_vel(pos,  v)
        vel = vel + v*one_4pi
      enddo

      ! calculate the influence of the wake panels
      do ie=1,size(wpan)
        call wpan(ie)%p%compute_vel(pos,  v)
        vel = vel + v*one_4pi
      enddo

      ! calculate the influence of the wake rings
      do ie=1,size(wrin)
        call wrin(ie)%p%compute_vel(pos,  v)
        vel = vel+ v*one_4pi
      enddo

      !calculate the influence of the end vortex
      do ie=1,size(wvort)
        call wvort(ie)%compute_vel(pos,  v)
        vel = vel+ v*one_4pi
      enddo

      !at last, add the free stream velocity
      vel = vel + variable_wind(pos, sim_param%time)
      octree%leaves(lv)%p%cell_parts(ip)%p%vel = vel

      !> Vortex stretching from other elements
      if(sim_param%use_vs .and. sim_param%vs_elems) then
        stretch_elem = 0.0_wp
        ! Influence of the body elements
        do ie=1,size(elem)
          call elem(ie)%p%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        ! Influence of the wake panels
        do ie=1,size(wpan)
          call wpan(ie)%p%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        ! Influence of the end vortex
        do ie=1,size(wvort)
          call wvort(ie)%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        stretch = stretch + stretch_elem
      endif

      !evolve the position in time
      if(sim_param%use_vs .or. sim_param%use_vd) then
        !evolve the intensity in time
        octree%leaves(lv)%p%cell_parts(ip)%p%stretch = stretch
        octree%leaves(lv)%p%cell_parts(ip)%p%stretch_alone = stretch_alone
        if(sim_param%use_divfilt) octree%leaves(lv)%p%cell_parts(ip)%p%rotu = rotu
      endif


    endif
    enddo
    deallocate(velv, stretchv, rotuv)
  enddo
!Don't ask me why, but without this pragma it is much faster!
!$omp end parallel do
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated leaves interactions in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)
  tend = dust_time()
  octree%t_lv = tend-tsta
  if(sim_param%use_dyn_layers) then
    if(real(octree%t_lv/octree%t_mp, wp) .gt. sim_param%LeavesTimeRatio .and. &
        octree%nlevels .lt. sim_param%NMaxOctreeLevels) then
      call add_layer(octree)
    endif
  endif

end subroutine apply_multipole

!----------------------------------------------------------------------

!> Apply the local expansion and calculate interaction when stretching and divfilt are used (optimized without if) 
subroutine apply_multipole_optimized(part, octree, elem, wpan, wrin, wvort)
  type(t_vortpart_p), intent(in), target :: part(:)
  type(t_octree), intent(inout)          :: octree
  type(t_pot_elem_p), intent(in)         :: elem(:)
  type(t_pot_elem_p), intent(in)         :: wpan(:)
  type(t_pot_elem_p), intent(in)         :: wrin(:)
  class(t_vortline), intent(in)          :: wvort(:)

  integer                     :: i, j, k, lv, ip, ipp, m, ie, iln
  real(wp)                    :: Rnorm2, vel(3), pos(3), v(3), stretch(3), stretch_alone(3), str(3), alpha(3), dir(3)
  real(wp), allocatable       :: velv(:,:), stretchv(:,:), rotuv(:,:)
  real(wp)                    :: grad(3,3)
  real(t_realtime)            :: tsta , tend
  real(wp)                    :: turbvisc, ave_ros
  real(wp)                    :: rotu(3), ru(3)

  real(wp) :: grad_elem(3,3) , stretch_elem(3)

  tsta = dust_time()
  !for all the leaves apply the local expansion and then local interactions
  t0 = dust_time()
  
!$omp parallel do private(lv, ip, ie, vel, pos, m, i, j, k, ipp, Rnorm2, &
!$omp& v, stretch, stretch_alone, str, grad, alpha, dir, velv, stretchv, ave_ros, &
!$omp& grad_elem, stretch_elem, iln, rotuv, rotu, ru) schedule(dynamic)
  do lv = 1, octree%nleaves
    !I am on a leaf, cycle on all the particles inside the leaf
    allocate(velv(3,octree%leaves(lv)%p%npart),&
              stretchv(3,octree%leaves(lv)%p%npart), &
              rotuv(3,octree%leaves(lv)%p%npart))

    ave_ros = 0.0_wp

!$omp simd private(vel, stretch, stretch_alone, str, grad, pos, alpha, dir, m) reduction(+:ave_ros)
    do ip = 1,octree%leaves(lv)%p%npart
    if(.not. octree%leaves(lv)%p%cell_parts(ip)%p%free) then
      vel = 0.0_wp
      stretch = 0.0_wp
      stretch_alone = 0.0_wp
      grad = 0.0_wp
      rotu = 0.0_wp
      pos = octree%leaves(lv)%p%cell_parts(ip)%p%cen
      alpha = octree%leaves(lv)%p%cell_parts(ip)%p%mag * &
              octree%leaves(lv)%p%cell_parts(ip)%p%dir
      dir = octree%leaves(lv)%p%cell_parts(ip)%p%dir

      !first apply the local multipole expansion
      do m = 1,size(octree%leaves(lv)%p%mp%b,2)
        !velocity
        vel = vel + &
              octree%leaves(lv)%p%mp%b(:,m)* &
              product((pos-octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
        !stretching: grad_ij = d u_j / d x_i (transpose) 
        grad = grad + octree%leaves(lv)%p%mp%c(:,:,m)* &
                product((pos-octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
      enddo
      velv(:,ip) = vel
      str = matmul(alpha, transpose(grad))
      stretch = stretch + str
      stretchv(:,ip) = stretch
      rotu(1) = grad(2,3) - grad(3,2)
      rotu(2) = grad(3,1) - grad(1,3)
      rotu(3) = grad(1,2) - grad(2,1)
      rotuv(:,ip) = rotu
      ave_ros = ave_ros + sum((0.5_wp*(grad+transpose(grad)))**2)
      endif
    enddo !ip in the leave
!$omp end simd

    if(sim_param%use_vd)  then
      if(octree%leaves(lv)%p%npart .gt. 0) then
        ave_ros = ave_ros/real(octree%leaves(lv)%p%npart,wp)
      endif
      !if(sim_param%use_tv) then
      !  turbvisc = (0.11_wp*sim_param%VortexRad)**2 * &
      !             sqrt(2.0_wp*ave_ros)
      !else
      !  turbvisc = 0.0_wp
      !endif

    endif

    do ip = 1,octree%leaves(lv)%p%npart
    if(.not. octree%leaves(lv)%p%cell_parts(ip)%p%free) then

      vel = velv(:,ip)
      stretch = stretchv(:,ip)
      stretch_alone = stretchv(:,ip)
      grad = 0.0_wp
      rotu  = rotuv(:,ip)
      pos = octree%leaves(lv)%p%cell_parts(ip)%p%cen
      alpha = octree%leaves(lv)%p%cell_parts(ip)%p%mag * &
              octree%leaves(lv)%p%cell_parts(ip)%p%dir
      dir = octree%leaves(lv)%p%cell_parts(ip)%p%dir
      !if(sim_param%use_tv) then
      !  !turbvisc = (0.11_wp*sim_param%VortexRad)**2 * &
      !  !           sqrt(2.0_wp*ave_ros)
      !  octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = &
      !                (0.11_wp*octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex)**2 * &
      !                sqrt(2.0_wp*ave_ros)
      !else
        octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = 0.0_wp
      !endif
      !octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc = ! turbvisc


      !then interact with all the neighbouring cell particles
      do k=-1,1; do j=-1,1; do i=-1,1
        if(associated(octree%leaves(lv)%p%neighbours(i,j,k)%p)) then
        if(octree%leaves(lv)%p%neighbours(i,j,k)%p%active) then
          do ipp = 1,octree%leaves(lv)%p%neighbours(i,j,k)%p%npart
            call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                  %compute_quantities(pos, alpha,  octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex , &
                                      v, str, ru)
            vel = vel +v*one_4pi
            stretch = stretch + str*one_4pi
            stretch_alone = stretch_alone + str*one_4pi
            rotu = rotu + ru*one_4pi
            !> diffusion 
            call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                  %compute_diffusion(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                  octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
            stretch = stretch +str*(2.0_wp*sim_param%nu_inf + octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)
          enddo
        endif
        endif
      enddo; enddo; enddo

      !interact with the leaf neighbours
      do iln = 1,size(octree%leaves(lv)%p%leaf_neigh)
        do ipp = 1,octree%leaves(lv)%p%leaf_neigh(iln)%p%npart
          call  octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
                  %compute_quantities(pos, alpha,  octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                                      v, str, ru)
          vel = vel + v*one_4pi
          stretch = stretch + str*one_4pi
          stretch_alone = stretch_alone + str*one_4pi
          rotu = rotu + ru*one_4pi
          !> diffusion 
          call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
                %compute_diffusion(pos, alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
          stretch = stretch + str*(2.0_wp*sim_param%nu_inf+octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)
        enddo
      enddo

      !finally interact with the particles inside the cell
      do ipp = 1,octree%leaves(lv)%p%npart
        if (ipp .ne. ip) then
          call octree%leaves(lv)%p%cell_parts(ipp)%p%&
                compute_quantities(pos, alpha,  octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex, &
                                  v, str, ru) 
          vel = vel +v*one_4pi
          stretch = stretch + str*one_4pi
          stretch_alone = stretch_alone + str*one_4pi
          rotu = rotu + ru*one_4pi
          !> diffusion 
          call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_diffusion(pos, &
                alpha, octree%leaves(lv)%p%cell_parts(ip)%p%r_Vortex,&
                octree%leaves(lv)%p%cell_parts(ip)%p%vol, str)
          stretch = stretch + str*(2.0_wp*sim_param%nu_inf+octree%leaves(lv)%p%cell_parts(ip)%p%turbvisc)
        endif
      enddo

      !Calculate the interactions with all the other elements
      !calculate the influence of the solid bodies
      do ie=1,size(elem)
        call elem(ie)%p%compute_vel(pos,  v)
        vel = vel + v*one_4pi
      enddo

      ! calculate the influence of the wake panels
      do ie=1,size(wpan)
        call wpan(ie)%p%compute_vel(pos,  v)
        vel = vel + v*one_4pi
      enddo

      ! calculate the influence of the wake rings
      do ie=1,size(wrin)
        call wrin(ie)%p%compute_vel(pos,  v)
        vel = vel+ v*one_4pi
      enddo

      !calculate the influence of the end vortex
      do ie=1,size(wvort)
        call wvort(ie)%compute_vel(pos,  v)
        vel = vel+ v*one_4pi
      enddo

      !at last, add the free stream velocity
      vel = vel + variable_wind(pos, sim_param%time)
      octree%leaves(lv)%p%cell_parts(ip)%p%vel = vel

      !> Vortex stretching from other elements
      if(sim_param%use_vs .and. sim_param%vs_elems) then
        stretch_elem = 0.0_wp
        ! Influence of the body elements
        do ie=1,size(elem)
          call elem(ie)%p%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        ! Influence of the wake panels
        do ie=1,size(wpan)
          call wpan(ie)%p%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        ! Influence of the end vortex
        do ie=1,size(wvort)
          call wvort(ie)%compute_grad(pos, grad_elem)
          stretch_elem = stretch_elem + &
              matmul(transpose(grad_elem),alpha)*one_4pi
        enddo
        stretch = stretch + stretch_elem
      endif

      octree%leaves(lv)%p%cell_parts(ip)%p%stretch = stretch
      octree%leaves(lv)%p%cell_parts(ip)%p%stretch_alone = stretch_alone
      octree%leaves(lv)%p%cell_parts(ip)%p%rotu = rotu
      

    endif
    enddo
    deallocate(velv, stretchv, rotuv)
  enddo
!Don't ask me why, but without this pragma it is much faster!
!$omp end parallel do
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated leaves interactions in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)
  tend = dust_time()
  octree%t_lv = tend-tsta
  if(sim_param%use_dyn_layers) then
    if(real(octree%t_lv/octree%t_mp, wp) .gt. sim_param%LeavesTimeRatio .and. &
        octree%nlevels .lt. sim_param%NMaxOctreeLevels) then
      call add_layer(octree)
    endif
  endif

end subroutine apply_multipole_optimized

!> Apply multipole to the panels
subroutine apply_multipole_panels(octree, elem)
 type(t_octree), intent(inout) :: octree
 type(t_pot_elem_p), intent(in) :: elem(:)

 integer :: i, j, k, lv, ip, ipp, m, ie, iln
 real(wp) :: Rnorm2, vel(3), pos(3), v(3), stretch(3), str(3), alpha(3), dir(3)
 real(wp), allocatable :: velv(:,:), stretchv(:,:)
 real(wp) :: grad(3,3)
 real(t_realtime) :: tsta
 real(wp) :: turbvisc, ave_ros

 character(len=*), parameter :: this_sub_name = 'apply_multipole_panels'


  tsta = dust_time()
  !for all the leaves apply the local expansion and then local interactions
  t0 = dust_time()
!$omp parallel do private(lv, ip, ie, vel, pos, m, i, j, k, ipp, Rnorm2, &
!$omp& v, stretch, str, grad, alpha, dir, velv, stretchv, ave_ros, turbvisc, iln) schedule(dynamic)
  do lv = 1, octree%nleaves
    !I am on a leaf, cycle on all the particles inside the leaf

    do ip = 1,octree%leaves(lv)%p%npan

      vel = 0.0_wp
      pos = octree%leaves(lv)%p%cell_pans(ip)%p%cen

      !first apply the local multipole expansion
      do m = 1,size(octree%leaves(lv)%p%mp%b,2)
        !velocity
        vel = vel + &
          octree%leaves(lv)%p%mp%b(:,m)* &
          product((pos-octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
      enddo

      !then interact with all the neighbouring cell particles
      do k=-1,1; do j=-1,1; do i=-1,1
        if(associated(octree%leaves(lv)%p%neighbours(i,j,k)%p)) then
        if(octree%leaves(lv)%p%neighbours(i,j,k)%p%active) then
          do ipp = 1,octree%leaves(lv)%p%neighbours(i,j,k)%p%npart

            call octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p&
                 %compute_vel(pos, v)
            vel = vel +v*one_4pi
          enddo
        endif
        endif
      enddo; enddo; enddo

      !interact with the leaf neighbours
      do iln = 1,size(octree%leaves(lv)%p%leaf_neigh)
       do ipp = 1,octree%leaves(lv)%p%leaf_neigh(iln)%p%npart

         call octree%leaves(lv)%p%leaf_neigh(iln)%p%cell_parts(ipp)%p&
              %compute_vel(pos,  v)
         vel = vel +v*one_4pi
       enddo
      enddo

      !finally interact with the particles inside the cell
      do ipp = 1,octree%leaves(lv)%p%npart
        call octree%leaves(lv)%p%cell_parts(ipp)%p%compute_vel(pos, &
                                                      v)
        vel = vel +v*one_4pi
      enddo

      octree%leaves(lv)%p%cell_pans(ip)%p%uvort = &
      octree%leaves(lv)%p%cell_pans(ip)%p%uvort + vel
    enddo
  enddo
!Don't ask me why, but without this pragma it is much faster!
!$omp end parallel do
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated panels interactions in: ' , t1 - t0,' s.'
  if(sim_param%debug_level.ge.5) call printout(msg)

end subroutine apply_multipole_panels
!----------------------------------------------------------------------

!! DISCONTINUED
!!!> Check the quantity of particles in a cell
!!!!
!!!! If there are too few particles, these are bound together and applied to
!!!! the parent
!!recursive subroutine check_cell_content(cell, leaves)
!! type(t_cell), intent(inout) :: cell
!! type(t_cell_p), allocatable, intent(inout) :: leaves(:)
!!
!! integer :: sib
!! logical :: enough
!!
!!
!!  if(cell%npart .ge. min_part_4_cell) then
!!    ! Enough particles, since we started from the bottom, it is a leaf, set
!!    ! it as leaf and cut out the descendant
!!    call set_leaf(cell,.true.)
!!    call push_ptr(leaves, cell)
!!
!!  else
!!    !Not enough particles in this cell
!!    enough = .false.
!!    !We have to check that all the siblings have not enough particles
!!    !oviously only if the cell has siblings, otherwise is top level and
!!    !it will be a leaf anyway
!!    if(associated(cell%parent%p)) then
!!      do sib = 1,8
!!        if(cell%parent%p%children(sib)%p%npart .ge. min_part_4_cell) then
!!          enough = .true.
!!        endif
!!      enddo
!!    else
!!      !it does not have a parent, it is top level, force to be a leaf
!!      enough = .true.
!!    endif
!!
!!    if (.not. enough) then
!!      !all the siblings have not enough particles,  all are going to be
!!      ! set inactive, push the pointers to particles to the parent
!!      do sib = 1,8
!!        !TODO: for some reasons the generic does not work here
!!        call push_ptr(cell%parent%p%cell_parts, cell%parent%p%children(sib)%p%cell_parts)
!!      enddo
!!      !escalate this call upwards
!!      call check_cell_content(cell%parent%p, leaves)
!!
!!    else
!!      !the present cell has not enough particles, however the siblings
!!      !do, so the present cell becomes a leaf anyway
!!      call set_leaf(cell,.true.)
!!      call push_ptr(leaves, cell)
!!
!!    endif
!!
!!
!!  endif
!!
!!end subroutine check_cell_content

!----------------------------------------------------------------------

!> Set a cell leaf
!!
!! A cell is set as leaf, and if necessary the multipole data is
!! allocated. In all cases the multipole data is reset.
!! Note: leaves are not a dynamic vector, are allocated as a whole and
!! then only the first one are employed. This means that to set a leaf
!! in the correct position an increasing index is necessary, and must be
!! passed with intent(inout)
recursive subroutine set_leaf(cell,leaves,il,pexp)
 type(t_cell), intent(inout), target :: cell
 type(t_cell_p), intent(inout) :: leaves(:)
 integer, intent(inout) :: il
 type(t_polyexp), intent(in) :: pexp

  !I is a leaf
  cell%leaf = .true.

!$omp critical
  il = il+1
  leaves(il)%p => cell
!$omp end critical

  if(.not. cell%active) then
    !If it was not active, set active and allocate data
    cell%active = .true.
    call cell%mp%init(pexp)
  else
    !If it was alerady active, just reset the data
    call cell%mp%reset
  endif

end subroutine set_leaf

!----------------------------------------------------------------------

!> Set a cell and all the upward branch branch
!!
!! Recursively set all the upward branch as branch, allocating
!! the multipole data when necessary, resetting it to zero in any case
recursive subroutine set_branch(cell,pexp)
 type(t_cell), intent(inout) :: cell
 type(t_polyexp), intent(in) :: pexp

  cell%branch = .true.

  if(.not. cell%active) then
    !If it was not active, set active and allocate data
    cell%active = .true.
    call cell%mp%init(pexp)
  else
    !If it was alerady active, just reset the data
    call cell%mp%reset
  endif

!  if(associated(cell%parent%p)) call set_branch(cell%parent%p, pexp)

end subroutine set_branch

!----------------------------------------------------------------------

!> Set a cell inactive
!!
!! Everything is set to false and the multipole data are deallocated
recursive subroutine set_inactive(cell)
 type(t_cell), intent(inout) :: cell

  cell%branch = .false.
  cell%leaf = .false.

  if(cell%active) then
    cell%active = .false.
    call cell%mp%destroy
  endif

end subroutine set_inactive

!----------------------------------------------------------------------

!> Set a cell near wall
!!
!! The cell is set near wall and then all the children are set near wall
recursive subroutine set_near_wall(cell)
 type(t_cell), intent(inout) :: cell

 integer :: child

  cell%near_wall = .true.

  do child = 1,8
    if(associated(cell%children(child)%p)) &
       call set_near_wall(cell%children(child)%p)
  enddo

end subroutine set_near_wall

!----------------------------------------------------------------------

! DISCONTINUED
!> Add a certain number of particles upward in the tree
!!recursive subroutine add_part_upwards(cell, npart)
!! type(t_cell), intent(inout) :: cell
!! integer, intent(in) :: npart
!!
!!  cell%npart = cell%npart + npart
!!  if(associated(cell%parent%p)) then
!!    call add_part_upwards(cell%parent%p, npart)
!!  endif
!!
!!end subroutine add_part_upwards

!----------------------------------------------------------------------

!> Push another pointer to a list of element pointers
subroutine push_cell_ptr_scalar(pointers, element)
 type(t_cell_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_cell), intent(in), target :: element

 type(t_cell_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l

  l = size(pointers)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+1))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  pointers(l+1)%p => element

end subroutine push_cell_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer list to a list of element pointers
subroutine push_cell_ptr_vec(pointers, elements)
 type(t_cell_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_cell_p), target, intent(in) :: elements(:)

 type(t_cell_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, m

  l = size(pointers)
  m = size(elements)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+m))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  do i = 1,m
    pointers(l+i)%p => elements(i)%p
  enddo

end subroutine push_cell_ptr_vec

!----------------------------------------------------------------------

!> Push another pointer to a list of particles pointers
subroutine push_part_ptr_scalar(pointers, particle, last_elem)
 type(t_vortpart_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_vortpart), intent(in), target :: particle
 integer, optional :: last_elem

 type(t_vortpart_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, ife

  l = size(pointers)
  !Get first empty element: receive from
  if(.not.present(last_elem)) then
    ife = 1
    do i=1,l
      if(.not.associated(pointers(i)%p)) exit
      ife = ife + 1
    enddo
  else
    ife = last_elem+1
  endif

  if (ife .gt. l .or. ife .le. 0) then !the array is full, allocate new block
    do i=1,l
      tmp_ptr(i)%p => pointers(i)%p
    enddo
    deallocate(pointers); allocate(pointers(l+allocation_block))
    do i=1,l
      pointers(i)%p => tmp_ptr(i)%p
    enddo
    pointers(l+1)%p => particle
    do i=l+2,size(pointers)
      nullify(pointers(i)%p)
    enddo

  else !insert the pointer in the first available empty space

    pointers(ife)%p => particle

  endif

end subroutine push_part_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer list to a list of particle pointers
subroutine push_part_ptr_vec(pointers, elements, last_elem)
 type(t_vortpart_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_vortpart_p), target, intent(in) :: elements(:)
 integer, optional :: last_elem

 type(t_vortpart_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, m, ife, nblck

  l = size(pointers)
  !Get first empty element: receive from
  if(.not.present(last_elem)) then
    ife = 1
    do i=1,l
      if(.not.associated(pointers(i)%p)) exit
      ife = ife + 1
    enddo
  else
    ife = last_elem+1
  endif

  m = size(elements)

  if (ife+m-1 .gt. l .or. ife .le. 0) then !the array is full, allocate new blocks
    do i=1,ife-1
      tmp_ptr(i)%p => pointers(i)%p
    enddo
    nblck = ceiling(real(m-(l-ife+1))/real(allocation_block))
    deallocate(pointers); allocate(pointers(l+nblck*allocation_block))
    if(ife .le. 0) ife = 1
    do i=1,(ife-1)
      pointers(i)%p => tmp_ptr(i)%p
    enddo
    do i = 1,m
      pointers(ife-1+i)%p => elements(i)%p
    enddo
    do i = ife+m,size(pointers)
      nullify(pointers(i)%p)
    enddo

  else !there is enough space just to fill the additional array

    do i = 1,m
      pointers(ife-1+i)%p => elements(i)%p
    enddo

  endif

end subroutine push_part_ptr_vec

!----------------------------------------------------------------------

subroutine pull_part_ptr_scalar(pointers, idx)
 type(t_vortpart_p), allocatable, target, intent(inout) :: pointers(:)
 integer, intent(in) :: idx

 integer :: i

 !simply shift all the arrray to remove the pointer in position idx
 !the number of elements should be decreased outside here

 !pointers(idx:size(pointers)-1) = pointers(idx+1:size(pointers))
 !nullify(pointers(size(pointers))%p)
  do i=idx,size(pointers)-1
    pointers(i)%p => pointers(i+1)%p
  enddo

end subroutine pull_part_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer to a list of panels pointers
subroutine push_pan_ptr_scalar(pointers, particle, last_elem)
 type(t_pot_elem_p), allocatable, target, intent(inout) :: pointers(:)
 class(c_pot_elem), intent(in), target :: particle
 integer, intent(in), optional :: last_elem

 type(t_pot_elem_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, ife

  l = size(pointers)
  !Get first empty element: receive from
  if(.not.present(last_elem)) then
    ife = 1
    do i=1,l
      if(.not.associated(pointers(i)%p)) exit
      ife = ife + 1
    enddo
  else
    ife = last_elem+1
  endif

  if (ife .gt. l) then !the array is full, allocate new block
    do i=1,l
      tmp_ptr(i)%p => pointers(i)%p
    enddo
    deallocate(pointers); allocate(pointers(l+allocation_block))
    do i=1,l
      pointers(i)%p => tmp_ptr(i)%p
    enddo
    pointers(l+1)%p => particle
    do i=l+2,size(pointers)
      nullify(pointers(i)%p)
    enddo

  else !insert the pointer in the first available empty space

    pointers(ife)%p => particle

  endif

end subroutine push_pan_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer list to a list of panels pointers
subroutine push_pan_ptr_vec(pointers, elements, last_elem)
 type(t_pot_elem_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_pot_elem_p), target, intent(in) :: elements(:)
 integer, optional :: last_elem

 type(t_pot_elem_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, m, ife, nblck

  l = size(pointers)
  !Get first empty element: receive from
  if(.not.present(last_elem)) then
    ife = 1
    do i=1,l
      if(.not.associated(pointers(i)%p)) exit
      ife = ife + 1
    enddo
  else
    ife = last_elem+1
  endif

  m = size(elements)

  if (ife+m-1 .gt. l) then !the array is full, allocate new blocks
    do i=1,ife-1
      tmp_ptr(i)%p => pointers(i)%p
    enddo
    nblck = ceiling(real(m-(l-ife+1))/real(allocation_block))
    deallocate(pointers); allocate(pointers(l+nblck*allocation_block))
    do i=1,(ife-1)
      pointers(i)%p => tmp_ptr(i)%p
    enddo
    do i = 1,m
      pointers(ife-1+i)%p => elements(i)%p
    enddo
    do i = ife+m,size(pointers)
      nullify(pointers(i)%p)
    enddo

  else !there is enough space just to fill the additional array

    do i = 1,m
      pointers(ife-1+i)%p => elements(i)%p
    enddo

  endif

end subroutine push_pan_ptr_vec

!----------------------------------------------------------------------

!> Reset the particles stored inside the cell
!!
!! The number of particles, the particles pointers and the brach/leaves
!! flags are re-set. Activity flag is kept since it marks whether the
!! multipole data are allocated or not, will be set during sorting
subroutine reset_cell(cell)
 type(t_cell), intent(inout) :: cell
 integer :: i

 !deallocate(cell%cell_parts)
 !allocate(cell%cell_parts(0))
 !cell%npart = 0
 !deallocate(cell%cell_pans)
 !allocate(cell%cell_pans(0))
 !cell%npan = 0
 !!cell%active = .true.
 !cell%branch = .false.
 !cell%leaf = .false.
 !cell%near_wall = .false.
 !cell%ave_vortmag = 0.0_wp
 !deallocate(cell%leaf_neigh)
 !allocate(cell%leaf_neigh(0))




  do concurrent (i=1:cell%npart)
    nullify(cell%cell_parts(i)%p)
  enddo
  cell%npart = 0

  do concurrent (i=1:cell%npan)
    nullify(cell%cell_pans(i)%p)
  enddo
  cell%npan = 0
  !cell%active = .true.
  cell%branch = .false.
  cell%leaf = .false.
  cell%near_wall = .false.
  cell%ave_vortmag = 0.0_wp
  deallocate(cell%leaf_neigh)
  allocate(cell%leaf_neigh(0))

end subroutine reset_cell
!----------------------------------------------------------------------

end module mod_octree
