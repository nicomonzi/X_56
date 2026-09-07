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


!> Module to treat input-output for lifting line components
!!
module mod_ll_io

use mod_param, only: &
  wp, max_char_len, nl, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption

use mod_math, only: &
  geoseries, geoseries_both
!----------------------------------------------------------------------

implicit none

public :: read_mesh_ll, read_xac_offset

private

character(len=max_char_len) :: msg
character(len=*), parameter :: this_mod_name = 'mod_ll_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

subroutine read_mesh_ll(mesh_file,ee,rr, y_fountain, &
                        airfoil_list_actual, nelem_span_list, &
                        i_airfoil_e , normalised_coord_e    , &
                        npoints_chord_tot , nelem_span_tot  , &
                        chord_p,theta_p,theta_e)

  character(len=*), intent(in)                            :: mesh_file
  integer  , allocatable, intent(out)                     :: ee(:,:)
  real(wp) , allocatable, intent(out)                     :: rr(:,:)
  character(len=max_char_len), allocatable , intent(out)  :: airfoil_list_actual(:)
  character(len=max_char_len), allocatable                :: airfoil_list(:)
  integer  , allocatable, intent(out)                     :: nelem_span_list(:)
  integer  , allocatable, intent(out)                     :: i_airfoil_e(:,:)
  real(wp) , allocatable, intent(out)                     :: normalised_coord_e(:,:)
  integer  ,              intent(out)                     :: npoints_chord_tot, nelem_span_tot
  real(wp) , allocatable, intent(out)                     :: chord_p(:),theta_p(:),theta_e(:)
  real(wp) , intent(out)                                  :: y_fountain 
  type(t_parse)                                           :: pmesh_prs

  logical                                                 :: twist_linear_interp
  integer                                                 :: nelem_chord, nelem_chord_tot 
  integer                                                 :: npoint_chord_tot, npoint_span_tot
  integer                                                 :: nRegions, nSections, nAirfoils 
  integer                                                 :: rr_size , ee_size , ispace
  integer                                                 :: iRegion, iSection, iAirfoil, iChord, iElement, iPoint
  real(wp)                                                :: ref_chord_fraction
  real(wp), allocatable                                   :: ref_point(:)
  !> data read from file
  ! sections ---
  real(wp), allocatable                                   :: chord_list(:), twist_list(:)
  ! regions  ---
  real(wp), allocatable                                   :: span_list(:), sweep_list(:), dihed_list(:)
  real(wp), allocatable                                   :: y_ref_list(:), r_in_list(:), r_ob_list(:) 
  character(len=max_char_len), allocatable                :: type_span_list(:)
  integer                                                 :: n_type_span
  character                                               :: ElType
  logical                                                 :: symmetry, mesh_flat
  !> Sections 1. 2.
  real(wp), allocatable                                   :: rrSection1(:,:) , rrSection2(:,:)
  real(wp)                                                :: dx_ref , dy_ref , dz_ref
  integer                                                 :: ista , iend , ich

  ! Linear interpolation of the twist angle
  real(wp)                                                :: w1, w2

  integer :: i_aero1, i_aero2, iSpan, i

  character(len=*), parameter :: this_sub_name = 'read_mesh_parametric'


  !Prepare all the parameters to be read in the file
  ! Global parameters
  call pmesh_prs%CreateStringOption('el_type', &
                'element type (temporary) p panel v vortex ring')
  call pmesh_prs%CreateLogicalOption('mesh_symmetry', &
                'symmetry yes/no','F' )
  call pmesh_prs%CreateLogicalOption('mesh_flat', &
                'flat mesh yes/no','F' )
  call pmesh_prs%CreateRealArrayOption('starting_point',&
                'Starting point (inboard TE), (x, y, z)', &
                '(/0.0, 0.0, 0.0/)',&
                multiple=.false.);
  call pmesh_prs%CreateRealOption('reference_chord_fraction',&
                'Reference chord fraction', &
                '0.0',&
                multiple=.false.);
  call pmesh_prs%CreateLogicalOption('twist_linear_interpolation',&
                'Linear interpolation of the twist angle, for the whole component',&
                'F',&
                multiple=.false.);
  ! TODO: do we really need to allow the user to define the parameter above?
  ! the reference chord fraction is a number in the range [0,1] that defines
  ! - where the reference point is located in the root chord
  ! - the line related to the sweep angle
  ! - the point around which the sections are rotated to impose twist


  ! Section parameters
  call pmesh_prs%CreateRealOption(     'chord', 'section chord', &
                multiple=.true.);
  call pmesh_prs%CreateRealOption(     'twist', 'section twist angle',&
                multiple=.true.);
  call pmesh_prs%CreateStringOption( 'airfoil_table', 'section airfoil',&
                multiple=.true.);

  ! Region parameters
  call pmesh_prs%CreateRealOption(    'span', 'region span',&
                multiple=.true.);
  call pmesh_prs%CreateRealOption(   'sweep', 'region sweep angle [degrees]',&
                multiple=.true.);
  call pmesh_prs%CreateRealOption(   'dihed', 'region dihedral angle [degrees]',&
                multiple=.true.);
  call pmesh_prs%CreateIntOption( 'nelem_span', 'number of span-wise elements in the region',&
                multiple=.true.);
  call pmesh_prs%CreateStringOption('type_span', 'type of span-wise division: &
                                    & uniform, cosine, cosineIB, cosineOB, equalarea', &
                multiple=.true.);

  call pmesh_prs%CreateRealOption( 'r', 'growth ratio of the elements at edge', &
                multiple=.true.);
  call pmesh_prs%CreateRealOption( 'r_ib', 'growth ratio of the elements inboard', &
                multiple=.true.);
  call pmesh_prs%CreateRealOption( 'r_ob', 'growth ratio of the elements at outboard', &
                multiple=.true.);
  call pmesh_prs%CreateRealOption( 'y_refinement', 'spanwise station to which the refinement start', &
                multiple=.true.); 
  call pmesh_prs%CreateRealOption( 'y_fountain', &
                                  'Upper y-coordinates beyond which the wake particles are not released', '0.0');
  !read the parameters
  call pmesh_prs%read_options(mesh_file, printout_val=.true.)

  nelem_chord     = 1          
  nelem_chord_tot = 1          
  ElType  = getstr(pmesh_prs,'el_type')
  symmetry= getlogical(pmesh_prs,'mesh_symmetry')
  mesh_flat= getlogical(pmesh_prs,'mesh_flat')
  y_fountain = getreal(pmesh_prs,'y_fountain') 

  if ( mesh_flat .and. trim(ElType) .ne. 'l' ) then
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &flat mesh option is available only for lifting line elements.')
  end if

  nSections = countoption(pmesh_prs,'chord')
  nRegions  = countoption(pmesh_prs,'span')

  write(msg,*) nl,'  Number of sections: ',nSections; call printout(msg)
  write(msg,*) ' Number of regions:  ',nRegions; call printout(msg)

  ! Check that nSections = nRegion + 1
  if ( nSections .ne. nRegions + 1 ) then
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of sections different from number of regions + 1.')
  end if

  ref_chord_fraction = getreal(pmesh_prs,'reference_chord_fraction')
  ref_point          = getrealarray(pmesh_prs,'starting_point',3)

  twist_linear_interp= getlogical(pmesh_prs,'twist_linear_interpolation')

  !Check the number of inputs
  if(countoption(pmesh_prs,'nelem_span') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of "nelem_span" different from number of regions.')
  if(countoption(pmesh_prs,'dihed') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of "dihed" different from number of regions.')
  if(countoption(pmesh_prs,'sweep') .ne. nRegions ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of "sweep" different from number of regions.')
  if(countoption(pmesh_prs,'twist') .ne. nSections ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of "twist" different from number of sections.')
  if(countoption(pmesh_prs,'airfoil_table') .ne. nSections ) &
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of "airfoil_table" different from number of sections.')

  ! Get total number of elements and initialize arrays
  allocate(nelem_span_list(nRegions));   nelem_span_list = 0
  allocate(      span_list(nRegions));         span_list = 0.0_wp
  allocate(     sweep_list(nRegions));        sweep_list = 0.0_wp
  allocate(     dihed_list(nRegions));        dihed_list = 0.0_wp
  allocate(     y_ref_list(nRegions));        y_ref_list = 0.0_wp
  allocate(      r_in_list(nRegions));         r_in_list = 0.0_wp
  allocate(      r_ob_list(nRegions));         r_ob_list = 0.0_wp 
  nelem_span_tot = 0
  do iRegion = 1,nRegions
    nelem_span_list(iRegion) = getint(pmesh_prs,'nelem_span')
    span_list(iRegion)  = getreal(pmesh_prs,'span' )
    sweep_list(iRegion) = getreal(pmesh_prs,'sweep')
    dihed_list(iRegion) = getreal(pmesh_prs,'dihed')
    nelem_span_tot = nelem_span_tot + nelem_span_list(iRegion)
  enddo

  ! type_span
  allocate( type_span_list(nRegions));!   type_span_list = ''
  n_type_span = countoption(pmesh_prs,'type_span')
  if ( n_type_span .eq. 0 ) then ! default is 'uniform'
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = 'uniform'
    end do
  else if ( n_type_span .eq. nRegions ) then
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = getstr(pmesh_prs,'type_span')
      !> check geometry series type 
      if (trim(type_span_list(iRegion)) .eq. 'geoseries') then 
        y_ref_list(iRegion)      = getreal(pmesh_prs,   'y_refinement')
        write(*,*) ' y_ref_list(',iRegion,') : ' , y_ref_list(iRegion)
        r_in_list(iRegion)       = getreal(pmesh_prs,   'r_ib')
        r_ob_list(iRegion)       = getreal(pmesh_prs,   'r_ob') 
      elseif (trim(type_span_list(iRegion)) .eq. 'geoseries_ib') then 
        r_in_list(iRegion)       = getreal(pmesh_prs,   'r_ib')
      elseif (trim(type_span_list(iRegion)) .eq. 'geoseries_ob') then
        r_ob_list(iRegion)       = getreal(pmesh_prs,   'r_ob')
      endif
    end do
  else
    write(*,*) ' mesh_file   : ' , trim(mesh_file)
    write(*,*) ' n_type_span : ' , n_type_span
    write(*,*) ' nRegions    : ' , nRegions
    call error(this_sub_name, this_mod_name, 'Inconsistent input: &
          &number of type_span different from number of regions.')
  end if

  allocate(chord_list  (nSections))  ; chord_list = 0.0_wp
  allocate(twist_list  (nSections))  ; twist_list = 0.0_wp
  allocate(airfoil_list(nSections))

  do iSection= 1,nSections
    chord_list(iSection)   = getreal(pmesh_prs,'chord')
    twist_list(iSection)   = getreal(pmesh_prs,'twist')
    airfoil_list(iSection) = getstr(pmesh_prs,'airfoil_table')
  enddo

  npoint_chord_tot = nelem_chord_tot + 1
  npoint_span_tot  = nelem_span_tot  + 1

  ee_size = nelem_span_tot
  rr_size = npoint_span_tot * npoint_chord_tot

  ! Directly build connectivity matrix
  allocate(ee(4,ee_size)) ; ee = 0
  do iChord = 1,nelem_chord_tot
    do iSpan = 1,nelem_span_tot
      iElement =  nelem_chord_tot*(iSpan-1) + iChord
      iPoint   = npoint_chord_tot*(iSpan-1) + iChord
      ee(1,iElement) = iPoint + npoint_chord_tot         ! iPoint
      ee(2,iElement) = iPoint                            ! iPoint + 1
      ee(3,iElement) = iPoint + 1                        ! iPoint + npoint_chord_tot + 1
      ee(4,iElement) = iPoint + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
    enddo
  enddo

  allocate(rr   (3,rr_size)) ; rr      = 0.0_wp
  allocate(chord_p(npoint_span_tot)) ; chord_p = 0.0_wp
  allocate(theta_p(npoint_span_tot)) ; theta_p = 0.0_wp

  allocate(i_airfoil_e       (2,nelem_span_tot)) ; i_airfoil_e        = 0
  allocate(normalised_coord_e(2,nelem_span_tot)) ; normalised_coord_e = 0.0_wp
  allocate(theta_e           (  nelem_span_tot)) ; theta_e            = 0.0_wp

! ! Initialize the span division to the maximum dimension
  allocate(rrSection1(3,npoint_chord_tot)) ; rrSection1 = 0.0_wp
  allocate(rrSection2(3,npoint_chord_tot)) ; rrSection2 = 0.0_wp

  ! Initialise dr_ref for the definition of the actual reference_point
  !  of each bay (by updating)
  dx_ref = 0.0_wp
  dy_ref = 0.0_wp
  dz_ref = 0.0_wp


  twist_list = twist_list * pi / 180
  ista = 1 ; iend = npoint_chord_tot
  ich = 1

  do iRegion = 1,nRegions

    if ( iRegion .gt. 1 ) then  ! first section = last section of the previous region
      rrSection1 = rrSection2
    else                        ! build points
      rrSection1(:,1) = ref_point

      ! Rotate the section around the reference line with the twist angle
      !
      ! For flat meshes the section is not rotated, but the
      ! normal/tangent vector are

      if ( mesh_flat ) then
        rrSection1(1,2) = rrSection1(1,1) + chord_list(iRegion)
        rrSection1(2,2) = rrSection1(2,1)
        rrSection1(3,2) = rrSection1(3,1)
      else
        rrSection1(1,2) = rrSection1(1,1) + chord_list(iRegion) * cos(twist_list(iRegion))
        rrSection1(2,2) = rrSection1(2,1)
        rrSection1(3,2) = rrSection1(3,1) - chord_list(iRegion) * sin(twist_list(iRegion))
      endif

      chord_p(ich) = chord_list(1)
      theta_p(ich) = twist_list(1)

      ! Update rr
      rr(:,ista:iend) = rrSection1

    end if

    ! Define the spacing type based on the input string
    if ( trim(type_span_list(iRegion)) .eq. 'uniform' ) then
      ispace = 0

    else if ( trim(type_span_list(iRegion)) .eq. 'cosine' ) then
      ispace = 1

    else if ( trim(type_span_list(iRegion)) .eq. 'cosineOB' ) then
      ispace = 2

    else if ( trim(type_span_list(iRegion)) .eq. 'cosineIB' ) then
      ispace = 3

    else if ( trim(type_span_list(iRegion)) .eq. 'equalarea' ) then
      ispace = 4

    else if ( trim(type_span_list(iRegion)) .eq. 'geoseries' ) then
      ispace = 5 

    else if ( trim(type_span_list(iRegion)) .eq. 'geoseries_ob' ) then
      ispace = 6 

    else if ( trim(type_span_list(iRegion)) .eq. 'geoseries_ib' ) then
      ispace = 7

    else
      write(*,*) ' mesh_file   : ' , trim(mesh_file)
      write(*,*) ' type_span_list(',iRegion,') : ' , &
                    trim(type_span_list(iRegion))
      call error(this_sub_name, this_mod_name, 'Inconsistent input: &
            & type_span must be equal to uniform, cosine, cosineIB,&
            & cosineOB, equalarea, geoseries, geoseries_ob, geoseries_ib.')
    end if

    !> save before update
    dx_ref =  span_list(iRegion) * tan( sweep_list(iRegion)* pi / 180.0_wp ) ! + dx_ref
    dy_ref =  span_list(iRegion)                                             ! + dy_ref
    dz_ref =  span_list(iRegion) * tan( dihed_list(iRegion)* pi / 180.0_wp ) ! + dz_ref

    ! Rotate the section around the reference line with the twist angle
    !
    ! For flat meshes the section is not rotated, but the
    ! normal/tangent vector are

    rrSection2(1,1) = rrSection1(1,1) + dx_ref
    rrSection2(2,1) = rrSection1(2,1) + dy_ref
    rrSection2(3,1) = rrSection1(3,1) + dz_ref

    if ( mesh_flat ) then
      rrSection2(1,2) = rrSection2(1,1) + chord_list(iRegion+1)
      rrSection2(2,2) = rrSection2(2,1)
      rrSection2(3,2) = rrSection2(3,1)
    else
      rrSection2(1,2) = rrSection2(1,1) + chord_list(iRegion+1) * cos(twist_list(iRegion+1))
      rrSection2(2,2) = rrSection2(2,1)
      rrSection2(3,2) = rrSection2(3,1) - chord_list(iRegion+1) * sin(twist_list(iRegion+1))
    endif

    ! Interpolation of the nodes of the region i (between sections i and i+1)
    do iSpan = 1 , nelem_span_list(iRegion)

      ista = iend + 1
      iend = iend + npoint_chord_tot
      ich = ich + 1

      w2 = spacing_weights ( ispace, iSpan, nelem_span_list(iRegion) , &
                            r_in_list(iRegion), r_ob_list(iRegion), y_ref_list(iRegion))
      w1 = 1.0_wp - w2

      chord_p(ich) = w1*chord_list(iRegion) + w2*chord_list(iRegion+1)
      theta_p(ich) = w1*twist_list(iRegion) + w2*twist_list(iRegion+1)


      ! if linear twist -> interpolate twist
      !        othrwise -> interpolate coordinates
      if ( twist_linear_interp .and. .not. mesh_flat ) then

        rr(1,ista) = rrSection1(1,1) + w2*dx_ref
        rr(2,ista) = rrSection1(2,1) + w2*dy_ref
        rr(3,ista) = rrSection1(3,1) + w2*dz_ref

        rr(1,iend) = rr(1,ista) + cos(theta_p(ich)) * chord_p(ich)
        rr(2,iend) = rr(2,ista)
        rr(3,iend) = rr(3,ista) - sin(theta_p(ich)) * chord_p(ich)

      else

        rr(:,ista:iend) = w1*rrSection1 + w2*rrSection2

      end if

      w2 = spacing_weights ( ispace, iSpan-1, nelem_span_list(iRegion), & 
                            r_in_list(iRegion), r_ob_list(iRegion), y_ref_list(iRegion) )
      w1 = 1.0_wp - w2

      ! For flat meshes the section is not rotated, but the
      ! normal/tangent vector are

      if ( mesh_flat ) then
        theta_e(ich-1) = 0.5_wp*(w1*twist_list(iRegion) + w2*twist_list(iRegion+1)) &
                        + 0.5_wp*theta_p(ich)
      endif

    end do

  end do


! Save the fields related to the definition of actual airfoils
! airfoil_list(i) may be equal to 'interp' --> no def. of airfoil on the i-th section
!                                              --> inteprolation
  ! Check airfoil_list()
  if ( trim(airfoil_list(1)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The first "airfoil"&
          & cannot be set as "interp".')
  end if
  if ( trim(airfoil_list(nSections)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The last "airfoil"&
          & cannot be set as "interp".')
  end if

  nAirfoils = 0
  do i = 1 , nSections
    if ( trim(airfoil_list(i)) .ne. 'interp' ) nAirfoils = nAirfoils + 1
  end do

  allocate(airfoil_list_actual(nAirfoils))
  iAirfoil = 0
  do i = 1 , nSections
    if ( trim(airfoil_list(i)) .ne. 'interp' ) then
      iAirfoil = iAirfoil + 1
      airfoil_list_actual(iAirfoil) = trim(airfoil_list(i))
      call check_file_exists(airfoil_list_actual(iAirfoil), this_sub_name, &
            this_mod_name)
    end if
  end do

  i_aero2 = 1 ; ista = 1 ; iAirfoil = 0
  do iAirfoil = 1 , nAirfoils - 1

    i_aero1 = i_aero2
    do i = i_aero1 + 1 , nSections
      if ( airfoil_list(i) .ne. 'interp' ) then
        i_aero2 = i
        exit
      end if
    end do

    iend = sum(nelem_span_list(1:i_aero2-1))

    do iSpan = ista , iend
      i_airfoil_e(1, iSpan) = iAirfoil
      i_airfoil_e(2, iSpan) = iAirfoil+1

      normalised_coord_e(2,iSpan) = ( rr(2,(iSpan+1)*2) - rr(2,ista*2) ) / &
                                    ( rr(2,iend *2+1) - rr(2,ista*2) )

      if ( iSpan .ne. ista ) then ! else = 0.0_wp
        normalised_coord_e(1,iSpan) = normalised_coord_e(2,iSpan-1)
      end if
    end do

    ista = iend + 1

  end do

! === new-2019-02-06 ===



  ! optional output ----
  npoints_chord_tot = npoint_chord_tot
  ! optional output ----


end subroutine read_mesh_ll

!----------------------------------------------------------------------

!> Function that computes the weights used to interpolate date between
!! two points.
!!
!! The interval defined by the extrema A and B, is discretized by N segments
!! and N+1 points such that
!!
!!    x_i = A + w_i*(B-A)     with i in (0,N)
!!
!! where `w_i` is the weight generated by this function for the i-th point.
!!
!! If a uniform spacing is used the distribution becomes
!!
!!    x_i = A + i/N*(B*A)
!!
!! such that x_0 = A and x_N = B
!!
!! Four options are currently available
!!
!!   ISPACE = 0    Uniform spacing
!!   ISPACE = 1    cosine spacing
!!   ISPACE = 2    left cosine spacing
!!   ISPACE = 3    right cosine spacing
!!   ISPACE = 4    equalarea spacing
!!   ISPACE = 5    geoseries 
!!   ISPACE = 6    geoseries_ob
!!   ISPACE = 7    geoseries_ib

function spacing_weights ( itype, i, n, r_ib, r_ob, y_refinement ) result(w)

  integer,      intent(in)  :: itype, i, n
  real(wp),     intent(in)  :: r_ib, r_ob, y_refinement
  real(wp)                  :: w
  real(wp),     allocatable :: division(:), divisionIB(:), divisionOB(:)

  select case (itype)

  case(1) ! cosine
    w = 0.5_wp - 0.5_wp * cos( real(i,wp)*pi/ real(n,wp) )

  case(2) ! left cosine 
    w = sin( 0.5_wp*real(i,wp)*pi/ real(n,wp) )

  case(3) ! right cosine
    w = 1.0_wp - cos( 0.5_wp*real(i,wp)*pi/ real(n,wp) )

  case(4) ! equalarea 
    w =  sqrt(real(i,wp)/real(n,wp)) 

  case(5) ! geoseries
    call geoseries_both(0.0_wp, 1.0_wp, n, y_refinement, r_ib, r_ob, division)
    w = division(i + 1)

  case(6) ! geoseries_ob
    call geoseries(0.0_wp, 1.0_wp, n, r_ob, divisionIB, divisionOB) 
    w = divisionOB(i + 1)

  case(7) ! geoseries_ib
    call geoseries(0.0_wp, 1.0_wp, n, r_ib, divisionIB, divisionOB) 
    w = divisionIB(i + 1)

  case default ! uniform
    w = real(i,wp) / real(n,wp)
  end select

end function spacing_weights

!----------------------------------------------------------------------
!> Set offset between xac and the structural coupling node (useless)
subroutine read_xac_offset( filen, rr, xac )
  character(len=*)     , intent(in)  :: filen
  real(wp)             , intent(in)  :: rr(:,:)
  real(wp), allocatable, intent(out) :: xac(:)

  real(wp), allocatable :: y_coord(:), yin(:), xacin(:)

  integer :: np, i, io, ny, j
  integer :: fid = 21

  !> Set y_coord array
  np = size(rr,2) / 2
  allocate(y_coord(np))
  do i = 1, np
    y_coord(i) = rr(2,2*i-1)
  end do
  !> Non-dimensionalization of y_coord
  y_coord = ( y_coord - y_coord(1) ) / ( y_coord(np) - y_coord(1) )

  !> Read offset_xac_file
  open(unit=fid, file=trim(filen))
  ny = 0
  do
    read(fid,*,iostat=io)
    if ( io .lt. 0 ) then; exit
    else                 ; ny = ny + 1
    end if
  end do
  close(fid)

  allocate(yin(ny), xacin(ny))
  open(unit=fid, file=trim(filen))
  do i = 1, ny
    read(fid,*) yin(i), xacin(i)
  end do
  close(fid)

  !> Non-dimensionalization of yin array
  yin = ( yin - yin(1) ) / ( yin(ny) - yin(1) )

  !> Interpolation of xacin in xac
  allocate(xac(np))
  xac(1)  = xacin(1)
  xac(np) = xacin(ny)
  do i = 2, np-1
    do j = 1, ny-1
      if ( ( y_coord(i)-yin(j)   .gt. 0.0_wp ) .and. &
            ( y_coord(i)-yin(j+1) .le. 0.0_wp ) ) then
        xac(i) = (xacin(j+1) - xacin(j)) / &
                  (yin(j+1) - yin(j)) * &
                  (y_coord(i) - yin(j)) + xacin(j)
      end if
    end do
  end do


end subroutine read_xac_offset

!----------------------------------------------------------------------

end module mod_ll_io
