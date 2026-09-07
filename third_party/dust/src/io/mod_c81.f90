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


! SUBROUTINEs:
!   read_c81_table
!   interp2d

!> Module to treat the c81 tables
module mod_c81

use mod_param, only: &
  wp, max_char_len , nl, pi

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, new_file_unit

implicit none

!-----------------------------------

! public and private

!-----------------------------------

type t_aero_2d_par
  real(wp) , allocatable :: par1(:)
  real(wp) , allocatable :: par2(:)
  real(wp) , allocatable :: cf(:,:)
end type t_aero_2d_par

!-----------------------------------

type t_aero_2d_tab

  !> value of the third parameter, Re
  real(wp)                          :: Re
  
  type(t_aero_2d_par) , allocatable :: coeff(:)
  real(wp) , allocatable            :: dclda(:,:)

  real(wp) , allocatable            :: clmax(:) , alcl0(:) , cdmin (:)
  real(wp) , allocatable            :: clstall_neg(:) , clstall_pos(:)
  real(wp) , allocatable            :: alstall_neg(:) , alstall_pos(:)

end type t_aero_2d_tab

!-----------------------------------

!> Tables containing aerodynamic coefficients

! Tabulated data for lifting line elements
type t_aero_tab

  !> airfoil file
  character(len=max_char_len) :: airfoil_file

  !> airfoil data id
  integer :: id

  !> tables cl(a,M) for each Re number
  type(t_aero_2d_tab) , allocatable :: aero_coeff(:)

end type t_aero_tab

!-----------------------------------

character(len=*), parameter :: this_mod_name='mod_c81'
character(len=max_char_len) :: msg

contains

!-----------------------------------
! Vahana input file style '../vahana_input/vahana4.inp'
! HP: the aerodynamic coefficients are given for the same values
! of al_i, M_i
subroutine read_c81_table ( filen , coeff )
  character(len=max_char_len) , intent(in) :: filen
  type(t_aero_tab) , intent(inout) :: coeff

  character(len=max_char_len) :: line , string
  integer :: cl1 , cl2 , cd1 , cd2 , cm1 , cm2 !  table_dim(2)
  integer :: cp1(3) , cp2(3)
  real(wp) :: dummy
  integer :: dummy_int , iblnk
  integer :: fid ,ierr, i1 , iRe , nRe , i_c
  real(wp) :: Re

  integer , parameter :: n_coeff = 3   ! cl , cd , cm

  ! Reynolds effect corrections ---
  real(wp) , parameter :: alMin = -180.0_wp , alMax = 180.0_wp
  real(wp) :: alcl0tmp, alcl0
  integer :: iAl , iMa
  real(wp) :: al1 , al2 , c1 , c2
  integer :: ind_cl0 , istallp , istallm, stall_found
  ! Reynolds effect corrections ---

  character(len=*), parameter :: this_sub_name = 'read_c81_table'

  call new_file_unit(fid,ierr)
  open(unit=fid,file=trim(adjustl(filen)))
  ! First tree lines containing unintelligilbe parameters
  read(fid,*) nRe , dummy_int , dummy_int ! 2 0 0
  read(fid,*) ! 0 1
  read(fid,*) ! 0.158 0.158

  allocate( coeff%aero_coeff(nRe) )

  do iRe = 1 , nRe

    ! For each Reynolds number (1 Reynolds number, 3 coefficients)
    read(fid,*) ! COMMENT#1
    read(fid,*) Re  , dummy
    read(fid,'(A)') line
      iblnk = index(line,' ')
      string = trim(adjustl(line(iblnk:)))
      read(string,'(6I2)') cl2 , cl1 , cd2 , cd1 , cm2 , cm1
      !                    NMa , Nal , NMa , Nal , NMa , Nal
    cp1 = (/ cl1 , cd1 , cm1 /)
    cp2 = (/ cl2 , cd2 , cm2 /)

    ! check that Reynolds number are defined in increasing order
    if ( iRe .gt. 1 ) then
      do i1 = 1 , iRe-1
        if ( coeff%aero_coeff(i1)%Re .ge. Re ) then
          write(msg,*) ' in file: ' , trim(adjustl(filen)) , ' Reynolds number used to define'// &
                      ' .c81 tables must be in increasing order, but ',nl, &
                      ' -> Re(', iRe , '): ' , Re, nl,&
                      ' -> Re(', i1  , '): ' , coeff%aero_coeff(i1)%Re
          call error(this_sub_name, this_mod_name, msg)
        end if
      end do
    end if

    ! Reynolds number
    coeff%aero_coeff(iRe)%Re = Re
    allocate(coeff%aero_coeff(iRe)%coeff(n_coeff))
    ! allocate coefficient structures
    do i_c = 1 , n_coeff
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%par1(cp1(i_c)))
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%par2(cp2(i_c)))
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%cf(cp1(i_c),cp2(i_c)))
    end do


    do i_c = 1 , n_coeff
      read(fid,*) coeff%aero_coeff(iRe)%coeff(i_c)%par2   ! Mach numbers
      do i1 = 1 , cp1(i_c)
        read(fid,*) coeff%aero_coeff(iRe)%coeff(i_c)%par1(i1) , &    ! alpha
                    coeff%aero_coeff(iRe)%coeff(i_c)%cf(i1,:)        ! adim.coeff
      end do
    end do

    ! === table containing the partial derivative dCL/dalpha, ===
    ! === to be used in iterative processes                   ===
    ! 1st order finite difference at the extreme values of alpha, 2nd order for inner alphas
    allocate(coeff%aero_coeff(iRe)%dclda(cl1,cl2)) ! coeff1: cl

    coeff%aero_coeff(iRe)%dclda(1,:) = &
      ( coeff%aero_coeff(iRe)%coeff(1)%cf(2,:) - coeff%aero_coeff(iRe)%coeff(1)%cf(1,:) ) / &
      ( coeff%aero_coeff(iRe)%coeff(1)%par1(2) - coeff%aero_coeff(iRe)%coeff(1)%par1(1) )
    coeff%aero_coeff(iRe)%dclda(cl1,:) = &
      ( coeff%aero_coeff(iRe)%coeff(1)%cf(cl1,:) - coeff%aero_coeff(iRe)%coeff(1)%cf(cl1-1,:) ) / &
      ( coeff%aero_coeff(iRe)%coeff(1)%par1(cl1) - coeff%aero_coeff(iRe)%coeff(1)%par1(cl1-1) )
    do i1 = 2 , cl1-1 ! loop over par1: alpha
      coeff%aero_coeff(iRe)%dclda(i1,:) = &
        ( coeff%aero_coeff(iRe)%coeff(1)%cf(i1+1,:) - coeff%aero_coeff(iRe)%coeff(1)%cf(i1-1,:) ) / &
        ( coeff%aero_coeff(iRe)%coeff(1)%par1(i1+1) - coeff%aero_coeff(iRe)%coeff(1)%par1(i1-1) )
    end do

    if ( sim_param%llReynoldsCorrections .or. sim_param%vl_correction) then
      ! === Parameters for corrections for Reynolds number effects ===
      ! Find clmax, alcl0, cdmin for each (Re,M)
      allocate(coeff%aero_coeff(iRe)%clmax(cl2)) ; coeff%aero_coeff(iRe)%clmax = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%alcl0(cl2)) ; coeff%aero_coeff(iRe)%alcl0 = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%cdmin(cd2)) ; coeff%aero_coeff(iRe)%cdmin = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%clstall_pos(cd2)) ; coeff%aero_coeff(iRe)%clstall_pos = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%clstall_neg(cd2)) ; coeff%aero_coeff(iRe)%clstall_neg = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%alstall_pos(cd2)) ; coeff%aero_coeff(iRe)%alstall_pos = -333.0_wp
      allocate(coeff%aero_coeff(iRe)%alstall_neg(cd2)) ; coeff%aero_coeff(iRe)%alstall_neg = -333.0_wp

      do iMa = 1 , cl2

        ! --- clmax --- cl = coeff(1)
        coeff%aero_coeff(iRe)%clmax(iMa) = maxval( coeff%aero_coeff(iRe)%coeff(1)%cf(:,iMa) )
        ! --- al(cl=0) ---
        alcl0 = alMin
        ind_cl0 = -333
        do iAl = 1 , cl1
          if ( ( coeff%aero_coeff(iRe)%coeff(1)%par1(iAl) .ge. alMin ) .and. &
                ( coeff%aero_coeff(iRe)%coeff(1)%par1(iAl) .lt. alMax ) ) then
            if ( coeff%aero_coeff(iRe)%coeff(1)%cf(iAl  ,iMa) * &
                coeff%aero_coeff(iRe)%coeff(1)%cf(iAl+1,iMa) .le. 1e-16_wp ) then

              al1 = coeff%aero_coeff(iRe)%coeff(1)%par1(iAl)
              al2 = coeff%aero_coeff(iRe)%coeff(1)%par1(iAl+1)
              c1  = coeff%aero_coeff(iRe)%coeff(1)%cf(iAl  ,iMa)
              c2  = coeff%aero_coeff(iRe)%coeff(1)%cf(iAl+1,iMa)

              alcl0tmp = al1 + (al2-al1) * (-c1)/(c2-c1)
              if(abs(alcl0tmp)<abs(alcl0)) then
                alcl0 = alcl0tmp
                ind_cl0 = iAl  ! index where cl changes sign: later used to find stall+ and stall-
              endif

            end if
          end if
        end do ! alpha
        coeff%aero_coeff(iRe)%alcl0(iMa) = alcl0

        ! --- cdmin --- cd = coeff(2)
        coeff%aero_coeff(iRe)%cdmin(iMa) = minval( coeff%aero_coeff(iRe)%coeff(2)%cf(:,iMa) )
        
        ! Throw an error if cl0 is not found TODO improve dealing with edge cases
        if (ind_cl0 .le. 0) call error(this_sub_name, this_mod_name, 'Error in finding CL = 0')

        ! --- positive and negative stall: alpha and cl ---
        istallp = ind_cl0 ; istallm = ind_cl0
        ! positive stall
        stall_found = 0

        ! ---
        do while ( ( istallp+1 .lt. cl1 ) .and. &
                    ( coeff%aero_coeff(iRe)%coeff(1)%cf(istallp+1,iMa) .gt. &
                      coeff%aero_coeff(iRe)%coeff(1)%cf(istallp  ,iMa) ) .and. &
                    ( stall_found .eq. 0 ) &
                    )
          istallp = istallp + 1

        end do
        ! negative stall
        stall_found = 0
        do while ( ( istallm-1 .gt. 1 ) .and. &
                  ( coeff%aero_coeff(iRe)%coeff(1)%cf(istallm-1,iMa) .lt. &
                    coeff%aero_coeff(iRe)%coeff(1)%cf(istallm  ,iMa) ) .and. &
                  ( stall_found .eq. 0 ) )
          istallm = istallm - 1

        end do

        coeff%aero_coeff(iRe)%alstall_pos(iMa) = coeff%aero_coeff(iRe)%coeff(1)%par1(istallp)
        coeff%aero_coeff(iRe)%alstall_neg(iMa) = coeff%aero_coeff(iRe)%coeff(1)%par1(istallm)
        coeff%aero_coeff(iRe)%clstall_pos(iMa) = coeff%aero_coeff(iRe)%coeff(1)%cf(istallp,iMa)
        coeff%aero_coeff(iRe)%clstall_neg(iMa) = coeff%aero_coeff(iRe)%coeff(1)%cf(istallm,iMa)

      end do ! Mach number

    end if 
      ! === Parameters for corrections for Reynolds number effects ===

  end do ! Reynolds number

  close(fid)

end subroutine read_c81_table

!-----------------------------------
!> Routine to compute aerodynamic coefficients,
!>  given the geometry ( airfoil_data, csi , airfoil_id) and
!>        the aerodynamic parameters( aero_par = (/ al , M , Re /)

! TODO: add some checks : if reyn1 .eq. reyn2 to avoid singularity
!       clear implementation
!       ...
subroutine interp_aero_coeff ( airfoil_data ,  csi , airfoil_id , &
                              aero_par_in, aero_coeff , &
                                                    dcl_da )
  type(t_aero_tab) , intent(in) :: airfoil_data(:)
  real(wp) , intent(in) :: csi
  integer  , intent(in) :: airfoil_id(2)
  real(wp) , intent(in) :: aero_par_in(3) ! (/al,M,Re/)
  real(wp) , allocatable , intent(out) :: aero_coeff(:)
  real(wp) :: aero_par(3) ! (/al,M,Re/)

  real(wp) , optional    , intent(out) :: dcl_da
  real(wp) :: dcl_da1,  dcl_da2
! newton cleaning
! real(wp)               , intent(out) :: dclda
! real(wp)               , intent(out) :: al0
! newton cleaning

  real(wp) :: al , mach , reyn
  real(wp) :: cf1(3) , cf2(3)
  real(wp) , allocatable :: coeff1(:) , coeff2(:)
  real(wp) , allocatable :: coeff_airfoil(:,:)
  real(wp) , allocatable ::dcl_da_airfoil(:)
  integer :: nRe, i_a, id_a

  real(wp) :: reyn1 , reyn2
  integer  :: irey

  ! Reynolds effect correction ----
  real(wp) :: n_fact , k_fact
  real(wp) :: aero_par_re(2)
  integer :: nmach , imach
  real(wp) :: machend , mach1 , mach2
  ! Reynolds effect correction ----

  ! dclda ----
  real(wp) :: al01
  ! dclda ----

  ! n factor for the corrections of aerodynamic coeffs, with (Re/Re_table)^n
  if ( sim_param%llReynoldsCorrections ) then
    n_fact = sim_param%llReynoldsCorrectionsNfact
  else ! not used -> initalised to 0.0_wp
    n_fact = 0.0_wp
  end if

  ! workaround to keep IN/OUT/INOUT old declarations
  aero_par = aero_par_in

  al   = aero_par(1)
  mach = aero_par(2)
  reyn = aero_par(3)

  ! al must be cyclic in [-180.0,180.0]
  al = -real(floor((al+180.0_wp)/360.0_wp),wp) * 360.0_wp + al
  aero_par(1) = al

  cf1 = 0.0_wp
  cf2 = 0.0_wp

  if ( .not. sim_param%llReynoldsCorrections ) then

    do i_a = 1 , 2

      id_a = airfoil_id(i_a)
      nRe = size(airfoil_data(id_a)%aero_coeff)

      ! Some checks ----
      if ( nRe .eq. 1 ) then ! .c81 defined just for one Reynolds number

        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(1)%coeff , &
                                                    aero_par(1:2) , coeff1 , &
                                                                  dcl_da1 )

        if ( .not. allocated(coeff_airfoil) ) then
          allocate(coeff_airfoil(2,size(coeff1)))
        end if
        if ( .not. allocated(dcl_da_airfoil) ) then
          allocate(dcl_da_airfoil(2))
        end if

        coeff_airfoil(i_a,:) = coeff1

        ! === dcl_da derivative ===
        dcl_da_airfoil(i_a) = dcl_da1

      else ! .c81 defined for more than one Reynolds number

        ! Some checks ----
        if ( reyn .lt. airfoil_data(id_a)%aero_coeff(1)%Re ) then
          reyn1 = airfoil_data(id_a)%aero_coeff(1)%Re
          reyn2 = airfoil_data(id_a)%aero_coeff(2)%Re
          irey  = 1
        else if ( reyn .gt. airfoil_data(id_a)%aero_coeff(nRe)%Re ) then
          reyn1 = airfoil_data(id_a)%aero_coeff(nRe-1)%Re
          reyn2 = airfoil_data(id_a)%aero_coeff(nRe)%Re
          irey  = nRe-1
        else

          irey  = 1
          reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
          do while ( ( reyn .ge. reyn1 ) .and. ( irey .lt. nRe ) )
            irey = irey + 1
            reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
          end do
          irey = irey - 1
          reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
          reyn2 = airfoil_data(id_a)%aero_coeff(irey+1)%Re

        end if

        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey  )%coeff , &
                                                        aero_par(1:2) , coeff1 , &
                                                                        dcl_da1 )
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey+1)%coeff , &
                                                        aero_par(1:2) , coeff2 , &
                                                                        dcl_da2 )

        if ( .not. allocated(coeff_airfoil) ) allocate(coeff_airfoil(2,size(coeff1)))
        coeff_airfoil(i_a,:) = ( coeff1 * ( reyn2 - reyn ) + coeff2 * ( reyn - reyn1 ) ) /  &
                              ( reyn2 - reyn1 )

        ! === dcl_da derivative ===
        if ( present( dcl_da ) ) then
          if ( .not. allocated(dcl_da_airfoil) ) allocate(dcl_da_airfoil(2))
          dcl_da_airfoil(i_a) = ( dcl_da1 * ( reyn2 - reyn ) + dcl_da2 * ( reyn - reyn1 ) ) /  &
                                ( reyn2 - reyn1 )
        end if

      end if

    end do

    allocate( aero_coeff(size(coeff_airfoil,2)) )
    aero_coeff = coeff_airfoil(1,:) * ( 1.0_wp-csi ) + coeff_airfoil(2,:) * csi

    if ( present( dcl_da ) ) then
      dcl_da = dcl_da_airfoil(1) * ( 1.0_wp-csi ) + dcl_da_airfoil(2) * csi
    end if

  else ! new: taking into account Reynolds effects

    do i_a = 1 , 2

      id_a = airfoil_id(i_a)
      nRe = size(airfoil_data(id_a)%aero_coeff)

      if ( reyn .le. airfoil_data(id_a)%aero_coeff( 1 )%Re ) then

        irey = 1
        k_fact = ( reyn / airfoil_data(id_a)%aero_coeff(irey)%Re ) ** n_fact

        ! aero_par taking into account the Reynolds effect:
        ! --- find al(cl=0), for the desired mach number ---
        nmach = size(airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2)
        imach = 1
        mach1   = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        machend = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(nmach)
        do while ( (  mach .ge. mach1 ) .and. &
                    ( imach .lt. nmach  ) )
          imach = imach + 1
          mach1 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        end do
        imach = imach - 1
        mach1 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        mach2 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach+1)

        al01 = airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach) + &
               (mach-mach1)/(mach2-mach1) * &
              ( airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach+1) - &
                airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach  ) )

        ! --- correct alpha ---
        aero_par_re    = aero_par(1:2)
        aero_par_re(1) = ( aero_par(1) - al01 ) / k_fact + al01
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey)%coeff , &
                                                        aero_par_re , coeff1 , &
                                                                      dcl_da1 )

        if ( .not. allocated( coeff_airfoil) )  allocate( coeff_airfoil(2,size(coeff1)))
        if ( .not. allocated(dcl_da_airfoil) )  allocate(dcl_da_airfoil(2))
        coeff_airfoil(i_a,1) = coeff1(1) * k_fact

        ! --- cd , cm ---
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey)%coeff , &
                                                      aero_par(1:2) , coeff1  )

        coeff1(2) = coeff1(2) / k_fact  ! cd correction
        coeff_airfoil(i_a,2:3) = coeff1(2:3)

        if ( present( dcl_da ) ) then
          dcl_da1 = dcl_da1 * k_fact     ! correct dcl_da derivative with k_fact
          dcl_da_airfoil(i_a) = dcl_da1
        end if

      else if ( reyn .ge. airfoil_data(id_a)%aero_coeff(nRe)%Re ) then
        ! --- Reynolds effects with semi-empirical laws ---
        irey = nRe
        k_fact = ( reyn / airfoil_data(id_a)%aero_coeff(irey)%Re ) ** n_fact

        ! TODO: write a subroutine to find the range of parameters
        ! aero_par taking into account the Reynolds effect:
        ! --- find al(cl=0), for the desired mach number ---
        nmach = size(airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2)
        imach = 1
        mach1   = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        machend = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(nmach)
        do while ( ( mach .ge. mach1 ) .and. &
                  ( imach .lt. nmach  ) )
          imach = imach + 1
          mach1 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        end do
        imach = imach - 1
        mach1 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach)
        mach2 = airfoil_data(id_a)%aero_coeff(irey)%coeff(1)%par2(imach+1)

        al01 = airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach) + &
               (mach-mach1)/(mach2-mach1) * &
              ( airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach+1) - &
                airfoil_data(id_a)%aero_coeff(irey)%alcl0(imach  ) )

        ! --- cl ---
        aero_par_re    = aero_par(1:2)
        aero_par_re(1) = ( aero_par(1) - al01 ) / k_fact + al01    ! correct alpha
        
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey)%coeff , &
                                                        aero_par_re , coeff1 , &
                                                                      dcl_da1 )

        if ( .not. allocated( coeff_airfoil) )  allocate(coeff_airfoil(2,size(coeff1)))
        if ( .not. allocated(dcl_da_airfoil) )  allocate(dcl_da_airfoil(2))
        coeff_airfoil(i_a,1) = coeff1(1) * k_fact

        ! --- cd , cm ---
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey)%coeff , &
                                                    aero_par(1:2) , coeff1    )

        coeff1(2) = coeff1(2) / k_fact  ! cd correction
        coeff_airfoil(i_a,2:3) = coeff1(2:3)

        if ( present( dcl_da ) ) then
          dcl_da1 = dcl_da1 * k_fact     ! correct dcl_da derivative with k_fact
          dcl_da_airfoil(i_a) = dcl_da1
        end if

      else
        ! --- linear interpolation ---
        ! find the smallest range of Re in tables including the desired Re
        irey  = 1
        reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
        do while ( ( reyn .ge. reyn1 ) .and. ( irey .lt. nRe ) )
          irey = irey + 1
          reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
        end do
        irey = irey - 1
        reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
        reyn2 = airfoil_data(id_a)%aero_coeff(irey+1)%Re

        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey  )%coeff , &
                                                        aero_par(1:2) , coeff1 , &
                                                                        dcl_da1 )
        call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey+1)%coeff , &
                                                        aero_par(1:2) , coeff2 , &
                                                                        dcl_da2 )

        if ( .not. allocated(coeff_airfoil) )  allocate(coeff_airfoil(2,size(coeff1)))
        coeff_airfoil(i_a,:) = ( coeff1 * ( reyn2 - reyn ) + coeff2 * ( reyn - reyn1 ) ) /  &
                                ( reyn2 - reyn1 )

        ! === dcl_da derivative ===
        if ( present( dcl_da ) ) then
          dcl_da_airfoil(i_a) = ( dcl_da1 * ( reyn2 - reyn ) + dcl_da2 * ( reyn - reyn1 ) ) /  &
                                ( reyn2 - reyn1 )
        end if

      end if

    end do

    allocate( aero_coeff(size(coeff_airfoil,2)) )
    aero_coeff = coeff_airfoil(1,:) * ( 1.0_wp-csi ) + coeff_airfoil(2,:) * csi

    if ( present( dcl_da ) ) then
      dcl_da = dcl_da_airfoil(1) * ( 1.0_wp-csi ) + dcl_da_airfoil(2) * csi
    end if

  end if


  ! deallocate
  deallocate(coeff_airfoil)
  if ( allocated(dcl_da_airfoil) ) deallocate(dcl_da_airfoil)

! ----
! write the routine to find the parameter range
!  contains
!
!  subroutine find_param_range()
!
!
!  end subroutine
! ----

end subroutine interp_aero_coeff

!-----------------------------------
!
subroutine interp2d_aero_coeff ( aero_coeff , x , c , dcl_da )
  type(t_aero_2d_par) , intent(in) :: aero_coeff(:)
  real(wp), intent(in)  :: x(:)     ! aero_par (al,M,Re)
  real(wp), allocatable , intent(out) :: c(:)
  real(wp), optional    , intent(out) :: dcl_da

  integer  :: n1 , n2 , nc
  real(wp) :: csi1 , csi2
  real(wp) :: phi1 , phi2 , phi3 , phi4

  integer :: ic , i1 , i2

  character(len=*), parameter :: this_sub_name='interp2d_aero_coeff'



  nc = size(aero_coeff)

  allocate(c(nc)) ; c = 0.0_wp

  ! only 2d parameters are allowed (al,M)
  if ( size(x) .ne. 2 ) then
    call error(this_sub_name, this_mod_name, 'Attempting to 2D interpolate&
    & data with more dimensions')
  end if

  ! Do it once ...
  ! Check dimensions ---------
  do ic = 1 , nc

    n1 = size(aero_coeff(ic)%par1)
    n2 = size(aero_coeff(ic)%par2)

    if ( (size(aero_coeff(ic)%cf,1).ne.n1) .or. &
        (size(aero_coeff(ic)%cf,2).ne.n2)) then
      write(msg,*) 'Error in the aerodynamic coefficients tables, force &
                &coefficient matrix size',shape(aero_coeff(ic)%cf), 'different&
                & from corresponding alpha and mach vectors sizes',n1,n2
      call error(this_sub_name, this_mod_name, msg)
    end if

    ! Check range of parameters the parameters are supposed to be defined in an
    ! increasing order
    if ( x(1) .lt. aero_coeff(ic)%par1(1) .or. &
        x(1) .gt. aero_coeff(ic)%par1(n1) ) then
        write(msg,*) 'Trying to interpolate aerodynamic coefficients at an &
        &angle ',x(1),' which is outside the table span of angles from ', &
        aero_coeff(ic)%par1(1),'to',aero_coeff(ic)%par1(n1)
        call error(this_sub_name, this_mod_name, msg)
    endif
    if ( x(2) .lt. aero_coeff(ic)%par2(1) .or. &
        x(2) .gt. aero_coeff(ic)%par2(n2) ) then
        write(msg,*) 'Trying to interpolate aerodynamic coefficients at a &
        &Mach number ',x(2),' which is outside the table span of Mach numbers &
        &from ', aero_coeff(ic)%par2(1),'to',aero_coeff(ic)%par2(n2)
        call error(this_sub_name, this_mod_name, msg)
    endif
    ! Check dimensions ---------

    i1 = 1
    do while ( (aero_coeff(ic)%par1(i1) .le. x(1)) )
      i1 = i1 + 1
    end do
    i2 = 1
    do while ( (aero_coeff(ic)%par2(i2) .le. x(2)) )
      i2 = i2 + 1
    end do

    ! par1: alpha , par2: Mach
    csi1 = 2.0_wp * ( x(1) - 0.5_wp*(aero_coeff(ic)%par1(i1-1)+aero_coeff(ic)%par1(i1)) ) / &
          (aero_coeff(ic)%par1(i1)-aero_coeff(ic)%par1(i1-1))

    csi2 = 2.0_wp * ( x(2) - 0.5_wp*(aero_coeff(ic)%par2(i2-1)+aero_coeff(ic)%par2(i2)) ) / &
          (aero_coeff(ic)%par2(i2)-aero_coeff(ic)%par2(i2-1))

    phi1 =   0.25_wp * ( 1.0_wp + csi1 ) * ( 1.0_wp + csi2 )
    phi2 =   0.25_wp * ( 1.0_wp - csi1 ) * ( 1.0_wp + csi2 )
    phi3 =   0.25_wp * ( 1.0_wp - csi1 ) * ( 1.0_wp - csi2 )
    phi4 =   0.25_wp * ( 1.0_wp + csi1 ) * ( 1.0_wp - csi2 )

    c(ic) = phi1 * aero_coeff(ic)%cf(i1  ,i2  ) + &
            phi2 * aero_coeff(ic)%cf(i1-1,i2  ) + &
            phi3 * aero_coeff(ic)%cf(i1-1,i2-1) + &
            phi4 * aero_coeff(ic)%cf(i1  ,i2-1)

    ! === Compute dcl_da derivative (rad), if required for Jacobian matrix ===
    if ( present(dcl_da) ) then
      if ( ic .eq. 1 ) then ! cl table -> compute dcl_da
        dcl_da = ( aero_coeff(ic)%cf( i1  ,i2 ) + aero_coeff(ic)%cf( i1  ,i2-1) &
                - aero_coeff(ic)%cf( i1-1,i2 ) - aero_coeff(ic)%cf( i1-1,i2-1) ) / &
                ( 2.0_wp * ( aero_coeff(ic)%par1( i1 ) - aero_coeff(ic)%par1( i1-1 ) ) * &
                            pi / 180.0_wp )
      end if
    end if

  end do


end subroutine interp2d_aero_coeff

!-----------------------------------

end module mod_c81
