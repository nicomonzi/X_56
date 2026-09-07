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
!!          Andrea Colli
!!          Alberto Savino
!!          Alessandro Cocco
!!=========================================================================

!> Module to treat the most simple input-output from ascii formatted data
!! files
!!
!! BUG: older versions of cgns expose a method called "error" and so dust
!! error cannot be employed.
module mod_cgns_io

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  !error,
  warning, info, printout, new_file_unit, dust_abort

use mod_stringtools, only: &
  IsInList, StripSpaces

use cgns, only: & 
  cgsize_t, CG_OPEN_F, CG_ERROR_EXIT_F, CG_ERROR_PRINT_F,         &
  CG_ZONE_TYPE_F, CG_ZONE_READ_F, CG_NZONES_F,                    &
  CG_NBASES_F, CG_BASE_READ_F, CG_NSECTIONS_F, CG_SECTION_READ_F

!----------------------------------------------------------------------

implicit none

!----------------------------------------------------------------------
! CGNS ENUM 
!----------------------------------------------------------------------
integer, parameter :: NODE     =  2
integer, parameter :: BAR_2    =  3
integer, parameter :: BAR_3    =  4
integer, parameter :: TRI_3    =  5
integer, parameter :: TRI_6    =  6
integer, parameter :: QUAD_4   =  7
integer, parameter :: QUAD_8   =  8
integer, parameter :: QUAD_9   =  9
integer, parameter :: TETRA_4  = 10
integer, parameter :: TETRA_10 = 11
integer, parameter :: PYRA_5   = 12
integer, parameter :: PYRA_14  = 13
integer, parameter :: PENTA_6  = 14
integer, parameter :: PENTA_15 = 15
integer, parameter :: PENTA_18 = 16
integer, parameter :: HEXA_8   = 17
integer, parameter :: HEXA_20  = 18
integer, parameter :: HEXA_27  = 19
integer, parameter :: MIXED    = 20
integer, parameter :: NGON_n  = 21

integer, parameter :: Integer    = 2
integer, parameter :: RealSingle = 3
integer, parameter :: RealDouble = 4
integer, parameter :: Character  = 5 

integer, parameter :: Structured   =  2
integer, parameter :: Unstructured =  3

integer, parameter :: ALL_OK         = 0
integer, parameter :: ERROR          = 1
integer, parameter :: NODE_NOT_FOUND = 2
integer, parameter :: INCORRECT_PATH = 3

integer, parameter :: MODE_READ   = 0
integer, parameter :: MODE_WRITE  = 1
integer, parameter :: MODE_CLOSED = 2
integer, parameter :: MODE_MODIFY = 3

! External Fortran subroutine interfaces to cgns library, 
! not included in the cgnslib_f.h header
external :: CG_COORD_READ_F, CG_POLY_ELEMENTS_READ_F

public :: read_mesh_cgns

private

character(len=*), parameter :: this_mod_name = 'mod_cgns_io'
character(len=max_char_len) :: msg

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!>Read a mesh in a cgns format employing the cgns API
!!
!! WARNING: this is still experimental and based on an incomplete reverse
!! engineering of cgns mesh files, it shall be extended to be more reliable
!! and general
subroutine read_mesh_cgns(mesh_file, sectionNamesUsed, ee, rr, ee_virtual, rr_virtual)

 character(len=*), intent(in) :: mesh_file
 character(len=*), intent(in) :: sectionNamesUsed(:)
 integer, allocatable, intent(out) :: ee(:,:), ee_virtual(:,:)
 real(wp) , allocatable, intent(out) :: rr(:,:), rr_virtual(:,:)



 character(len=max_char_len), allocatable :: sectionNames(:)
 logical, allocatable :: selectedSection(:)
 integer :: index_file , ier  , nbase , ibase
 character(len=32) :: basename , zonename , sectionname, coordname(3)
 character(len=72) :: eltype
 integer :: icelldim , ndim , nzone , izone , id , isec , nsec
 integer :: ieltype ,  nbelem, nnod
 integer ::  nNodesUsed, posNode
 integer(cgsize_t) :: nNodes, nelem, idl
 integer ::  iprec , parent_flag
 integer(cgsize_t) :: nelem_zone
 integer :: nSectionsUsed, posSection
 logical :: sectionFound

 integer(cgsize_t) :: isize(3), nstart , nend , range_min, range_max


 real(wp) , allocatable :: coordinateList(:)
 integer, allocatable :: nodeMap(:), nodeList(:)

 integer, pointer     :: nmixed(:)
 integer(cgsize_t), pointer     :: pdata(:)
 integer(cgsize_t), allocatable :: elemcg(:,:)

 integer(cgsize_t), allocatable :: estart(:) , eend(:)
 integer, allocatable :: etype(:)
! integer, allocatable :: ee_cgns(:)

 integer :: i1 , i, ielem_section, iNode
 integer(cgsize_t) :: ielem

 integer(cgsize_t), allocatable :: ConnectOffset(:)

 character(len=*), parameter :: this_sub_name = 'read_mesh_cgns'

  ! Name of the arrays containint the grid coordinates
  coordname(1) = 'CoordinateX'
  coordname(2) = 'CoordinateY'
  coordname(3) = 'CoordinateZ'

  nSectionsUsed = size(sectionNamesUsed)

  call printout(nl//' Loading geometry from cgns file: '//trim(mesh_file))
  call CG_OPEN_F(trim(mesh_file), MODE_READ, index_file, ier)
  if ( ier .ne. ALL_OK ) then
    call printout(' CGNS error when reading cgns file!')
    call CG_ERROR_EXIT_F()
  end if

! Number of bases

  call CG_NBASES_F(INDEX_FILE, nbase, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
  if (nbase /= 1) then
    write(msg,'(A,I2)') ' Number of bases in the cgns file: ', nbase
    call printout(trim(msg))
    !call error(this_sub_name, this_mod_name,'More than one base in the &
    !  &CGNS file, not supported')
    call printout('More than one base in the CGNS file, not supported')
    call dust_abort()
  end if
  ibase = 1

! Base data

  call CG_BASE_READ_F(INDEX_FILE, ibase, basename, icelldim, ndim, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

!  Special case for 2d: if icelldim=2 we assume a pure 2d mesh with
!  coordinates in x and y only (not a surface mesh in 3d)

! Check.
!So far it makes work a mesh done from surface cad (not from solid cad)

!  if ((icelldim == 2) .and. (ndim > icelldim)) ndim = icelldim

!  Number of zones. Assume a zone is a region. Loop over the zones

  call CG_NZONES_F(INDEX_FILE, ibase, nzone, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

  if ( nzone .ne. 1 ) then
    write(msg,'(A,I2)') ' Number of zones in the cgns base: ', nzone
    call printout(trim(msg))
    !call error(this_sub_name, this_mod_name,'More than one zon in the &
    !  &CGNS file, not supported')
    call printout('More than one zone in the CGNS file, not supported')
    call dust_abort()
  end if

  do izone = 1,nzone

! *** edge format ***
!   call APPEND_CHILD_N(pmsh,'region',overwrite=.false.)
!   preg => NTH_CHILD(pmsh,izone,'region')
!   call APPEND_CHILD_L0(preg,'region_name','volume_elements')

!  Read zone data

    call CG_ZONE_READ_F(INDEX_FILE, ibase, izone, zonename, isize, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

!  Check if unstructured zone. todo

    call CG_ZONE_TYPE_F(INDEX_FILE, ibase, izone, id, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
    if ( id /= Unstructured) then
      !call error(this_sub_name, this_mod_name,'No unstructured grid was &
      !  &found in the CGNS file')
      call printout('No unstructured grid was found in the CGNS file')
      call dust_abort()
    end if

    nNodes = isize(1)

    call printout('CGNS zone recap:'//nl)
    write(msg,'(A,I9,A)') ' Number of nodes:      ', isize(1),nl
    call printout(trim(msg))
    write(msg,'(A,I9,A)') ' Number of cells:      ', isize(2),nl
    call printout(trim(msg))
    write(msg,'(A,I2,A)') ' Coordinate dimension: ', ndim,nl
    call printout(trim(msg))
    write(msg,'(A,I2,A)') ' Mesh elements dimension: ', icelldim,nl
    call printout(trim(msg))

#if (DUST_PRECISION==1)
    iprec = RealSingle
#elif (DUST_PRECISION==2)  
    iprec = RealDouble
#endif 

!  Find number of sections
!
    call CG_NSECTIONS_F(INDEX_FILE, ibase, izone, nsec, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
    write(msg,'(A,I3)') ' Number of sections in the zone: ', nsec
    call printout(trim(msg))

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! find the total number of elements and section names
    allocate(estart(nsec)) ; estart = 0
    allocate(eend(nsec))   ; eend   = 0
    allocate(etype(nsec))  ; etype = 0
    allocate(sectionNames(nsec));
    allocate(selectedSection(nsec)); selectedSection = .false.
    do isec = 1 , nsec

      call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                             ieltype, nstart, nend, nbelem, parent_flag, ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
      nelem = nend - nstart + 1

      eend(isec)   = nend
      estart(isec) = nstart
      etype(isec) = ieltype

      call StripSpaces(sectionname)
      sectionNames(isec) = sectionname

      write(msg,'(A,I3,A,I10,2A)') ' Available section number:', isec,&
                           '   Number of elements: ', nelem, &
                           '   Name: ', trim(sectionname)//nl
      call printout(trim(msg))

    end do


    ! Get number of elements considering the selection of sections
    ! some checks for cell elements (ignoring point and lines)
    if (nSectionsUsed > 0) then
      nelem_zone = 0
      write(msg,'(A,I3)') ' Number of sections to be loaded: ', nSectionsUsed
      call printout(trim(msg))

      do isec = 1, nSectionsUsed
        call printout(' Section '//trim(sectionNamesUsed(isec))//&
          &' will be loaded'//nl)

        sectionFound = IsInList(sectionNamesUsed(isec), sectionNames, posSection)
        if (sectionFound) then

          if (etype(isec) .lt. TRI_3) then
            call printout('ERROR selected section with &
            &name '//trim(sectionNamesUsed(isec))//' does not contain cells')
            call dust_abort()
          endif
            
          selectedSection(posSection) = .true.
          nelem_zone = nelem_zone + eend(posSection) - estart(posSection) + 1
        else
          call printout('ERROR during CGNS section loading: section with &
            &name '//trim(sectionNamesUsed(isec))//' not found'//nl//&
            & '        see above for the list of available sections')
            call dust_abort()
        endif
      enddo
    else
      call printout('All found CGNS sections will be loaded')
      ! All sections will be selected
      nelem_zone = 0
      do isec = 1, nsec
        if (etype(isec) .ge. TRI_3) then
          selectedSection(isec) = .true.
          nelem_zone = nelem_zone + eend(isec) - estart(isec) + 1
        else
          call printout('Section with &
            &name '//trim(sectionNames(isec))//' skipped')
        endif
      end do
      if (nelem_zone .eq. 0) then
        call printout('ERROR during CGNS section loading: no suitable cell &
         &elements found in any section')
            call dust_abort()
      endif
    endif

    write(msg,'(A,I0,A)') 'Total number of loaded elements: ', nelem_zone,nl
    call printout(trim(msg))

    allocate(ee(4,nelem_zone)) ; ee = 0
    allocate(ee_virtual(4,nelem_zone)) ; ee_virtual = 0
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Save connectivity
    ielem = 0
    do isec = 1, nsec
      if (selectedSection(isec)) then
        call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                               ieltype, nstart, nend, nbelem, parent_flag, ier)
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
        nelem = nend - nstart + 1

        ! Set dust element name.
        select case (ieltype)
        case(NODE)
          eltype = 'node'
          nnod = 1
        case(BAR_2)
          eltype = 'bar2'
          nnod = 2
        case(TRI_3)
          eltype = 'tria3'
          nnod = 3
        case(QUAD_4)
          eltype = 'quad4'
          nnod = 4
        case(MIXED)
          eltype = 'mixed'
          nnod = 9
        case default
          ! only TRI and QUAD are supported
          ! node and BAR should have been already excluded
          write(msg,*) ' error: unsupported CGNS element type: ', ieltype, &
                     'NONE '!, trim(ElementTypeName(ieltype))
          call printout(trim(msg))
          call dust_abort()
        end select


        !  Read cgns element data
        if (eltype == 'mixed') then
          allocate(elemcg(nnod*nelem,1))
        else
          allocate(elemcg(nnod, nelem))
        end if
        elemcg = 0_cgsize_t

        if (parent_flag==1) then
          allocate(pdata(nelem*4))
        else
          allocate(pdata(3))
        end if
        pdata = 0_cgsize_t

        allocate(ConnectOffset(nelem+1))
        ConnectOffset = 0_cgsize_t
        call CG_POLY_ELEMENTS_READ_F(INDEX_FILE, ibase, izone, isec, &
                                elemcg, ConnectOffset, pdata, ier)
        if ( ier /= ALL_OK ) call CG_ERROR_PRINT_F()
        call check_cgns_overflow(elemcg)

        deallocate(ConnectOffset)

        !  Store data, mixed elements get special treatment
        if (eltype == 'mixed') then
          allocate(nmixed(mixed))
          nmixed = 0
          ielem_section = 0
          I = 1
          do while (ielem_section < nelem)
            ielem = ielem + 1
            ielem_section = ielem_section + 1
            select case (elemcg(I,1))
            case(NODE)
              nnod = 1
            case(BAR_2)
              nnod = 2
            case(TRI_3)
              nnod = 3
            case(QUAD_4)
              nnod = 4
            case(TETRA_4)
              nnod = 4
            case(PYRA_5) ! TODO are these useful?
              nnod = 5
            case(PENTA_6)
              nnod = 6
            case(HEXA_8)
              nnod = 8
            case default
              write(msg,*) ' error: unsupported element type: ', elemcg(I,1), &
                          'NONE'!, trim(ElementTypeName(elemcg(I,1)))
              call printout(trim(msg))
              call dust_abort()
            end select
            ee(1:nnod,ielem) = int(elemcg(I+1:I+nnod,1),4)
            nmixed(elemcg(I,1)) = nmixed(elemcg(I,1)) + 1
            I = I + nnod + 1
          end do
          deallocate(nmixed)

        else
          do idl = 1, eend(isec) - estart(isec) + 1
            ielem = ielem + 1
            ee(1:nnod,ielem) = int(elemcg(:,idl),4)
          end do
        end if


        deallocate(elemcg)

      endif
    end do ! Loop on sections


    ! If a selection of elements is performed performs also a selection on nodes
    if (nSectionsUsed > 0) then
      nNodesUsed = 0
      allocate(nodeMap(nNodes)); nodeMap = 0
      allocate(nodeList(nNodes)); nodeList = 0

      ! Select nodes used
      do ielem = 1, nelem_zone
        do i1 = 1, 4
          if (ee(i1,ielem) > 0) then
            posNode = ee(i1,ielem)
            if (nodeMap(posNode) .eq. 0) then
              nNodesUsed = nNodesUsed + 1;
              nodeMap(posNode) = nNodesUsed
              nodeList(nNodesUsed) = posNode
            endif
            ee(i1,ielem) = nodeMap(posNode)
          endif
        enddo
      enddo

      ! Get nodes
      allocate(coordinateList(nNodes))
      allocate(rr(ndim,nNodesUsed))
      do i = 1, ndim
        coordinateList = 0.0_wp
        range_min = 1
        range_max = isize(1)
        call CG_COORD_READ_F(INDEX_FILE, ibase, izone, coordname(i), &
                             iprec, range_min, range_max, coordinateList, ier); 
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

        do iNode = 1, nNodesUsed
          rr(i,iNode) = coordinateList(nodeList(iNode))
        enddo

      enddo

    ! No selection of sections performed, read all nodes
    else
      allocate(coordinateList(nNodes))
      allocate(rr(ndim,nNodes))
      do i = 1, ndim
        coordinateList = 0.0_wp;

        !TODO: passing rr(1,:) is not a contiguous memory location and
        !needs a temporary, and rises a warning. Consider either using
        !an explicit temporary or a transposition of rr
        range_min = 1
        range_max = isize(1)
        call CG_COORD_READ_F(INDEX_FILE, ibase, izone, coordname(i), &
                             iprec, range_min, range_max, coordinateList, ier);
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

        rr(i,:) = coordinateList
      enddo
      

    endif

  end do ! Loop on zones
  
  ! Virtual mesh
  allocate(rr_virtual(size(rr,1), size(rr, 2))); rr_virtual = 0.0_wp
  rr_virtual = rr
  ee_virtual = ee 
  
  deallocate(estart, eend, etype)

end subroutine read_mesh_cgns

!----------------------------------------------------------------------

!>Check the cgns  arrays for overflow
!!
!! Newst cgns versions support 8 Byte indexing. We surely do not.
!! We handle all that arrives from CGNS in the native integer that the current
!! library uses, but at the end we must cast some of the received data to
!! our native int4 types. So we need to check for overflow before doing that.
subroutine check_cgns_overflow(cgns_array)
 integer(cgsize_t), intent(in) :: cgns_array(:,:)

 integer :: prot_int

 if(any(cgns_array .gt. huge(prot_int)) .or. &
  any(cgns_array .lt. -huge(prot_int))) then

    call printout('ERROR while importing CGNS data: the file contains data &
     &bigger than a 4 Byte integer which is used internally for indexing in &
     &DUST. Either a really gargantuan mesh is being loaded (and there is &
     & no reason to use it in DUST) or there is some issue with the CGNS file')
    call dust_abort()

 endif

end subroutine check_cgns_overflow

!----------------------------------------------------------------------

end module mod_cgns_io
