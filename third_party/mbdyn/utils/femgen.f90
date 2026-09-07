!===============================================================================
! FEMGEN - Finite Element Modal Generator for MBDyn
!===============================================================================
!
! PURPOSE
!   Converts modal analysis results from MSC Nastran SOL 103 (normal modes)
!   or SOL 106 (nonlinear static + modal) into a modal-data file (.fem)
!   suitable for MBDyn's modal joint element.
!
! INPUT FILES
!   mbdyn_modal.mat : Nastran OUTPUT4 binary file containing matrices:
!                    - MHH    : Modal mass matrix (n_modes x n_modes)
!                    - KHH    : Modal stiffness matrix (n_modes x n_modes)
!                    - LUMPMS : Lumped mass diagonal (6*n_nodes x 1)
!
!   <name>.op2 : Standard Nastran OP2 binary file containing:
!               - GEOM1  : Grid coordinates (requires PARAM,OGEOM,YES)
!               - OUG1   : Eigenvector mode shapes
!               - OUG1   : Static displacements (for SOL 106)
!
! OUTPUT FILE
!   <name>.fem : MBDyn modal data file with record groups:
!                1. Header (revision, counts)
!                2. Node list
!                3. Initial modal displacements (optional)
!                4. Initial modal velocities (optional)
!                5-7. Nodal X, Y, Z coordinates
!                8. Mode shapes
!                9. Modal mass matrix
!                10. Modal stiffness matrix
!                11. Lumped mass diagonal
!
! PUBLIC INTERFACE
!   integer(c_int) function femgen(outname,
!                                    is_modal_displacement,
!                                    is_modal_velocity,
!                                    mass_direction) bind(C)
!     outname : Output filename (C string, max 72 chars)
!     is_modal_displacement    : Include initial modal displacements (0=no, 1=yes)
!     is_modal_velocity        : Include initial modal velocities (0=no, 1=yes)
!     mass_direction           : Lumped mass treatment (-1=auto, 0=none, 1-3=use component)
!     Returns                  : 0 on success, non-zero on failure
!
! AUTHORS
!   Giuseppe Quaranta: <giuseppe.quaranta@polimi.it>
!   Alessandro Cocco : <alessandro.cocco@polimi.it>
!
!===============================================================================

module femgen_mod
   use, intrinsic :: iso_fortran_env, only: &
      int32, real32, real64
   use, intrinsic :: iso_c_binding, only: &
      c_int, c_int32_t, c_char, c_null_char, c_float, c_double

   implicit none
   private

   public :: femgen_c_wrapper

   !---------------------------------------------------------------------------
   ! Module Constants
   !---------------------------------------------------------------------------
   !> Tolerance for snapping small values to zero in output
   real(real64), parameter :: FEMGEN_ZERO_TOL = 1.0e-16_real64

   !> Maximum number of FEM nodes supported
   integer, parameter :: MAX_NODES = 20000

   !> Maximum words per OUTPUT4 matrix record
   integer, parameter :: OUTPUT4_MAX_RECORD_WORDS = 2000000

   !> Maximum words per OP2 data block
   integer, parameter :: OP2_MAX_BLOCK_WORDS = 5000000

   !> Maximum length of input/output filenames (including null terminator)
   integer, parameter :: MAX_NAME_LENGTH = 256

   !> OP2 table name constants (little-endian 4-byte words)
   integer(int32), parameter :: TNAME_GEOM1_W1 = int(z'4D4F4547', int32)  ! 'GEOM'
   integer(int32), parameter :: TNAME_GEOM1_W2 = int(z'20202031', int32)  ! '1   '
   integer(int32), parameter :: TNAME_OUG1_W1 = int(z'3147554F', int32)  ! 'OUG1'

   !---------------------------------------------------------------------------
   ! Derived Types
   !---------------------------------------------------------------------------

   !> Container for complete modal model data
   !>
   !> Holds all information extracted from Nastran files needed to generate
   !> the MBDyn .fem output file.
   type :: femgen_model
      integer :: n_modes = 0                       !< Number of vibration modes
      integer :: n_nodes = 0                       !< Number of FEM nodes
      integer :: phix_nrow = 0                     !< PHIX (eigenvectors) row dimension
      integer :: phix_ncol = 0                     !< PHIX (eigenvectors) column dimension
      integer(int32), allocatable :: node_ids(:)   !< External node ID labels [n_nodes]
      integer(int32), allocatable :: internal_ids(:) !< Internal NASTRAN IDs [n_nodes]
      real(real64), allocatable :: xg(:)           !< Global X coordinates [n_nodes]
      real(real64), allocatable :: yg(:)           !< Global Y coordinates [n_nodes]
      real(real64), allocatable :: zg(:)           !< Global Z coordinates [n_nodes]
      real(real64), allocatable :: mhh(:, :)        !< Modal mass matrix [n_modes, n_modes]
      real(real64), allocatable :: khh(:, :)        !< Modal stiffness matrix [n_modes, n_modes]
      real(real64), allocatable :: lump_vec(:)     !< Lumped mass diagonal [6*n_nodes]
      real(real64), allocatable :: phix(:, :)       !< Mode shape matrix [6*n_nodes, n_modes]
   end type femgen_model

   !> Nastran OP2 binary file reader
   !>
   !> Provides buffered reading of Nastran OP2 format files, handling
   !> the complex record structure with variable-length data blocks.
   !> Note: OP2 is the standard Nastran binary output format containing
   !> geometry tables (GEOM1, GEOM2, etc.) and results (OUG1, OES1, etc.)
   type :: op2_reader
      integer :: file_unit = -1                 !< Fortran unit number for input file
      integer :: log_unit = 6                   !< Unit for diagnostic messages (stdout)
      logical :: record_active = .false.        !< True if within a logical record
      logical :: buffer_valid = .false.         !< True if buffer contains unread data
      integer :: buffer_pos = 0                 !< Current read position in buffer
      integer :: buffer_len = 0                 !< Number of valid words in buffer
      integer(int32), allocatable :: buffer(:)  !< Word buffer for block data
      logical :: silent = .false.               !< Suppress error messages (for EOF drain)
   contains
      procedure :: init => op2_init
      procedure :: read_block_header => op2_read_block_header
      procedure :: read_words => op2_read_words
      procedure :: skip_record => op2_skip_record
   end type op2_reader

contains

   !===========================================================================
   ! C BINDING WRAPPER
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Main entry point for C/C++ callers
   !>
   !> Converts Nastran modal data to MBDyn format.
   !> Reads <input_name>.op2 for geometry/modes and <input_name>.mat for matrices.
   !>
   !> @param[in] input_name                Nastran output basename (C string)
   !> @param[in] output_name               Output .fem filename (empty = use input_name)
   !> @param[in] is_modal_displacement     Include initial modal displacements (0/1)
   !> @param[in] is_modal_velocity         Include initial modal velocities (0/1)
   !> @param[in] mass_direction            Lumped mass component selector (-1=auto, 0-3)
   !> @return                              0=success, non-zero=error code
   !---------------------------------------------------------------------------
   integer(c_int) function femgen_c_wrapper(input_name_c, &
                                            output_name_c, &
                                            is_modal_displacement, &
                                            is_modal_velocity, &
                                            mass_direction) bind(C, name='femgen')
      character(kind=c_char), intent(in) :: input_name_c(MAX_NAME_LENGTH)
      character(kind=c_char), intent(in) :: output_name_c(MAX_NAME_LENGTH)
      integer(c_int32_t), value, intent(in) :: is_modal_displacement
      integer(c_int32_t), value, intent(in) :: is_modal_velocity
      integer(c_int32_t), value, intent(in) :: mass_direction

      type(femgen_model)                   :: model
      integer                              :: log_unit, idx, status
      integer                              :: input_len, output_len
      integer                              :: opt_init_disp, opt_init_vel, opt_mass_idx
      character(len=MAX_NAME_LENGTH - 1)     :: input_name, output_name
      character(len=:), allocatable        :: op2_path, mat_path, output_path
      logical                              :: mat_exists

      femgen_c_wrapper = 0_c_int
      log_unit = 6
      input_name = ' '
      output_name = ' '

      ! Convert C strings to Fortran
      do idx = 1, MAX_NAME_LENGTH - 1
         if (input_name_c(idx) == c_null_char) exit
         input_name(idx:idx) = input_name_c(idx)
      end do

      do idx = 1, MAX_NAME_LENGTH - 1
         if (output_name_c(idx) == c_null_char) exit
         output_name(idx:idx) = output_name_c(idx)
      end do

      write (log_unit, '(A)') 'FEMGEN: translation of Binary NASTRAN modal data'
      write (log_unit, '(A)') '        for use in MBDyn''s modal joint element'
      write (log_unit, '(A)') ' '

      input_len = len_trim(input_name)
      if (input_len <= 0) then
         write (log_unit, '(A)') 'ERROR: No input name provided'
         femgen_c_wrapper = 1_c_int
         return
      end if

      ! Strip .op2 or .mat extension if user provided the full filename
      if (input_len > 4) then
         if (lower1(input_name(input_len - 3:input_len - 3)) == '.' .and. &
             lower1(input_name(input_len - 2:input_len - 2)) == 'o' .and. &
             lower1(input_name(input_len - 1:input_len - 1)) == 'p' .and. &
             lower1(input_name(input_len:input_len)) == '2') then
            input_len = input_len - 4  ! Strip .op2 extension
         else if (lower1(input_name(input_len - 3:input_len - 3)) == '.' .and. &
                  lower1(input_name(input_len - 2:input_len - 2)) == 'm' .and. &
                  lower1(input_name(input_len - 1:input_len - 1)) == 'a' .and. &
                  lower1(input_name(input_len:input_len)) == 't') then
            input_len = input_len - 4  ! Strip .mat extension
         end if
      end if

      ! Construct OP2 path (always <model>.op2)
      op2_path = trim(input_name(1:input_len))//'.op2'

      ! Find MAT file: try <model>.mat first, then fallback to mbdyn_modal.mat
      mat_path = trim(input_name(1:input_len))//'.mat'
      inquire (file=mat_path, exist=mat_exists)
      if (.not. mat_exists) then
         mat_path = 'mbdyn_modal.mat'
         inquire (file=mat_path, exist=mat_exists)
         if (.not. mat_exists) then
            write (log_unit, '(A)') 'ERROR: Cannot find .mat file'
            write (log_unit, '(A,A)') '  Tried: ', trim(input_name(1:input_len))//'.mat'
            write (log_unit, '(A)') '  Tried: mbdyn_modal.mat'
            femgen_c_wrapper = 1_c_int
            return
         end if
      end if

      write (log_unit, '(A,A)') ' Input OP2 file: ', trim(op2_path)
      write (log_unit, '(A,A)') ' Input MAT file: ', trim(mat_path)

      ! Determine output filename
      output_len = len_trim(output_name)
      if (output_len > 0) then
         ! Use user-specified output name
         output_path = trim(output_name(1:output_len))
      else
         ! Default: use input basename
         output_path = trim(input_name(1:input_len))
      end if

      ! Add .fem extension if needed
      if (.not. ends_with_fem(output_path)) then
         output_path = output_path//'.fem'
      end if

      ! Read Nastran input files
      call femgen_read(mat_path, op2_path, model, log_unit, status)
      if (status /= 0) then
         femgen_c_wrapper = int(status, c_int)
         return
      end if

      ! Write MBDyn output file
      write (log_unit, '(A,A)') 'Output to file: ', trim(output_path)
      opt_init_disp = int(is_modal_displacement)
      opt_init_vel = int(is_modal_velocity)
      opt_mass_idx = int(mass_direction)
      call femgen_write(output_path, model, opt_init_disp, opt_init_vel, opt_mass_idx, status)
      if (status /= 0) then
         femgen_c_wrapper = int(status, c_int)
         return
      end if
   end function femgen_c_wrapper

   !===========================================================================
   ! MAIN FILE READER
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Read modal data from Nastran OP2 and OUTPUT4 files
   !>
   !> Unified reading for both SOL 103 (linear modal) and SOL 106 (nonlinear
   !> static + modal). Reads geometry from GEOM1 table and mode shapes from OUG1.
   !> For SOL 106, also reads static displacements and computes deformed coords.
   !>
   !> @param[in]    mat_path   Path to OUTPUT4 file (.mat - matrices)
   !> @param[in]    tab_path   Path to OP2 file (geometry and mode shapes)
   !> @param[inout] model      Model structure to populate
   !> @param[in]    log_unit   Fortran unit for diagnostic output
   !> @param[out]   status     0=success, non-zero=error
   !---------------------------------------------------------------------------
   subroutine femgen_read(mat_path, tab_path, model, log_unit, status)
      character(len=*), intent(in)    :: mat_path
      character(len=*), intent(in)    :: tab_path
      type(femgen_model), intent(inout) :: model
      integer, intent(in)    :: log_unit
      integer, intent(out)   :: status

      type(op2_reader) :: reader
      real(real64), allocatable :: temp_matrix(:, :), lump_matrix(:, :)
      real(real64), allocatable :: static_disp_x(:), static_disp_y(:), static_disp_z(:)
      real :: timestamp(7)
      integer(int32), allocatable :: word_buffer(:)
      integer(int32) :: file_label(2), table_name(2)
      character(len=8) :: matrix_name

      integer :: op2_unit, out4_unit, io_ret, num_modes, words_read, ret_code
      integer :: node_idx, mode_idx, row_count, col_count
      integer :: lump_rows, lump_cols, read_status
      integer :: approach_code, oug1_table_count, total_tables
      logical :: op2_is_open, out4_is_open
      logical :: has_geom1, has_static_disp, has_eigenvectors
      logical :: done_scanning, is_eof
      status = 0
      op2_is_open = .false.
      out4_is_open = .false.
      has_geom1 = .false.
      has_static_disp = .false.
      has_eigenvectors = .false.
      op2_unit = 11
      out4_unit = 12

      ! Open OUTPUT4 (matrix) file
      open (out4_unit, file=mat_path, status='old', form='unformatted', iostat=io_ret)
      if (io_ret /= 0) then
         write (log_unit, *) 'FEMGEN error while opening input file ', trim(mat_path), ', IRET=', io_ret
         call read_fail_and_cleanup(2, 'FEMGEN failed to open OUTPUT4 file')
         return
      end if
      out4_is_open = .true.

      ! Open OP2 (results/geometry) file
      open (op2_unit, file=tab_path, status='old', form='unformatted', iostat=io_ret)
      if (io_ret /= 0) then
         write (log_unit, *) 'FEMGEN error while opening input file ', trim(tab_path), ', IRET=', io_ret
         call read_fail_and_cleanup(2, 'FEMGEN failed to open OP2 file')
         return
      end if
      op2_is_open = .true.

      ! Read MHH (modal mass matrix)
      call output4_read_matrix(out4_unit, matrix_name, temp_matrix, row_count, col_count, read_status)
      if (read_status /= 0) then
         call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
         return
      end if
      if (matrix_name(1:3) /= 'MHH') then
         call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
         return
      end if

      num_modes = row_count
      model%n_modes = num_modes
      allocate (model%mhh(num_modes, num_modes))
      model%mhh = temp_matrix(1:num_modes, 1:num_modes)
      deallocate (temp_matrix)

      ! Read KHH (modal stiffness matrix)
      call output4_read_matrix(out4_unit, matrix_name, temp_matrix, row_count, col_count, read_status)
      if (read_status /= 0) then
         call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
         return
      end if
      if (matrix_name(1:3) /= 'KHH') then
         call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
         return
      end if

      allocate (model%khh(num_modes, num_modes))
      model%khh = 0.0_real64
      ! Handle diagonal eigenvalue format (ncol=5)
      if (row_count == num_modes .and. col_count == 5) then
         do node_idx = 1, num_modes
            do mode_idx = 1, num_modes
               model%khh(mode_idx, node_idx) = model%mhh(mode_idx, node_idx)*temp_matrix(node_idx, 1)
            end do
         end do
         ! Symmetrize
         do mode_idx = 1, num_modes
            do node_idx = mode_idx + 1, num_modes
               model%khh(mode_idx, node_idx) = 0.5_real64*(model%khh(mode_idx, node_idx) + model%khh(node_idx, mode_idx))
               model%khh(node_idx, mode_idx) = model%khh(mode_idx, node_idx)
            end do
         end do
      else
         model%khh(1:row_count, 1:col_count) = temp_matrix(1:row_count, 1:col_count)
      end if
      deallocate (temp_matrix)

      ! Read LUMPMS (lumped mass diagonal)
      call output4_read_matrix(out4_unit, matrix_name, lump_matrix, lump_rows, lump_cols, read_status)
      if (read_status /= 0) then
         call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
         return
      end if
      if (matrix_name(1:6) /= 'LUMPMS') then
         call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
         return
      end if

      allocate (model%lump_vec(lump_rows))
      model%lump_vec = lump_matrix(:, 1)
      deallocate (lump_matrix)

      ! Prepare word buffer for OP2 reading
      allocate (word_buffer(12*MAX_NODES))
      word_buffer = 0_int32

      ! Initialize OP2 reader and read first table
      call reader%init(op2_unit, log_unit, file_label, read_status)
      if (read_status /= 0) then
         call read_fail_and_cleanup(3, 'Error initializing OP2 reader')
         return
      end if
      call reader%read_block_header(table_name, timestamp, read_status)
      if (read_status /= 0) then
         call read_fail_and_cleanup(3, 'Error reading OP2 table header')
         return
      end if

      ! =========================================================================
      ! Table Scanning: Unified OP2 reading for both SOL 103 and SOL 106
      ! OP2 contains: GEOM1 (original coords) + OUG1 (eigenvectors or static+eigenvectors)
      ! =========================================================================

      oug1_table_count = 0

      ! Scan for GEOM1 (original coordinates) and OUG1 tables
      do while (table_name(1) /= TNAME_OUG1_W1)

         ! Check for GEOM1 specifically (both words must match)
         if (table_name(1) == TNAME_GEOM1_W1 .and. table_name(2) == TNAME_GEOM1_W2) then
            if (.not. has_geom1) then
               call read_geom1_grids(reader, word_buffer, model%n_nodes, model%node_ids, &
                                     model%xg, model%yg, model%zg, read_status)
               if (read_status == 0 .and. model%n_nodes > 0) then
                  has_geom1 = .true.
                  allocate (model%internal_ids(model%n_nodes))
                  model%internal_ids = model%node_ids
               end if
            end if
         else
            ! Skip this table
            call op2_skip_table(reader, word_buffer, ret_code)
            if (ret_code == -99) then
               call read_fail_and_cleanup(3, 'Error: EOF while scanning tables')
               return
            end if
         end if

         call reader%read_block_header(table_name, timestamp, read_status)
         if (read_status /= 0) then
            call read_fail_and_cleanup(3, 'Error reading OP2 table header')
            return
         end if
      end do

      ! Verify we have geometry from GEOM1
      if (.not. has_geom1) then
         call read_fail_and_cleanup(3, 'Error: GEOM1 table not found in OP2 file - ensure PARAM,OGEOM,YES is set')
         return
      end if

      ! Scan all OUG1 tables to get eigenvectors and optionally static displacements
      if (table_name(1) == TNAME_OUG1_W1) then
         done_scanning = .false.
         oug1_table_count = 0
         total_tables = 0
         is_eof = .false.

         do while (.not. done_scanning .and. .not. is_eof)
            total_tables = total_tables + 1

            if (table_name(1) == TNAME_OUG1_W1) then
               oug1_table_count = oug1_table_count + 1

               ! Read IDENT record to check approach_code
               call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
               if (ret_code == -99) then
                  is_eof = .true.
                  exit
               end if

               if (ret_code /= 2) then
                  approach_code = word_buffer(1)

                  if (approach_code == 102 .or. approach_code == 101) then
                     ! Static displacement OUG1 - read and store (keep last one)

                     if (allocated(static_disp_x)) deallocate (static_disp_x)
                     if (allocated(static_disp_y)) deallocate (static_disp_y)
                     if (allocated(static_disp_z)) deallocate (static_disp_z)

                     call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
                     if (ret_code == -99) then
                        is_eof = .true.
                        exit
                     end if

                     call parse_oug1_static_data(word_buffer, words_read, model%node_ids, &
                                                 static_disp_x, static_disp_y, static_disp_z)
                     has_static_disp = .true.

                  else if (approach_code == 22) then
                     ! Eigenvector OUG1 - read mode shapes
                     call read_oug1_eigenvectors(reader, word_buffer, model, ret_code)
                     if (ret_code == -99) then
                        is_eof = .true.
                        exit
                     end if
                     has_eigenvectors = .true.
                  end if
               end if

               ! Skip to end of this OUG1 table
               do while (ret_code /= 2 .and. ret_code /= -99)
                  call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
               end do
               if (ret_code == -99) then
                  is_eof = .true.
                  exit
               end if

            else
               ! Not OUG1 - skip this table
               call op2_skip_table(reader, word_buffer, ret_code)
               if (ret_code == -99) then
                  is_eof = .true.
                  exit
               end if
            end if

            ! Read next table header (silence reader as EOF is expected)
            reader%silent = .true.
            call reader%read_block_header(table_name, timestamp, read_status)
            reader%silent = .false.
            if (read_status /= 0) then
               done_scanning = .true.
            end if
         end do

         ! Apply static displacements to get deformed geometry
         if (has_static_disp) then
            model%xg = model%xg + static_disp_x
            model%yg = model%yg + static_disp_y
            model%zg = model%zg + static_disp_z
            deallocate (static_disp_x, static_disp_y, static_disp_z)
         end if

         ! Check if we have all data from OP2 scan
         if (allocated(model%node_ids) .and. allocated(model%phix)) then
            call cleanup_local()
            return
         end if
      else
         call read_fail_and_cleanup(3, 'Error: OUG1 table not found in OP2 file')
         return
      end if

      ! All processing complete
      call cleanup_local()
      return

   contains

      subroutine read_fail_and_cleanup(error_code, error_msg)
         integer, intent(in) :: error_code
         character(len=*), intent(in) :: error_msg
         status = error_code
         if (len_trim(error_msg) > 0) write (log_unit, '(A)') trim(error_msg)
         call cleanup_local()
      end subroutine read_fail_and_cleanup

      subroutine cleanup_local()
         if (allocated(temp_matrix)) deallocate (temp_matrix)
         if (allocated(lump_matrix)) deallocate (lump_matrix)
         if (allocated(word_buffer)) deallocate (word_buffer)
         if (allocated(static_disp_x)) deallocate (static_disp_x)
         if (allocated(static_disp_y)) deallocate (static_disp_y)
         if (allocated(static_disp_z)) deallocate (static_disp_z)
         if (op2_is_open) then; close (op2_unit); op2_is_open = .false.; end if
         if (out4_is_open) then; close (out4_unit); out4_is_open = .false.; end if
      end subroutine cleanup_local

   end subroutine femgen_read

   !===========================================================================
   ! MAIN FILE WRITER
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Write MBDyn modal data file (.fem)
   !>
   !> Generates a formatted text file containing all modal data in the format
   !> expected by MBDyn's modal joint element.
   !>
   !> @param[in]    output_path       Output file path
   !> @param[inout] model             Model data (may be modified for mass fix)
   !> @param[in]    write_init_disp   Include initial modal displacements (0/1)
   !> @param[in]    write_init_vel    Include initial modal velocities (0/1)
   !> @param[in]    mass_component    Lumped mass treatment option
   !> @param[out]   status            0=success, non-zero=error
   !---------------------------------------------------------------------------
   subroutine femgen_write(output_path, model, write_init_disp, write_init_vel, mass_component, status)
      character(len=*), intent(in)    :: output_path
      type(femgen_model), intent(inout) :: model
      integer, intent(in)    :: write_init_disp, write_init_vel, mass_component
      integer, intent(out)   :: status

      integer :: file_unit, io_ret, node_idx, mode_idx, dof_idx, row_end
      real(real64) :: mass_value

      status = 0
      file_unit = 10

      open (file_unit, file=output_path, iostat=io_ret)
      if (io_ret /= 0) then
         write (*, *) 'FEMGEN error opening output file, IRET=', io_ret
         status = 2
         return
      end if

      ! Record Group 1: Header
      write (file_unit, '(A24)') '** MBDyn MODAL DATA FILE'
      write (file_unit, '(A18)') '** NODE SET "ALL" '
      write (file_unit, *) ' '; write (file_unit, *) ' '
      write (file_unit, '(A25)') '** RECORD GROUP 1, HEADER'
      write (file_unit, '(A31, A39)') '**   REVISION,  NODE,  NORMAL, ', 'ATTACHMENT, CONSTRAINT, REJECTED MODES.'
      write (file_unit, '(6X, A4, 5I10)') 'REV0', model%n_nodes, model%n_modes, 0, 0, 0

      ! Record Group 2: Node list
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A34,A9)') '** RECORD GROUP 2, FINITE ELEMENT ', 'NODE LIST'
      do node_idx = 1, model%n_nodes, 6
         row_end = min(node_idx + 5, model%n_nodes)
         write (file_unit, '(6I10)') (int(model%node_ids(mode_idx)), mode_idx=node_idx, row_end)
      end do

      ! Record Group 3: Initial modal displacements (optional)
      if (write_init_disp /= 0) then
         write (file_unit, '(A2)') '**'
         write (file_unit, '(A33, A13)') '** RECORD GROUP 3, INITIAL MODAL ', 'DISPLACEMENTS'
         write (file_unit, '(500(1X,1PE17.10))') (0.0*mode_idx, mode_idx=1, model%n_modes)
      end if

      ! Record Group 4: Initial modal velocities (optional)
      if (write_init_vel /= 0) then
         write (file_unit, '(A2)') '**'
         write (file_unit, '(A32,A10)') '** RECORD GROUP 4, INITIAL MODAL ', 'VELOCITIES'
         write (file_unit, '(500(1X,1PE17.10))') (0.0*mode_idx, mode_idx=1, model%n_modes)
      end if

      ! Record Groups 5-7: Nodal coordinates
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A38)') '** RECORD GROUP 5, NODAL X COORDINATES'
      do node_idx = 1, model%n_nodes
         write (file_unit, '(E17.10)') snap_zero(real(model%xg(node_idx), kind=real64))
      end do
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A38)') '** RECORD GROUP 6, NODAL Y COORDINATES'
      do node_idx = 1, model%n_nodes
         write (file_unit, '(E17.10)') snap_zero(real(model%yg(node_idx), kind=real64))
      end do
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A38)') '** RECORD GROUP 7, NODAL Z COORDINATES'
      do node_idx = 1, model%n_nodes
         write (file_unit, '(E17.10)') snap_zero(real(model%zg(node_idx), kind=real64))
      end do

      ! Record Group 8: Mode shapes
      ! Since OUG1 is required, phix always contains full data for all nodes
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A30)') '** RECORD GROUP 8, MODE SHAPES'

      ! Verify PHIX dimensions (should be 6*n_nodes rows, n_modes cols from OUG1)
      if (model%phix_nrow /= 6*model%n_nodes .or. model%phix_ncol /= model%n_modes) then
         write (*, *) 'PHIX dimensions mismatch: expected', 6*model%n_nodes, 'x', model%n_modes, &
            ' got', model%phix_nrow, 'x', model%phix_ncol
         status = 3
         close (file_unit)
         return
      end if

      ! Write all mode shapes (PHIX stored as rows=DOFs, cols=modes)
      do mode_idx = 1, model%n_modes
         write (file_unit, '(A26,I2)') '**    NORMAL MODE SHAPE # ', mode_idx
         do node_idx = 1, model%n_nodes
            write (file_unit, '(6(1X,1PE17.10))') &
               (snap_zero(model%phix(6*(node_idx - 1) + dof_idx, mode_idx)), dof_idx=1, 6)
         end do
      end do

      ! Record Group 9: Modal mass matrix
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A36)') '** RECORD GROUP 9, MODAL MASS MATRIX'
      do mode_idx = 1, model%n_modes
         write (file_unit, '(500(1X,1PE17.10))') (snap_zero(model%mhh(mode_idx, node_idx)), node_idx=1, model%n_modes)
      end do

      ! Record Group 10: Modal stiffness matrix
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A42)') '** RECORD GROUP 10, MODAL STIFFNESS MATRIX'
      do mode_idx = 1, model%n_modes
         write (file_unit, '(500(1X,1PE17.10))') (snap_zero(model%khh(mode_idx, node_idx)), node_idx=1, model%n_modes)
      end do

      ! Process lumped mass based on mass_component option
      if (mass_component == -1) then
         ! Auto-detect: check consistency across translational DOFs
         do node_idx = 0, size(model%lump_vec) - 6, 6
            if (abs(model%lump_vec(node_idx + 1) - model%lump_vec(node_idx + 2)) > FEMGEN_ZERO_TOL .or. &
                abs(model%lump_vec(node_idx + 1) - model%lump_vec(node_idx + 3)) > FEMGEN_ZERO_TOL) then
               mass_value = model%lump_vec(node_idx + 1)
               if (abs(mass_value) <= FEMGEN_ZERO_TOL) mass_value = model%lump_vec(node_idx + 2)
               if (abs(mass_value) <= FEMGEN_ZERO_TOL) mass_value = model%lump_vec(node_idx + 3)
               do dof_idx = 1, 3
                  if (abs(model%lump_vec(node_idx + dof_idx)) >= FEMGEN_ZERO_TOL .and. &
                      abs(model%lump_vec(node_idx + dof_idx) - mass_value) > FEMGEN_ZERO_TOL) then
                     write (*, *) 'Inconsistent lumped mass values for node #', node_idx/6 + 1
                     status = 3
                     close (file_unit)
                     return
                  end if
               end do
               do dof_idx = 1, 3
                  model%lump_vec(node_idx + dof_idx) = snap_zero(mass_value)
               end do
            end if
         end do
      else if (mass_component /= 0) then
         ! Use specified component (1, 2, or 3) for all translational DOFs
         do node_idx = 0, size(model%lump_vec) - 6, 6
            mass_value = model%lump_vec(node_idx + mass_component)
            do dof_idx = 1, 3
               model%lump_vec(node_idx + dof_idx) = snap_zero(mass_value)
            end do
         end do
      end if

      ! Record Group 11: Lumped mass diagonal
      write (file_unit, '(A2)') '**'
      write (file_unit, '(A38,A12)') '** RECORD GROUP 11, DIAGONAL OF LUMPED', ' MASS MATRIX'
      do node_idx = 0, size(model%lump_vec) - 1, 6
         write (file_unit, '(500(1X,1PE17.10))') (snap_zero(model%lump_vec(node_idx + dof_idx)), dof_idx=1, 6)
      end do
      write (file_unit, '(A14)') '** END OF FILE'
      close (file_unit)
   end subroutine femgen_write

   !===========================================================================
   ! OUTPUT4 MATRIX READER
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Read a matrix from Nastran OUTPUT4 binary file
   !>
   !> Parses a single matrix block from an OUTPUT4 file, handling both dense
   !> and sparse column formats, single and double precision data.
   !>
   !> @param[in]  file_unit    Fortran unit number (file must be open)
   !> @param[out] matrix_name  8-character matrix name from header
   !> @param[out] matrix       Allocated matrix data (nrow x ncol)
   !> @param[out] nrow         Number of rows
   !> @param[out] ncol         Number of columns
   !> @param[out] error_code   0=success, 1=bad type, 2=bad dims, 3=overflow, 4=I/O error
   !---------------------------------------------------------------------------
   subroutine output4_read_matrix(file_unit, matrix_name, matrix, nrow, ncol, error_code)
      integer, intent(in) :: file_unit
      character(len=8), intent(out) :: matrix_name
      real(real64), allocatable, intent(out) :: matrix(:, :)
      integer, intent(out) :: nrow, ncol
      integer, intent(out) :: error_code

      integer(int32) :: form_code, data_type       !< Matrix format and precision type
      integer(int32) :: col_idx, start_row         !< Current column and starting row
      integer(int32) :: num_words                  !< Words in current record
      integer        :: io_status, col_num, word_idx
      integer(int32), allocatable :: word_buffer(:)

      error_code = 0
      allocate (word_buffer(0))

      ! Read matrix header: ncol, nrow, form, type, name
      read (file_unit, iostat=io_status) ncol, nrow, form_code, data_type, matrix_name
      if (io_status /= 0) then
         call fail_with_code(4)
         return
      end if

      ! Validate dimensions
      if (ncol <= 0 .or. nrow <= 0) then
         call fail_with_code(2)
         return
      end if

      ! Check data type (1=single, 2=double)
      if (data_type > 2) then
         rewind (file_unit)
         call fail_with_code(1)
         return
      end if

      allocate (matrix(nrow, ncol))
      matrix = 0.0_real64

      ! Read column data
      do col_num = 1, ncol
         read (file_unit, iostat=io_status) col_idx, start_row, num_words
         if (io_status /= 0) then
            call fail_with_code(4)
            return
         end if
         backspace (file_unit)

         if (num_words < 0 .or. num_words > OUTPUT4_MAX_RECORD_WORDS) then
            call fail_with_code(3)
            return
         end if

         call ensure_buffer_size(int(num_words))
         read (file_unit, iostat=io_status) col_idx, start_row, num_words, (word_buffer(word_idx), word_idx=1, num_words)
         if (io_status /= 0) then
            call fail_with_code(4)
            return
         end if

         if (col_idx > ncol) exit

         if (start_row == 0) then
            call fill_sparse_column(nrow, data_type, col_idx, word_buffer(1:num_words), matrix)
         else
            call fill_dense_column(nrow, data_type, col_idx, start_row, word_buffer(1:num_words), matrix)
         end if
      end do

      ! Check for trailer record
      read (file_unit, iostat=io_status) col_idx, start_row, num_words
      if (io_status == 0) then
         if (col_idx <= ncol) backspace (file_unit)
      end if

      call cleanup()

   contains
      !> Set error code and cleanup
      subroutine fail_with_code(code)
         integer, intent(in) :: code
         error_code = code
         call cleanup()
      end subroutine fail_with_code

      !> Deallocate temporary buffer
      subroutine cleanup()
         if (allocated(word_buffer)) deallocate (word_buffer)
      end subroutine cleanup

      !> Ensure buffer has capacity for given number of words
      subroutine ensure_buffer_size(needed)
         integer, intent(in) :: needed
         integer(int32), allocatable :: temp(:)
         integer :: new_size
         if (needed <= 0) return
         if (allocated(word_buffer)) then
            if (size(word_buffer) >= needed) return
         end if
         if (.not. allocated(word_buffer)) then
            new_size = needed
         else if (size(word_buffer) <= 0) then
            new_size = needed
         else
            new_size = max(needed, size(word_buffer)*2)
         end if
         allocate (temp(new_size))
         call move_alloc(temp, word_buffer)
      end subroutine ensure_buffer_size

      !---------------------------------------------------------------------------
      !> Fill matrix column with dense (contiguous) data
      !>
      !> @param[in]    nrows      Total matrix rows
      !> @param[in]    ntype      Data type (1=single, 2=double)
      !> @param[in]    icol       Column index
      !> @param[in]    irow0      Starting row index
      !> @param[in]    w          Raw word data
      !> @param[inout] a          Target matrix
      !---------------------------------------------------------------------------
      subroutine fill_dense_column(nrows, ntype, icol, irow0, w, a)
         integer, intent(in) :: nrows
         integer(int32), intent(in) :: ntype
         integer(int32), intent(in) :: icol, irow0
         integer(int32), intent(in) :: w(:)
         real(real64), intent(inout) :: a(:, :)
         integer :: val_idx, row, num_values
         if (ntype == 1) then
            do val_idx = 1, size(w)
               row = int(irow0) + val_idx - 1
               if (row >= 1 .and. row <= nrows) then
                  a(row, int(icol)) = real(word_to_real32(w(val_idx)), real64)
               end if
            end do
         else
            num_values = size(w)/2
            do val_idx = 1, num_values
               row = int(irow0) + val_idx - 1
               if (row >= 1 .and. row <= nrows) then
                  a(row, int(icol)) = words_to_real64(w(2*val_idx - 1), w(2*val_idx))
               end if
            end do
         end if
      end subroutine fill_dense_column

      !---------------------------------------------------------------------------
      !> Fill matrix column with sparse (string-encoded) data
      !>
      !> Sparse format uses "strings" where each starts with a header word:
      !> bits 31-16 = length L, bits 15-0 = starting row. Data follows header.
      !>
      !> @param[in]    nrows      Total matrix rows
      !> @param[in]    ntype      Data type (1=single, 2=double)
      !> @param[in]    icol       Column index
      !> @param[in]    w          Raw word data with string headers
      !> @param[inout] a          Target matrix
      !---------------------------------------------------------------------------
      subroutine fill_sparse_column(nrows, ntype, icol, w, a)
         integer, intent(in) :: nrows
         integer(int32), intent(in) :: ntype
         integer(int32), intent(in) :: icol
         integer(int32), intent(in) :: w(:)
         real(real64), intent(inout) :: a(:, :)
         integer :: string_pos, val_idx, row
         integer(int32) :: string_len, str_start_row, data_words
         string_pos = 1
         do
            if (string_pos > size(w)) exit
            ! Decode string header: length in upper 16 bits, start row in lower 16
            string_len = iand(ishft(w(string_pos), -16), int(z'FFFF', int32))
            str_start_row = iand(w(string_pos), int(z'FFFF', int32))
            data_words = string_len - 1
            if (ntype == 1) then
               do val_idx = 1, int(data_words)
                  row = int(str_start_row) + val_idx - 1
                  if (row >= 1 .and. row <= nrows) then
                     a(row, int(icol)) = real(word_to_real32(w(string_pos + val_idx)), real64)
                  end if
               end do
            else
               do val_idx = 1, int(data_words), 2
                  row = int(str_start_row) + (val_idx - 1)/2
                  if (row >= 1 .and. row <= nrows) then
                     a(row, int(icol)) = words_to_real64(w(string_pos + val_idx), w(string_pos + val_idx + 1))
                  end if
               end do
            end if
            string_pos = string_pos + int(string_len)
            if (string_pos >= size(w)) exit
         end do
      end subroutine fill_sparse_column
   end subroutine output4_read_matrix

   !===========================================================================
   ! OP2 FILE PARSER - PUBLIC INTERFACE
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Parse a Nastran OP2 file and extract geometry and modal data
   !>
   !> This is the main public interface for OP2 file parsing. It reads:
   !>   - GEOM1 table: Grid point coordinates (requires PARAM,OGEOM,YES in Nastran)
   !>   - OUG1 tables: Eigenvectors (mode shapes) and optional static displacements
   !>
   !> For SOL 103 (normal modes): Only eigenvectors are extracted
   !> For SOL 106 (nonlinear static + modal): Static displacements are applied
   !>   to coordinates to get deformed geometry, then eigenvectors are extracted
   !>
   !> OP2 File Structure:
   !>   - Header: Date, timestamp, file label
   !>   - Tables: Each table has header (name, timestamp) followed by data records
   !>   - Table names are 8 characters (2 x 4-byte words)
   !>   - Data records contain variable-length blocks with continuation markers
   !>
   !> @param[in]    op2_path      Path to OP2 file
   !> @param[in]    n_modes       Expected number of modes (from modal mass matrix)
   !> @param[out]   n_nodes       Number of grid points found
   !> @param[out]   node_ids      External node IDs [n_nodes]
   !> @param[out]   xg, yg, zg    Grid coordinates [n_nodes] (deformed if SOL 106)
   !> @param[out]   phix          Mode shapes [6*n_nodes, n_modes]
   !> @param[in]    log_unit      Unit for diagnostic output (default: 6)
   !> @param[out]   status        0=success, non-zero=error code
   !>
   !> Error codes:
   !>   1 = File open error
   !>   2 = OP2 header parse error
   !>   3 = GEOM1 table not found (check PARAM,OGEOM,YES)
   !>   4 = OUG1 table not found (eigenvectors missing)
   !>   5 = I/O error during parsing
   !---------------------------------------------------------------------------
   subroutine op2_parse_file(op2_path, n_modes, n_nodes, node_ids, xg, yg, zg, phix, log_unit, status)
      character(len=*), intent(in) :: op2_path
      integer, intent(in) :: n_modes
      integer, intent(out) :: n_nodes
      integer(int32), allocatable, intent(out) :: node_ids(:)
      real(real64), allocatable, intent(out) :: xg(:), yg(:), zg(:)
      real(real64), allocatable, intent(out) :: phix(:, :)
      integer, intent(in), optional :: log_unit
      integer, intent(out) :: status

      type(op2_reader) :: reader
      integer(int32), allocatable :: word_buffer(:)
      integer(int32) :: file_label(2), table_name(2)
      real :: timestamp(7)
      real(real64), allocatable :: static_disp_x(:), static_disp_y(:), static_disp_z(:)

      integer :: op2_unit, io_ret, log_u, read_status, ret_code, words_read
      integer :: approach_code, oug1_count
      logical :: op2_open, has_geom1, has_eigenvectors, has_static_disp
      logical :: done_scanning, is_eof

      ! Initialize
      status = 0
      op2_unit = 20
      log_u = 6
      if (present(log_unit)) log_u = log_unit
      op2_open = .false.
      has_geom1 = .false.
      has_eigenvectors = .false.
      has_static_disp = .false.
      n_nodes = 0

      ! Open OP2 file
      open (op2_unit, file=op2_path, status='old', form='unformatted', iostat=io_ret)
      if (io_ret /= 0) then
         write (log_u, '(A,A)') 'OP2 ERROR: Cannot open file ', trim(op2_path)
         status = 1
         return
      end if
      op2_open = .true.

      ! Allocate work buffer
      allocate (word_buffer(12*MAX_NODES))
      word_buffer = 0_int32

      ! Initialize OP2 reader
      call reader%init(op2_unit, log_u, file_label, read_status)
      if (read_status /= 0) then
         write (log_u, '(A)') 'OP2 ERROR: Invalid OP2 header format'
         status = 2
         call cleanup_op2_parse()
         return
      end if

      ! Read first table header
      call reader%read_block_header(table_name, timestamp, read_status)
      if (read_status /= 0) then
         write (log_u, '(A)') 'OP2 ERROR: Cannot read first table header'
         status = 2
         call cleanup_op2_parse()
         return
      end if

      !=========================================================================
      ! Phase 1: Scan for GEOM1 table (grid coordinates)
      !=========================================================================
      do while (table_name(1) /= TNAME_OUG1_W1)
         if (table_name(1) == TNAME_GEOM1_W1 .and. table_name(2) == TNAME_GEOM1_W2) then
            ! Found GEOM1 - read grid coordinates
            write (log_u, '(A)') ' OP2: Reading GEOM1 (grid coordinates)'
            call read_geom1_grids(reader, word_buffer, n_nodes, node_ids, xg, yg, zg, read_status)
            if (read_status == 0 .and. n_nodes > 0) then
               has_geom1 = .true.
            end if
         else
            ! Skip non-GEOM1 table
            call op2_skip_table(reader, word_buffer, ret_code)
            if (ret_code == -99) then
               status = 5
               call cleanup_op2_parse()
               return
            end if
         end if

         ! Read next table header
         call reader%read_block_header(table_name, timestamp, read_status)
         if (read_status /= 0) then
            write (log_u, '(A)') 'OP2 ERROR: Unexpected end while scanning for OUG1'
            status = 4
            call cleanup_op2_parse()
            return
         end if
      end do

      ! Verify GEOM1 was found
      if (.not. has_geom1) then
         write (log_u, '(A)') 'OP2 ERROR: GEOM1 table not found'
         write (log_u, '(A)') '           Ensure PARAM,OGEOM,YES is set in Nastran input'
         status = 3
         call cleanup_op2_parse()
         return
      end if

      !=========================================================================
      ! Phase 2: Process OUG1 tables (eigenvectors and static displacements)
      !=========================================================================

      ! Allocate mode shape matrix
      allocate (phix(6*n_nodes, n_modes))
      phix = 0.0_real64

      oug1_count = 0
      done_scanning = .false.
      is_eof = .false.

      do while (.not. done_scanning .and. .not. is_eof)
         if (table_name(1) == TNAME_OUG1_W1) then
            oug1_count = oug1_count + 1

            ! Read IDENT record to determine OUG1 type
            call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
            if (ret_code == -99) then
               is_eof = .true.
               exit
            end if

            if (ret_code /= 2) then
               approach_code = word_buffer(1)

               select case (approach_code)
               case (101, 102)  ! Static displacement
                  write (log_u, '(A,I0)') ' OP2: Reading static displacement OUG1 #', oug1_count

                  ! Re-read static displacement data
                  if (allocated(static_disp_x)) deallocate (static_disp_x)
                  if (allocated(static_disp_y)) deallocate (static_disp_y)
                  if (allocated(static_disp_z)) deallocate (static_disp_z)

                  call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
                  if (ret_code == -99) then
                     is_eof = .true.
                     exit
                  end if

                  call parse_oug1_static_data(word_buffer, words_read, node_ids, &
                                              static_disp_x, static_disp_y, static_disp_z)
                  has_static_disp = .true.

               case (22)  ! Eigenvectors
                  call parse_oug1_eigenvectors(reader, word_buffer, n_nodes, n_modes, &
                                               node_ids, phix, ret_code)
                  if (ret_code == -99) then
                     is_eof = .true.
                     exit
                  end if
                  has_eigenvectors = .true.

               case default
                  ! Unknown OUG1 type - skip it
                  write (log_u, '(A,I0,A,I0)') ' OP2: Skipping OUG1 #', oug1_count, ' (approach_code=', approach_code, ')'
               end select
            end if

            ! Skip to end of this OUG1 table
            do while (ret_code /= 2 .and. ret_code /= -99)
               call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
            end do
            if (ret_code == -99) then
               is_eof = .true.
               exit
            end if

         else
            ! Not OUG1 - skip this table
            call op2_skip_table(reader, word_buffer, ret_code)
            if (ret_code == -99) then
               is_eof = .true.
               exit
            end if
         end if

         ! Read next table header
         call reader%read_block_header(table_name, timestamp, read_status)
         if (read_status /= 0) then
            done_scanning = .true.
         end if
      end do

      !=========================================================================
      ! Phase 3: Apply static displacements if present (SOL 106)
      !=========================================================================
      if (has_static_disp) then
         xg = xg + static_disp_x
         yg = yg + static_disp_y
         zg = zg + static_disp_z
         write (log_u, '(A)') ' OP2: Applied static displacements to get deformed coordinates'
         deallocate (static_disp_x, static_disp_y, static_disp_z)
      end if

      ! Verify eigenvectors were found
      if (.not. has_eigenvectors) then
         write (log_u, '(A)') 'OP2 ERROR: No eigenvector data found in OUG1 tables'
         status = 4
         call cleanup_op2_parse()
         return
      end if

      write (log_u, '(A,I0,A,I0,A)') ' OP2: Successfully parsed ', n_nodes, ' nodes, ', n_modes, ' modes'
      call cleanup_op2_parse()

   contains

      subroutine cleanup_op2_parse()
         if (allocated(word_buffer)) deallocate (word_buffer)
         if (allocated(static_disp_x)) deallocate (static_disp_x)
         if (allocated(static_disp_y)) deallocate (static_disp_y)
         if (allocated(static_disp_z)) deallocate (static_disp_z)
         if (op2_open) close (op2_unit)
      end subroutine cleanup_op2_parse

   end subroutine op2_parse_file

   ! ---------------------------------------------------------------------------
   ! OP2 Reader Logic
   ! ---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   !> Initialize OP2 reader and parse file header
   !>
   !> Validates the OP2 file format by reading the standard header sequence
   !> (date, timestamp, label) and prepares for table reading.
   !> OP2 is the standard Nastran binary output format.
   !>
   !> @param[in]  file_unit   Fortran unit number (file must be open)
   !> @param[in]  log_unit    Optional unit for diagnostic output (default: 6)
   !> @param[out] label       2-word file label from header
   !> @param[out] status      0=success, non-zero=error
   !---------------------------------------------------------------------------
   subroutine op2_init(this, file_unit, log_unit, label, status)
      class(op2_reader), intent(inout) :: this
      integer, intent(in) :: file_unit
      integer, intent(in), optional :: log_unit
      integer(int32), intent(out) :: label(2)
      integer, intent(out) :: status

      integer(int32) :: key, month, day, year, timestamp(7)
      integer :: idx, io_status
      logical :: has_full_header

      status = 0
      this%file_unit = file_unit
      if (present(log_unit)) this%log_unit = log_unit
      if (.not. allocated(this%buffer)) allocate (this%buffer(200000))
      rewind (this%file_unit)

      ! 1. Header Key
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return

      ! Check if file has full header (key=3) or simplified header (key=2)
      has_full_header = (key == 3_int32)

      if (.not. has_full_header) then
         ! Simplified format: rewind and skip full header parsing
         rewind (this%file_unit)
         label(1) = 0_int32
         label(2) = 0_int32
         this%record_active = .false.
         this%buffer_valid = .false.
         this%buffer_pos = 0
         this%buffer_len = 0
         return
      end if

      ! 2. Date
      read (this%file_unit, iostat=io_status) month, day, year
      call check_read_error(io_status)
      if (status /= 0) return

      ! 3. Info Key
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 7_int32)
      if (status /= 0) return

      ! 4. Timestamp
      read (this%file_unit, iostat=io_status) (timestamp(idx), idx=1, 7)
      call check_read_error(io_status)
      if (status /= 0) return

      ! 5. Key=2
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 2_int32)
      if (status /= 0) return

      ! 6. Label
      read (this%file_unit, iostat=io_status) label
      call check_read_error(io_status)
      if (status /= 0) return

      ! 7. Key=-1
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, -1_int32)
      if (status /= 0) return

      ! 8. Key=0
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 0_int32)
      if (status /= 0) return

      this%record_active = .false.
      this%buffer_valid = .false.
      this%buffer_pos = 0
      this%buffer_len = 0

   contains
      subroutine check_read_error(io_stat)
         integer, intent(in) :: io_stat
         if (io_stat /= 0) then
            call fail_op2(this, 'IOPEN READ ERROR', 0_int32, status)
         end if
      end subroutine check_read_error

      subroutine check_key_value(val, expected)
         integer(int32), intent(in) :: val, expected
         if (val /= expected) then
            call fail_op2(this, 'IOPEN BAD KEY', val, status)
         end if
      end subroutine check_key_value

   end subroutine op2_init

   !---------------------------------------------------------------------------
   !> Read OP2 table block header
   !>
   !> Parses the header sequence that precedes each data table, extracting
   !> the table name and timestamp information.
   !>
   !> @param[out] table_name  2-word table identifier (e.g., 'GPL ', 'BGPDT')
   !> @param[out] timestamp   7-element timestamp array
   !> @param[out] status      0=success, non-zero=error
   !---------------------------------------------------------------------------
   subroutine op2_read_block_header(this, table_name, timestamp, status)
      class(op2_reader), intent(inout) :: this
      integer(int32), intent(out)  :: table_name(2)
      real, intent(out)  :: timestamp(7)
      integer, intent(out)  :: status

      integer(int32) :: key, header_info(2)
      integer :: io_status

      status = 0

      ! 1. Start Marker
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 2_int32)
      if (status /= 0) return

      ! 2. Table Name
      read (this%file_unit, iostat=io_status) table_name
      call check_read_error(io_status)
      if (status /= 0) return

      ! 3. Separator
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, -1_int32)
      if (status /= 0) return

      ! 4. Info Marker
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 7_int32)
      if (status /= 0) return

      ! 5. Timestamp
      read (this%file_unit, iostat=io_status) timestamp
      call check_read_error(io_status)
      if (status /= 0) return

      ! 6. Separator (some formats use -2, others use 3 - try to detect)
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      if (key /= -2_int32 .and. key /= 3_int32) then
         call fail_op2(this, 'IHEADR BAD KEY after timestamp', key, status)
         return
      end if

      ! 7. Block Index - path depends on previous key
      if (key == 3_int32) then
         ! When key=3, next comes -2, then 1, then 0
         read (this%file_unit, iostat=io_status) key
         call check_read_error(io_status)
         if (status /= 0) return
         call check_key_value(key, -2_int32)
         if (status /= 0) return
      end if

      ! Now read block index (happens in both paths)
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 1_int32)
      if (status /= 0) return

      ! 8. Separator
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, 0_int32)
      if (status /= 0) return

      ! 9. Length Check (this will be the next key regardless of path taken)
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      if (key < 2_int32) then
         call fail_op2(this, 'IHEADR BAD KEY', key, status)
         return
      end if

      ! 10. Header Data
      read (this%file_unit, iostat=io_status) header_info
      call check_read_error(io_status)
      if (status /= 0) return

      ! 11. End Marker
      read (this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, -3_int32)
      if (status /= 0) return

      this%record_active = .false.
      this%buffer_valid = .false.
      this%buffer_pos = 0
      this%buffer_len = 0

   contains
      subroutine check_read_error(io_stat)
         integer, intent(in) :: io_stat
         if (io_stat /= 0) then
            call fail_op2(this, 'IHEADR READ ERROR', 0_int32, status)
         end if
      end subroutine check_read_error

      subroutine check_key_value(val, expected)
         integer(int32), intent(in) :: val, expected
         if (val /= expected) then
            call fail_op2(this, 'IHEADR BAD KEY', val, status)
         end if
      end subroutine check_key_value

   end subroutine op2_read_block_header

   !---------------------------------------------------------------------------
   !> Read words from OP2 data record
   !>
   !> Reads requested number of 32-bit words from the current data record,
   !> handling multi-block records transparently. Supports skipping and
   !> partial reads.
   !>
   !> @param[out]    dest_array     Destination array for read words
   !> @param[in]     words_requested Number of words to read (negative=skip)
   !> @param[in]     finish_record  If non-zero, skip to end of record after read
   !> @param[out]    words_read     Actual number of words transferred
   !> @param[out]    return_code    0=more data, 1=end of record, 2=end of table, -99=error
   !---------------------------------------------------------------------------
   subroutine op2_read_words(this, dest_array, words_requested, finish_record, words_read, return_code)
      class(op2_reader), intent(inout) :: this
      integer(int32), intent(inout) :: dest_array(:)
      integer, intent(in) :: words_requested
      integer, intent(in) :: finish_record
      integer, intent(out) :: words_read, return_code

      integer :: words_remaining, dest_pos, available, idx, io_status
      integer(int32) :: key

      words_remaining = words_requested
      dest_pos = 0
      return_code = 0
      words_read = 0

      do
         ! Start new record if needed
         if (.not. this%record_active) then
            read (this%file_unit, iostat=io_status) key
            call check_io(io_status, 'IREAD READ ERROR')
            if (return_code /= 0) return

            ! Check for end-of-table marker (key=0) - indicates all tables consumed
            if (key == 0_int32) then
               words_read = 0
               return_code = 2  ! End of table
               return
            end if

            call check_key(key, 1_int32)
            if (return_code /= 0) return

            read (this%file_unit, iostat=io_status) key
            call check_io(io_status, 'IREAD READ ERROR')
            if (return_code /= 0) return

            call check_key(key, 0_int32)
            if (return_code /= 0) return

            this%record_active = .true.
         end if

         ! Refill buffer if empty
         if (.not. this%buffer_valid) then
            read (this%file_unit, iostat=io_status) key
            call check_io_recover(io_status, 'IREAD READ ERROR')
            if (return_code /= 0) return

            if (key < 0_int32) then
               words_read = dest_pos
               return_code = 1
               this%record_active = .false.
               return
            end if
            if (key == 0_int32) then
               words_read = dest_pos
               return_code = 2
               this%record_active = .false.
               return
            end if
            if (key > int(OP2_MAX_BLOCK_WORDS, int32)) then
               write (this%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', key
               return_code = -99
               this%record_active = .false.
               return
            end if
            call ensure_buffer_capacity(this, int(key))
            read (this%file_unit, iostat=io_status) this%buffer(1:int(key))
            call check_io_recover(io_status, 'IREAD READ ERROR')
            if (return_code /= 0) return

            this%buffer_len = int(key)
            this%buffer_pos = 0
            this%buffer_valid = .true.
         end if

         ! Process request
         if (words_remaining < 0) then
            call skip_within_record(this, -words_remaining)
            words_remaining = 0
         else if (words_remaining == 0) then
            exit
         else
            available = this%buffer_len - this%buffer_pos
            if (available >= words_remaining) then
               do idx = 1, words_remaining
                  dest_array(dest_pos + idx) = this%buffer(this%buffer_pos + idx)
               end do
               dest_pos = dest_pos + words_remaining
               this%buffer_pos = this%buffer_pos + words_remaining
               words_remaining = 0
            else
               do idx = 1, available
                  dest_array(dest_pos + idx) = this%buffer(this%buffer_pos + idx)
               end do
               dest_pos = dest_pos + available
               words_remaining = words_remaining - available
               this%buffer_valid = .false.
               this%buffer_pos = 0
               this%buffer_len = 0
               cycle
            end if
         end if
         if (words_remaining == 0) exit
      end do

      if (finish_record /= 0) then
         call skip_to_end_of_record(this)
         words_read = dest_pos
         return_code = 1
         this%record_active = .false.
         return
      end if

      words_read = dest_pos
      if (this%buffer_valid) then
         return_code = 0
      else
         return_code = 1
         this%record_active = .false.
      end if

   contains
      subroutine check_io(io_stat, msg)
         integer, intent(in) :: io_stat
         character(len=*), intent(in) :: msg
         if (io_stat /= 0) then
            if (.not. this%silent) write (this%log_unit, '(A)') msg
            return_code = -99
         end if
      end subroutine check_io

      subroutine check_io_recover(io_stat, msg)
         integer, intent(in) :: io_stat
         character(len=*), intent(in) :: msg
         if (io_stat /= 0) then
            if (.not. this%silent) write (this%log_unit, '(A)') msg
            return_code = -99
            this%record_active = .false.
         end if
      end subroutine check_io_recover

      subroutine check_key(val, expected)
         integer(int32), intent(in) :: val, expected
         if (val /= expected) then
            if (.not. this%silent) write (this%log_unit, '(A,1X,I0)') 'IREAD BAD KEY:', val
            return_code = -99
         end if
      end subroutine check_key

      subroutine skip_within_record(rdr, num_skip)
         class(op2_reader), intent(inout) :: rdr
         integer, intent(in) :: num_skip
         integer :: remaining, io_stat
         integer(int32) :: blk_key
         remaining = num_skip
         do while (remaining > 0)
            if (.not. rdr%buffer_valid) then
               read (rdr%file_unit, iostat=io_stat) blk_key
               call check_io_recover(io_stat, 'IREAD READ ERROR')
               if (return_code /= 0) then
                  rdr%buffer_valid = .false.
                  return
               end if

               if (blk_key <= 0_int32) then
                  rdr%record_active = .false.
                  rdr%buffer_valid = .false.
                  return
               end if
               if (blk_key > int(OP2_MAX_BLOCK_WORDS, int32)) then
                  write (rdr%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', blk_key
                  return_code = -99
                  rdr%record_active = .false.
                  rdr%buffer_valid = .false.
                  return
               end if
               call ensure_buffer_capacity(rdr, int(blk_key))
               read (rdr%file_unit, iostat=io_stat) rdr%buffer(1:int(blk_key))
               call check_io_recover(io_stat, 'IREAD READ ERROR')
               if (return_code /= 0) then
                  rdr%buffer_valid = .false.
                  return
               end if
               rdr%buffer_len = int(blk_key)
               rdr%buffer_pos = 0
               rdr%buffer_valid = .true.
            end if
            if (rdr%buffer_len - rdr%buffer_pos >= remaining) then
               rdr%buffer_pos = rdr%buffer_pos + remaining
               remaining = 0
            else
               remaining = remaining - (rdr%buffer_len - rdr%buffer_pos)
               rdr%buffer_valid = .false.
               rdr%buffer_pos = 0
               rdr%buffer_len = 0
            end if
         end do
      end subroutine skip_within_record

      subroutine skip_to_end_of_record(rdr)
         class(op2_reader), intent(inout) :: rdr
         integer(int32) :: blk_key
         integer :: io_stat
         rdr%buffer_valid = .false.
         rdr%buffer_pos = 0
         rdr%buffer_len = 0
         do
            read (rdr%file_unit, iostat=io_stat) blk_key
            call check_io_recover(io_stat, 'IREAD READ ERROR')
            if (return_code /= 0) return

            if (blk_key <= 0_int32) exit
            if (blk_key > int(OP2_MAX_BLOCK_WORDS, int32)) then
               write (rdr%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', blk_key
               return_code = -99
               rdr%record_active = .false.
               return
            end if
            call ensure_buffer_capacity(rdr, int(blk_key))
            read (rdr%file_unit, iostat=io_stat) rdr%buffer(1:int(blk_key))
            call check_io_recover(io_stat, 'IREAD READ ERROR')
            if (return_code /= 0) return
         end do
      end subroutine skip_to_end_of_record
   end subroutine op2_read_words

   !---------------------------------------------------------------------------
   !> Parse static displacement data from OUG1 DATA record
   !>
   !> OUG1 DATA format: [node_id*10+device(i4), type_code(i4), t1(f32), t2(f32),
   !> t3(f32), r1(f32), r2(f32), r3(f32)] per node (8 words = 32 bytes)
   !>
   !> @param[in]  word_buffer  Raw data from OUG1 DATA record
   !> @param[in]  words_read   Number of words in buffer
   !> @param[in]  node_ids     Reference node ID array for ordering
   !> @param[out] disp_x,y,z   Displacement arrays (allocated here)
   !> @param[in]  log_unit     Unit for diagnostic output
   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   !> Read eigenvector mode shapes from OUG1 table (approach_code=22)
   !>
   !> Wrapper for parse_oug1_eigenvectors that works with femgen_model type.
   !> Allocates model%phix if needed and delegates parsing to the shared function.
   !>
   !> @param[inout] reader      OP2 reader (positioned after first IDENT)
   !> @param[inout] word_buffer Work buffer for reading
   !> @param[inout] model       Model structure (phix will be allocated/filled)
   !> @param[out]   ret_code    Return code (0=success, -99=EOF)
   !---------------------------------------------------------------------------
   subroutine read_oug1_eigenvectors(reader, word_buffer, model, ret_code)
      type(op2_reader), intent(inout) :: reader
      integer(int32), intent(inout) :: word_buffer(:)
      type(femgen_model), intent(inout) :: model
      integer, intent(out) :: ret_code

      ! Allocate phix if not already allocated
      if (.not. allocated(model%phix)) then
         allocate (model%phix(6*model%n_nodes, model%n_modes))
         model%phix = 0.0_real64
         model%phix_nrow = 6*model%n_nodes
         model%phix_ncol = model%n_modes
      end if

      ! Delegate to shared parsing function
      call parse_oug1_eigenvectors(reader, word_buffer, model%n_nodes, model%n_modes, &
                                   model%node_ids, model%phix, ret_code)
   end subroutine read_oug1_eigenvectors

   subroutine parse_oug1_static_data(word_buffer, words_read, node_ids, disp_x, disp_y, disp_z)
      integer(int32), intent(in) :: word_buffer(:)
      integer, intent(in) :: words_read
      integer(int32), intent(in) :: node_ids(:)
      real(real64), allocatable, intent(out) :: disp_x(:), disp_y(:), disp_z(:)

      integer :: n_nodes, i, ext_id, node_idx, words_per_node, iword
      real(real32) :: t1_32, t2_32, t3_32

      n_nodes = size(node_ids)
      allocate (disp_x(n_nodes), disp_y(n_nodes), disp_z(n_nodes))
      disp_x = 0.0_real64
      disp_y = 0.0_real64
      disp_z = 0.0_real64

      ! Determine format: 8 words (real32) or 14 words (real64) per node
      if (mod(words_read, 8) == 0) then
         words_per_node = 8
      else if (mod(words_read, 14) == 0) then
         words_per_node = 14
      else
         words_per_node = 8  ! Default guess
      end if

      ! Parse displacement data
      iword = 1
      do while (iword + words_per_node - 1 <= words_read)
         ! Extract node ID (encoded as ext_id * 10 + device_code)
         ext_id = word_buffer(iword)/10

         ! Find matching node index
         node_idx = 0
         do i = 1, n_nodes
            if (node_ids(i) == ext_id) then
               node_idx = i
               exit
            end if
         end do

         if (node_idx > 0) then
            if (words_per_node == 8) then
               ! Real32 format: [id, type, t1, t2, t3, r1, r2, r3]
               t1_32 = word_to_real32(word_buffer(iword + 2))
               t2_32 = word_to_real32(word_buffer(iword + 3))
               t3_32 = word_to_real32(word_buffer(iword + 4))
               disp_x(node_idx) = real(t1_32, real64)
               disp_y(node_idx) = real(t2_32, real64)
               disp_z(node_idx) = real(t3_32, real64)
            else
               ! Real64 format: [id, type, t1(2), t2(2), t3(2), r1(2), r2(2), r3(2)]
               disp_x(node_idx) = words_to_real64(word_buffer(iword + 2), word_buffer(iword + 3))
               disp_y(node_idx) = words_to_real64(word_buffer(iword + 4), word_buffer(iword + 5))
               disp_z(node_idx) = words_to_real64(word_buffer(iword + 6), word_buffer(iword + 7))
            end if
         end if

         iword = iword + words_per_node
      end do
   end subroutine parse_oug1_static_data

   !---------------------------------------------------------------------------
   !> Read GEOM1 table to extract original grid coordinates
   !>
   !> GEOM1 contains GRID card data with record type 4501.
   !> Each GRID: [id(i4), cp(i4), x(f32), y(f32), z(f32), cd(i4), ps(i4), seid(i4)]
   !> Total 8 words per grid point (32 bytes)
   !>
   !> @param[inout] reader      OUTPUT2 reader object
   !> @param[in]    word_buffer Work buffer for reading
   !> @param[out]   n_nodes     Number of grid points found
   !> @param[out]   node_ids    Array of node IDs (allocated here)
   !> @param[out]   xg,yg,zg    Coordinate arrays (allocated here)
   !> @param[out]   status      0=success, non-zero=error
   !---------------------------------------------------------------------------
   subroutine read_geom1_grids(reader, word_buffer, n_nodes, node_ids, xg, yg, zg, status)
      type(op2_reader), intent(inout) :: reader
      integer(int32), intent(inout) :: word_buffer(:)
      integer, intent(out) :: n_nodes
      integer(int32), allocatable, intent(out) :: node_ids(:)
      real(real64), allocatable, intent(out) :: xg(:), yg(:), zg(:)
      integer, intent(out) :: status

      integer :: words_read, ret_code, record_type, i, grid_count
      integer(int32), allocatable :: temp_ids(:)
      real(real64), allocatable :: temp_x(:), temp_y(:), temp_z(:)
      integer :: max_grids, start_offset

      status = 0
      max_grids = size(word_buffer)/8  ! 8 words per GRID
      allocate (temp_ids(max_grids), temp_x(max_grids), temp_y(max_grids), temp_z(max_grids))
      grid_count = 0

      ! GEOM1 structure:
      !   IDENT record
      !   Data records with format: [record_type(i4), extra(i4), count(i4), data...]
      !   For GRID (type 4501): data is [nid, cp, x, y, z, cd, ps, seid] * count

      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)  ! IDENT

      ! Read data records looking for record type 4501 (GRID)
      do
         call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
         if (ret_code == 2) exit  ! End of table
         if (ret_code == -99) then
            status = 1
            return
         end if

         ! First word is record type (4501 for GRID)
         record_type = word_buffer(1)

         if (record_type == 4501) then
            ! GRID record format: [type=4501, subtype, count, then GRID data...]
            ! GRID data starts at word 4 (1-indexed): [nid, cp, x, y, z, cd, ps, seid]
            ! Each grid is 8 words
            start_offset = 4  ! Skip type, subtype, count
            do i = start_offset, words_read - 8 + 1, 8
               if (word_buffer(i) == -1) exit  ! End marker
               if (word_buffer(i) <= 0) cycle  ! Skip invalid entries
               grid_count = grid_count + 1
               if (grid_count > max_grids) exit

               temp_ids(grid_count) = word_buffer(i)      ! Node ID
               temp_x(grid_count) = real(word_to_real32(word_buffer(i + 2)), real64)  ! X
               temp_y(grid_count) = real(word_to_real32(word_buffer(i + 3)), real64)  ! Y
               temp_z(grid_count) = real(word_to_real32(word_buffer(i + 4)), real64)  ! Z
            end do
         end if
      end do

      ! Copy to output arrays
      n_nodes = grid_count
      if (grid_count > 0) then
         allocate (node_ids(n_nodes), xg(n_nodes), yg(n_nodes), zg(n_nodes))
         node_ids = temp_ids(1:n_nodes)
         xg = temp_x(1:n_nodes)
         yg = temp_y(1:n_nodes)
         zg = temp_z(1:n_nodes)
      end if

      deallocate (temp_ids, temp_x, temp_y, temp_z)
   end subroutine read_geom1_grids

   !---------------------------------------------------------------------------
   !> Parse eigenvector mode shapes from OUG1 data into array
   !>
   !> Core parsing logic shared by both femgen_read (model-based) and
   !> op2_parse_file (array-based) interfaces.
   !>
   !> @param[inout] reader      OP2 reader (positioned after first IDENT)
   !> @param[inout] word_buffer Work buffer for reading
   !> @param[in]    n_nodes     Number of nodes
   !> @param[in]    n_modes     Number of modes to read
   !> @param[in]    node_ids    Node ID array for matching
   !> @param[inout] phix        Mode shape matrix [6*n_nodes, n_modes]
   !> @param[out]   ret_code    Return code (0=success, -99=EOF)
   !---------------------------------------------------------------------------
   subroutine parse_oug1_eigenvectors(reader, word_buffer, n_nodes, n_modes, node_ids, phix, ret_code)
      type(op2_reader), intent(inout) :: reader
      integer(int32), intent(inout) :: word_buffer(:)
      integer, intent(in) :: n_nodes, n_modes
      integer(int32), intent(in) :: node_ids(:)
      real(real64), intent(inout) :: phix(:, :)
      integer, intent(out) :: ret_code

      integer :: words_read, mode_idx, i
      integer :: words_per_node, iword, ext_id, node_idx
      real(real32) :: t1, t2, t3, r1, r2, r3

      mode_idx = 1

      do while (mode_idx <= n_modes)
         call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
         if (ret_code == -99) return
         if (ret_code == 2) exit

         ! Skip IDENT records (146 words)
         if (words_read == 146) then
            call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
            if (ret_code == -99) return
            if (ret_code == 2) exit
         end if

         ! Parse DATA record (8 words per node: id*10+dev, type, t1-t3, r1-r3)
         words_per_node = 8
         iword = 1
         do while (iword + words_per_node - 1 <= words_read)
            ext_id = word_buffer(iword)/10

            node_idx = 0
            do i = 1, n_nodes
               if (node_ids(i) == ext_id) then
                  node_idx = i
                  exit
               end if
            end do

            if (node_idx > 0) then
               t1 = word_to_real32(word_buffer(iword + 2))
               t2 = word_to_real32(word_buffer(iword + 3))
               t3 = word_to_real32(word_buffer(iword + 4))
               r1 = word_to_real32(word_buffer(iword + 5))
               r2 = word_to_real32(word_buffer(iword + 6))
               r3 = word_to_real32(word_buffer(iword + 7))

               phix(6*(node_idx - 1) + 1, mode_idx) = real(t1, real64)
               phix(6*(node_idx - 1) + 2, mode_idx) = real(t2, real64)
               phix(6*(node_idx - 1) + 3, mode_idx) = real(t3, real64)
               phix(6*(node_idx - 1) + 4, mode_idx) = real(r1, real64)
               phix(6*(node_idx - 1) + 5, mode_idx) = real(r2, real64)
               phix(6*(node_idx - 1) + 6, mode_idx) = real(r3, real64)
            end if

            iword = iword + words_per_node
         end do

         mode_idx = mode_idx + 1
      end do

      ret_code = 0
   end subroutine parse_oug1_eigenvectors

   !===========================================================================
   ! UTILITY FUNCTIONS
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Convert 32-bit integer word to single-precision real
   !>
   !> Reinterprets the bit pattern of a 32-bit integer as IEEE single-precision
   !> floating point, used for reading OUTPUT4 single-precision data.
   !>
   !> @param[in] word  32-bit integer containing IEEE float bits
   !> @return          Equivalent single-precision real value
   !---------------------------------------------------------------------------
   pure function word_to_real32(word) result(x)
      integer(int32), intent(in) :: word
      real(real32) :: x
      x = transfer(word, x)
   end function word_to_real32

   !---------------------------------------------------------------------------
   !> Convert two 32-bit integer words to double-precision real
   !>
   !> Reinterprets a pair of 32-bit integers as IEEE double-precision
   !> floating point, used for reading OUTPUT4 double-precision data.
   !>
   !> @param[in] word1  Low 32 bits of IEEE double
   !> @param[in] word2  High 32 bits of IEEE double
   !> @return           Equivalent double-precision real value
   !---------------------------------------------------------------------------
   pure function words_to_real64(word1, word2) result(x)
      integer(int32), intent(in) :: word1, word2
      real(real64) :: x
      integer(int32) :: w(2)
      w(1) = word1
      w(2) = word2
      x = transfer(w, x)
   end function words_to_real64

   !---------------------------------------------------------------------------
   !> Snap very small values to exactly zero
   !>
   !> Values with magnitude below FEMGEN_ZERO_TOL are set to zero,
   !> improving output readability and avoiding numerical noise.
   !>
   !> @param[in] x  Input value
   !> @return       Zero if |x| <= tolerance, otherwise x unchanged
   !---------------------------------------------------------------------------
   pure elemental function snap_zero(x) result(y)
      real(real64), intent(in) :: x
      real(real64) :: y
      if (abs(x) <= FEMGEN_ZERO_TOL) then
         y = 0.0_real64
      else
         y = x
      end if
   end function snap_zero

   !===========================================================================
   ! SHARED HELPER SUBROUTINES
   !===========================================================================

   !---------------------------------------------------------------------------
   !> Skip an entire OP2 table by reading until end-of-table marker
   !>
   !> @param[inout] reader      OP2 reader object
   !> @param[inout] word_buffer Work buffer for reading
   !> @param[out]   ret_code    Return code (2=end of table, -99=EOF/error)
   !---------------------------------------------------------------------------
   subroutine op2_skip_table(reader, word_buffer, ret_code)
      type(op2_reader), intent(inout) :: reader
      integer(int32), intent(inout) :: word_buffer(:)
      integer, intent(out) :: ret_code
      integer :: words_read
      reader%silent = .true.
      do
         call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
         if (ret_code == 2 .or. ret_code == -99) exit
      end do
      reader%silent = .false.
   end subroutine op2_skip_table

   !---------------------------------------------------------------------------
   !> Skip current OP2 record
   !>
   !> Advances file position past remaining data in current record.
   !---------------------------------------------------------------------------
   subroutine op2_skip_record(this)
      class(op2_reader), intent(inout) :: this
      integer(int32) :: dummy(1)
      integer :: num_read, ret_code
      ret_code = 0
      do while (ret_code == 0)
         call this%read_words(dummy, 0, 1, num_read, ret_code)
      end do
   end subroutine op2_skip_record

   !---------------------------------------------------------------------------
   !> Report OP2 read error
   !---------------------------------------------------------------------------
   subroutine fail_op2(this, msg, key, status)
      class(op2_reader), intent(in)  :: this
      character(len=*), intent(in)  :: msg
      integer(int32), intent(in)  :: key
      integer, intent(out) :: status

      if (.not. this%silent) then
         if (key /= 0_int32) then
            write (this%log_unit, '(A,1X,I0)') trim(msg)//':', key
         else
            write (this%log_unit, '(A)') trim(msg)
         end if
      end if
      status = 1
   end subroutine fail_op2

   !---------------------------------------------------------------------------
   !> Ensure OP2 reader buffer has capacity for given number of words
   !---------------------------------------------------------------------------
   subroutine ensure_buffer_capacity(this, needed)
      class(op2_reader), intent(inout) :: this
      integer, intent(in)    :: needed

      integer(int32), allocatable :: temp(:)
      integer :: current_size, new_size
      integer, parameter :: INITIAL_SIZE = 200000

      if (needed <= 0) return
      if (.not. allocated(this%buffer)) then
         allocate (this%buffer(max(needed, INITIAL_SIZE)))
         return
      end if

      current_size = size(this%buffer)
      if (current_size >= needed) return
      new_size = max(needed, current_size*2)
      allocate (temp(new_size))
      temp(1:current_size) = this%buffer
      call move_alloc(temp, this%buffer)
   end subroutine ensure_buffer_capacity

   !---------------------------------------------------------------------------
   !> Convert single character to lowercase
   !>
   !> @param[in] c  Single character (any case)
   !> @return       Lowercase equivalent (A-Z -> a-z), others unchanged
   !---------------------------------------------------------------------------
   pure function lower1(c) result(lc)
      character(len=1), intent(in) :: c
      character(len=1) :: lc
      integer :: ia
      ia = iachar(c)
      if (ia >= iachar('A') .and. ia <= iachar('Z')) then
         lc = achar(ia + (iachar('a') - iachar('A')))
      else
         lc = c
      end if
   end function lower1

   !---------------------------------------------------------------------------
   !> Check if filename ends with .fem extension (case-insensitive)
   !>
   !> @param[in] name  Filename string
   !> @return          True if name ends with .fem/.FEM/.Fem etc.
   !---------------------------------------------------------------------------
   pure function ends_with_fem(name) result(ok)
      character(len=*), intent(in) :: name
      logical :: ok
      integer :: n
      ok = .false.
      n = len_trim(name)
      if (n < 4) return
      ok = (lower1(name(n - 3:n - 3)) == '.') .and. &
           (lower1(name(n - 2:n - 2)) == 'f') .and. &
           (lower1(name(n - 1:n - 1)) == 'e') .and. &
           (lower1(name(n:n)) == 'm')
   end function ends_with_fem

end module femgen_mod
