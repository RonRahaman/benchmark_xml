program main

  use read_xml_primitives
  use write_xml_primitives
  use xmlparse
  use global
  use constants
  use error
  use dict_header,      only: DictIntInt, ElemKeyValueCI, dict_clear_ii
  use string,           only: lower_case, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with
  use list_header,      only: ListChar, ListReal
  use geometry_header,  only: Cell, Surface, Lattice

  implicit none

  type(DictIntInt) :: cells_in_univ_dict ! used to count how many cells each
                                         ! universe contains

  character(len=50) :: outfile = 'benchmark_xml.txt'
  character(len=2) :: path_input = './'
  double precision :: time
  integer :: fsize, istat

  open(unit=500, file='geometry.xml', status='old')
  inquire( 500, size=fsize)
  close(unit=500)

  print *, "***calling read_xml_file_geometry (varstring) with filesize = ", fsize
  call benchmark_xml_fortran_new(time)
  print *, 'read_xml_file_geometry (fox) completed in (seconds):', time


  open(unit=600, file=trim(outfile), status='unknown', position='append')
  write(600, '(a20, i20, f20.10)') 'v3', fsize, time
  close(unit=600)


  contains


    subroutine benchmark_xml_fortran_new(time)
      double precision, intent(out) :: time
      integer :: cstart, cend, crate
      call system_clock(cstart, crate)
      call read_geometry_xml_fortran_new()
      call system_clock(cend, crate)
      time = dble(cend-cstart)/crate 
      if (allocated(cells)) deallocate(cells)
      if (allocated(lattices)) deallocate(lattices)
      if (allocated(universes)) deallocate(universes)
      if (allocated(surfaces)) deallocate(surfaces)
      if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)
      call dict_clear_ii(cell_dict)
      call dict_clear_ii(universe_dict)
      call dict_clear_ii(lattice_dict)
      call dict_clear_ii(surface_dict)
      call dict_clear_ii(cells_in_univ_dict)

    end subroutine benchmark_xml_fortran_new

  subroutine read_geometry_xml_fortran_new()

    use xml_data_geometry_t

    integer :: i, j, k, m
    integer :: n
    integer :: n_x, n_y, n_z
    integer :: universe_num
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    integer :: mid
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    type(Cell),    pointer :: c => null()
    type(Surface), pointer :: s => null()
    type(Lattice), pointer :: lat => null()

    ! Display output message
    message = "Reading geometry XML file..."
    call write_message(5)

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      message = "Geometry XML file '" // trim(filename) // "' does not exist!"
      call fatal_error()
    end if

    ! Parse geometry.xml file
    call read_xml_file_geometry_t(filename)

    ! Get number of <cell> tags
    n_cells = size(cell_)

    ! Check for no cells
    if (n_cells == 0) then
      message = "No cells found in geometry.xml!"
      call fatal_error()
    end if

    ! Allocate cells array
    allocate(cells(n_cells))

    if (check_overlaps) then
      allocate(overlap_check_cnt(n_cells))
      overlap_check_cnt = 0
    end if

    n_universes = 0
    do i = 1, n_cells
      c => cells(i)

      ! Copy data into cells
      c % id       = cell_(i) % id
      c % universe = cell_(i) % universe
      c % fill     = cell_(i) % fill

      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has_key(c % id)) then
        message = "Two or more cells use the same unique ID: " // to_str(c % id)
        call fatal_error()
      end if

      ! Read material
      word = cell_(i) % material
      call lower_case(word)
      select case(word)
      case ('void')
        c % material = MATERIAL_VOID

      case ('')
        ! This case is called if no material was specified
        c % material = 0

      case default
        c % material = int(str_to_int(word), 4)

        ! Check for error
        if (c % material == ERROR_INT) then
          message = "Invalid material specified on cell " // to_str(c % id)
          call fatal_error()
        end if
      end select

      ! Check to make sure that either material or fill was specified
      if (c % material == NONE .and. c % fill == NONE) then
        message = "Neither material nor fill was specified for cell " // & 
             trim(to_str(c % id))
        call fatal_error()
      end if

      ! Check to make sure that both material and fill haven't been
      ! specified simultaneously
      if (c % material /= NONE .and. c % fill /= NONE) then
        message = "Cannot specify material and fill simultaneously"
        call fatal_error()
      end if

      ! Check to make sure that surfaces were specified
      if (.not. associated(cell_(i) % surfaces)) then
        message = "No surfaces specified for cell " // &
             trim(to_str(c % id))
        call fatal_error()
      end if

      ! Allocate array for surfaces and copy
      n = size(cell_(i) % surfaces)
      c % n_surfaces = n
      allocate(c % surfaces(n))
      c % surfaces = cell_(i) % surfaces

      ! Rotation matrix
      if (associated(cell_(i) % rotation)) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a rotation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of rotation parameters
        n = size(cell_(i) % rotation)
        if (n /= 3) then
          message = "Incorrect number of rotation parameters on cell " // &
               to_str(c % id)
          call fatal_error()
        end if

        ! Copy rotation angles in x,y,z directions
        phi   = -cell_(i) % rotation(1) * PI/180.0
        theta = -cell_(i) % rotation(2) * PI/180.0
        psi   = -cell_(i) % rotation(3) * PI/180.0

        ! Calculate rotation matrix based on angles given
        allocate(c % rotation(3,3))
        c % rotation = reshape((/ &
             cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta), &
             -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), &
             cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), &
             sin(phi)*cos(theta), &
             sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi), &
             -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi), &
             cos(phi)*cos(theta) /), (/ 3,3 /))
      end if

      ! Translation vector
      if (associated(cell_(i) % translation)) then
        ! Translations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a translation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of translation parameters
        n = size(cell_(i) % translation)
        if (n /= 3) then
          message = "Incorrect number of translation parameters on cell " &
               // to_str(c % id)
          call fatal_error()
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        c % translation = cell_(i) % translation
      end if

      ! Add cell to dictionary
      call cell_dict % add_key(c % id, i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      universe_num = cell_(i) % universe
      if (.not. cells_in_univ_dict % has_key(universe_num)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % add_key(universe_num, n_universes)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get_key(universe_num)
      end if
      call cells_in_univ_dict % add_key(universe_num, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! This variable is used to check whether at least one boundary condition was
    ! applied to a surface
    boundary_exists = .false.

    ! Get number of <surface> tags
    n_surfaces = size(surface_)

    ! Check for no surfaces
    if (n_surfaces == 0) then
      message = "No surfaces found in geometry.xml!"
      call fatal_error()
    end if

    ! Allocate cells array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      s => surfaces(i)

      ! Copy data into cells
      s % id = surface_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (surface_dict % has_key(s % id)) then
        message = "Two or more surfaces use the same unique ID: " // &
             to_str(s % id)
        call fatal_error()
      end if

      ! Copy and interpret surface type
      word = surface_(i) % type
      call lower_case(word)
      select case(trim(word))
      case ('x-plane')
        s % type = SURF_PX
        coeffs_reqd  = 1
      case ('y-plane')
        s % type = SURF_PY
        coeffs_reqd  = 1
      case ('z-plane')
        s % type = SURF_PZ
        coeffs_reqd  = 1
      case ('plane')
        s % type = SURF_PLANE
        coeffs_reqd  = 4
      case ('x-cylinder')
        s % type = SURF_CYL_X
        coeffs_reqd  = 3
      case ('y-cylinder')
        s % type = SURF_CYL_Y
        coeffs_reqd  = 3
      case ('z-cylinder')
        s % type = SURF_CYL_Z
        coeffs_reqd  = 3
      case ('sphere')
        s % type = SURF_SPHERE
        coeffs_reqd  = 4
      case ('x-cone')
        s % type = SURF_CONE_X
        coeffs_reqd  = 4
      case ('y-cone')
        s % type = SURF_CONE_Y
        coeffs_reqd  = 4
      case ('z-cone')
        s % type = SURF_CONE_Z
        coeffs_reqd  = 4
      case default
        message = "Invalid surface type: " // trim(surface_(i) % type)
        call fatal_error()
      end select

      ! Check to make sure that the proper number of coefficients
      ! have been specified for the given type of surface. Then copy
      ! surface coordinates.

      n = size(surface_(i) % coeffs)
      if (n < coeffs_reqd) then
        message = "Not enough coefficients specified for surface: " // & 
             trim(to_str(s % id))
        call fatal_error()
      elseif (n > coeffs_reqd) then
        message = "Too many coefficients specified for surface: " // &
             trim(to_str(s % id))
        call fatal_error()
      else
        allocate(s % coeffs(n))
        s % coeffs = surface_(i) % coeffs
      end if

      ! Boundary conditions
      word = surface_(i) % boundary
      call lower_case(word)
      select case (trim(word))
      case ('transmission', 'transmit', '')
        s % bc = BC_TRANSMIT
      case ('vacuum')
        s % bc = BC_VACUUM
        boundary_exists = .true.
      case ('reflective', 'reflect', 'reflecting')
        s % bc = BC_REFLECT
        boundary_exists = .true.
      case default
        message = "Unknown boundary condition '" // trim(word) // &
             "' specified on surface " // trim(to_str(s % id))
        call fatal_error()
      end select

      ! Add surface to dictionary
      call surface_dict % add_key(s % id, i)

    end do

    ! Check to make sure a boundary condition was applied to at least one
    ! surface
    if (.not. boundary_exists) then
      message = "No boundary conditions were applied to any surfaces!"
      call fatal_error()
    end if

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    ! Allocate lattices array
    n_lattices = size(lattice_)
    allocate(lattices(n_lattices))

    do i = 1, n_lattices
      lat => lattices(i)

      ! ID of lattice
      lat % id = lattice_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        message = "Two or more lattices use the same unique ID: " // &
             to_str(lat % id)
        call fatal_error()
      end if

      ! Read lattice type
      word = lattice_(i) % type
      call lower_case(word)
      select case (trim(word))
      case ('rect', 'rectangle', 'rectangular')
        lat % type = LATTICE_RECT
      case ('hex', 'hexagon', 'hexagonal')
        lat % type = LATTICE_HEX
      case default
        message = "Invalid lattice type: " // trim(lattice_(i) % type)
        call fatal_error()
      end select

      ! Read number of lattice cells in each dimension
      n = size(lattice_(i) % dimension)
      if (n /= 2 .and. n /= 3) then
        message = "Lattice must be two or three dimensions."
        call fatal_error()
      end if

      lat % n_dimension = n
      allocate(lat % dimension(n))
      lat % dimension = lattice_(i) % dimension

      ! Read lattice lower-left location
      if (size(lattice_(i) % dimension) /= size(lattice_(i) % lower_left)) then
        message = "Number of entries on <lower_left> must be the same as &
             &the number of entries on <dimension>."
        call fatal_error()
      end if

      allocate(lat % lower_left(n))
      lat % lower_left = lattice_(i) % lower_left

      ! Read lattice widths
      if (size(lattice_(i) % width) /= size(lattice_(i) % lower_left)) then
        message = "Number of entries on <width> must be the same as &
             &the number of entries on <lower_left>."
        call fatal_error()
      end if

      allocate(lat % width(n))
      lat % width = lattice_(i) % width

      ! Copy number of dimensions
      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if
      allocate(lat % universes(n_x, n_y, n_z))

      ! Check that number of universes matches size
      if (size(lattice_(i) % universes) /= n_x*n_y*n_z) then
        message = "Number of universes on <universes> does not match size of &
             &lattice " // trim(to_str(lat % id)) // "."
        call fatal_error()
      end if

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = lattice_(i) % &
                 universes(j + n_x*k + n_x*n_y*(m-1))
          end do
        end do
      end do

      ! Read material for area outside lattice
      mid = lattice_(i) % outside
      if (mid == 0 .or. mid == MATERIAL_VOID) then
        lat % outside = MATERIAL_VOID
      else
        lat % outside = mid
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, i)

    end do

  end subroutine read_geometry_xml_fortran_new


    subroutine write_message(level)

      integer, optional :: level ! verbosity level

      integer :: i_start    ! starting position
      integer :: i_end      ! ending position
      integer :: line_wrap  ! length of line
      integer :: length     ! length of message
      integer :: last_space ! index of last space (relative to start)

      ! Set length of line
      line_wrap = 80

      ! Only allow master to print to screen
      !if (.not. master .and. present(level)) return

      if (.not. present(level) .or. level <= verbosity) then
        ! Determine length of message
        length = len_trim(message)

        i_start = 0
        do
          if (length - i_start < line_wrap - 1) then
            ! Remainder of message will fit on line
            write(*, fmt='(1X,A)') message(i_start+1:length)
            exit

          else
            ! Determine last space in current line
            last_space = index(message(i_start+1:i_start+line_wrap), &
                 ' ', BACK=.true.)
            if (last_space == 0) then 
              i_end = min(length + 1, i_start+line_wrap) - 1
              write(*, fmt='(1X,A)') message(i_start+1:i_end)
            else
              i_end = i_start + last_space
              write(*, fmt='(1X,A)') message(i_start+1:i_end-1)
            end if

            ! Write up to last space

            ! Advance starting position
            i_start = i_end
            if (i_start > length) exit
          end if
        end do
      end if

    end subroutine write_message

    subroutine expand_natural_element(name, xs, density, list_names, &
         list_density)

      character(*),   intent(in)    :: name
      character(*),   intent(in)    :: xs
      real(8),        intent(in)    :: density
      type(ListChar), intent(inout) :: list_names
      type(ListReal), intent(inout) :: list_density

      character(2) :: element_name

      element_name = name(1:2)
      call lower_case(element_name)

      select case (element_name)
      case ('h')
        call list_names % append('1001.' // xs)
        call list_density % append(density * 0.999885_8)
        call list_names % append('1002.' // xs)
        call list_density % append(density * 0.000115_8)

      case ('he')
        call list_names % append('2003.' // xs)
        call list_density % append(density * 0.00000134_8)
        call list_names % append('2004.' // xs)
        call list_density % append(density * 0.99999866_8)

      case ('li')
        call list_names % append('3006.' // xs)
        call list_density % append(density * 0.0759_8)
        call list_names % append('3007.' // xs)
        call list_density % append(density * 0.9241_8)

      case ('be')
        call list_names % append('4009.' // xs)
        call list_density % append(density)

      case ('b')
        call list_names % append('5010.' // xs)
        call list_density % append(density * 0.199_8)
        call list_names % append('5011.' // xs)
        call list_density % append(density * 0.801_8)

      case ('c')
        ! The evaluation of Carbon in ENDF/B-VII.1 and JEFF 3.1.2 is a natural
        ! element, i.e. it's not possible to split into C-12 and C-13.
        call list_names % append('6000.' // xs)
        call list_density % append(density)

      case ('n')
        call list_names % append('7014.' // xs)
        call list_density % append(density * 0.99636_8)
        call list_names % append('7015.' // xs)
        call list_density % append(density * 0.00364_8)

      case ('o')
        ! O-18 does not exist in ENDF/B-VII.1 or JEFF 3.1.2 so its 0.205% has been
        ! added to O-16. The isotopic abundance for O-16 is ordinarily 99.757%.
        call list_names % append('8016.' // xs)
        call list_density % append(density * 0.99962_8)
        call list_names % append('8017.' // xs)
        call list_density % append(density * 0.00038_8)

      case ('f')
        call list_names % append('9019.' // xs)
        call list_density % append(density)

      case ('ne')
        call list_names % append('10020.' // xs)
        call list_density % append(density * 0.9048_8)
        call list_names % append('10021.' // xs)
        call list_density % append(density * 0.0027_8)
        call list_names % append('10022.' // xs)
        call list_density % append(density * 0.0925_8)

      case ('na')
        call list_names % append('11023.' // xs)
        call list_density % append(density)

      case ('mg')
        call list_names % append('12024.' // xs)
        call list_density % append(density * 0.7899_8)
        call list_names % append('12025.' // xs)
        call list_density % append(density * 0.1000_8)
        call list_names % append('12026.' // xs)
        call list_density % append(density * 0.1101_8)

      case ('al')
        call list_names % append('13027.' // xs)
        call list_density % append(density)

      case ('si')
        call list_names % append('14028.' // xs)
        call list_density % append(density * 0.92223_8)
        call list_names % append('14029.' // xs)
        call list_density % append(density * 0.04685_8)
        call list_names % append('14030.' // xs)
        call list_density % append(density * 0.03092_8)

      case ('p')
        call list_names % append('15031.' // xs)
        call list_density % append(density)

      case ('s')
        call list_names % append('16032.' // xs)
        call list_density % append(density * 0.9499_8)
        call list_names % append('16033.' // xs)
        call list_density % append(density * 0.0075_8)
        call list_names % append('16034.' // xs)
        call list_density % append(density * 0.0425_8)
        call list_names % append('16036.' // xs)
        call list_density % append(density * 0.0001_8)

      case ('cl')
        call list_names % append('17035.' // xs)
        call list_density % append(density * 0.7576_8)
        call list_names % append('17037.' // xs)
        call list_density % append(density * 0.2424_8)

      case ('ar')
        call list_names % append('18036.' // xs)
        call list_density % append(density * 0.003336_8)
        call list_names % append('18038.' // xs)
        call list_density % append(density * 0.000629_8)
        call list_names % append('18040.' // xs)
        call list_density % append(density * 0.996035_8)

      case ('k')
        call list_names % append('19039.' // xs)
        call list_density % append(density * 0.932581_8)
        call list_names % append('19040.' // xs)
        call list_density % append(density * 0.000117_8)
        call list_names % append('19041.' // xs)
        call list_density % append(density * 0.067302_8)

      case ('ca')
        call list_names % append('20040.' // xs)
        call list_density % append(density * 0.96941_8)
        call list_names % append('20042.' // xs)
        call list_density % append(density * 0.00647_8)
        call list_names % append('20043.' // xs)
        call list_density % append(density * 0.00135_8)
        call list_names % append('20044.' // xs)
        call list_density % append(density * 0.02086_8)
        call list_names % append('20046.' // xs)
        call list_density % append(density * 0.00004_8)
        call list_names % append('20048.' // xs)
        call list_density % append(density * 0.00187_8)

      case ('sc')
        call list_names % append('21045.' // xs)
        call list_density % append(density)

      case ('ti')
        call list_names % append('22046.' // xs)
        call list_density % append(density * 0.0825_8)
        call list_names % append('22047.' // xs)
        call list_density % append(density * 0.0744_8)
        call list_names % append('22048.' // xs)
        call list_density % append(density * 0.7372_8)
        call list_names % append('22049.' // xs)
        call list_density % append(density * 0.0541_8)
        call list_names % append('22050.' // xs)
        call list_density % append(density * 0.0518_8)

      case ('v')
        ! The evaluation of Vanadium in ENDF/B-VII.1 and JEFF 3.1.2 is a natural
        ! element. The IUPAC isotopic composition specifies the following
        ! breakdown which is not used:
        !   V-50 =  0.250%
        !   V-51 = 99.750%
        call list_names % append('23000.' // xs)
        call list_density % append(density)

      case ('cr')
        call list_names % append('24050.' // xs)
        call list_density % append(density * 0.04345_8)
        call list_names % append('24052.' // xs)
        call list_density % append(density * 0.83789_8)
        call list_names % append('24053.' // xs)
        call list_density % append(density * 0.09501_8)
        call list_names % append('24054.' // xs)
        call list_density % append(density * 0.02365_8)

      case ('mn')
        call list_names % append('25055.' // xs)
        call list_density % append(density)

      case ('fe')
        call list_names % append('26054.' // xs)
        call list_density % append(density * 0.05845_8)
        call list_names % append('26056.' // xs)
        call list_density % append(density * 0.91754_8)
        call list_names % append('26057.' // xs)
        call list_density % append(density * 0.02119_8)
        call list_names % append('26058.' // xs)
        call list_density % append(density * 0.00282_8)

      case ('co')
        call list_names % append('27059.' // xs)
        call list_density % append(density)

      case ('ni')
        call list_names % append('28058.' // xs)
        call list_density % append(density * 0.68077_8)
        call list_names % append('28060.' // xs)
        call list_density % append(density * 0.26223_8)
        call list_names % append('28061.' // xs)
        call list_density % append(density * 0.011399_8)
        call list_names % append('28062.' // xs)
        call list_density % append(density * 0.036346_8)
        call list_names % append('28064.' // xs)
        call list_density % append(density * 0.009255_8)

      case ('cu')
        call list_names % append('29063.' // xs)
        call list_density % append(density * 0.6915_8)
        call list_names % append('29065.' // xs)
        call list_density % append(density * 0.3085_8)

      case ('zn')
        ! The evaluation of Zinc in ENDF/B-VII.1 is a natural element. The IUPAC
        ! isotopic composition specifies the following breakdown which is not used
        ! here:
        !   Zn-64 = 48.63%
        !   Zn-66 = 27.90%
        !   Zn-67 =  4.10%
        !   Zn-68 = 18.75%
        !   Zn-70 =  0.62%
        call list_names % append('30000.' // xs)
        call list_density % append(density)

      case ('ga')
        ! JEFF 3.1.2 does not have evaluations for Ga-69 and Ga-71, only for
        ! natural Gallium, so this may cause problems.
        call list_names % append('31069.' // xs)
        call list_density % append(density * 0.60108_8)
        call list_names % append('31071.' // xs)
        call list_density % append(density * 0.39892_8)

      case ('ge')
        call list_names % append('32070.' // xs)
        call list_density % append(density * 0.2057_8)
        call list_names % append('32072.' // xs)
        call list_density % append(density * 0.2745_8)
        call list_names % append('32073.' // xs)
        call list_density % append(density * 0.0775_8)
        call list_names % append('32074.' // xs)
        call list_density % append(density * 0.3650_8)
        call list_names % append('32076.' // xs)
        call list_density % append(density * 0.0773_8)

      case ('as')
        call list_names % append('33075.' // xs)
        call list_density % append(density)

      case ('se')
        call list_names % append('34074.' // xs)
        call list_density % append(density * 0.0089_8)
        call list_names % append('34076.' // xs)
        call list_density % append(density * 0.0937_8)
        call list_names % append('34077.' // xs)
        call list_density % append(density * 0.0763_8)
        call list_names % append('34078.' // xs)
        call list_density % append(density * 0.2377_8)
        call list_names % append('34080.' // xs)
        call list_density % append(density * 0.4961_8)
        call list_names % append('34082.' // xs)
        call list_density % append(density * 0.0873_8)

      case ('br')
        call list_names % append('35079.' // xs)
        call list_density % append(density * 0.5069_8)
        call list_names % append('35081.' // xs)
        call list_density % append(density * 0.4931_8)

      case ('kr')
        call list_names % append('36078.' // xs)
        call list_density % append(density * 0.00355_8)
        call list_names % append('36080.' // xs)
        call list_density % append(density * 0.02286_8)
        call list_names % append('36082.' // xs)
        call list_density % append(density * 0.11593_8)
        call list_names % append('36083.' // xs)
        call list_density % append(density * 0.11500_8)
        call list_names % append('36084.' // xs)
        call list_density % append(density * 0.56987_8)
        call list_names % append('36086.' // xs)
        call list_density % append(density * 0.17279_8)

      case ('rb')
        call list_names % append('37085.' // xs)
        call list_density % append(density * 0.7217_8)
        call list_names % append('37087.' // xs)
        call list_density % append(density * 0.2783_8)

      case ('sr')
        call list_names % append('38084.' // xs)
        call list_density % append(density * 0.0056_8)
        call list_names % append('38086.' // xs)
        call list_density % append(density * 0.0986_8)
        call list_names % append('38087.' // xs)
        call list_density % append(density * 0.0700_8)
        call list_names % append('38088.' // xs)
        call list_density % append(density * 0.8258_8)

      case ('y')
        call list_names % append('39089.' // xs)
        call list_density % append(density)

      case ('zr')
        call list_names % append('40090.' // xs)
        call list_density % append(density * 0.5145_8)
        call list_names % append('40091.' // xs)
        call list_density % append(density * 0.1122_8)
        call list_names % append('40092.' // xs)
        call list_density % append(density * 0.1715_8)
        call list_names % append('40094.' // xs)
        call list_density % append(density * 0.1738_8)
        call list_names % append('40096.' // xs)
        call list_density % append(density * 0.0280_8)

      case ('nb')
        call list_names % append('41093.' // xs)
        call list_density % append(density)

      case ('mo')
        call list_names % append('42092.' // xs)
        call list_density % append(density * 0.1453_8)
        call list_names % append('42094.' // xs)
        call list_density % append(density * 0.0915_8)
        call list_names % append('42095.' // xs)
        call list_density % append(density * 0.1584_8)
        call list_names % append('42096.' // xs)
        call list_density % append(density * 0.1667_8)
        call list_names % append('42097.' // xs)
        call list_density % append(density * 0.0960_8)
        call list_names % append('42098.' // xs)
        call list_density % append(density * 0.2439_8)
        call list_names % append('42100.' // xs)
        call list_density % append(density * 0.0982_8)

      case ('ru')
        call list_names % append('44096.' // xs)
        call list_density % append(density * 0.0554_8)
        call list_names % append('44098.' // xs)
        call list_density % append(density * 0.0187_8)
        call list_names % append('44099.' // xs)
        call list_density % append(density * 0.1276_8)
        call list_names % append('44100.' // xs)
        call list_density % append(density * 0.1260_8)
        call list_names % append('44101.' // xs)
        call list_density % append(density * 0.1706_8)
        call list_names % append('44102.' // xs)
        call list_density % append(density * 0.3155_8)
        call list_names % append('44104.' // xs)
        call list_density % append(density * 0.1862_8)

      case ('rh')
        call list_names % append('45103.' // xs)
        call list_density % append(density)

      case ('pd')
        call list_names % append('46102.' // xs)
        call list_density % append(density * 0.0102_8)
        call list_names % append('46104.' // xs)
        call list_density % append(density * 0.1114_8)
        call list_names % append('46105.' // xs)
        call list_density % append(density * 0.2233_8)
        call list_names % append('46106.' // xs)
        call list_density % append(density * 0.2733_8)
        call list_names % append('46108.' // xs)
        call list_density % append(density * 0.2646_8)
        call list_names % append('46110.' // xs)
        call list_density % append(density * 0.1172_8)

      case ('ag')
        call list_names % append('47107.' // xs)
        call list_density % append(density * 0.51839_8)
        call list_names % append('47109.' // xs)
        call list_density % append(density * 0.48161_8)

      case ('cd')
        call list_names % append('48106.' // xs)
        call list_density % append(density * 0.0125_8)
        call list_names % append('48108.' // xs)
        call list_density % append(density * 0.0089_8)
        call list_names % append('48110.' // xs)
        call list_density % append(density * 0.1249_8)
        call list_names % append('48111.' // xs)
        call list_density % append(density * 0.1280_8)
        call list_names % append('48112.' // xs)
        call list_density % append(density * 0.2413_8)
        call list_names % append('48113.' // xs)
        call list_density % append(density * 0.1222_8)
        call list_names % append('48114.' // xs)
        call list_density % append(density * 0.2873_8)
        call list_names % append('48116.' // xs)
        call list_density % append(density * 0.0749_8)

      case ('in')
        call list_names % append('49113.' // xs)
        call list_density % append(density * 0.0429_8)
        call list_names % append('49115.' // xs)
        call list_density % append(density * 0.9571_8)

      case ('sn')
        call list_names % append('50112.' // xs)
        call list_density % append(density * 0.0097_8)
        call list_names % append('50114.' // xs)
        call list_density % append(density * 0.0066_8)
        call list_names % append('50115.' // xs)
        call list_density % append(density * 0.0034_8)
        call list_names % append('50116.' // xs)
        call list_density % append(density * 0.1454_8)
        call list_names % append('50117.' // xs)
        call list_density % append(density * 0.0768_8)
        call list_names % append('50118.' // xs)
        call list_density % append(density * 0.2422_8)
        call list_names % append('50119.' // xs)
        call list_density % append(density * 0.0859_8)
        call list_names % append('50120.' // xs)
        call list_density % append(density * 0.3258_8)
        call list_names % append('50122.' // xs)
        call list_density % append(density * 0.0463_8)
        call list_names % append('50124.' // xs)
        call list_density % append(density * 0.0579_8)

      case ('sb')
        call list_names % append('51121.' // xs)
        call list_density % append(density * 0.5721_8)
        call list_names % append('51123.' // xs)
        call list_density % append(density * 0.4279_8)

      case ('te')
        call list_names % append('52120.' // xs)
        call list_density % append(density * 0.0009_8)
        call list_names % append('52122.' // xs)
        call list_density % append(density * 0.0255_8)
        call list_names % append('52123.' // xs)
        call list_density % append(density * 0.0089_8)
        call list_names % append('52124.' // xs)
        call list_density % append(density * 0.0474_8)
        call list_names % append('52125.' // xs)
        call list_density % append(density * 0.0707_8)
        call list_names % append('52126.' // xs)
        call list_density % append(density * 0.1884_8)
        call list_names % append('52128.' // xs)
        call list_density % append(density * 0.3174_8)
        call list_names % append('52130.' // xs)
        call list_density % append(density * 0.3408_8)

      case ('i')
        call list_names % append('53127.' // xs)
        call list_density % append(density)

      case ('xe')
        call list_names % append('54124.' // xs)
        call list_density % append(density * 0.000952_8)
        call list_names % append('54126.' // xs)
        call list_density % append(density * 0.000890_8)
        call list_names % append('54128.' // xs)
        call list_density % append(density * 0.019102_8)
        call list_names % append('54129.' // xs)
        call list_density % append(density * 0.264006_8)
        call list_names % append('54130.' // xs)
        call list_density % append(density * 0.040710_8)
        call list_names % append('54131.' // xs)
        call list_density % append(density * 0.212324_8)
        call list_names % append('54132.' // xs)
        call list_density % append(density * 0.269086_8)
        call list_names % append('54134.' // xs)
        call list_density % append(density * 0.104357_8)
        call list_names % append('54136.' // xs)
        call list_density % append(density * 0.088573_8)

      case ('cs')
        call list_names % append('55133.' // xs)
        call list_density % append(density)

      case ('ba')
        call list_names % append('56130.' // xs)
        call list_density % append(density * 0.00106_8)
        call list_names % append('56132.' // xs)
        call list_density % append(density * 0.00101_8)
        call list_names % append('56134.' // xs)
        call list_density % append(density * 0.02417_8)
        call list_names % append('56135.' // xs)
        call list_density % append(density * 0.06592_8)
        call list_names % append('56136.' // xs)
        call list_density % append(density * 0.07854_8)
        call list_names % append('56137.' // xs)
        call list_density % append(density * 0.11232_8)
        call list_names % append('56138.' // xs)
        call list_density % append(density * 0.71698_8)

      case ('la')
        call list_names % append('57138.' // xs)
        call list_density % append(density * 0.0008881_8)
        call list_names % append('57139.' // xs)
        call list_density % append(density * 0.9991119_8)

      case ('ce')
        call list_names % append('58136.' // xs)
        call list_density % append(density * 0.00185_8)
        call list_names % append('58138.' // xs)
        call list_density % append(density * 0.00251_8)
        call list_names % append('58140.' // xs)
        call list_density % append(density * 0.88450_8)
        call list_names % append('58142.' // xs)
        call list_density % append(density * 0.11114_8)

      case ('pr')
        call list_names % append('59141.' // xs)
        call list_density % append(density)

      case ('nd')
        call list_names % append('60142.' // xs)
        call list_density % append(density * 0.27152_8)
        call list_names % append('60143.' // xs)
        call list_density % append(density * 0.12174_8)
        call list_names % append('60144.' // xs)
        call list_density % append(density * 0.23798_8)
        call list_names % append('60145.' // xs)
        call list_density % append(density * 0.08293_8)
        call list_names % append('60146.' // xs)
        call list_density % append(density * 0.17189_8)
        call list_names % append('60148.' // xs)
        call list_density % append(density * 0.05756_8)
        call list_names % append('60150.' // xs)
        call list_density % append(density * 0.05638_8)

      case ('sm')
        call list_names % append('62144.' // xs)
        call list_density % append(density * 0.0307_8)
        call list_names % append('62147.' // xs)
        call list_density % append(density * 0.1499_8)
        call list_names % append('62148.' // xs)
        call list_density % append(density * 0.1124_8)
        call list_names % append('62149.' // xs)
        call list_density % append(density * 0.1382_8)
        call list_names % append('62150.' // xs)
        call list_density % append(density * 0.0738_8)
        call list_names % append('62152.' // xs)
        call list_density % append(density * 0.2675_8)
        call list_names % append('62154.' // xs)
        call list_density % append(density * 0.2275_8)

      case ('eu')
        call list_names % append('63151.' // xs)
        call list_density % append(density * 0.4781_8)
        call list_names % append('63153.' // xs)
        call list_density % append(density * 0.5219_8)

      case ('gd')
        call list_names % append('64152.' // xs)
        call list_density % append(density * 0.0020_8)
        call list_names % append('64154.' // xs)
        call list_density % append(density * 0.0218_8)
        call list_names % append('64155.' // xs)
        call list_density % append(density * 0.1480_8)
        call list_names % append('64156.' // xs)
        call list_density % append(density * 0.2047_8)
        call list_names % append('64157.' // xs)
        call list_density % append(density * 0.1565_8)
        call list_names % append('64158.' // xs)
        call list_density % append(density * 0.2484_8)
        call list_names % append('64160.' // xs)
        call list_density % append(density * 0.2186_8)

      case ('tb')
        call list_names % append('65159.' // xs)
        call list_density % append(density)

      case ('dy')
        call list_names % append('66156.' // xs)
        call list_density % append(density * 0.00056_8)
        call list_names % append('66158.' // xs)
        call list_density % append(density * 0.00095_8)
        call list_names % append('66160.' // xs)
        call list_density % append(density * 0.02329_8)
        call list_names % append('66161.' // xs)
        call list_density % append(density * 0.18889_8)
        call list_names % append('66162.' // xs)
        call list_density % append(density * 0.25475_8)
        call list_names % append('66163.' // xs)
        call list_density % append(density * 0.24896_8)
        call list_names % append('66164.' // xs)
        call list_density % append(density * 0.28260_8)

      case ('ho')
        call list_names % append('67165.' // xs)
        call list_density % append(density)

      case ('er')
        call list_names % append('68162.' // xs)
        call list_density % append(density * 0.00139_8)
        call list_names % append('68164.' // xs)
        call list_density % append(density * 0.01601_8)
        call list_names % append('68166.' // xs)
        call list_density % append(density * 0.33503_8)
        call list_names % append('68167.' // xs)
        call list_density % append(density * 0.22869_8)
        call list_names % append('68168.' // xs)
        call list_density % append(density * 0.26978_8)
        call list_names % append('68170.' // xs)
        call list_density % append(density * 0.14910_8)

      case ('tm')
        call list_names % append('69169.' // xs)
        call list_density % append(density)

      case ('yb')
        call list_names % append('70168.' // xs)
        call list_density % append(density * 0.00123_8)
        call list_names % append('70170.' // xs)
        call list_density % append(density * 0.02982_8)
        call list_names % append('70171.' // xs)
        call list_density % append(density * 0.1409_8)
        call list_names % append('70172.' // xs)
        call list_density % append(density * 0.2168_8)
        call list_names % append('70173.' // xs)
        call list_density % append(density * 0.16103_8)
        call list_names % append('70174.' // xs)
        call list_density % append(density * 0.32026_8)
        call list_names % append('70176.' // xs)
        call list_density % append(density * 0.12996_8)

      case ('lu')
        call list_names % append('71175.' // xs)
        call list_density % append(density * 0.97401_8)
        call list_names % append('71176.' // xs)
        call list_density % append(density * 0.02599_8)

      case ('hf')
        call list_names % append('72174.' // xs)
        call list_density % append(density * 0.0016_8)
        call list_names % append('72176.' // xs)
        call list_density % append(density * 0.0526_8)
        call list_names % append('72177.' // xs)
        call list_density % append(density * 0.1860_8)
        call list_names % append('72178.' // xs)
        call list_density % append(density * 0.2728_8)
        call list_names % append('72179.' // xs)
        call list_density % append(density * 0.1362_8)
        call list_names % append('72180.' // xs)
        call list_density % append(density * 0.3508_8)

      case ('ta')
        call list_names % append('73180.' // xs)
        call list_density % append(density * 0.0001201_8)
        call list_names % append('73181.' // xs)
        call list_density % append(density * 0.9998799_8)

      case ('w')
        ! ENDF/B-VII.0 does not have W-180 so this may cause problems. However, it
        ! has been added as of ENDF/B-VII.1
        call list_names % append('74180.' // xs)
        call list_density % append(density * 0.0012_8)
        call list_names % append('74182.' // xs)
        call list_density % append(density * 0.2650_8)
        call list_names % append('74183.' // xs)
        call list_density % append(density * 0.1431_8)
        call list_names % append('74184.' // xs)
        call list_density % append(density * 0.3064_8)
        call list_names % append('74186.' // xs)
        call list_density % append(density * 0.2843_8)

      case ('re')
        call list_names % append('75185.' // xs)
        call list_density % append(density * 0.3740_8)
        call list_names % append('75187.' // xs)
        call list_density % append(density * 0.6260_8)

      case ('os')
        call list_names % append('76184.' // xs)
        call list_density % append(density * 0.0002_8)
        call list_names % append('76186.' // xs)
        call list_density % append(density * 0.0159_8)
        call list_names % append('76187.' // xs)
        call list_density % append(density * 0.0196_8)
        call list_names % append('76188.' // xs)
        call list_density % append(density * 0.1324_8)
        call list_names % append('76189.' // xs)
        call list_density % append(density * 0.1615_8)
        call list_names % append('76190.' // xs)
        call list_density % append(density * 0.2626_8)
        call list_names % append('76192.' // xs)
        call list_density % append(density * 0.4078_8)

      case ('ir')
        call list_names % append('77191.' // xs)
        call list_density % append(density * 0.373_8)
        call list_names % append('77193.' // xs)
        call list_density % append(density * 0.627_8)

      case ('pt')
        call list_names % append('78190.' // xs)
        call list_density % append(density * 0.00012_8)
        call list_names % append('78192.' // xs)
        call list_density % append(density * 0.00782_8)
        call list_names % append('78194.' // xs)
        call list_density % append(density * 0.3286_8)
        call list_names % append('78195.' // xs)
        call list_density % append(density * 0.3378_8)
        call list_names % append('78196.' // xs)
        call list_density % append(density * 0.2521_8)
        call list_names % append('78198.' // xs)
        call list_density % append(density * 0.07356_8)

      case ('au')
        call list_names % append('79197.' // xs)
        call list_density % append(density)

      case ('hg')
        call list_names % append('80196.' // xs)
        call list_density % append(density * 0.0015_8)
        call list_names % append('80198.' // xs)
        call list_density % append(density * 0.0997_8)
        call list_names % append('80199.' // xs)
        call list_density % append(density * 0.1687_8)
        call list_names % append('80200.' // xs)
        call list_density % append(density * 0.2310_8)
        call list_names % append('80201.' // xs)
        call list_density % append(density * 0.1318_8)
        call list_names % append('80202.' // xs)
        call list_density % append(density * 0.2986_8)
        call list_names % append('80204.' // xs)
        call list_density % append(density * 0.0687_8)

      case ('tl')
        call list_names % append('81203.' // xs)
        call list_density % append(density * 0.2952_8)
        call list_names % append('81205.' // xs)
        call list_density % append(density * 0.7048_8)

      case ('pb')
        call list_names % append('82204.' // xs)
        call list_density % append(density * 0.014_8)
        call list_names % append('82206.' // xs)
        call list_density % append(density * 0.241_8)
        call list_names % append('82207.' // xs)
        call list_density % append(density * 0.221_8)
        call list_names % append('82208.' // xs)
        call list_density % append(density * 0.524_8)

      case ('bi')
        call list_names % append('83209.' // xs)
        call list_density % append(density)

      case ('th')
        call list_names % append('90232.' // xs)
        call list_density % append(density)

      case ('pa')
        call list_names % append('91231.' // xs)
        call list_density % append(density)

      case ('u')
        call list_names % append('92234.' // xs)
        call list_density % append(density * 0.000054_8)
        call list_names % append('92235.' // xs)
        call list_density % append(density * 0.007204_8)
        call list_names % append('92238.' // xs)
        call list_density % append(density * 0.992742_8)

      case default
        message = "Cannot expand element: " // name
        call fatal_error()

      end select

    end subroutine expand_natural_element





end program main
