program main

use read_xml_primitives
use write_xml_primitives
use xmlparse
use global
use constants
use error
use dict_header,      only: DictIntInt, ElemKeyValueCI, dict_clear_ii, &
  dict_clear_ci
use string,           only: lower_case, to_str, str_to_int, str_to_real, &
  starts_with, ends_with
use list_header,      only: ListChar, ListReal
use geometry_header,  only: Cell, Surface, Lattice
use tally_header,     only: TallyObject, TallyFilter
use tally_initialize, only: add_tallies
use xml_data_geometry_t

implicit none

type(DictIntInt) :: cells_in_univ_dict ! used to count how many cells each
! universe contains

character(len=50) :: outfile = 'benchmark.txt'
character(len=50) :: mode
character(len=2) :: path_input = './'
double precision :: time
integer :: fsize, istat

call get_command_argument(1, mode)

select case(trim(mode))
case ('geometry')

  open(unit=500, file='geometry.xml', status='old')
  inquire( 500, size=fsize)
  close(unit=500)

  print *, "***calling read_xml_file_geometry (v4) with filesize = ", fsize
  call benchmark_geometry(time)
  print *, 'read_xml_file_geometry (v4) completed in (seconds):', time


  open(unit=600, file='geometry_'//trim(outfile), status='unknown', position='append')
  write(600, '(a20, i20, f20.10)') 'v4', fsize, time
  close(unit=600)

case ('tallies')

  open(unit=500, file='tallies.xml', status='old')
  inquire( 500, size=fsize)
  close(unit=500)

  print *, "***calling read_xml_file_tallies (v4) with filesize = ", fsize
  call benchmark_tallies(time)
  print *, 'read_xml_file_tallies (v4) completed in (seconds):', time

  open(unit=600, file='tallies_'//trim(outfile), status='unknown', position='append')
  write(600, '(a20, i20, f20.10)') 'v4', fsize, time
  close(unit=600)

case default
  print *, 'Error: must specify mode for benchmark (geometry, tallies)'
end select


contains


subroutine benchmark_geometry(time)
  double precision, intent(out) :: time
  integer :: cstart, cend, crate

  call system_clock(cstart, crate)

  call read_geometry_xml()

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

end subroutine benchmark_geometry

subroutine benchmark_tallies(time)
  double precision, intent(out) :: time
  integer :: cstart, cend, crate

  call system_clock(cstart, crate)

  call read_tallies_xml()

  call system_clock(cend, crate)
  time = dble(cend-cstart)/crate 

  if (allocated(meshes)) deallocate(meshes)
  if (allocated(tallies)) deallocate(tallies)
  !if (associated(user_tallies)) deallocate(user_tallies)
  !if (associated(cmfd_tallies)) deallocate(cmfd_tallies)
  call dict_clear_ii(mesh_dict)
  call dict_clear_ii(tally_dict)
  call dict_clear_ci(nuclide_dict)

end subroutine benchmark_tallies

subroutine read_tallies_xml

  use xml_data_tallies_t

  integer :: i             ! loop over user-specified tallies
  integer :: j             ! loop over words
  integer :: k             ! another loop index
  integer :: l             ! another loop index
  integer :: id            ! user-specified identifier
  integer :: i_mesh        ! index in meshes array
  integer :: n             ! size of arrays in mesh specification
  integer :: n_words       ! number of words read
  integer :: n_filters     ! number of filters
  integer :: n_new         ! number of new scores to add based on Pn tally
  integer :: n_scores      ! number of tot scores after adjusting for Pn tally
  integer :: n_order       ! Scattering order requested
  integer :: n_order_pos   ! Position of Scattering order in score name string
  integer :: MT            ! user-specified MT for score
  logical :: file_exists   ! does tallies.xml file exist?
  character(MAX_LINE_LEN) :: filename
  character(MAX_WORD_LEN) :: word
  character(MAX_WORD_LEN) :: score_name
  type(ElemKeyValueCI), pointer :: pair_list => null()
  type(TallyObject),    pointer :: t => null()
  type(StructuredMesh), pointer :: m => null()
  type(TallyFilter), allocatable :: filters(:) ! temporary filters

  ! Check if tallies.xml exists
  filename = trim(path_input) // "tallies.xml"
  inquire(FILE=filename, EXIST=file_exists)
  if (.not. file_exists) then
    ! Since a tallies.xml file is optional, no error is issued here
    return
  end if

  ! Display output message
  message = "Reading tallies XML file..."
  call write_message(5)

  ! Parse tallies.xml file
  call read_xml_file_tallies_t(filename)

  ! ==========================================================================
  ! DETERMINE SIZE OF ARRAYS AND ALLOCATE

  ! Check for user meshes
  if (.not. associated(mesh_)) then
    n_user_meshes = 0
  else
    n_user_meshes = size(mesh_)
    if (cmfd_run) then
      n_meshes = n_user_meshes + n_cmfd_meshes
    else
      n_meshes = n_user_meshes
    end if
  end if

  ! Allocate mesh array
  if (n_meshes > 0) allocate(meshes(n_meshes))

  ! Check for user tallies
  if (.not. associated(tally_)) then
    n_user_tallies = 0
    message = "No tallies present in tallies.xml file!"
    call warning()
  else
    n_user_tallies = size(tally_)
  end if

  ! Allocate tally array
  if (n_user_tallies > 0) then
    call add_tallies("user", n_user_tallies)
  end if

  ! Check for <assume_separate> setting
  call lower_case(separate_)
  if (separate_ == 'true' .or. separate_ == '1') assume_separate = .true.

  ! ==========================================================================
  ! READ MESH DATA

  do i = 1, n_user_meshes
    m => meshes(i)

    ! copy mesh id
    m % id = mesh_(i) % id

    ! Check to make sure 'id' hasn't been used
    if (mesh_dict % has_key(m % id)) then
      message = "Two or more meshes use the same unique ID: " // &
        to_str(m % id)
      call fatal_error()
    end if

    ! Read mesh type
    call lower_case(mesh_(i) % type)
    select case (mesh_(i) % type)
    case ('rect', 'rectangle', 'rectangular')
      m % type = LATTICE_RECT
    case ('hex', 'hexagon', 'hexagonal')
      m % type = LATTICE_HEX
    case default
      message = "Invalid mesh type: " // trim(mesh_(i) % type)
      call fatal_error()
    end select

    ! Determine number of dimensions for mesh
    n = size(mesh_(i) % dimension)
    if (n /= 2 .and. n /= 3) then
      message = "Mesh must be two or three dimensions."
      call fatal_error()
    end if
    m % n_dimension = n

    ! Allocate attribute arrays
    allocate(m % dimension(n))
    allocate(m % lower_left(n))
    allocate(m % width(n))
    allocate(m % upper_right(n))

    ! Check that dimensions are all greater than zero
    if (any(mesh_(i) % dimension <= 0)) then
      message = "All entries on the <dimension> element for a tally mesh &
        &must be positive."
      call fatal_error()
    end if

    ! Read dimensions in each direction
    m % dimension = mesh_(i) % dimension

    ! Read mesh lower-left corner location
    if (m % n_dimension /= size(mesh_(i) % lower_left)) then
      message = "Number of entries on <lower_left> must be the same as &
        &the number of entries on <dimension>."
      call fatal_error()
    end if
    m % lower_left = mesh_(i) % lower_left

    ! Make sure either upper-right or width was specified
    if (associated(mesh_(i) % upper_right) .and. &
      associated(mesh_(i) % width)) then
    message = "Cannot specify both <upper_right> and <width> on a &
      &tally mesh."
    call fatal_error()
  end if

  ! Make sure either upper-right or width was specified
  if (.not. associated(mesh_(i) % upper_right) .and. &
    .not. associated(mesh_(i) % width)) then
  message = "Must specify either <upper_right> and <width> on a &
    &tally mesh."
  call fatal_error()
end if

if (associated(mesh_(i) % width)) then
  ! Check to ensure width has same dimensions
  if (size(mesh_(i) % width) /= size(mesh_(i) % lower_left)) then
    message = "Number of entries on <width> must be the same as the &
      &number of entries on <lower_left>."
    call fatal_error()
  end if

  ! Check for negative widths
  if (any(mesh_(i) % width < ZERO)) then
    message = "Cannot have a negative <width> on a tally mesh."
    call fatal_error()
  end if

  ! Set width and upper right coordinate
  m % width = mesh_(i) % width
  m % upper_right = m % lower_left + m % dimension * m % width

  elseif (associated(mesh_(i) % upper_right)) then
  ! Check to ensure width has same dimensions
  if (size(mesh_(i) % upper_right) /= size(mesh_(i) % lower_left)) then
    message = "Number of entries on <upper_right> must be the same as &
      &the number of entries on <lower_left>."
    call fatal_error()
  end if

  ! Check that upper-right is above lower-left
  if (any(mesh_(i) % upper_right < mesh_(i) % lower_left)) then
    message = "The <upper_right> coordinates must be greater than the &
      &<lower_left> coordinates on a tally mesh."
    call fatal_error()
  end if

  ! Set width and upper right coordinate
  m % upper_right = mesh_(i) % upper_right
  m % width = (m % upper_right - m % lower_left) / m % dimension
end if

! Set volume fraction
m % volume_frac = ONE/real(product(m % dimension),8)

! Add mesh to dictionary
call mesh_dict % add_key(m % id, i)
    end do

    ! ==========================================================================
    ! READ TALLY DATA

    READ_TALLIES: do i = 1, n_user_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Set tally type to volume by default
      t % type = TALLY_VOLUME

      ! It's desirable to use a track-length esimator for tallies since
      ! generally more events will score to the tally, reducing the
      ! variance. However, for tallies that require information on
      ! post-collision parameters (e.g. tally with an energyout filter) the
      ! analog esimator must be used.

      t % estimator = ESTIMATOR_TRACKLENGTH

      ! Copy material id
      t % id = tally_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (tally_dict % has_key(t % id)) then
        message = "Two or more tallies use the same unique ID: " // &
          to_str(t % id)
        call fatal_error()
      end if

      ! Copy tally label
      t % label = tally_(i) % label

      ! =======================================================================
      ! READ DATA FOR FILTERS

      ! In older versions, tally filters were specified with a <filters>
      ! element followed by sub-elements <cell>, <mesh>, etc. This checks for
      ! the old format and if it is present, raises an error

      if (size(tally_(i) % filters) > 0) then
        message = "Tally filters should be specified with multiple <filter> &
          &elements. Did you forget to change your <filters> element?"
        call fatal_error()
      end if

      if (associated(tally_(i) % filter)) then
        ! Determine number of filters
        n_filters = size(tally_(i) % filter)

        ! Allocate filters array
        t % n_filters = n_filters
        allocate(t % filters(n_filters))

        READ_FILTERS: do j = 1, n_filters
          ! Convert filter type to lower case
          call lower_case(tally_(i) % filter(j) % type)

          ! Determine number of bins
          n_words = size(tally_(i) % filter(j) % bins)

          ! Determine type of filter
          select case (tally_(i) % filter(j) % type)
          case ('cell')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('cellborn')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELLBORN

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('material')
            ! Set type of filter
            t % filters(j) % type = FILTER_MATERIAL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('universe')
            ! Set type of filter
            t % filters(j) % type = FILTER_UNIVERSE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('surface')
            message = "Surface filter is not yet supported!"
            call fatal_error()

            ! Set type of filter
            t % filters(j) % type = FILTER_SURFACE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('mesh')
            ! Set type of filter
            t % filters(j) % type = FILTER_MESH

            ! Check to make sure multiple meshes weren't given
            if (n_words /= 1) then
              message = "Can only have one mesh filter specified."
              call fatal_error()
            end if

            ! Determine id of mesh
            id = int(str_to_int(tally_(i) % filter(j) % bins(1)),4)

            ! Get pointer to mesh
            if (mesh_dict % has_key(id)) then
              i_mesh = mesh_dict % get_key(id)
              m => meshes(i_mesh)
            else
              message = "Could not find mesh " // trim(to_str(id)) // &
                " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if

            ! Determine number of bins -- this is assuming that the tally is
            ! a volume tally and not a surface current tally. If it is a
            ! surface current tally, the number of bins will get reset later
            t % filters(j) % n_bins = product(m % dimension)

            ! Allocate and store index of mesh
            allocate(t % filters(j) % int_bins(1))
            t % filters(j) % int_bins(1) = i_mesh

          case ('energy')
            ! Set type of filter
            t % filters(j) % type = FILTER_ENERGYIN

            ! Set number of bins
            t % filters(j) % n_bins = n_words - 1

            ! Allocate and store bins
            allocate(t % filters(j) % real_bins(n_words))
            do k = 1, n_words
              t % filters(j) % real_bins(k) = str_to_real(&
                tally_(i) % filter(j) % bins(k))
            end do

          case ('energyout')
            ! Set type of filter
            t % filters(j) % type = FILTER_ENERGYOUT

            ! Set number of bins
            t % filters(j) % n_bins = n_words - 1

            ! Allocate and store bins
            allocate(t % filters(j) % real_bins(n_words))
            do k = 1, n_words
              t % filters(j) % real_bins(k) = str_to_real(&
                tally_(i) % filter(j) % bins(k))
            end do

            ! Set to analog estimator
            t % estimator = ESTIMATOR_ANALOG

          case default
            ! Specified tally filter is invalid, raise error
            message = "Unknown filter type '" // trim(tally_(i) % &
              filter(j) % type) // "' on tally " // &
              trim(to_str(t % id)) // "."
            call fatal_error()

          end select

          ! Set find_filter, e.g. if filter(3) has type FILTER_CELL, then
          ! find_filter(FILTER_CELL) would be set to 3.

          t % find_filter(t % filters(j) % type) = j

        end do READ_FILTERS

        ! Check that both cell and surface weren't specified
        if (t % find_filter(FILTER_CELL) > 0 .and. &
          t % find_filter(FILTER_SURFACE) > 0) then
        message = "Cannot specify both cell and surface filters for tally " &
          // trim(to_str(t % id))
        call fatal_error()
      end if

    else
      ! No filters were specified
      t % n_filters = 0
    end if

    ! =======================================================================
    ! READ DATA FOR NUCLIDES

    if (associated(tally_(i) % nuclides)) then
      if (tally_(i) % nuclides(1) == 'all') then
        ! Handle special case <nuclides>all</nuclides>
        allocate(t % nuclide_bins(n_nuclides_total + 1))

        ! Set bins to 1, 2, 3, ..., n_nuclides_total, -1
        t % nuclide_bins(1:n_nuclides_total) = &
          (/ (j, j=1, n_nuclides_total) /)
        t % nuclide_bins(n_nuclides_total + 1) = -1

        ! Set number of nuclide bins
        t % n_nuclide_bins = n_nuclides_total + 1

        ! Set flag so we can treat this case specially
        t % all_nuclides = .true.
      else
        ! Any other case, e.g. <nuclides>U-235 Pu-239</nuclides>
        n_words = size(tally_(i) % nuclides)
        allocate(t % nuclide_bins(n_words))
        do j = 1, n_words
          ! Check if total material was specified
          if (tally_(i) % nuclides(j) == 'total') then
            t % nuclide_bins(j) = -1
            cycle
          end if

          ! Check if xs specifier was given
          if (ends_with(tally_(i) % nuclides(j), 'c')) then
            word = tally_(i) % nuclides(j)
          else
            if (default_xs == '') then
              ! No default cross section specified, search through nuclides
              pair_list => nuclide_dict % keys()
              do while (associated(pair_list))
                if (starts_with(pair_list % key, &
                  tally_(i) % nuclides(j))) then
                word = pair_list % key(1:150)
                exit
              end if

              ! Advance to next
              pair_list => pair_list % next
            end do

            ! Check if no nuclide was found
            if (.not. associated(pair_list)) then
              message = "Could not find the nuclide " // trim(&
                tally_(i) % nuclides(j)) // " specified in tally " &
                // trim(to_str(t % id)) // " in any material."
              call fatal_error()
            end if
            deallocate(pair_list)
          else
            ! Set nuclide to default xs
            word = trim(tally_(i) % nuclides(j)) // "." // default_xs
          end if
        end if

        ! Check to make sure nuclide specified is in problem
        if (.not. nuclide_dict % has_key(word)) then
          message = "The nuclide " // trim(word) // " from tally " // &
            trim(to_str(t % id)) // " is not present in any material."
          call fatal_error()
        end if

        ! Set bin to index in nuclides array
        t % nuclide_bins(j) = nuclide_dict % get_key(word)
      end do

      ! Set number of nuclide bins
      t % n_nuclide_bins = n_words
    end if

  else
    ! No <nuclides> were specified -- create only one bin will be added
    ! for the total material.
    allocate(t % nuclide_bins(1))
    t % nuclide_bins(1) = -1
    t % n_nuclide_bins = 1
  end if

  ! =======================================================================
  ! READ DATA FOR SCORES

  if (associated(tally_(i) % scores)) then
    ! Loop through scores and determine if a scatter-p# input was used
    ! to allow for proper pre-allocating of t % score_bins
    ! This scheme allows multiple scatter-p# to be requested by the user
    ! if so desired
    n_words = size(tally_(i) % scores)
    n_new = 0
    do j = 1, n_words
      call lower_case(tally_(i) % scores(j))
      ! Find if scores(j) is of the form 'scatter-p'
      ! If so, get the number and do a select case on that.
      score_name = tally_(i) % scores(j)
      if (starts_with(score_name,'scatter-p')) then
        n_order_pos = scan(score_name,'0123456789')
        n_order = int(str_to_int( &
          score_name(n_order_pos:(len_trim(score_name)))),4)
        if (n_order > SCATT_ORDER_MAX) then
          ! Throw a warning. Set to the maximum number.
          ! The above scheme will essentially take the absolute value
          message = "Invalid scattering order of " // trim(to_str(n_order)) // &
            " requested. Setting to the maximum permissible value, " // &
            trim(to_str(SCATT_ORDER_MAX))
          call warning()
          n_order = SCATT_ORDER_MAX
          tally_(i) % scores(j) = SCATT_ORDER_MAX_PNSTR
        end if
        n_new = n_new + n_order
      end if
    end do
    n_scores = n_words + n_new

    ! Allocate accordingly
    allocate(t % score_bins(n_scores))
    allocate(t % scatt_order(n_scores))
    t % scatt_order = 0
    j = 0
    do l = 1, n_words
      j = j + 1
      ! Get the input string in scores(l) but if scatter-n or scatter-pn
      ! then strip off the n, and store it as an integer to be used later
      ! Peform the select case on this modified (number removed) string
      score_name = tally_(i) % scores(l)
      if (starts_with(score_name,'scatter-p')) then
        n_order_pos = scan(score_name,'0123456789')
        n_order = int(str_to_int( &
          score_name(n_order_pos:(len_trim(score_name)))),4)
        if (n_order > SCATT_ORDER_MAX) then
          ! Throw a warning. Set to the maximum number.
          ! The above scheme will essentially take the absolute value
          message = "Invalid scattering order of " // trim(to_str(n_order)) // &
            " requested. Setting to the maximum permissible value, " // &
            trim(to_str(SCATT_ORDER_MAX))
          call warning()
          n_order = SCATT_ORDER_MAX
        end if
        score_name = "scatter-pn"
      else if (starts_with(score_name,'scatter-')) then
        n_order_pos = scan(score_name,'0123456789')
        n_order = int(str_to_int( &
          score_name(n_order_pos:(len_trim(score_name)))),4)
        if (n_order > SCATT_ORDER_MAX) then
          ! Throw a warning. Set to the maximum number.
          ! The above scheme will essentially take the absolute value
          message = "Invalid scattering order of " // trim(to_str(n_order)) // &
            " requested. Setting to the maximum permissible value, " // &
            trim(to_str(SCATT_ORDER_MAX))
          call warning()
          n_order = SCATT_ORDER_MAX
        end if
        score_name = "scatter-n"
      end if

      select case (score_name)
      case ('flux')
        ! Prohibit user from tallying flux for an individual nuclide
        if (.not. (t % n_nuclide_bins == 1 .and. &
          t % nuclide_bins(1) == -1)) then
        message = "Cannot tally flux for an individual nuclide."
        call fatal_error()
      end if

      t % score_bins(j) = SCORE_FLUX
      if (t % find_filter(FILTER_ENERGYOUT) > 0) then
        message = "Cannot tally flux with an outgoing energy filter."
        call fatal_error()
      end if
    case ('total')
      t % score_bins(j) = SCORE_TOTAL
      if (t % find_filter(FILTER_ENERGYOUT) > 0) then
        message = "Cannot tally total reaction rate with an &
          &outgoing energy filter."
        call fatal_error()
      end if
    case ('scatter')
      t % score_bins(j) = SCORE_SCATTER

    case ('nu-scatter')
      t % score_bins(j) = SCORE_NU_SCATTER

      ! Set tally estimator to analog
      t % estimator = ESTIMATOR_ANALOG
    case ('scatter-n')
      if (n_order == 0) then
        t % score_bins(j) = SCORE_SCATTER
      else
        t % score_bins(j) = SCORE_SCATTER_N
        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG
      end if
      t % scatt_order(j) = n_order

    case ('scatter-pn')
      t % estimator = ESTIMATOR_ANALOG
      ! Setup P0:Pn
      t % score_bins(j : j + n_order) = SCORE_SCATTER_PN
      t % scatt_order(j : j + n_order) = n_order
      j = j + n_order

    case('transport')
      t % score_bins(j) = SCORE_TRANSPORT

      ! Set tally estimator to analog
      t % estimator = ESTIMATOR_ANALOG
    case ('diffusion')
      message = "Diffusion score no longer supported for tallies, & 
        &please remove"
      call fatal_error()
    case ('n1n')
      t % score_bins(j) = SCORE_N_1N

      ! Set tally estimator to analog
      t % estimator = ESTIMATOR_ANALOG
    case ('n2n')
      t % score_bins(j) = N_2N

    case ('n3n')
      t % score_bins(j) = N_3N

    case ('n4n')
      t % score_bins(j) = N_4N

    case ('absorption')
      t % score_bins(j) = SCORE_ABSORPTION
      if (t % find_filter(FILTER_ENERGYOUT) > 0) then
        message = "Cannot tally absorption rate with an outgoing &
          &energy filter."
        call fatal_error()
      end if
    case ('fission')
      t % score_bins(j) = SCORE_FISSION
      if (t % find_filter(FILTER_ENERGYOUT) > 0) then
        message = "Cannot tally fission rate with an outgoing &
          &energy filter."
        call fatal_error()
      end if
    case ('nu-fission')
      t % score_bins(j) = SCORE_NU_FISSION
      if (t % find_filter(FILTER_ENERGYOUT) > 0) then
        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG
      end if
    case ('kappa-fission')
      t % score_bins(j) = SCORE_KAPPA_FISSION
    case ('current')
      t % score_bins(j) = SCORE_CURRENT
      t % type = TALLY_SURFACE_CURRENT

      ! Check to make sure that current is the only desired response
      ! for this tally
      if (n_words > 1) then
        message = "Cannot tally other scoring functions in the same &
          &tally as surface currents. Separate other scoring &
          &functions into a distinct tally."
        call fatal_error()
      end if

      ! Since the number of bins for the mesh filter was already set
      ! assuming it was a volume tally, we need to adjust the number
      ! of bins

      ! Get index of mesh filter
      k = t % find_filter(FILTER_MESH)

      ! Get pointer to mesh
      i_mesh = t % filters(k) % int_bins(1)
      m => meshes(i_mesh)

      ! We need to increase the dimension by one since we also need
      ! currents coming into and out of the boundary mesh cells.
      t % filters(k) % n_bins = product(m % dimension + 1)

      ! Copy filters to temporary array
      allocate(filters(t % n_filters + 1))
      filters(1:t % n_filters) = t % filters

      ! Move allocation back -- filters becomes deallocated during
      ! this call
      call move_alloc(FROM=filters, TO=t%filters)

      ! Add surface filter
      t % n_filters = t % n_filters + 1
      t % filters(t % n_filters) % type = FILTER_SURFACE
      t % filters(t % n_filters) % n_bins = 2 * m % n_dimension
      allocate(t % filters(t % n_filters) % int_bins(&
        2 * m % n_dimension))
      if (m % n_dimension == 2) then
        t % filters(t % n_filters) % int_bins = (/ IN_RIGHT, &
          OUT_RIGHT, IN_FRONT, OUT_FRONT /)
        elseif (m % n_dimension == 3) then
        t % filters(t % n_filters) % int_bins = (/ IN_RIGHT, &
          OUT_RIGHT, IN_FRONT, OUT_FRONT, IN_TOP, OUT_TOP /)
      end if
      t % find_filter(FILTER_SURFACE) = t % n_filters

    case ('events')
      t % score_bins(j) = SCORE_EVENTS

    case default
      ! Assume that user has specified an MT number
      MT = int(str_to_int(score_name))

      if (MT /= ERROR_INT) then
        ! Specified score was an integer
        if (MT > 1) then
          t % score_bins(j) = MT
        else
          message = "Invalid MT on <scores>: " // &
            trim(tally_(i) % scores(j))
          call fatal_error()
        end if

      else
        ! Specified score was not an integer
        message = "Unknown scoring function: " // &
          trim(tally_(i) % scores(j))
        call fatal_error()
      end if

    end select
  end do
  t % n_score_bins = n_scores
  t % n_user_score_bins = n_words
else
  message = "No <scores> specified on tally " // trim(to_str(t % id)) &
    // "."
  call fatal_error()
end if

! =======================================================================
! SET TALLY ESTIMATOR

! Check if user specified estimator
if (len_trim(tally_(i) % estimator) > 0) then
  select case(tally_(i) % estimator)
  case ('analog')
    t % estimator = ESTIMATOR_ANALOG

  case ('tracklength', 'track-length', 'pathlength', 'path-length')
    ! If the estimator was set to an analog estimator, this means the
    ! tally needs post-collision information
    if (t % estimator == ESTIMATOR_ANALOG) then
      message = "Cannot use track-length estimator for tally " &
        // to_str(t % id)
      call fatal_error()
    end if

    ! Set estimator to track-length estimator
    t % estimator = ESTIMATOR_TRACKLENGTH

  case default
    message = "Invalid estimator '" // trim(tally_(i) % estimator) &
      // "' on tally " // to_str(t % id)
    call fatal_error()
  end select
  end if

  ! Add tally to dictionary
  call tally_dict % add_key(t % id, i)

end do READ_TALLIES

  end subroutine read_tallies_xml
  subroutine read_geometry_xml

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

  end subroutine read_geometry_xml


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
