module global
  use constants
  use dict_header,      only: DictCharInt, DictIntInt
  use geometry_header,  only: Cell, Universe, Lattice, Surface
  use tally_header,     only: TallyObject, TallyMap, TallyResult
  use mesh_header,      only: StructuredMesh

  implicit none

  type(Cell),      allocatable, target :: cells(:)
  type(Universe),  allocatable, target :: universes(:)
  type(Lattice),   allocatable, target :: lattices(:)
  type(Surface),   allocatable, target :: surfaces(:)
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict
  type(DictIntInt) :: surface_dict
  character(MAX_LINE_LEN) :: message
  character(3):: default_xs
  integer :: verbosity = 7
  integer :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials
  integer :: n_plots     ! # of plots
  integer :: n_nuclides_total ! Number of nuclide cross section tables
  integer :: n_sab_tables     ! Number of S(a,b) thermal scattering tables
  integer :: run_mode = NONE

  ! Flag for enabling cell overlap checking during transport
  logical                  :: check_overlaps = .false.
  integer(8), allocatable  :: overlap_check_cnt(:)

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  type(DictIntInt) :: mesh_dict
  type(DictIntInt) :: tally_dict
  type(DictCharInt) :: nuclide_dict

  logical :: cmfd_run = .false.

  type(StructuredMesh), allocatable, target :: meshes(:)
  type(TallyObject),    allocatable, target :: tallies(:)

  ! Pointers for different tallies
  type(TallyObject), pointer :: user_tallies(:) => null()
  type(TallyObject), pointer :: cmfd_tallies(:) => null()

  ! Starting index (minus 1) in tallies for each tally group
  integer :: i_user_tallies = -1
  integer :: i_cmfd_tallies = -1

  ! Global tallies
  !   1) collision estimate of k-eff
  !   2) track-length estimate of k-eff
  !   3) leakage fraction

  type(TallyResult), target :: global_tallies(N_GLOBAL_TALLIES)

  ! Tally map structure
  type(TallyMap), allocatable :: tally_maps(:)

  integer :: n_meshes       = 0 ! # of structured meshes
  integer :: n_user_meshes  = 0 ! # of structured user meshes
  integer :: n_tallies      = 0 ! # of tallies
  integer :: n_user_tallies = 0 ! # of user tallies

  ! Normalization for statistics
  integer :: n_realizations = 0 ! # of independent realizations
  real(8) :: total_weight       ! total starting particle weight in realization

  ! Flag for turning tallies on
  logical :: tallies_on = .false.
  logical :: active_batches = .false.

  ! Assume all tallies are spatially distinct
  logical :: assume_separate = .false.

  ! Use confidence intervals for results instead of standard deviations
  logical :: confidence_intervals = .false.

  integer :: n_cmfd_meshes              = 1 ! # of structured meshes
  integer :: n_cmfd_tallies             = 3 ! # of user-defined tallies
end module global
