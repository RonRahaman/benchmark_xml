module global
  use constants
  use dict_header,      only: DictCharInt, DictIntInt
  use geometry_header,  only: Cell, Universe, Lattice, Surface

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
end module global
