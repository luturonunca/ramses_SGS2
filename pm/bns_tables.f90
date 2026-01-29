module bns_tables
  use amr_parameters, only: dp
  implicit none

  integer, parameter :: nbns_m = 1
  integer, parameter :: nbns_z = 1
  real(dp), allocatable, dimension(:,:) :: bns_prob
  real(dp), allocatable, dimension(:) :: bns_mass
  real(dp), allocatable, dimension(:) :: bns_z
  logical :: bns_tables_ready = .false.

contains

  subroutine init_bns_tables()
    if (bns_tables_ready) return
    allocate(bns_prob(nbns_m, nbns_z))
    allocate(bns_mass(nbns_m))
    allocate(bns_z(nbns_z))
    bns_prob = 1.0d0
    bns_mass = 0.0d0
    bns_z = 0.0d0
    bns_tables_ready = .true.
  end subroutine init_bns_tables

  real(dp) function bns_prob_of(mass, z)
    real(dp), intent(in) :: mass, z
    if (.not. bns_tables_ready) call init_bns_tables()
    bns_prob_of = 1.0d0
  end function bns_prob_of

end module bns_tables
