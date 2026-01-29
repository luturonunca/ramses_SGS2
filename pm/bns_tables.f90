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

  subroutine bns_draw(mass, z, p_bns, m1, m2, vk1, vk2, t_sn2, t_merge)
    real(dp), intent(in) :: mass, z
    real(dp), intent(out) :: p_bns, m1, m2, vk1, vk2, t_sn2, t_merge
    if (.not. bns_tables_ready) call init_bns_tables()
    p_bns = 1.0d0
    m1 = 9.0d0
    m2 = 9.0d0
    vk1 = 200.0d0
    vk2 = 200.0d0
    t_sn2 = 5.0d0
    t_merge = 5.0d0
  end subroutine bns_draw

  real(dp) function bns_m2_of(mass, z)
    real(dp), intent(in) :: mass, z
    if (.not. bns_tables_ready) call init_bns_tables()
    bns_m2_of = 9.0d0
  end function bns_m2_of

end module bns_tables
