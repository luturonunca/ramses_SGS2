module bns_tables
  use amr_commons, only: myid
  use amr_parameters, only: dp, bns_efficiency_table
  implicit none

  real(dp), allocatable, dimension(:) :: bns_prob
  real(dp), allocatable, dimension(:) :: bns_z
  logical :: bns_tables_ready = .false.
  integer, parameter :: bns_table_lu = 17

contains

  subroutine init_bns_tables()
    integer :: nrows, i, ios
    real(dp) :: zval, pval
    logical :: ok
    character(len=256) :: line
    if (bns_tables_ready) return
    inquire(file=trim(bns_efficiency_table), exist=ok)
    if(.not. ok)then
       if(myid==1)write(*,*)'BNS efficiency table not found: ', trim(bns_efficiency_table)
       call clean_stop
    endif
    open(unit=bns_table_lu,file=trim(bns_efficiency_table),status='old',form='formatted')
    nrows = 0
    do
       read(bns_table_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       if(line(1:1) == '#') cycle
       nrows = nrows + 1
    end do
    if(nrows <= 0)then
       if(myid==1)write(*,*)'BNS efficiency table is empty: ', trim(bns_efficiency_table)
       call clean_stop
    endif
    rewind(bns_table_lu)
    allocate(bns_prob(nrows))
    allocate(bns_z(nrows))
    i = 0
    do
       read(bns_table_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       if(line(1:1) == '#') cycle
       read(line,*,iostat=ios) zval, pval
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad line in BNS efficiency table: ', trim(line)
          call clean_stop
       endif
       i = i + 1
       bns_z(i) = zval
       bns_prob(i) = pval
    end do
    close(bns_table_lu)
    bns_tables_ready = .true.
  end subroutine init_bns_tables

  subroutine bns_draw(mass, z, p_bns, m1, m2, vk1, vk2, t_sn2, t_merge)
    real(dp), intent(in) :: mass, z
    real(dp), intent(out) :: p_bns, m1, m2, vk1, vk2, t_sn2, t_merge
    real(dp) :: p_per_msun
    if (.not. bns_tables_ready) call init_bns_tables()
    call bns_prob_of_z(z, p_per_msun)
    p_bns = max(0.0d0, min(1.0d0, p_per_msun * mass))
    m1 = 9.0d0
    m2 = 9.0d0
    vk1 = 200.0d0
    vk2 = 200.0d0
    t_sn2 = 30.0d0
    t_merge = 20.0d0
  end subroutine bns_draw

  real(dp) function bns_m2_of(mass, z)
    real(dp), intent(in) :: mass, z
    if (.not. bns_tables_ready) call init_bns_tables()
    bns_m2_of = 9.0d0
  end function bns_m2_of

  subroutine bns_prob_of_z(z, p_per_msun)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: p_per_msun
    integer :: i, n
    real(dp) :: w
    n = size(bns_z)
    if(z <= bns_z(1))then
       p_per_msun = bns_prob(1)
    else if(z >= bns_z(n))then
       p_per_msun = bns_prob(n)
    else
       p_per_msun = bns_prob(n)
       do i = 1, n-1
          if(z >= bns_z(i) .and. z < bns_z(i+1))then
             w = (z - bns_z(i)) / (bns_z(i+1) - bns_z(i))
             p_per_msun = (1.0d0 - w) * bns_prob(i) + w * bns_prob(i+1)
             exit
          endif
       end do
    endif
  end subroutine bns_prob_of_z

end module bns_tables
