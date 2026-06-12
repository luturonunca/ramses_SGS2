module bns_tables
  use amr_commons, only: myid
  use amr_parameters, only: dp, bns_table_dir, fixed_bns_vars, &
       & bns_t_sn2, bns_t_merge, bns_v_kick1, bns_v_kick2, bns_efficiency, bns_merger_frac
  use pm_commons, only: localseed
  use random, only: ranf
  implicit none

  real(dp), allocatable, dimension(:) :: bns_prob
  real(dp), allocatable, dimension(:) :: bns_frac_merge
  real(dp), allocatable, dimension(:) :: bns_z
  character(len=32), allocatable, dimension(:) :: bns_z_label
  logical :: bns_tables_ready = .false.
  integer, parameter :: bns_table_lu = 17
  integer, parameter :: bns_hist_lu = 18
  integer :: current_z_index = -1
  real(dp), allocatable, dimension(:) :: vk1_bin_min, vk1_bin_max, vk1_cdf
  real(dp), allocatable, dimension(:) :: vk2_bin_min, vk2_bin_max, vk2_cdf
  real(dp), allocatable, dimension(:) :: tsn2_bin_min, tsn2_bin_max, tsn2_cdf
  real(dp), allocatable, dimension(:) :: tmerge_bin_min, tmerge_bin_max, tmerge_cdf

contains

  subroutine init_bns_tables()
    integer :: nrows, i, ios
    real(dp) :: zval, pval, fval
    logical :: ok
    character(len=256) :: line, table_path
    character(len=32) :: ztxt, ptxt, ftxt
    integer :: dir_len
    if (bns_tables_ready) return
    dir_len = len_trim(bns_table_dir)
    if(dir_len <= 0)then
       if(myid==1)write(*,*)'BNS table directory is empty'
       call clean_stop
    endif
    if(bns_table_dir(dir_len:dir_len) == '/')then
       table_path = trim(bns_table_dir)//'bns_efficiency_Z.dat'
    else
       table_path = trim(bns_table_dir)//'/bns_efficiency_Z.dat'
    endif
    inquire(file=trim(table_path), exist=ok)
    if(.not. ok)then
       if(myid==1)write(*,*)'BNS efficiency table not found: ', trim(table_path)
       call clean_stop
    endif
    open(unit=bns_table_lu,file=trim(table_path),status='old',form='formatted')
    nrows = 0
    do
       read(bns_table_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       if(line(1:1) == '#') cycle
       nrows = nrows + 1
    end do
    if(nrows <= 0)then
       if(myid==1)write(*,*)'BNS efficiency table is empty: ', trim(table_path)
       call clean_stop
    endif
    rewind(bns_table_lu)
    allocate(bns_prob(nrows))
    allocate(bns_frac_merge(nrows))
    allocate(bns_z(nrows))
    allocate(bns_z_label(nrows))
    i = 0
    do
       read(bns_table_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       if(line(1:1) == '#') cycle
       read(line,*,iostat=ios) ztxt, ptxt, ftxt
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad line in BNS efficiency table: ', trim(line)
          call clean_stop
       endif
       read(ztxt,*,iostat=ios) zval
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad metallicity in BNS efficiency table: ', trim(line)
          call clean_stop
       endif
       read(ptxt,*,iostat=ios) pval
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad efficiency in BNS efficiency table: ', trim(line)
          call clean_stop
       endif
       read(ftxt,*,iostat=ios) fval
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad merger fraction in BNS efficiency table: ', trim(line)
          call clean_stop
       endif
       i = i + 1
       bns_z(i) = zval
       bns_prob(i) = pval
       bns_frac_merge(i) = fval
       bns_z_label(i) = ztxt
    end do
    close(bns_table_lu)
    bns_tables_ready = .true.
  end subroutine init_bns_tables

  subroutine bns_draw(mass, z, p_bns, m1, m2, vk1, vk2, t_sn2, t_merge)
    real(dp), intent(in) :: mass, z
    real(dp), intent(out) :: p_bns, m1, m2, vk1, vk2, t_sn2, t_merge
    real(dp) :: p_per_msun, frac_merge, u
    integer :: iz
    if(fixed_bns_vars) then
       p_bns   = max(0.0d0, min(1.0d0, bns_efficiency * mass))
       vk1     = bns_v_kick1
       vk2     = bns_v_kick2
       t_sn2   = bns_t_sn2
       call ranf(localseed, u)
       if(u <= bns_merger_frac) then
          t_merge = bns_t_merge
       else
          t_merge = 1.0d30
       end if
       m1 = 9.0d0
       m2 = 9.0d0
       return
    end if
    if (.not. bns_tables_ready) call init_bns_tables()
    call bns_prob_of_z(z, p_per_msun)
    p_bns = max(0.0d0, min(1.0d0, p_per_msun * mass))
    call bns_find_z_index(z, iz)
    call bns_load_histograms(iz)
    call draw_from_hist(vk1_bin_min, vk1_bin_max, vk1_cdf, vk1)
    call draw_from_hist(vk2_bin_min, vk2_bin_max, vk2_cdf, vk2)
    call draw_from_hist(tsn2_bin_min, tsn2_bin_max, tsn2_cdf, t_sn2)
    call bns_merge_frac_of_z(z, frac_merge)
    call ranf(localseed, u)
    if(u <= frac_merge) then
       call draw_from_hist(tmerge_bin_min, tmerge_bin_max, tmerge_cdf, t_merge)
    else
       t_merge = 1.0d30
    end if
    m1 = 9.0d0
    m2 = 9.0d0
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

  subroutine bns_merge_frac_of_z(z, frac)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: frac
    integer :: i, n
    real(dp) :: w
    n = size(bns_z)
    if(z <= bns_z(1))then
       frac = bns_frac_merge(1)
    else if(z >= bns_z(n))then
       frac = bns_frac_merge(n)
    else
       frac = bns_frac_merge(n)
       do i = 1, n-1
          if(z >= bns_z(i) .and. z < bns_z(i+1))then
             w = (z - bns_z(i)) / (bns_z(i+1) - bns_z(i))
             frac = (1.0d0 - w) * bns_frac_merge(i) + w * bns_frac_merge(i+1)
             exit
          endif
       end do
    endif
  end subroutine bns_merge_frac_of_z

  subroutine bns_find_z_index(z, iz)
    real(dp), intent(in) :: z
    integer, intent(out) :: iz
    integer :: i, n
    n = size(bns_z)
    if(z <= bns_z(1))then
       iz = 1
    else if(z >= bns_z(n))then
       iz = n
    else
       iz = n
       do i = 1, n-1
          if(z >= bns_z(i) .and. z < bns_z(i+1))then
             iz = i
             exit
          endif
       end do
    endif
  end subroutine bns_find_z_index

  subroutine bns_load_histograms(iz)
    integer, intent(in) :: iz
    character(len=256) :: base_dir, path
    integer :: dir_len
    if(iz == current_z_index) return
    if(allocated(vk1_bin_min)) deallocate(vk1_bin_min, vk1_bin_max, vk1_cdf)
    if(allocated(vk2_bin_min)) deallocate(vk2_bin_min, vk2_bin_max, vk2_cdf)
    if(allocated(tsn2_bin_min)) deallocate(tsn2_bin_min, tsn2_bin_max, tsn2_cdf)
    if(allocated(tmerge_bin_min)) deallocate(tmerge_bin_min, tmerge_bin_max, tmerge_cdf)
    dir_len = len_trim(bns_table_dir)
    if(dir_len <= 0)then
       if(myid==1)write(*,*)'BNS table directory is empty'
       call clean_stop
    endif
    if(bns_table_dir(dir_len:dir_len) == '/')then
       base_dir = trim(bns_table_dir)//'chabrier_Z_'//trim(bns_z_label(iz))
    else
       base_dir = trim(bns_table_dir)//'/chabrier_Z_'//trim(bns_z_label(iz))
    endif
    path = trim(base_dir)//'/hist_v_kick_sn1.csv'
    call read_hist_file(path, vk1_bin_min, vk1_bin_max, vk1_cdf)
    path = trim(base_dir)//'/hist_v_kick_sn2.csv'
    call read_hist_file(path, vk2_bin_min, vk2_bin_max, vk2_cdf)
    path = trim(base_dir)//'/hist_t_sn_2_log_yr.csv'
    call read_hist_file(path, tsn2_bin_min, tsn2_bin_max, tsn2_cdf, 1.0d-6)
    path = trim(base_dir)//'/hist_t_merge_log_yr.csv'
    call read_hist_file(path, tmerge_bin_min, tmerge_bin_max, tmerge_cdf, 1.0d-6)
    current_z_index = iz
  end subroutine bns_load_histograms

  subroutine read_hist_file(filename, bin_min, bin_max, cdf, scale)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: bin_min(:), bin_max(:), cdf(:)
    real(dp), intent(in), optional :: scale
    integer :: nrows, i, ios
    real(dp) :: vmin, vmax, prob, sum_prob
    character(len=256) :: line
    logical :: ok
    real(dp) :: scale_val
    inquire(file=trim(filename), exist=ok)
    if(.not. ok)then
       if(myid==1)write(*,*)'BNS histogram file not found: ', trim(filename)
       call clean_stop
    endif
    open(unit=bns_hist_lu,file=trim(filename),status='old',form='formatted')
    read(bns_hist_lu,'(A)',iostat=ios) line
    if(ios /= 0)then
       if(myid==1)write(*,*)'Cannot read BNS histogram header: ', trim(filename)
       call clean_stop
    endif
    nrows = 0
    do
       read(bns_hist_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       nrows = nrows + 1
    end do
    if(nrows <= 0)then
       if(myid==1)write(*,*)'BNS histogram file is empty: ', trim(filename)
       call clean_stop
    endif
    rewind(bns_hist_lu)
    read(bns_hist_lu,'(A)',iostat=ios) line
    allocate(bin_min(nrows), bin_max(nrows), cdf(nrows))
    sum_prob = 0.0d0
    scale_val = 1.0d0
    if(present(scale)) scale_val = scale
    i = 0
    do
       read(bns_hist_lu,'(A)',iostat=ios) line
       if(ios /= 0) exit
       if(len_trim(line) == 0) cycle
       read(line,*,iostat=ios) vmin, vmax, prob
       if(ios /= 0)then
          if(myid==1)write(*,*)'Bad line in BNS histogram file: ', trim(line)
          call clean_stop
       endif
       i = i + 1
       bin_min(i) = vmin * scale_val
       bin_max(i) = vmax * scale_val
       cdf(i) = prob
       sum_prob = sum_prob + prob
    end do
    close(bns_hist_lu)
    if(sum_prob <= 0.0d0)then
       if(myid==1)write(*,*)'BNS histogram has non-positive sum: ', trim(filename)
       call clean_stop
    endif
    do i = 1, nrows
       cdf(i) = cdf(i) / sum_prob
       if(i > 1) cdf(i) = cdf(i) + cdf(i-1)
    end do
    cdf(nrows) = 1.0d0
  end subroutine read_hist_file

  subroutine draw_from_hist(bin_min, bin_max, cdf, value)
    real(dp), intent(in) :: bin_min(:), bin_max(:), cdf(:)
    real(dp), intent(out) :: value
    real(dp) :: u
    integer :: i, n
    n = size(cdf)
    call ranf(localseed, u)
    do i = 1, n
       if(u <= cdf(i))then
          call ranf(localseed, u)
          value = bin_min(i) + (bin_max(i) - bin_min(i)) * u
          return
       endif
    end do
    value = bin_max(n)
  end subroutine draw_from_hist

end module bns_tables
