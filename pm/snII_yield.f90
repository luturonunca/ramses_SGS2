module snII_yield
  use amr_commons, only: myid
  use amr_parameters, only: dp
  implicit none

  integer,parameter::nz_snII=7
  real(dp),dimension(nz_snII)::log_Zgrid_snII
  real(dp),dimension(nz_snII)::log_fFe_snII
  real(dp),dimension(nz_snII)::log_fMg_snII
  logical  :: snII_yield_ready = .false.

contains

  subroutine init_snII_yield()
    !----------------------------------------------------------------------------
    !  Limongi & Chieffi (2018) SNII yields with Prantzos et al. (2018) initial
    !  distribution of rotational velocities and Karakas (2010) for AGB stars,
    !  Chabrier (2005) IMF-averaged over [8,30] Msun (snyield_model=1 in
    !  Dubois et al. 2024, Sect. 2.4 / mechanical_fine.f90:SNII_yield).
    !  Mg is artificially doubled there to match observed [Mg/Si] (Limongi &
    !  Chieffi underproduce Mg relative to Si/solar).
    !  Trimmed to just the Fe and Mg mass fractions of the total SNII ejecta,
    !  tabulated on a metallicity grid and log-log interpolated; no other
    !  elements, no dust.
    !----------------------------------------------------------------------------
    if(snII_yield_ready) return

    log_Zgrid_snII = (/ &
     &-4.8712778d0,-3.8712778d0,-2.8712778d0,-2.4712777d0,-2.1712778d0,-1.8712777d0,-1.5712777d0 /)
    log_fFe_snII   = (/ &
     &-2.3678155d0,-2.3513636d0,-2.3487234d0,-2.3716861d0,-2.3260777d0,-2.2553820d0,-2.2553820d0 /)
    log_fMg_snII   = (/ &
     &-2.2677263d0,-2.2441224d0,-2.2427171d0,-2.2612530d0,-2.2412496d0,-2.2056181d0,-2.2056181d0 /)

    if(myid==1) then
       write(*,*) 'SNII Fe/Mg yields (Limongi & Chieffi 2018, Dubois et al. 2024):'
       write(*,*) '   Z grid (log10) = ',log_Zgrid_snII
       write(*,*) '   f_Fe (log10)   = ',log_fFe_snII
       write(*,*) '   f_Mg (log10)   = ',log_fMg_snII
    endif

    snII_yield_ready = .true.
  end subroutine init_snII_yield

  subroutine snII_yield_FeMg(zp_star, fFe_z, fMg_z)
    !----------------------------------------------------------------------------
    ! Returns the Fe and Mg mass fractions of the total SNII ejecta mass for a
    ! star of metallicity zp_star, by log-log interpolation on the tabulated
    ! grid above (no extrapolation beyond the grid edges).
    !----------------------------------------------------------------------------
    real(dp),intent(in) ::zp_star
    real(dp),intent(out)::fFe_z,fMg_z
    real(dp)::log_Zstar,fz
    integer ::izg

    call init_snII_yield()

    log_Zstar = log10(max(zp_star,1d-10))

    izg=1
    do while(izg<nz_snII-1 .and. log_Zstar>log_Zgrid_snII(izg+1))
       izg=izg+1
    end do

    fz = (log_Zgrid_snII(izg+1)-log_Zstar)/(log_Zgrid_snII(izg+1)-log_Zgrid_snII(izg))
    fz = min(max(fz,0d0),1d0)

    fFe_z = 10d0**(log_fFe_snII(izg)*fz + log_fFe_snII(izg+1)*(1d0-fz))
    fMg_z = 10d0**(log_fMg_snII(izg)*fz + log_fMg_snII(izg+1)*(1d0-fz))

  end subroutine snII_yield_FeMg

end module snII_yield
