module snIa_yield
  use amr_commons, only: myid
  use amr_parameters, only: dp
  implicit none

  real(dp) :: mejecta_Ia
  real(dp) :: f_Fe_Ia
  real(dp) :: f_Mg_Ia
  logical  :: snIa_yield_ready = .false.

contains

  subroutine init_snIa_yield()
    !----------------------------------------------------------------------------
    !  Iwamoto et al. (1999) W70 (carbon-deflagration model)
    !            (updated from Nomoto et al. 1997)
    !            (better fit with Tycho observation)
    !  Trimmed to just the total ejecta mass and the Fe/Mg mass fractions
    !  needed by snIa_fine; no C/N/O/Si/S, no dust.
    !----------------------------------------------------------------------------
    !       12C ,     13C ,   14N ,    15N ,    16O ,    17O ,    18O ,    19F ,
    !       20Ne,     21Ne,   22Ne,    23Na,    24Mg,    25Mg,    26Mg,    27Al,
    !       28Si,     29Si,   30Si,    31P ,    32S ,    33S ,    34S ,    36S ,
    !       35Cl,     37Cl,   36Ar,    38Ar,    40Ar,    39K ,    41K ,    40Ca,
    !       42Ca,     43Ca,   44Ca,    46Ca,    48Ca,    45Sc,    46Ti,    47Ti,
    !       48Ti,     49Ti,   50Ti,    50V ,    51V ,    50Cr,    52Cr,    53Cr,
    !       54Cr,     55Mn,   54Fe,    56Fe,    57Fe,    58Fe,    59Co,    58Ni,
    !       60Ni,     61Ni,   62Ni,    64Ni,    63Cu,    65Cu,    64Zn,    66Zn,
    !       67Zn,     68Zn
    !----------------------------------------------------------------------------
    real(dp) :: yield_snIa(1:66)

    if(snIa_yield_ready) return

    yield_snIa = (/ &  ! Msun per SN
     &5.08E-02,1.56E-09,3.31E-08,4.13E-07,1.33E-01,3.33E-10,2.69E-10,1.37E-10,&
     &2.29E-03,2.81E-08,2.15E-08,1.41E-05,1.58E-02,1.64E-07,1.87E-07,1.13E-04,&
     &1.42E-01,5.79E-05,7.12E-05,9.12E-05,9.14E-02,6.07E-05,1.74E-05,3.41E-11,&
     &1.06E-05,5.56E-06,1.91E-02,6.60E-07,3.42E-12,1.67E-06,4.83E-07,1.81E-02,&
     &1.06E-08,6.17E-08,1.38E-05,1.01E-09,2.47E-09,3.85E-08,3.49E-07,4.08E-07,&
     &3.13E-04,2.94E-06,1.04E-04,1.22E-08,4.27E-05,6.65E-05,7.73E-03,5.66E-04,&
     &9.04E-04,6.66E-03,7.30E-02,6.80E-01,1.92E-02,2.96E-03,9.68E-04,8.34E-02,&
     &1.47E-02,2.15E-04,1.85E-03,1.65E-05,3.00E-06,8.33E-07,7.01E-05,6.26E-06,&
     &7.28E-09,1.13E-08/)

    mejecta_Ia = sum(yield_snIa)
    f_Fe_Ia    = sum(yield_snIa(51:54)) / mejecta_Ia  ! 54Fe,56Fe,57Fe,58Fe
    f_Mg_Ia    = sum(yield_snIa(13:15)) / mejecta_Ia  ! 24Mg,25Mg,26Mg

    if(myid==1) then
       write(*,*) 'SNIa yield (Iwamoto et al. 1999, W70):'
       write(*,*) '   mejecta_Ia = ',mejecta_Ia,' Msun'
       write(*,*) '   f_Fe_Ia    = ',f_Fe_Ia
       write(*,*) '   f_Mg_Ia    = ',f_Mg_Ia
    endif

    snIa_yield_ready = .true.
  end subroutine init_snIa_yield

end module snIa_yield
