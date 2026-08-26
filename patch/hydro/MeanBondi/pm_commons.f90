module pm_commons

  use amr_parameters
  use pm_parameters
  use random

  implicit none

  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)    ::msink,xmsink
  real(dp),allocatable,dimension(:)    ::msink_new,msink_all
  real(dp),allocatable,dimension(:)    ::msmbh,msmbh_new,msmbh_all
  real(dp),allocatable,dimension(:)    ::oksink_new,oksink_all
  real(dp),allocatable,dimension(:)    ::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)    ::dMsink_overdt,dMBHoverdt
  real(dp),allocatable,dimension(:)    ::dMBHoverdt_fraction
  real(dp),allocatable,dimension(:)    ::dMsmbh_overdt,dMBHoverdt_smbh,dMBHoverdt_fraction_smbh
  real(dp),allocatable,dimension(:)    ::rho_gas,volume_gas,eps_sink,dMtorque_sink
  real(dp),allocatable,dimension(:,:)  ::vel_gas
  real(dp),allocatable,dimension(:)    ::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)    ::wden,weth,wvol,wdiv,wden_new,weth_new,wvol_new,wdiv_new
  real(dp),allocatable,dimension(:)    ::wfrac, wfvol, wfrac_new, wfvol_new
  real(dp),allocatable,dimension(:)    ::wv2, wc2, wv2_new, wc2_new
  real(dp),allocatable,dimension(:)    ::wvr2, wvphi2, wvr2_new, wvphi2_new
  real(dp),allocatable,dimension(:)    ::wcold_w, whot_w, wcold_w_new, whot_w_new
  real(dp),allocatable,dimension(:)    ::wcold_rho, whot_rho, wcold_rho_new, whot_rho_new
  real(dp),allocatable,dimension(:)    ::wcold_mass, wcold_mass_new
  real(dp),allocatable,dimension(:)    ::wtotal_mass, wtotal_mass_new
  real(dp),allocatable,dimension(:)    ::whot_cs2, whot_v2, whot_cs2_new, whot_v2_new
  real(dp),allocatable,dimension(:)    ::wcold_vphi2, wcold_cs2, wcold_vphi2_new, wcold_cs2_new
  real(dp),allocatable,dimension(:)    ::dMtorque2_sink, dMbondi2_sink
  real(dp),allocatable,dimension(:)    ::dMtorque_rot_sink, dMbondi_norot_sink, dMtorque_star_sink
  real(dp),allocatable,dimension(:)    ::wstar_mass, wstar_rot_mass
  real(dp),allocatable,dimension(:)    ::wstar_mass_new, wstar_rot_mass_new
  real(dp),allocatable,dimension(:)    ::sigma2_coll_sink
  real(dp),allocatable,dimension(:)    ::wsigma2_coll, wsigma2_coll_w
  real(dp),allocatable,dimension(:)    ::wsigma2_coll_new, wsigma2_coll_w_new
  real(dp),allocatable,dimension(:)    ::wrot_mass, wrot_mass_new
  real(dp),allocatable,dimension(:)    ::wdc_cold_mass, wdc_cold_mass_new
  real(dp),allocatable,dimension(:)    ::wdc_infall_mass, wdc_infall_mass_new ! same as wdc_cold_mass, additionally masked by S_inf_loc*S_rinf_loc when use_infall_mass
  real(dp),allocatable,dimension(:)    ::wdc_cold_j2mass, wdc_cold_j2mass_new
  real(dp),allocatable,dimension(:)    ::wdc_infall_j2mass, wdc_infall_j2mass_new ! same as wdc_cold_j2mass, additionally masked by S_inf_loc*S_rinf_loc when use_infall_mass, so eps_dc's j2_cold_mean matches the population M_cold_dc is actually drawn from instead of the unmasked one
  real(dp),allocatable,dimension(:)    ::wdc_tot_mass,  wdc_tot_mass_new
  real(dp),allocatable,dimension(:)    ::wdc_hot_w,     wdc_hot_w_new
  real(dp),allocatable,dimension(:)    ::wdc_hot_rho,   wdc_hot_rho_new
  real(dp),allocatable,dimension(:)    ::wdc_hot_cs2,   wdc_hot_cs2_new
  real(dp),allocatable,dimension(:)    ::wdc_hot_v2,    wdc_hot_v2_new
  real(dp),allocatable,dimension(:)    ::wdc_cold_w,    wdc_cold_w_new
  real(dp),allocatable,dimension(:)    ::wdc_cold_rho,  wdc_cold_rho_new
  real(dp),allocatable,dimension(:)    ::wff_part_mass, wff_part_mass_new ! star+DM mass in cloud for tff_include_particles
  ! Per-level storage for the freefall reservoir sums above: collect_acczone_avg(ilevel) is
  ! called once per level and resets/rebuilds wdc_*_new from scratch each time, so without a
  ! level dimension only the most-recently-processed level's contribution survives to
  ! compute_accretion_rate instead of the sum over the whole accretion zone (mirrors
  ! weighted_density(isink,ilevel) below).
  real(dp),allocatable,dimension(:,:)  ::wdc_cold_mass_lvl, wdc_cold_j2mass_lvl, wdc_tot_mass_lvl
  real(dp),allocatable,dimension(:,:)  ::wdc_infall_mass_lvl, wdc_infall_j2mass_lvl
  ! j2_crit from the previous sink update, per sink: the per-cell infall mask (S_inf_loc) needs
  ! j2_crit = G*(M_enc_dc+M_bh_dc)*R0_dc, but that is only known after compute_accretion_rate has
  ! summed wdc_*_lvl over every level -- strictly after the per-level collection loop that would
  ! build wdc_infall_mass has already run. Lagging by one sink update breaks that circularity,
  ! the same way c2sink/sigma2sink/r2sink below are computed once and then reused.
  real(dp),allocatable,dimension(:)    ::j2_crit_sink
  ! r_inf^2 = G*M_bh_dc/v_bondi_turb^2 from the previous sink update. Unlike j2_crit_sink, this
  ! one has no gas-mass term (point-mass only -- the standard sphere-of-influence definition,
  ! avoiding a feedback loop where a wider gate pulls in more enclosed/infall mass, which would
  ! then widen r_inf further), so v_bondi_turb/M_bh_dc are actually known before the per-cell
  ! collection loop runs. It is still lagged by one sink update purely because r2_inf_sink is
  ! only written inside compute_accretion_rate, which runs after this step's collection loop
  ! that needs the previous value -- a call-order constraint, not a mass circularity. v_bondi_turb
  ! adds the subgrid turbulent dispersion sigma2sink to v_bondi (c_s^2+v_rel^2+sigma_turb^2,
  ! Krumholz & McKee 2005), so this radius shrinks with local turbulence rather than only
  ! bulk/thermal motion. Used as a second,
  ! independent multiplicative mask (S_rinf_loc) restricting wdc_infall_mass to gas within the
  ! sink's actual gravitational influence radius, on top of the angular-momentum mask (S_inf_loc)
  ! -- the accretion-zone radius R0_dc=ir_cloud*dx_min is resolution-set and can be much larger
  ! than where the sink's gravity actually dominates the ambient gas motion.
  real(dp),allocatable,dimension(:)    ::r2_inf_sink
  real(dp),allocatable,dimension(:,:)  ::wdc_hot_w_lvl, wdc_hot_rho_lvl, wdc_hot_cs2_lvl, wdc_hot_v2_lvl
  real(dp),allocatable,dimension(:,:)  ::wff_part_mass_lvl
  ! Spatial hash (cell list) over sink positions, rebuilt each call to
  ! collect_acczone_avg, so collect_sigma_coll_np only checks sinks whose
  ! cloud-sized bin is adjacent to a given particle instead of every sink.
  integer,allocatable,dimension(:)     ::sink_hash_head, sink_hash_next
  integer,allocatable,dimension(:)     ::sink_bin_ix, sink_bin_iy, sink_bin_iz
  integer::nsink_hash=0
  real(dp),allocatable,dimension(:)    ::rho_cold_sink, rho_hot_sink
  real(dp),allocatable,dimension(:)    ::dMdc_cold_sink, dMdc_hot_sink
  real(dp),allocatable,dimension(:)    ::wnorot_w, wnorot_rho, wnorot_cs2, wnorot_v2
  real(dp),allocatable,dimension(:)    ::wnorot_w_new, wnorot_rho_new, wnorot_cs2_new, wnorot_v2_new
  real(dp),allocatable,dimension(:)    ::v2sink, c2sink, r2sink
  real(dp),allocatable,dimension(:)    ::wsigma2, sigma2sink,wsigma2_new
  real(dp),allocatable,dimension(:,:)  ::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)  ::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)  ::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)  ::lsink,lsink_new,lsink_all
  real(dp),allocatable,dimension(:,:)  ::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:)    ::graddescent_over_dt
  real(dp),allocatable,dimension(:,:)  ::xsink_graddescent
  real(dp),allocatable,dimension(:,:)  ::weighted_density,weighted_volume,weighted_ethermal,weighted_divergence
  real(dp), allocatable :: weighted_fraction(:,:)
  real(dp), allocatable :: weighted_fraction_weight(:,:)  
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:)    ::rho_sink_tff
  real(dp),allocatable,dimension(:)    ::msum_overlap
  integer,allocatable,dimension(:)     ::idsink,idsink_new,idsink_old,idsink_all
  logical,allocatable,dimension(:)     ::ok_blast_agn,ok_blast_agn_all
  logical,allocatable,dimension(:)     ::direct_force_sink
  logical,allocatable,dimension(:)     ::new_born,new_born_all,new_born_new
  integer,allocatable,dimension(:)     ::idsink_sort
  integer::ncloud_sink,ncloud_sink_massive
  integer::nindsink=0
  integer::sinkint_level=0         ! maximum level currently active is where the global sink variables are updated
  real(dp)::ssoft                  ! sink softening lenght in code units

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)  ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)  ::vp       ! Velocities
  real(dp),allocatable,dimension(:)    ::mp,mp0   ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)    ::ptcl_phi ! Potential of particle added by AP for output purposes
#endif
  real(dp),allocatable,dimension(:)    ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)    ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)    ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)    ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)    ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)    ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)    ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)    ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1

  ! Particle types
  integer, parameter   :: NFAMILIES=5
  integer(1),parameter :: FAM_DM=1, FAM_STAR=2, FAM_CLOUD=3, FAM_DEBRIS=4, FAM_OTHER=5, FAM_UNDEF=127
  integer(1),parameter :: FAM_TRACER_GAS=0
  integer(1),parameter :: FAM_TRACER_DM=-1, FAM_TRACER_STAR=-2, FAM_TRACER_CLOUD=-3, FAM_TRACER_DEBRIS=-4, FAM_TRACER_OTHER=-5

  ! Customize here for particle tags within particle types (e.g. different kind of stars).
  ! Note that the type should be integer(1) (1 byte integers) for memory concerns.
  ! Also don't forget to create a function is_<type>_<tag>. See the wiki for a more complete example.
  ! By default, the tag is always 0.

  ! Particle keys for outputing. They should match the above particle
  ! types, except for 'under' family
  character(len=13), dimension(-NFAMILIES:NFAMILIES), parameter :: particle_family_keys = (/ &
       ' other_tracer', 'debris_tracer', ' cloud_tracer', '  star_tracer', ' other_tracer', &
       '   gas_tracer', &
       '           DM', '         star', '        cloud', '       debris', '        other'/)

  type(part_t), allocatable, dimension(:) :: typep  ! Particle type array

contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  elemental logical pure function is_DM(typep)
    type(part_t), intent(in) :: typep
    is_DM = typep%family == FAM_DM
  end function is_DM

  elemental logical pure function is_star(typep)
    type(part_t), intent(in) :: typep
    is_star = typep%family == FAM_STAR
  end function is_star

  elemental logical pure function is_cloud(typep)
    type(part_t), intent(in) :: typep
    is_cloud = typep%family == FAM_CLOUD
  end function is_cloud

  elemental logical pure function is_debris(typep)
    type(part_t), intent(in) :: typep
    is_debris = typep%family == FAM_DEBRIS
  end function is_debris

  elemental logical pure function is_tracer(typep)
    type(part_t), intent(in) :: typep
    is_tracer = typep%family <= 0
  end function is_tracer

  elemental logical pure function is_not_tracer(typep)
    type(part_t), intent(in) :: typep
    is_not_tracer = typep%family > 0
  end function is_not_tracer

  elemental logical pure function is_not_DM(typep)
    type(part_t), intent(in) :: typep
    is_not_DM = typep%family /= FAM_DM
  end function is_not_DM
  

  elemental function part2int (part)
    ! Convert a particle into an integer
    ! This saves some space e.g. when communicating
    integer :: part2int
    type(part_t), intent(in) :: part

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    part2int = (int(part%family) + a) * b + (int(part%tag) + a)
  end function part2int

  elemental function int2part(index)
    ! Convert from an index to particle type
    type(part_t) :: int2part
    integer, intent(in) :: index

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    int2part%family = int(index / b - a, 1)
    int2part%tag = int(mod(index, b) - a, 1)
  end function int2part

  function props2type(idpii, tpii, mpii)
    use amr_commons
    use pm_parameters, only : part_t

    ! Converts from "old" ramses to "new" ramses
    !
    ! Here's the match, add yours here for backward compatibility purposes
    ! DM     tpii == 0
    ! stars  tpii != 0 and idpii > 0
    ! sinks  tpii != 0 and idpii < 0
    !
    ! This is mostly for support of GRAFFIC I/O.
    ! The reason we use idpii instead of idp is to prevent name clashes
    real(dp), intent(in) :: tpii, mpii
    integer, intent(in)  :: idpii

    type(part_t) :: props2type

    if (tpii == 0) then
       props2type%family = FAM_DM
    else if (idpii > 0) then
       props2type%family = FAM_STAR
    else if (idpii < 0) then
       props2type%family = FAM_CLOUD
    else if (mpii == 0) then
       props2type%family = FAM_TRACER_GAS
    end if
    props2type%tag = 0
  end function props2type
end module pm_commons
