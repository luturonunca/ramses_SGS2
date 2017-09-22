subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,idim,nboundary_true=0
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::ek_bound

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
#if NENER>0
       & ,prad_region &
#endif
#if NVAR>NDIM+2+NENER
       & ,var_region &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
#if NENER>0
       & ,gamma_rad &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type,sink_refine
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
#if NENER>0
       & ,prad_bound &
#endif
#if NVAR>NDIM+2+NENER
       & ,var_bound &
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow
  namelist/physics_params/omega_b,cooling,haardt_madau,metal,isothermal &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star,jeans_ncells &
       & ,eta_sn,eta_ssn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc,kappa_IR &
       & ,J21,a_spec,z_ave,z_reion,ind_rsink,delayed_cooling,T2max &
       & ,self_shielding,smbh,agn,momentum_feedback &
       & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff,t_diss,t_sne &
       & ,sf_virial,sf_trelax,sf_tdiss,sf_model,sf_log_properties,sf_imf &
       & ,mass_star_max,mass_sne_min,sf_compressive
#ifdef grackle
   namelist/grackle_params/use_grackle,grackle_with_radiative_cooling,grackle_primordial_chemistry,grackle_metal_cooling &
       & ,grackle_UVbackground,grackle_cmb_temperature_floor,grackle_h2_on_dust,grackle_photoelectric_heating &
       & ,grackle_use_volumetric_heating_rate,grackle_use_specific_heating_rate,grackle_three_body_rate,grackle_cie_cooling &
       & ,grackle_h2_optical_depth_approximation,grackle_ih2co,grackle_ipiht,grackle_NumberOfTemperatureBins,grackle_CaseBRecombination &
       & ,grackle_Compton_xray_heating,grackle_LWbackground_sawtooth_suppression,grackle_NumberOfDustTemperatureBins,grackle_use_radiative_transfer &
       & ,grackle_radiative_transfer_coupled_rate_solver,grackle_radiative_transfer_intermediate_step,grackle_radiative_transfer_hydrogen_only &
       & ,grackle_self_shielding_method,grackle_Gamma,grackle_photoelectric_heating_rate,grackle_HydrogenFractionByMass &
       & ,grackle_DeuteriumToHydrogenRatio,grackle_SolarMetalFractionByMass,grackle_TemperatureStart,grackle_TemperatureEnd &
       & ,grackle_DustTemperatureStart,grackle_DustTemperatureEnd,grackle_LWbackground_intensity,grackle_UVbackground_redshift_on &
       & ,grackle_UVbackground_redshift_off,grackle_UVbackground_redshift_fullon,grackle_UVbackground_redshift_drop &
       & ,grackle_cloudy_electron_fraction_factor,grackle_data_file
#endif

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue
#ifdef grackle
  rewind(1)
  read(1,NML=grackle_params)
#endif
#ifdef ATON
  if(aton)call read_radiation_params(1)
#endif

  !--------------------------------------------------
  ! Check for dm only cosmo run
  !--------------------------------------------------
  if(.not.hydro)then
     omega_b = 0.0D0
  endif

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(t_star>0)then
     star=.true.
     pic=.true.
  else if(eps_star>0)then
     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
  if(metal.and.nvar<(ndim+3))then
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+3'
     if(myid==1)write(*,*)'Modify hydro_parameters.f90 and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
  if(nvar<(ndim+2+nener))then
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check ind_rsink
  !--------------------------------------------------
  if(ind_rsink<=0.0d0)then
     if(myid==1)write(*,*)'Error in the namelist'
     if(myid==1)write(*,*)'Check ind_rsink'
     nml_ok=.false.
  end if

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
  inener=ndim+3 ! MUST BE THIS VALUE !!!
  imetal=nener+ndim+3
  idelay=imetal
  if(metal)idelay=imetal+1
  ivirial1=idelay
  ivirial2=idelay
  if(delayed_cooling)then
     ivirial1=idelay+1
     ivirial2=idelay+1
  endif
  if(sf_virial)then
     if(sf_compressive) ivirial2=ivirial1+1
  endif
  ixion=ivirial2
  if(sf_virial)ixion=ivirial2+1
  ichem=ixion
  if(aton)ichem=ixion+1
  if(myid==1.and.hydro.and.(nvar>ndim+2)) then
     write(*,'(A50)')"__________________________________________________"
     write(*,*) 'Hydro var indices:'
#if NENER>0
     write(*,*) '   inener   = ',inener
#endif
     if(metal)           write(*,*) '   imetal   = ',imetal
     if(delayed_cooling) write(*,*) '   idelay   = ',idelay
     if(sf_virial)then
        write(*,*) '   ivirial1 = ',ivirial1
        if(sf_compressive) write(*,*) '   ivirial2 = ',ivirial2
     endif
     if(aton)            write(*,*) '   ixion    = ',ixion
#ifdef RT
     if(rt) write(*,*) '   iIons    = ',ichem
#endif
     write(*,'(A50)')"__________________________________________________"
  endif

  ! Last variable is ichem

end subroutine read_hydro_params

