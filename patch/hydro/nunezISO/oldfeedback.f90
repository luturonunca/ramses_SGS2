!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose) write(*,111)ilevel

  ! Gather star particles only.

#if NDIM==3
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).lt.0.and.tp(ipart).ne.0d0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(idp(ipart).lt.0.and.tp(ipart).ne.0d0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

#endif

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each debris particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! uold.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::t0,ESN,mejecta,zloss
  real(dp)::dx,dx_loc,scale,vol_loc,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  !integer,dimension(1:nvector),save::igrid_son,ind_son
  !integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Massive star lifetime from Myr to code units
  t0=t_sne*1d6*(365.*24.*3600.)/scale_t

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(M_SNII*2d33)/scale_v**2*eff_sfbk

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do

  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)


  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sn2'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=floor(x(j,idim),kind=4)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Compute individual time steps                                               
  do j=1,np
     if(ok(j))then
        if(levelp(ind_part(j))>=ilevel)then
           dteff(j)=dtnew(levelp(ind_part(j)))
        else
           dteff(j)=dtold(levelp(ind_part(j)))
        end if
     endif
  end do

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     if(ok(j))then
        mloss(j)=0d0
        mzloss(j)=0d0
        ethermal(j)=0d0
     endif
  end do

  do j=1,np
     if(ok(j))then
        birth_time=tp(ind_part(j))
        ! Make sure that we don't count feedback twice
        if( idp(ind_part(j)).le.0 .and. tp(ind_part(j)).ne.0d0 .and. &
             & tp(ind_part(j)).lt.(t-t0) )then
           if(use_initial_mass)then
              mejecta=eta_sn*mp0(ind_part(j))
           else
              mejecta=eta_sn*mp(ind_part(j))
           endif
           mloss(j)=mloss(j)+mejecta/vol_loc
           ! Thermal energy
           ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc
           ! Metallicity
           if(metal)then
              zloss=yield+(1d0-yield)*zp(ind_part(j))
              mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc
           endif
           ! Reduce star particle mass
           mp(ind_part(j))=mp(ind_part(j))-mejecta
           idp(ind_part(j))=-idp(ind_part(j))
        endif
     end if
  end do

  ! Update hydro variables due to feedback
  do j=1,np
     if(ok(j))then
        ! Specific kinetic energy of the star
        ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
             &          +vp(ind_part(j),2)**2 &
             &          +vp(ind_part(j),3)**2)
        ! Update hydro variable in NGP cell
        uold(indp(j),1)=uold(indp(j),1)+mloss(j)
        uold(indp(j),2)=uold(indp(j),2)+mloss(j)*vp(ind_part(j),1)
        uold(indp(j),3)=uold(indp(j),3)+mloss(j)*vp(ind_part(j),2)
        uold(indp(j),4)=uold(indp(j),4)+mloss(j)*vp(ind_part(j),3)
        uold(indp(j),5)=uold(indp(j),5)+mloss(j)*ekinetic(j)+ethermal(j)
     endif
  end do
  if(metal)then
     do j=1,np
        if(ok(j))then
           uold(indp(j),imetal)=uold(indp(j),imetal)+mzloss(j)
        endif
     end do
  endif

#endif
  
end subroutine feedbk
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine checks SN events in cells where a
  ! star particle has been spawned.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_tot,info,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  integer ,dimension(:),allocatable::indSN,iSN_myid
  real(dp),dimension(:),allocatable::mSN,ZSN,vol_gas,ekBlast
  real(dp),dimension(:),allocatable::mloadSN,ZloadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,vloadSN
  integer ,dimension(:),allocatable::itemp
  real(dp),dimension(:),allocatable::mSN_tot,ZSN_tot
  real(dp),dimension(:,:),allocatable::xSN_tot,vSN_tot
  integer::isort

  if(.not. hydro)return
  if(ndim.ne.3)return

  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Time delay for SN explosion from Myr to code units
  t0=t_sne*(1d6*365.*24.*3600.)/scale_t

  !------------------------------------------------------
  ! Gather star particles eligible for a SN event
  !------------------------------------------------------

  nSN_tot=0
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0        
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        nSN_tot=nSN_tot+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  enddo

  nSN_icpu=0
  nSN_icpu(myid)=nSN_tot  
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif
  nSN_tot=sum(nSN_icpu(1:ncpu))

  if(myid==1.and.nSN_tot>0)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  if (nSN_tot .eq. 0) return
 
  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN_tot(1:nSN_tot,1:3),vSN_tot(1:nSN_tot,1:3),mSN_tot(1:nSN_tot))
  allocate(ZSN_tot(1:nSN_tot),itemp(1:nSN_tot))

  xSN_tot=0.;vSN_tot=0.;mSN_tot=0.;ZSN_tot=0.
  !------------------------------------------------------
  ! Give position and mass of the star to the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid        
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then
                 iSN=iSN+1
                 xSN_tot(iSN,1)=xp(ipart,1)
                 xSN_tot(iSN,2)=xp(ipart,2)
                 xSN_tot(iSN,3)=xp(ipart,3)
                 vSN_tot(iSN,1)=vp(ipart,1)
                 vSN_tot(iSN,2)=vp(ipart,2)
                 vSN_tot(iSN,3)=vp(ipart,3)
                 if(use_initial_mass)then
                    mSN_tot(iSN)  =eta_sn*mp0(ipart)
                 else
                    mSN_tot(iSN)  =eta_sn*mp(ipart)
                 endif
                 if(metal)ZSN_tot(iSN)=zp(ipart)
                 ! Remove the mass ejected by the SN
                 mp(ipart) =mp(ipart)-eta_sn*mp(ipart)
                 idp(ipart)=-idp(ipart)
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
     
        igrid=next(igrid)   ! Go to next grid
     end do
  enddo

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN_tot,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN_tot,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN_tot,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN_tot,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN_tot=xSN_all
  vSN_tot=vSN_all
  mSN_tot=mSN_all
  ZSN_tot=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,ZSN_all)
#endif

  call getSNonmyid(itemp,nSN,xSN_tot,nSN_tot)

  ! Allocate the arrays for the position and the mass of the SN
  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN))
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0

  do iSN=1,nSN
     isort=itemp(iSN)
     iSN_myid(iSN)=isort
     xSN(iSN,1)=xSN_tot(isort,1)
     xSN(iSN,2)=xSN_tot(isort,2)
     xSN(iSN,3)=xSN_tot(isort,3)
     vSN(iSN,1)=vSN_tot(isort,1)
     vSN(iSN,2)=vSN_tot(isort,2)
     vSN(iSN,3)=vSN_tot(isort,3)
     mSN(iSN)  =mSN_tot(isort)
     ZSN(iSN)  =ZSN_tot(isort)
  enddo
  deallocate(xSN_tot,vSN_tot,mSN_tot,ZSN_tot,itemp)

  allocate(vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN),indSN(1:nSN))
  allocate(mloadSN(1:nSN),ZloadSN(1:nSN),vloadSN(1:nSN,1:3))

  ! Compute the grid discretization effects
  call average_SN(xSN,vSN,vol_gas,dq,ekBlast,indSN,nSN,nSN_tot,&
                &iSN_myid,mSN,mloadSN,ZSN,ZloadSN,vloadSN,.true.)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,vloadSN)

  deallocate(xSN,vSN,mSN,ZSN,iSN_myid)
  deallocate(indSN,vol_gas,dq,ekBlast)
  deallocate(mloadSN,ZloadSN,vloadSN)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vSN,vol_gas,dq,ekBlast,ind_blast,nSN,nSN_tot,iSN_myid,&
                    &mSN,mloadSN,ZSN,ZloadSN,vloadSN,global_search)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  ! and do the mass loading process
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,nSN_tot,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,d,u,v,w,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::eint,ekk,mload,Zload
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::mSN,mloadSN,ZSN,ZloadSN,vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq,u2Blast,vloadSN
  real(dp)::pxmi,pxma,pymi,pyma,pzmi,pzma,rlocmax
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN_tot)::vol_gas_mpi,vol_gas_all
  real(dp),dimension(1:nSN_tot)::mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all
  real(dp),dimension(1:nSN_tot,1:3)::dq_mpi,u2Blast_mpi,vloadSN_mpi
  real(dp),dimension(1:nSN_tot,1:3)::dq_all,u2Blast_all,vloadSN_all
#endif
  logical ,dimension(1:nvector),save::ok
  integer ,dimension(1:nSN)::iSN_myid
  integer::ind_SN
  logical::global_search

  if(verbose)write(*,*)'Entering average_SN'
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
 
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! To reduce time spent to search for cells
  pxmi=minval(xSN(:,1))
  pxma=maxval(xSN(:,1))
  pymi=minval(xSN(:,2))
  pyma=maxval(xSN(:,2))
  pzmi=minval(xSN(:,3))
  pzma=maxval(xSN(:,3))

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1;mloadSN=0.0;ZloadSN=0.0;vloadSN=0.0

  ! Loop over levels - determine dq first
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim

     ! compute the rough estimate of the region where stars are located 
     rlocmax=MAX(dx_loc/2d0,rmax)

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rlocmax).or.(x.gt.pxma+rlocmax).or.&
                   &(y.lt.pymi-rlocmax).or.(y.gt.pyma+rlocmax).or.&
                   &(z.lt.pzmi-rlocmax).or.(z.gt.pzma+rlocmax)) then 
                    ok(i)=.false.
                 endif
              endif
           enddo

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale   !
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale   !
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale   !
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels


#ifndef WITHOUTMPI
  if(global_search)then
     vol_gas_mpi=0d0; dq_mpi=0d0
     vol_gas_all=0d0; dq_all=0d0
     
     ! Put the nSN size arrays into nSN_tot size arrays to synchronize processors
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        vol_gas_mpi(ind_SN)=vol_gas(iSN)
        dq_mpi     (ind_SN,1)=dq     (iSN,1)
        dq_mpi     (ind_SN,2)=dq     (iSN,2)
        dq_mpi     (ind_SN,3)=dq     (iSN,3)
     enddo
     call MPI_ALLREDUCE(vol_gas_mpi,vol_gas_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dq_mpi     ,dq_all     ,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     
     vol_gas_mpi=vol_gas_all
     dq_mpi     =dq_all
     
     ! Put the nSN_tot size arrays into nSN size arrays
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        vol_gas(iSN)=vol_gas_mpi(ind_SN)
        dq     (iSN,1)=dq_mpi     (ind_SN,1)
        dq     (iSN,2)=dq_mpi     (ind_SN,2)
        dq     (iSN,3)=dq_mpi     (ind_SN,3)
     enddo
  endif
#endif

  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
     endif
  enddo

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim

     ! compute the rough estimate of the region where stars are located 
     rlocmax=MAX(dx_loc/2d0,rmax)

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rlocmax).or.(x.gt.pxma+rlocmax).or.&
                   &(y.lt.pymi-rlocmax).or.(y.gt.pyma+rlocmax).or.&
                   &(z.lt.pzmi-rlocmax).or.(z.gt.pzma+rlocmax)) then 
                    ok(i)=.false.
                 endif
              endif
           enddo

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale   !
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale   !
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale   !
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax - dq(iSN,1)
                       v=dyy/rmax - dq(iSN,2)
                       w=dzz/rmax - dq(iSN,3)
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                       d=uold(ind_blast(iSN),1)
                       u=uold(ind_blast(iSN),2)/d
                       v=uold(ind_blast(iSN),3)/d
                       w=uold(ind_blast(iSN),4)/d
                       ekk=0.5d0*d*(u*u+v*v+w*w)
                       eint=uold(ind_blast(iSN),5)-ekk
                       ! Mass loading factor of the Sedov explosion
                       ! Ensure that no more that 90% of the gas content is removed
                       mload=min(f_w*mSN(iSN),0.90d0*d*vol_loc)
                       mloadSN(iSN)=mSN(iSN)+mload


                       ! Update gas mass and metal content in the cell
                       if(metal)then
                          Zload=uold(ind_blast(iSN),imetal)/d
                          ZloadSN(iSN)=( mload*Zload + yield*(1-zSN(iSN))*mSN(iSN) + zSN(iSN)*mSN(iSN)) / mloadSN(iSN)
                          uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal)-Zload*mload/vol_loc
                       endif
                       d=uold(ind_blast(iSN),1)-mload/vol_loc

                       uold(ind_blast(iSN),1)=d
                       uold(ind_blast(iSN),2)=d*u
                       uold(ind_blast(iSN),3)=d*v
                       uold(ind_blast(iSN),4)=d*w
                       uold(ind_blast(iSN),5)=eint+0.5d0*d*(u*u+v*v+w*w)

                       vloadSN(iSN,1)=(mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(iSN)
                       vloadSN(iSN,2)=(mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(iSN)
                       vloadSN(iSN,3)=(mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels


#ifndef WITHOUTMPI
  if(global_search)then
     u2Blast_mpi=0d0; mloadSN_mpi=0d0; ZloadSN_mpi=0d0; vloadSN_mpi=0d0
     ! Put the nSN size arrays into nSN_tot size arrays to synchronize processors
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        mloadSN_mpi(ind_SN)=mloadSN(iSN)
        ZloadSN_mpi(ind_SN)=ZloadSN(iSN)
        vloadSN_mpi(ind_SN,1)=vloadSN(iSN,1)
        vloadSN_mpi(ind_SN,2)=vloadSN(iSN,2)
        vloadSN_mpi(ind_SN,3)=vloadSN(iSN,3)
        u2Blast_mpi(ind_SN,1)=u2Blast(iSN,1)
        u2Blast_mpi(ind_SN,2)=u2Blast(iSN,2)
        u2Blast_mpi(ind_SN,3)=u2Blast(iSN,3)
     enddo
     call MPI_ALLREDUCE(mloadSN_mpi,mloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ZloadSN_mpi,ZloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vloadSN_mpi,vloadSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(u2Blast_mpi,u2Blast_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     mloadSN_mpi=mloadSN_all
     ZloadSN_mpi=ZloadSN_all
     vloadSN_mpi=vloadSN_all
     u2Blast_mpi=u2Blast_all
     ! Put the nSN_tot size arrays into nSN size arrays
     do iSN=1,nSN
        ind_SN=iSN_myid(iSN)
        mloadSN(iSN)=mloadSN_mpi(ind_SN)
        ZloadSN(iSN)=ZloadSN_mpi(ind_SN)
        vloadSN(iSN,1)=vloadSN_mpi(ind_SN,1)
        vloadSN(iSN,2)=vloadSN_mpi(ind_SN,2)
        vloadSN(iSN,3)=vloadSN_mpi(ind_SN,3)
        u2Blast(iSN,1)=u2Blast_mpi(ind_SN,1)
        u2Blast(iSN,2)=u2Blast_mpi(ind_SN,2)
        u2Blast(iSN,3)=u2Blast_mpi(ind_SN,3)
     enddo
  endif
#endif


  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)
        v2=u2Blast(iSN,2)
        w2=u2Blast(iSN,3)
        ekBlast(iSN)=0.5d0*(u2+v2+w2)
        if (ekBlast(iSN).le.0.01) ekBlast(iSN)=0d0
      endif
  end do
  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,d_gas,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,p_gas,vol_gas,uSedov,ekBlast,mloadSN,ZloadSN
  real(dp),dimension(1:nSN,1:3)::xSN,dq,vloadSN
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Ejecta specific energy (accounting for dilution)
  ESN=(1d51/(M_SNII*2d33))/scale_v**2*eff_sfbk

  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        d_gas=mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas (iSN)=d_gas*ESN
           uSedov(iSN)=0d0
        else
           p_gas (iSN)=(1d0-f_ek)*d_gas*ESN
           uSedov(iSN)=sqrt(f_ek*mSN(iSN)*ESN/ekBlast(iSN)/mloadSN(iSN))
        endif
     else
        p_gas(iSN)=mSN(iSN)*ESN/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       d_gas=mloadSN(iSN)/vol_gas(iSN)
                       ! Compute the density and the metal density of the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_gas*ZloadSN(iSN)

                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vloadSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vloadSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vloadSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        d_gas=mloadSN(iSN)/ekBlast(iSN)
        u=vloadSN(iSN,1)
        v=vloadSN(iSN,2)
        w=vloadSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_gas*ZloadSN(iSN)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getSNonmyid(iSN_myid,nSN_myid,xSN,nSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid
  integer::nSN_myid,ii,nSN
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,iSN,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:nSN,1:3)::xSN 

 
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
  do iSN=1,nSN

     cpu_read=.false.
     ! Compute boundaries for the SN cube of influence
     xxmin=(xSN(iSN,1)-drSN)/scale ; xxmax=(xSN(iSN,1)+drSN)/scale
     yymin=(xSN(iSN,2)-drSN)/scale ; yymax=(xSN(iSN,2)+drSN)/scale
     zzmi=(xSN(iSN,3)-drSN)/scale ; zzmax=(xSN(iSN,3)+drSN)/scale

     if(TRIM(ordering).eq.'hilbert')then
        
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
      (bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif
        
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+oneqdp)*dkey
        end do
        
        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if
     
     ! Create the index array for SN in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           iSN_myid(ii)=iSN
        endif
     enddo

  enddo
  
  ! Number of SN in processor myid
  nSN_myid=ii

end subroutine getSNonmyid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine mechanical_feedback_cell
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(1:ncpu)::nSNcomm_icpu_mpi
  real(dp),dimension(:),allocatable::mSN_mpi,ZSN_mpi
  real(dp),dimension(:,:),allocatable::xSN_mpi,vSN_mpi
  integer, dimension(:),allocatable::lSN_mpi
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine gathers SN events on a cell-by-cell basis
  ! Here 'local' variable means things that do not require MPI comm
  ! Taysun Kimm
  !----------------------------------------------------------------------
  ! local constants
  integer::igrid,npart1,ipart,jpart,next_part
  integer::nSNp,nSNp_tot,nSNc,nSNc_tot,nSNcomm,nSNcomm_tot,nSN_glo,nSN_loc
  integer::ilevel,ind,ivar,nx_loc,iSN,nsn_star,nSN
  integer,dimension(1:ncpu)::nSNcomm_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun,t0
  real(dp)::scale,dx,dx_loc,skip_loc(3),current_time
  integer ,dimension(:),allocatable::iSN_myid
  real(dp),dimension(:),allocatable::mSN_loc,ZSN_loc
  real(dp),dimension(:,:),allocatable::xSN_loc,vSN_loc
  real(dp),dimension(:),allocatable::mSN,ZSN
  real(dp),dimension(:,:),allocatable::xSN,vSN
  integer ,dimension(:),allocatable::itemp,lSN_loc,lSN_glo,lSN
  real(dp),dimension(:),allocatable::mSN_glo,ZSN_glo
  real(dp),dimension(:,:),allocatable::xSN_glo,vSN_glo
  real(dp),dimension(1:twotondim,1:3)::xc,p_ej
  real(dp),dimension(1:twotondim)::m_ej,Z_ej
  real(dp),dimension(1:ndim)::x0,x
  real(dp)::ttsta,ttend,mejecta
  integer::ind_grid,ind_cell,ngrid,ind_son,idim,iSN_loc,iSN_glo
  integer::isort,iskip,i,ix,iy,iz,info,nok_sn,ncpu_read
  logical::ok,done_star

  if(.not. hydro)return
  if(ndim.ne.3)return

#ifdef RT
  if(static) return
#endif

#ifndef WITHOUTMPI
  ttsta = MPI_WTIME(info)
#endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = (boxlen*scale_l)**3*scale_d/1.989d33 

  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_sne*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_sne*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif


  ! update the particle lists on ilevel > levelmin
  do ilevel=levelmin,nlevelmax
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! gather particle from the grid
  do ilevel=nlevelmax-1,1,-1
     call merge_tree_fine(ilevel)
  end do

  ! Scatter particle to the grid
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! gather particle
  do ilevel=nlevelmax,levelmin,-1
     call merge_tree_fine(ilevel)
  end do

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)


  !------------------------------------------------------
  ! Count cells eligible for a SN event
  !------------------------------------------------------
  nSNc     = 0  ! the total number of cells that will launch SN in myid
  nSNp     = 0  ! the total number of SN particles in myid
  nSNcomm  = 0  ! the total number of cells that need communications
  nSNc_tot = 0  ! total(nSN) for the entire cpus
  nSNp_tot = 0  ! total(nSNp) for the entire cpus
  do ilevel=levelmin,nlevelmax
     dx=0.5D0**ilevel
     dx_loc=dx*scale

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ngrid=active(ilevel)%ngrid
     do igrid=1,ngrid
        ind_grid=active(ilevel)%igrid(igrid)
        npart1=numbp(ind_grid) ! number of particles in the grid
        ok=.true.
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(ind_grid,idim) - dx
           end do

           m_ej=0d0 
           ipart=headp(ind_grid)
           do jpart=1,npart1
              next_part=nextp(ipart)
              ok=.false.
              if(sn2_real_delay)then
                 !if tp is younger than t_sne
                 if(idp(ipart).le.0.and.tp(ipart).ge.(current_time-t0))then
                    call get_number_of_sn2  (tp(ipart), zp(ipart), &
                                   & mp0(ipart)*scale_msun,mp(ipart)*scale_msun,nsn_star,done_star)
                    if(nsn_star>0)ok=.true.
                 endif
              else ! single SN event
                 !if tp is older than t_sne
                 if(idp(ipart).le.0.and.tp(ipart).le.(current_time-t0))then
                    ok=.true.
                 endif
              endif

              if(ok)then
                 ! find the cell index to get the position of it
                 ind_son=1  ! double-checked
                 do idim=1,ndim
                    i = int((xp(ipart,idim)/scale+skip_loc(idim)-x0(idim))/dx)
                    ind_son=ind_son+i*2**(idim-1)
                 end do
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+ind_grid
                 if(son(ind_cell)==0)then  ! leaf cell
                    nSNp=nSNp+1
                    m_ej(ind_son)=1.
                 endif
              endif

              ipart=next_part
           enddo

           nok_sn = count(m_ej>0d0)
           nSNc = nSNc + nok_sn ! add the number of cells that contain SN(e)

           if (nok_sn>0) then ! check how many SNe need MPI comm
              do ind=1,twotondim
                 if (m_ej(ind)>0.5) then
                    x(1) = (xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale 
                    x(2) = (xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale 
                    x(3) = (xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale 
                    call checkSNboundary(x,ncpu_read,ilevel)
                    if(ncpu_read>1) nSNcomm = nSNcomm + 1
                 endif
              enddo
           endif


        endif
     enddo
  enddo
  
  nSNp_tot=nSNp
  nSNc_tot=nSNc
  nSNcomm_tot=0 
  nSNcomm_icpu=0  ! This array is for cells that need MPI communications
  nSNcomm_icpu(myid)=nSNcomm
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  nSNcomm_icpu_mpi=0
  nSNp_tot=0
  call MPI_ALLREDUCE(nSNc    ,nSNc_tot    ,1,   MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nSNp    ,nSNp_tot    ,1,   MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nSNcomm_icpu,nSNcomm_icpu_mpi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSNcomm_icpu=nSNcomm_icpu_mpi
#endif
  nSNcomm_tot=sum(nSNcomm_icpu(1:ncpu))  ! the total number of cells that contain SNe for all cpus

  if(myid==1.and.nSNc_tot>0)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN/SNc/SNcomm to explode=',nSNp_tot,nSNc_tot,nSNcomm_tot
     write(*,*)'-----------------------------------------------'
  endif

  if (nSNc_tot .eq. 0) return


  ! local SN cells except the ones that will go to "global". This wouldn't require MPI
  nSN_loc = nSNc - nSNcomm
  nSN_glo = nSNcomm_tot

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN_loc(1:nSN_loc,1:3),vSN_loc(1:nSN_loc,1:3))
  allocate(mSN_loc(1:nSN_loc),ZSN_loc(1:nSN_loc),lSN_loc(1:nSN_loc))
  xSN_loc=0.;vSN_loc=0.;mSN_loc=0.;ZSN_loc=0.; lSN_loc=0

  allocate(xSN_glo(1:nSN_glo,1:3),vSN_glo(1:nSN_glo,1:3))
  allocate(mSN_glo(1:nSN_glo),ZSN_glo(1:nSN_glo),lSN_glo(1:nSN_glo))
  xSN_glo=0.;vSN_glo=0.;mSN_glo=0.;ZSN_glo=0.; lSN_glo=0

  iSN_loc=0
  if(myid==1) then
     iSN_glo=0
  else
     iSN_glo=sum(nSNcomm_icpu(1:myid-1))
  endif


  !---------------------------------------------------------------------
  ! Gather ejecta m, p, Z, etc., for local/global SN events, separately
  ! here 'global' means that it requires cells belonging to other cpus 
  !---------------------------------------------------------------------
  do ilevel=levelmin,nlevelmax
     dx=0.5D0**ilevel
     dx_loc=dx*scale

     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ngrid=active(ilevel)%ngrid
     do igrid=1,ngrid
        ind_grid=active(ilevel)%igrid(igrid)
        npart1=numbp(ind_grid) ! number of particles in the grid
        ok=.true.
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(ind_grid,idim) - dx
           end do

           m_ej=0d0; p_ej=0d0; Z_ej=0d0
           ipart=headp(ind_grid)
           do jpart=1,npart1
              next_part=nextp(ipart)
              ok=.false.
              if(sn2_real_delay)then
                 !if tp is younger than t_sne
                 if(idp(ipart).le.0.and.tp(ipart).ge.(current_time-t0))then
                    call get_number_of_sn2  (tp(ipart), zp(ipart), &
                                   & mp0(ipart)*scale_msun,mp(ipart)*scale_msun,nsn_star,done_star)
                    if(nsn_star>0)ok=.true.
                 endif
              else ! single SN event
                 !if tp is older than t_sne
                 if(idp(ipart).le.0.and.tp(ipart).le.(current_time-t0))then
                    ok=.true.
                 endif
              endif

              if(ok)then
                 ! find the cell index to get the position of it
                 ind_son=1  ! double-checked
                 do idim=1,ndim
                    i = int((xp(ipart,idim)/scale+skip_loc(idim)-x0(idim))/dx)
                    ind_son=ind_son+i*2**(idim-1)
                 end do
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+ind_grid
                 if(son(ind_cell)==0)then
                    nSNp=nSNp+1
                    if(sn2_real_delay)then
                       mejecta = M_SNII/scale_msun*nsn_star
                    else
                       if(use_initial_mass)then
                          mejecta=eta_sn*mp0(ipart)
                       else
                          mejecta=eta_sn*mp(ipart)
                       endif
                    endif                         
                    m_ej(ind_son) = m_ej(ind_son) + mejecta
                    p_ej(ind_son,1) = p_ej(ind_son,1) + mejecta*vp(ipart,1)
                    p_ej(ind_son,2) = p_ej(ind_son,2) + mejecta*vp(ipart,2)
                    p_ej(ind_son,3) = p_ej(ind_son,3) + mejecta*vp(ipart,3)
                    if(metal) Z_ej(ind_son) = Z_ej(ind_son) + mejecta*zp(ipart)
                    ! Remove the mass ejected by the SN
                    mp(ipart)  = mp(ipart) - mejecta
                    if(.not.sn2_real_delay) idp(ipart) = -idp(ipart)
                 endif
              endif
              ipart=next_part
           enddo

           nok_sn = count(m_ej>0d0)
           nSNc = nSNc + nok_sn ! add the number of cells that contain SN(e)

           if (nok_sn>0) then ! check how many SNe need MPI comm
              do ind=1,twotondim
                 if (m_ej(ind)>0d0) then
                    x(1) = (xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale 
                    x(2) = (xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale 
                    x(3) = (xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale 
                    call checkSNboundary(x,ncpu_read,ilevel)
                    if(ncpu_read>1)then  ! SNe that requires cells in other cpus
                       iSN_glo = iSN_glo + 1
                       xSN_glo(iSN_glo,1) = x(1)
                       xSN_glo(iSN_glo,2) = x(2)
                       xSN_glo(iSN_glo,3) = x(3)
                       mSN_glo(iSN_glo  ) = m_ej(ind)
                       vSN_glo(iSN_glo,1) = p_ej(ind,1)/m_ej(ind)
                       vSN_glo(iSN_glo,2) = p_ej(ind,2)/m_ej(ind)
                       vSN_glo(iSN_glo,3) = p_ej(ind,3)/m_ej(ind)
                       if(metal)ZSN_glo(iSN_glo)=Z_ej(ind)/m_ej(ind)
                       lSN_glo(iSN_glo  ) = ilevel
                    else                 ! SNe that are happy with cells in myid
                       iSN_loc = iSN_loc + 1
                       xSN_loc(iSN_loc,1) = x(1)
                       xSN_loc(iSN_loc,2) = x(2)
                       xSN_loc(iSN_loc,3) = x(3)
                       mSN_loc(iSN_loc  ) = m_ej(ind)
                       vSN_loc(iSN_loc,1) = p_ej(ind,1)/m_ej(ind)
                       vSN_loc(iSN_loc,2) = p_ej(ind,2)/m_ej(ind)
                       vSN_loc(iSN_loc,3) = p_ej(ind,3)/m_ej(ind)
                       if(metal)ZSN_loc(iSN_loc)=Z_ej(ind)/m_ej(ind)
                       lSN_loc(iSN_loc  ) = ilevel
                    endif
                 endif
              enddo
           endif

        endif
     enddo
  enddo
  


#ifndef WITHOUTMPI
  !---------------------------------------------------------------------
  !  Part-1 release: global SNe
  !---------------------------------------------------------------------
  if (nSN_glo>0) then

     allocate(xSN_mpi(1:nSN_glo,3),vSN_mpi(1:nSN_glo,3))
     allocate(mSN_mpi(1:nSN_glo),ZSN_mpi(1:nSN_glo),lSN_mpi(1:nSN_glo))
     xSN_mpi=0d0;vSN_mpi=0d0;mSN_mpi=0d0;ZSN_mpi=0d0;lSN_mpi=0
     call MPI_ALLREDUCE(xSN_glo,xSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vSN_glo,vSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mSN_glo,mSN_mpi,nSN_glo  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(lSN_glo,lSN_mpi,nSN_glo  ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
     if(metal)call MPI_ALLREDUCE(ZSN_glo,ZSN_mpi,nSN_glo,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     xSN_glo=xSN_mpi
     vSN_glo=vSN_mpi
     mSN_glo=mSN_mpi
     lSN_glo=lSN_mpi
     if(metal)ZSN_glo=ZSN_mpi
     deallocate(xSN_mpi,vSN_mpi,mSN_mpi,ZSN_mpi,lSN_mpi)

     allocate(itemp(1:nSN_glo))
     itemp=0
     call getSNonmyid2(itemp,nSN,xSN_glo,nSN_glo,lSN_glo)

     allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN),lSN(1:nSN))
     xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0; lSN=0

     do iSN=1,nSN
        isort = itemp(iSN)
        iSN_myid(iSN) = isort
        xSN(iSN,1) = xSN_glo(isort,1)
        xSN(iSN,2) = xSN_glo(isort,2)
        xSN(iSN,3) = xSN_glo(isort,3)
        vSN(iSN,1) = vSN_glo(isort,1)
        vSN(iSN,2) = vSN_glo(isort,2)
        vSN(iSN,3) = vSN_glo(isort,3)
        mSN(iSN  ) = mSN_glo(isort)
        ZSN(iSN  ) = ZSN_glo(isort)
        lSN(iSN  ) = lSN_glo(isort)
     end do

     ! Inject momentum
     call inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_glo,.true.)

     deallocate(xSN,vSN,mSN,ZSN,iSN_myid,itemp,lSN)

  endif
#endif
  deallocate(xSN_glo,vSN_glo,mSN_glo,ZSN_glo,lSN_glo)

  !---------------------------------------------------------------------
  !  Part-2 release: local SNe
  !---------------------------------------------------------------------
  allocate(itemp(1:nSN_loc))
  itemp=0
  call getSNonmyid2(itemp,nSN,xSN_loc,nSN_loc,lSN_loc)

  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN),lSN(1:nSN))
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0; lSN=0

  do iSN=1,nSN
     isort = itemp(iSN)
     iSN_myid(iSN) = isort
     xSN(iSN,1) = xSN_loc(isort,1)
     xSN(iSN,2) = xSN_loc(isort,2)
     xSN(iSN,3) = xSN_loc(isort,3)
     vSN(iSN,1) = vSN_loc(isort,1)
     vSN(iSN,2) = vSN_loc(isort,2)
     vSN(iSN,3) = vSN_loc(isort,3)
     mSN(iSN  ) = mSN_loc(isort)
     ZSN(iSN  ) = ZSN_loc(isort)
     lSN(iSN  ) = lSN_loc(isort)
  end do

  ! Inject momentum
  call inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_loc,.false.)

  deallocate(xSN,vSN,mSN,ZSN,iSN_myid,itemp,lSN)
  deallocate(xSN_loc,vSN_loc,mSN_loc,ZSN_loc)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo


#ifndef WITHOUTMPI
  ttend = MPI_WTIME(info)
  if(myid.eq.1)then
     write(*,*) 'Time elapsed in mechanical_feedback_cell [s]', sngl(ttend-ttsta)
  endif
#endif

end subroutine mechanical_feedback_cell
!################################################################
!################################################################
!################################################################
!################################################################
subroutine inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_glo,global_search)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(dp),dimension(1:nSN_glo,1:3)::vloadSN_mpi
  real(dp),dimension(1:nSN_glo,1:3)::dq_mpi
  real(dp),dimension(1:nSN_glo)::mloadSN_mpi,ZloadSN_mpi
#endif
  integer,parameter::n_neigh = 48 ! fixed ; DO NOT CHANGE
  integer::ilevel,iSN,nSN,nSN_glo,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache,ind_SN,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,drr,f_adjacent
  real(dp)::ekk,eint,eint2,u,v,w,d,vol_nbor,dload,pload(3),dloadloc
  real(dp)::scale,dx_min,dx_loc,vol_loc,rloose2,dx_sn,dd_sn
  real(dp)::scale_nh,scale_t2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp)::pxmi,pxma,pymi,pyma,pzmi,pzma,rmaxloc,dr_cell
  real(dp)::mload,vload,Zload,ESN,f_w_cell,navg,Tk,Tk0,Tshock
  real(dp)::f_esn,f_esn2,f_w_crit,num_sn,Zdepen,f_leftover,f_load
  real(dp),dimension(1:3)::skip_loc,x0
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,ZSN
  real(dp),dimension(1:nSN,1:3)::xSN,vSN
  integer, dimension(1:nSN)::lSN,ind_blast
  real(dp),dimension(1:nSN_glo)::mloadSN,ZloadSN
  real(dp),dimension(1:nSN_glo,1:3)::vloadSN,dq
  integer, dimension(1:nSN_glo)::iSN_myid
  logical::global_search
  logical,dimension(1:nvector),save::ok
!--------------------------------------------------------------
  integer::ista,iend,jsta,jend,ksta,kend,n48,icello
  integer::nm_sncell
  real(dp)::dm_ejecta
  real(dp),dimension(1:n_neigh,1:nSN)::rho_nbor,Z_nbor ! rho for neigboring cells
  integer ,dimension(1:n_neigh,1:nSN)::icell_nbor ! icell to speed up
  integer ,dimension(1:n_neigh,1:nSN)::lv_nbor
  integer ,dimension(1:n_neigh)::idx48
  real(dp),dimension(1:n_neigh,1:3)::v48
  real(dp)::u0,v0,w0
!--------------------------------------------------------------
! TKNOTE
!--------------------------------------------------------------
! There are recurrent numbers like '48' '56'
! The structure one can suppose is 4*4*4 of which the central 2*2*2 is a SN cell.
! The total number of cells in this structure is 64, but we decide to neglect 
! 8 vertices, making the total number to be 56. On top of that, 
! what we inject momentum is only for neibouring cells, making the relevant cells 
! to be 48. Note that we deposit nm_sncell/56 of the SN mass + mass in the SN cell. This is desirable 
! for numerical stability. If there are many stars with different ages in a cell,
! it will redistribute mass very efficiently, in which case the cell is gonna run 
! out of gas. 
 
  nm_sncell=4 ! mass fraction to be left over
  f_leftover=dble(nm_sncell)/dble(n_neigh+nm_sncell)  ! M_leftover in the SN cell = f_leftover*(m_cell(iSN) + m_ejecta)
  f_load    =1d0 - f_leftover
  
  ! mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**nlevelmax
  f_adjacent = 1.7 ! the longest distance=1.658 if we neglect 8 edge cells  

  ! conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  scale_msun=(scale_d*scale_l**3d0)/2d33 

  ! In the original code, M_SN=10Msun is used. If you want to stick to this,
  ! eff_sfbk=0.657510, eta_sn=0.31675331 (Chabrier)
  ! This does not mean that the actual efficiency is 66%....
  ! Here I change M_SNII to make things more intuitive.
  ESN=1d51*eff_sfbk
  f_esn  = 0.676 ! Blondin+(98) at t=trad (applies only for high f_w_cell)

  ! To reduce time spent to search for cells
  pxmi=minval(xSN(:,1))
  pxma=maxval(xSN(:,1))
  pymi=minval(xSN(:,2))
  pyma=maxval(xSN(:,2))
  pzmi=minval(xSN(:,3))
  pzma=maxval(xSN(:,3))

  ind_blast = -1

  mloadSN =0d0; ZloadSN=0d0; vloadSN=0d0; rho_nbor=0d0; Z_nbor=0d0; icell_nbor=-1


  !---------------------------------------------------------------------
  ! Gather rho,v,Z,icell for SN + adjacent cell
  !---------------------------------------------------------------------
  ! loop over levels
  do ilevel=levelmin,nlevelmax
     ! computing local volume (important for averaging hydro quantities) 
     dx=0.5d0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     rmaxloc = 2*dx_loc 
     rloose2 = (f_adjacent*dx_loc*2)**2d0

     ! cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx
        xc(ind,2)=(dble(iy)-0.5d0)*dx
        xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do

     ! loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=min(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rmaxloc).or.(x.gt.pxma+rmaxloc).or.&
                   &(y.lt.pymi-rmaxloc).or.(y.gt.pyma+rmaxloc).or.&
                   &(z.lt.pzmi-rmaxloc).or.(z.gt.pzma+rmaxloc)) then
                    ok(i)=.false.
                 endif
              endif
           enddo
 
           do i=1,ngrid
              if(ok(i))then
                 ! get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ind_SN = iSN_myid(iSN)
                    dxx = x - xSN(iSN,1)
                    dyy = y - xSN(iSN,2)
                    dzz = z - xSN(iSN,3)
                    drr = dxx*dxx + dyy*dyy + dzz*dzz
                    if(drr<=rloose2*2)then
                       dx_sn = 0.5d0**lSN(iSN)*scale    ! cell size at which iSN exist
                       dd_sn = (dx_loc+dx_sn)/2d0 + dx_sn/1d5 ! maximum separation at this level
                       dr_cell = max(abs(dxx),abs(dyy),abs(dzz))
                       if(dr_cell<dd_sn*2) then !if adjacent
                          x0(1) = xSN(iSN,1) - dx_sn ! left-bottom pos of the SN region (not cell)
                          x0(2) = xSN(iSN,2) - dx_sn
                          x0(3) = xSN(iSN,3) - dx_sn
                          
                          ista  = idnint(((x-dx_loc/2d0) - x0(1))/dx_sn*2d0)
                          iend  = idnint(((x+dx_loc/2d0) - x0(1))/dx_sn*2d0) -1
                          jsta  = idnint(((y-dx_loc/2d0) - x0(2))/dx_sn*2d0)
                          jend  = idnint(((y+dx_loc/2d0) - x0(2))/dx_sn*2d0) -1
                          ksta  = idnint(((z-dx_loc/2d0) - x0(3))/dx_sn*2d0)
                          kend  = idnint(((z+dx_loc/2d0) - x0(3))/dx_sn*2d0) -1

                          ! calculate special indices designed for momentum
                          ! injection. if this cell is twice larger than 
                          ! the SN  cell, the number of indices for this cell is
                          ! gonna be 4.
                          call get_idx48(ista,iend,jsta,jend,ksta,kend,idx48,n48)

                          ! remember surrounding densities
                          rho_nbor(idx48(1:n48),iSN)   =uold(ind_cell(i),1) 
                          if(metal)then
                             Z_nbor(idx48(1:n48),iSN)  =uold(ind_cell(i),imetal)/uold(ind_cell(i),1)
                          endif
                          icell_nbor(idx48(1:n48),iSN) =ind_cell(i)
                          lv_nbor(idx48(1:n48),iSN)    =ilevel
                          if(dr_cell<dx_loc/2d0)then  ! if the SN cell, not adjacent cells 
                             ind_blast(iSN)=ind_cell(i)
                             d = uold(ind_blast(iSN),1)
                             mload = f_load*d*vol_loc !min(f_w*mSN(iSN),0.90d0*d*vol_loc)
                             mloadSN(ind_SN) = mSN(iSN)*f_load + mload
                             if(metal)then
                                Zload = uold(ind_blast(iSN),imetal)/d
                                ZloadSN(ind_SN)=(mload*Zload + (yield+(1d0-yield)*ZSN(iSN))*mSN(iSN)*f_load)/mloadSN(ind_SN)
                             endif
                             u = uold(ind_blast(iSN),2)/d
                             v = uold(ind_blast(iSN),3)/d
                             w = uold(ind_blast(iSN),4)/d
                             vloadSN(ind_SN,1)=(f_load*mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(ind_SN)
                             vloadSN(ind_SN,2)=(f_load*mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(ind_SN)
                             vloadSN(ind_SN,3)=(f_load*mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(ind_SN)
                          endif
                       endif
                        
                    end if 
                 end do
              end if
           end do

        end do ! loop over cells
     end do ! loop over grids
  end do ! loop over levels


#ifndef WITHOUTMPI
  if(global_search)then
     mloadSN_mpi=0d0; vloadSN_mpi=0d0; ZloadSN_mpi=0d0
     call MPI_ALLREDUCE(mloadSN,mloadSN_mpi,nSN_glo  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vloadSN,vloadSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     if(metal)call MPI_ALLREDUCE(ZloadSN,ZloadSN_mpi,nSN_glo,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     mloadSN = mloadSN_mpi
     vloadSN = vloadSN_mpi
     ZloadSN = ZloadSN_mpi
  endif
#endif


  !---------------------------------------------------------------------
  ! Calculate the local mass loading and enforce momentum conservation
  !---------------------------------------------------------------------
  dq=0d0 ; v48=0d0
  call get_v48(v48)

  do iSN=1,nSN
     ind_SN = iSN_myid(iSN)
     dx=0.5d0**lSN(iSN)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim   
     mload=mloadSN(ind_SN)
     if(metal)Zload=ZloadSN(ind_SN)


     do i=1,n_neigh
        if(icell_nbor(i,iSN)>0)then
           ! idx48 is measured for sub volume of vol_loc/8
           dm_ejecta=f_load*mSN(iSN)/dble(n_neigh)
           f_w_cell = (mload/dble(n_neigh) + rho_nbor(i,iSN)*vol_loc/8d0) / dm_ejecta - 1d0 
           vol_nbor = (0.5d0**lv_nbor(i,iSN)*scale)**ndim
           dloadloc = mload/dble(n_neigh)/vol_loc  ! double-checked
           navg     = (dloadloc+rho_nbor(i,iSN))/2d0*scale_nH
           num_sn   = (mSN(iSN)/(M_SNII/scale_msun))
           if(metal)then
              Zdepen = (dloadloc*Zload+rho_nbor(i,iSN)*Z_nbor(i,iSN))/(dloadloc+rho_nbor(i,iSN))
              Zdepen = (max(0.01,Zdepen/0.02))**(-0.28) ! From Thornton et al. (1999)
           else
              Zdepen = 1d0
           endif 
           !From Blondin et al. (1998)
           f_w_crit = max(0d0,9d2/(f_esn*M_SNII)*num_sn**(-2d0/17d0)*navg**(-4d0/17d0)*Zdepen-1d0)
           if(f_w_cell.ge.f_w_crit)then !radiative phase (momentum-conserving phase)
              vload    = dsqrt(2d0*f_esn*ESN*(1+f_w_crit)/(M_SNII*2d33))/scale_v/(1+f_w_cell)/f_load
           else
              f_esn2   = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload    = dsqrt(2d0*f_esn2*ESN/(1+f_w_cell)/(M_SNII*2d33))/scale_v/f_load
           endif

           ! dq: net residual momentum for ind_SN
           dq(ind_SN,1) = dq(ind_SN,1) + ((1+f_w_cell)*dm_ejecta)*vload*v48(i,1) 
           dq(ind_SN,2) = dq(ind_SN,2) + ((1+f_w_cell)*dm_ejecta)*vload*v48(i,2)
           dq(ind_SN,3) = dq(ind_SN,3) + ((1+f_w_cell)*dm_ejecta)*vload*v48(i,3)
        endif
     enddo

     ! remove gas from SN cells
     if(ind_blast(iSN)>0)then  ! if SN cells are inside the 'myid' domain
         d=uold(ind_blast(iSN),1)
         u=uold(ind_blast(iSN),2)/d
         v=uold(ind_blast(iSN),3)/d
         w=uold(ind_blast(iSN),4)/d
         ekk=0.5d0*d*(u*u+v*v+w*w)
         eint=uold(ind_blast(iSN),5)-ekk
         if(metal)Zload=uold(ind_blast(iSN),imetal)/d

         ! remove 48/56*mcell gas
         d=d-(mload-mSN(iSN)*f_load)/vol_loc  
         uold(ind_blast(iSN),1)=d
         uold(ind_blast(iSN),2)=d*u
         uold(ind_blast(iSN),3)=d*v
         uold(ind_blast(iSN),4)=d*w
         uold(ind_blast(iSN),5)=eint+0.5*d*(u*u + v*v + w*w)

         ! add f_leftover*SN ejecta to this cell 
         d = mSN(iSN)*f_leftover/vol_loc
         u = vSN(iSN,1)
         v = vSN(iSN,2)
         w = vSN(iSN,3)
         uold(ind_blast(iSN),1)=uold(ind_blast(iSN),1)+d
         uold(ind_blast(iSN),2)=uold(ind_blast(iSN),2)+d*u
         uold(ind_blast(iSN),3)=uold(ind_blast(iSN),3)+d*v
         uold(ind_blast(iSN),4)=uold(ind_blast(iSN),4)+d*w
         uold(ind_blast(iSN),5)=uold(ind_blast(iSN),5)+0.5*d*(u*u+v*v+w*w) 
         if(metal)then
            uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal) &
                                      & - Zload*(mload-mSN(iSN)*f_load)/vol_loc &
                                      & + (yield + (1d0-yield)*ZSN(iSN))*mSN(iSN)*f_leftover/vol_loc
         endif
     endif
  enddo


#ifndef WITHOUTMPI
  if(global_search) then
     dq_mpi=0d0
     call MPI_ALLREDUCE(dq,dq_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     dq=dq_mpi
  endif
#endif
 

  !---------------------------------------------------------------------
  ! Inject momenta from SNe
  !---------------------------------------------------------------------
  do iSN=1,nSN
     ind_SN= iSN_myid(iSN)
     dx=0.5d0**lSN(iSN)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim   
     mload=mloadSN(ind_SN)
     do i=1,n_neigh
        if(icell_nbor(i,iSN)>0)then
           ! idx48 is measured for sub volume of vol_loc/8
           dm_ejecta= f_load*mSN(iSN)/dble(n_neigh)
           f_w_cell = (mload/dble(n_neigh) + rho_nbor(i,iSN)*vol_loc/8d0) / dm_ejecta - 1d0
           dloadloc = mload/dble(n_neigh)/vol_loc
           navg     = (dloadloc+rho_nbor(i,iSN))/2d0*scale_nH
           num_sn   = mSN(iSN)/(M_SNII/scale_msun)
           if(metal)then
              Zdepen = (dloadloc*Zload+rho_nbor(i,iSN)*Z_nbor(i,iSN))/(dloadloc+rho_nbor(i,iSN))
              Zdepen = (max(0.01,Zdepen/0.02))**(-0.28)
           else
              Zdepen = 1d0
           endif 

           f_w_crit = max(0d0,9d2/(f_esn*M_SNII)*num_sn**(-2d0/17d0)*navg**(-4d0/17d0)*Zdepen-1d0)
           if(f_w_cell.ge.f_w_crit)then !radiative phase
              vload    = dsqrt(2d0*f_esn*ESN*(1+f_w_crit)/(M_SNII*2d33))/scale_v/(1d0+f_w_cell)/f_load
           else
              f_esn2   = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload    = dsqrt(2d0*f_esn2*ESN/(1+f_w_cell)/(M_SNII*2d33))/scale_v/f_load
           endif

           icello=icell_nbor(i,iSN)
           d=uold(icello,1)      
           u=uold(icello,2)/d 
           v=uold(icello,3)/d 
           w=uold(icello,4)/d
           ekk=0.5d0*d*(u*u+v*v+w*w)
           eint=uold(icello,5)-ekk

           u0=u;v0=v;w0=w

           vol_nbor=(0.5d0**lv_nbor(i,iSN)*scale)**ndim
           ! momentum of the original star and cell
           if(mechanical_conserve_p)then
              !cons: maximum(v) can be > 4000 km/s if the SN explodes in a region
              !      where density varies abruptly
              pload(1) = ((mload*vloadSN(ind_SN,1)-dq(ind_SN,1))/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,1))/vol_nbor ! double-checked
              pload(2) = ((mload*vloadSN(ind_SN,2)-dq(ind_SN,2))/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,2))/vol_nbor
              pload(3) = ((mload*vloadSN(ind_SN,3)-dq(ind_SN,3))/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,3))/vol_nbor
           else
              pload(1) = (mload*vloadSN(ind_SN,1)/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,1))/vol_nbor ! double-checked
              pload(2) = (mload*vloadSN(ind_SN,2)/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,2))/vol_nbor
              pload(3) = (mload*vloadSN(ind_SN,3)/dble(n_neigh) + (1+f_w_cell)*dm_ejecta*vload*v48(i,3))/vol_nbor
           endif

           dload=mload/dble(n_neigh)/vol_nbor
           uold(icello,1)=d+dload
           uold(icello,2)=uold(icello,2)+pload(1)
           uold(icello,3)=uold(icello,3)+pload(2)
           uold(icello,4)=uold(icello,4)+pload(3)
 
           ! energy conservation (if momentum vector cancels out, it naturally becomes thermal)
           u=pload(1)/dload 
           v=pload(2)/dload 
           w=pload(3)/dload
           uold(icello,5)=uold(icello,5)+0.5*dload*(u*u+v*v+w*w)
           u=uold(icello,2)/(d+dload) 
           v=uold(icello,3)/(d+dload) 
           w=uold(icello,4)/(d+dload)
           eint2 = uold(icello,5)- 0.5*(d+dload)*(u*u+v*v+w*w)

           Tk  = eint2/(d+dload)*scale_T2
           if(T2maxSN>0.and.Tk>T2maxSN)then
              eint2 = (d+dload)*T2maxSN/scale_T2
              uold(icello,5) = 0.5*(d+dload)*(u*u+v*v+w*w) + eint2
           endif

           if(metal)then
               uold(icello,imetal)=uold(icello,imetal)+dload*ZloadSN(ind_SN)
           endif
        endif
     end do
  end do
 
  !if(.not.global_search) call clean_stop

end subroutine inject_momentum_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_idx48(ista,iend,jsta,jend,ksta,kend,idx48,n48)
   implicit none
   integer,parameter::n_neigh=48
   integer::ista,iend,jsta,jend,ksta,kend,i,j,k,kk,n48
   integer, dimension(1:n_neigh)::idx48
   integer, dimension(1:4,1:4,1:4)::indall
   logical::error

   n48=0
   idx48=0
   error=.false.
   if(ista<0)ista=0
   if(jsta<0)jsta=0
   if(ksta<0)ksta=0
   if(iend>4)error=.true.
   if(jend>4)error=.true.
   if(kend>4)error=.true.
   if(iend>3)iend=3
   if(jend>3)jend=3
   if(kend>3)kend=3


   indall=reshape((/&
         -1, 1, 2,-1,   3, 4, 5, 6,   7, 8, 9,10,  -1,11,12,-1,&
         13,14,15,16,  17,-1,-1,18,  19,-1,-1,20,  21,22,23,24,&
         25,26,27,28,  29,-1,-1,30,  31,-1,-1,32,  33,34,35,36,&
         -1,37,38,-1,  39,40,41,42,  43,44,45,46, -1,47,48,-1/),(/4,4,4/))
   
   do k=ksta,kend
   do j=jsta,jend
   do i=ista,iend
      kk = indall(i+1,j+1,k+1)
      if(kk>0)then
         n48=n48+1
         idx48(n48)=kk
      endif
   end do
   end do
   end do 
    
end subroutine get_idx48
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_v48(v48)
   use amr_commons
   implicit none
   integer,parameter::n_neigh=48 
   real(dp)::v48(1:n_neigh,3),x,y,z,r
   integer::i,j,k,ind,indall(1:64)
   logical::ok

   ind = 0 
   indall=-1
   ! centre: [0,0,0]
   do k=1,4
   do j=1,4
   do i=1,4
      ok = .true.
      if((i==1.or.i==4).and.&
         (j==1.or.j==4).and.&
         (k==1.or.k==4)) ok=.false. ! edge
      if((i==2.or.i==3).and.&
         (j==2.or.j==3).and.&
         (k==2.or.k==3)) ok=.false. ! centre
      if(ok)then
         ind=ind+1
         x = (i-1)+0.5d0 - 2 
         y = (j-1)+0.5d0 - 2 
         z = (k-1)+0.5d0 - 2 
         r = dsqrt(dble(x*x+y*y+z*z))
         v48(ind,1) = x/r 
         v48(ind,2) = y/r 
         v48(ind,3) = z/r 
         indall(i+(j-1)*4+(k-1)*4*4) = ind 
      endif
   end do
   end do
   end do
end subroutine get_v48
!################################################################
!################################################################
!################################################################
!################################################################

subroutine checkSNboundary(xSN,ncpu_read,mylevel)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::mylevel,ncpu_read
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,nx_loc,ilevel,lmax,bit_length,maxdom
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_loc,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:3)::xSN 

 
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=scale*0.5D0**mylevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(1.5d0*dx_loc*scale_l/aexp,rbubble*3.08d18) ! adaptive radius
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  cpu_read=.false.
  ! Compute boundaries for the SN cube of influence
  xxmin=(xSN(1)-drSN)/scale ; xxmax=(xSN(1)+drSN)/scale
  yymin=(xSN(2)-drSN)/scale ; yymax=(xSN(2)+drSN)/scale
  zzmin=(xSN(3)-drSN)/scale ; zzmax=(xSN(3)+drSN)/scale

  if(TRIM(ordering).eq.'hilbert')then
        
     dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
     do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xxmin*dble(maxdom))
        imax=imin+1
        jmin=int(yymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zzmin*dble(maxdom))
        kmax=kmin+1
     endif
     
     dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+oneqdp)*dkey
     end do
     
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if
  
end subroutine checkSNboundary
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getSNonmyid2(iSN_myid,nSN_myid,xSN,nSN,lvSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid,lvSN
  integer::nSN_myid,ii,nSN
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,iSN,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:nSN,1:3)::xSN 

 
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  !drSN=2d0*MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  !drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
  do iSN=1,nSN

     cpu_read=.false.
 
     drSN = dx_min*3*2d0**(nlevelmax-lvSN(iSN))
     ! Compute boundaries for the SN cube of influence
     xxmin=(xSN(iSN,1)-drSN)/scale ; xxmax=(xSN(iSN,1)+drSN)/scale
     yymin=(xSN(iSN,2)-drSN)/scale ; yymax=(xSN(iSN,2)+drSN)/scale
     zzmin=(xSN(iSN,3)-drSN)/scale ; zzmax=(xSN(iSN,3)+drSN)/scale

     if(TRIM(ordering).eq.'hilbert')then
        
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif
        
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+oneqdp)*dkey
        end do
        
        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if
     
     ! Create the index array for SN in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           iSN_myid(ii)=iSN
        endif
     enddo

  enddo
  
  ! Number of SN in processor myid
  nSN_myid=ii

end subroutine getSNonmyid2
!################################################################
!################################################################
!################################################################
!################################################################

