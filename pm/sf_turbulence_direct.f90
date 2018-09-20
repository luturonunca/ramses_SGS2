  ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
  ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor.
  ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation
  ! from neighbouring cell values and differentiate.
  ! Get neighbor cells if they exist, otherwise use straight injection from local cell
  ncell = 1 ! we just want the neighbors of that cell
  ind_cell2(1) = ind_cell(i)
  call getnbor(ind_cell2,ind_nbor,ncell,ilevel)
  d1           = uold(ind_nbor(1,1),1) ; d2 = uold(ind_nbor(1,2),1) ; d3 = uold(ind_nbor(1,3),1)
  d4           = uold(ind_nbor(1,4),1) ; d5 = uold(ind_nbor(1,5),1) ; d6 = uold(ind_nbor(1,6),1)
  sigma2       = 0d0 ; sigma2_comp = 0d0 ; sigma2_sole = 0d0
  trgv         = 0d0 ; divv = 0d0 ; curlva = 0d0 ; curlvb = 0d0 ; curlvc = 0d0
  flong        = 0d0
  !!!!!!!!!!!!!!!!!!
  ! Divergence terms
  !!!!!!!!!!!!!!!!!!
  ul        = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
  ur        = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
  if(sf_model.le.2) then
     fl     = (d2*f(ind_nbor(1,2),1)    + d*f(ind_cell(i),1))/(d2+d)
     fr     = (d1*f(ind_nbor(1,1),1)    + d*f(ind_cell(i),1))/(d1+d)
     flong  = flong+max((d2+d)/2*ul*fl-(d1+d)/2*ur*fr,0d0)
  endif
  sigma2_comp = sigma2_comp + (ur-ul)**2
  divv      = divv + (ur-ul)
  ul        = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
  ur        = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
  if(sf_model.le.2) then
     fl     = (d4*f(ind_nbor(1,4),2)    + d*f(ind_cell(i),2))/(d4+d)
     fr     = (d3*f(ind_nbor(1,3),2)    + d*f(ind_cell(i),2))/(d3+d)
     flong  = flong+max((d4+d)/2*ul*fl-(d3+d)/2*ur*fr,0d0)
  endif
  sigma2_comp = sigma2_comp + (ur-ul)**2
  divv      = divv + (ur-ul)
  ul        = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
  ur        = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
  if(sf_model.le.2) then
     fl     = (d6*f(ind_nbor(1,6),3)    + d*f(ind_cell(i),3))/(d6+d)
     fr     = (d5*f(ind_nbor(1,5),3)    + d*f(ind_cell(i),3))/(d5+d)
     flong  = flong+max((d6+d)/2*ul*fl-(d5+d)/2*ur*fr,0d0)
  endif
  sigma2_comp = sigma2_comp + (ur-ul)**2
  divv      = divv + (ur-ul)
  ftot      = flong
  !!!!!!!!!!!!
  ! Curl terms
  !!!!!!!!!!!!
  ul        = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
  ur        = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
  if(sf_model.le.2) then
     fl     = (d6*f(ind_nbor(1,6),2)    + d*f(ind_cell(i),2))/(d6+d)
     fr     = (d5*f(ind_nbor(1,5),2)    + d*f(ind_cell(i),2))/(d5+d)
     ftot   = ftot+abs((d6+d)/2*ul*fl-(d5+d)/2*ur*fr)
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlva    = curlva-(ur-ul)
  ul        = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
  ur        = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
  if(sf_model.le.2) then
     fl     = (d4*f(ind_nbor(1,4),3)    + d*f(ind_cell(i),3))/(d4+d)
     fr     = (d3*f(ind_nbor(1,3),3)    + d*f(ind_cell(i),3))/(d3+d)
     ftot   = ftot+abs((d4+d)/2*ul*fl-(d3+d)/2*ur*fr)
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlva    = (curlva + (ur-ul))
  ul        = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
  ur        = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
  if(sf_model.le.2) then
     fl     = (d6*f(ind_nbor(1,6),1)    + d*f(ind_cell(i),1))/(d6+d)
     fr     = (d5*f(ind_nbor(1,5),1)    + d*f(ind_cell(i),1))/(d5+d)
     ftot   = ftot+abs((d6+d)/2*ul*fl-(d5+d)/2*ur*fr)
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlvb    = curlvb+(ur-ul)
  ul        = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
  ur        = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
  if(sf_model.le.2) then
     fl     = (d2*f(ind_nbor(1,2),3)    + d*f(ind_cell(i),3))/(d2+d)
     fr     = (d1*f(ind_nbor(1,1),3)    + d*f(ind_cell(i),3))/(d1+d)
     ftot   = ftot+abs((d2+d)/2*ul*fl-(d1+d)/2*ur*fr)
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlvb    = (curlvb - (ur-ul))
  ul        = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
  ur        = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
  if(sf_model.le.2) then
     fl     = (d4*f(ind_nbor(1,4),1)    + d*f(ind_cell(i),1))/(d4+d)
     fr     = (d3*f(ind_nbor(1,3),1)    + d*f(ind_cell(i),1))/(d3+d)
     ftot   = ftot+abs((d4+d)/2*ul*fl-(d3+d)/2*ur*fr)
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlvc    = curlvc-(ur-ul)
  ul        = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
  ur        = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
  if(sf_model.le.2) then
     fl     = (d2*f(ind_nbor(1,2),2)    + d*f(ind_cell(i),2))/(d2+d)
     fr     = (d1*f(ind_nbor(1,1),2)    + d*f(ind_cell(i),2))/(d1+d)
     ftot   = ftot+abs((d2+d)/2*ul*fl-(d1+d)/2*ur*fr)
     pcomp  = flong/ftot
  endif
  sigma2_sole = sigma2_sole + (ur-ul)**2
  curlvc    = (curlvc + (ur-ul))
  sigma2    = sigma2_comp+sigma2_sole
  ! Trace of gradient velocity tensor
  trgv      = sigma2/dx_loc**2
  ! Velocity vector divergence
  divv      = divv/dx_loc
  ! Velocity vector curl
  curlv     = (curlva+curlvb+curlvc)/dx_loc
  divv2     = divv**2
  curlv2    = curlv**2
  if(sf_save_sigma2) then
     if(sf_compressive)then
        uold(ind_cell(i),ivirial1) = sigma2_comp
        uold(ind_cell(i),ivirial2) = sigma2_sole
     else
        uold(ind_cell(i),ivirial1) = sigma2
     endif
  end if
endif
