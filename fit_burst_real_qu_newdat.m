try
  mpi_init;
end
more off
format short g

if ~exist('dat')
  dat_raw=read_npy('filtered.npy');
  %dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:));
  dat_q=squeeze(dat_raw(:,2,:));
  dat_u=squeeze(dat_raw(:,3,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=sum(isnan(dat))==0;
  dat=dat(:,good_chan);
  dat_q=-1*dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  freq_use=freq(good_chan);

  %crap=read_chains('chains/chain_tt_kiyo_test.txt',0.2,'');crap=crap(:,3:end);
  crap2=load('chains/chain_tt_real_test3.txt_1');crap2=crap2(:,3:end);
  crap=load('chains/chain_pol_kiyo_rmpow_testcov.txt_1');crap=crap(round(0.2*end):end,3:end);
  ii=(crap(:,end)>1.9)&(crap(:,end)<2.1);
  crap=crap(ii,1:end-1);
  crap(:,4)=crap(:,4)*mean(crap2(:,4))/4;
  best_guess=mean(crap);
  best_guess(5)=mean(crap2(:,5));
  mycov=cov(crap);

   
end

myopts.nstep=500000;
myopts.outroot='chains/chain_qu_real_newdat_test.txt';
[pp,ll,nvec]=mcmc_burst_real_qu(dat_q,dat_u,myopts,freq_use,best_guess,0.5*mycov,dt);

return

qptr=float2ptr(-1*dat_q);
uptr=float2ptr(dat_u);
model=float2ptr(dat_q);
nvec=mean(dat_q.^2);nvec=1./nvec.^2;

get_burst_real_chisq_qu_cached_c(best_guess,qptr,uptr,freq_use,length(dat_q),nvec,dt,model);

return


tvec=dt*(-10:0.5:0)';
chivec=0*tvec;
model=floatcomplex2ptr(dataft);
dataptr=floatcomplex2ptr(dataft);

for j=1:length(tvec),
  gg=best_guess;gg(5)=gg(5)+tvec(j);
  chivec(j)=get_burst_chisq_cached_c(gg,dataptr,freq_use,length(dat),nvec,dt,imin,model);
end

gg0=best_guess;gg0(5)=gg0(5)-5*dt;
avec=(0:0.5:6)';chivec=0*avec;
for j=1:length(avec),
  gg=gg0;gg(4)=avec(j);
  chivec(j)=get_burst_chisq_cached_c(gg,dataptr,freq_use,length(dat),nvec,dt,imin,model);
end



[tvec chivec-min(chivec)]
return

myopts.nstep=300000;
myopts.outroot='chains/chain_tt_kiyo.txt';
[pp,ll,datft,nvec]=mcmc_burst_intensity(dat,myopts,freq_use,best_guess,0.25*mycov*0.5,dt,2.0);

return
