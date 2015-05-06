try
  mpi_init;
end

if ~exist('dat')
  dat_raw=read_npy('/scratch2/p/pen/kiyo/burst_data/filtered.npy');
  %dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:));
  dat_q=-1*squeeze(dat_raw(:,2,:));
  dat_u=squeeze(dat_raw(:,3,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=(sum(isnan(dat))==0)&(sum(isnan(dat_u))==0);
  dat=dat(:,good_chan);
  dat_q=dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  freq_use=freq(good_chan);
end

scat=0.05;
alpha=-1;
dm=625;
fwhm=3e-3;
t0=23577*dt;
amp=1.5;


crap=read_chains('chains/chain_pol_kiyo_test.txt',0.2,'');crap=crap(:,3:end);
best_guess=mean(crap);
mycov=cov(crap);

myopts.nstep=300000;
myopts.outroot='chains/chain_pol_kiyo_wcov.txt';

[pp,ll,qft,uft,nvec,imin]=mcmc_burst_pol2(dat_q,dat_u,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);

return

tvec=dt*(-20:0.5:10)';
chivec=0*tvec;
model=floatcomplex2ptr(qft);
qptr=floatcomplex2ptr(qft);
uptr=floatcomplex2ptr(uft);


phivec=(0:10:360)'*pi/180;
chivec=0*phivec;

for j=1:length(tvec),
  gg=best_guess;gg(8)=phivec(j);
  chivec(j)=get_burst_chisq_qu_cached_c(gg,qptr,uptr,freq_use,length(dat),nvec,dt,imin,model);
end
return


for j=1:length(tvec),
  gg=best_guess;gg(5)=gg(5)+tvec(j);
  chivec(j)=get_burst_chisq_qu_cached_c(gg,qptr,uptr,freq_use,length(dat),nvec,dt,imin,model);
end
return
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
