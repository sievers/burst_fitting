try
  mpi_init;
end

if ~exist('dat')
  dat_raw=read_npy('/scratch2/p/pen/kiyo/burst_data/filtered.npy');
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
  dat_q=dat_u(:,good_chan);
  dat_q=dat_u(:,good_chan);
  freq_use=freq(good_chan);
end

scat=0.05;
alpha=-1;
dm=625;
fwhm=3e-3;
t0=23577*dt;
amp=1.5;


crap=read_chains('chains/chain_ttfit.txt',0.2,'');crap=crap(:,2:end);
crap(:,5)=crap(:,5)-5*dt;  %appears Kiyo shifted some samples
crap(:,4)=crap(:,4);
%best_guess=mean(crap);
best_guess=[0.00196876661786207   0.00167518624882759     -8.03136105572414           3.914862908      24.1407924641379      623.111481648275];
mycov=cov(crap);



myopts.nstep=300000;
myopts.outroot='chains/chain_tt_kiyo.txt';
[pp,ll,dataft,nvec,imin]=mcmc_burst_intensity(dat,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);
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
