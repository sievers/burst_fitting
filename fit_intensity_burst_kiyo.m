try
  mpi_init;
end

if ~exist('dat')
  dat_raw=read_npy('/scratch2/p/pen/kiyo/burst_data/filtered.npy');
  dat_raw=read_npy('calibrated.npy');
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
best_guess=mean(crap);
mycov=cov(crap);



myopts.nstep=300000;
myopts.outroot='chains/chain_tt_kiyo.txt';
[pp,ll,datft,nvec]=mcmc_burst_intensity(dat,myopts,freq_use,best_guess,0.25*mycov*0.5,dt,2.0);

return
