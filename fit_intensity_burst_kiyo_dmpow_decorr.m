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
  crap=read_chains('chains/chain_tt_kiyo_test.txt',0.2,'');crap=crap(:,3:end);
end

scat=0.05;
alpha=-1;
dm=625;
fwhm=3e-3;
t0=23577*dt;
amp=1.5;


%crap=read_chains('chains/chain_ttfit.txt',0.2,'');crap=crap(:,2:end);
%crap(:,5)=crap(:,5)-5*dt;  %appears Kiyo shifted some samples
%crap(:,4)=crap(:,4);
%best_guess=mean(crap);
%best_guess=[0.00196876661786207   0.00167518624882759     -8.03136105572414           3.914862908      24.1407924641379      623.111481648275];

best_guess=mean(crap);
mycov=cov(crap);
best_guess(end+1)=-2;
mycov(end+1,end+1)=3e-4^2;


myopts.nstep=3e5;
myopts.outroot='chains/chain_tt_kiyo_dmpow_decorr.txt';
myopts.func=@get_burst_chisq_dmpow_cached_c;
[pp,ll,dataft,nvec,imin]=mcmc_burst_intensity(dat,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);

return

model=floatcomplex2ptr(dataft);
dataptr=floatcomplex2ptr(dataft);

chi_ref=get_burst_chisq_dmpow_cached_c(best_guess,dataptr,freq_use,length(dat),nvec,dt,imin,model);
gg=best_guess;gg(4)=0;
chi0=get_burst_chisq_dmpow_cached_c(gg,dataptr,freq_use,length(dat),nvec,dt,imin,model);


if 0
  %pvec=(-2.01:1e-3:-1.99)';
  pvec=(-0.001:1e-4:0.001)';
  chivec=0*pvec;
  
  for j=1:length(pvec),
    gg=best_guess;gg(end)=gg(end)+pvec(j);
    chivec(j)=get_burst_chisq_dmpow_cached_c(gg,dataptr,freq_use,length(dat),nvec,dt,imin,model);
  end
  [pvec chivec-chi_ref]
end

deps=2e-4;
delta_t=(-2e-3:1e-4:2e-3)';
chivec=0*delta_t;
for j=1:length(delta_t),
  gg=best_guess;gg(end)=gg(end)+deps;gg(5)=gg(5)+delta_t(j);
  chivec(j)=get_burst_chisq_dmpow_cached_c(gg,dataptr,freq_use,length(dat),nvec,dt,imin,model);
end
disp([delta_t chivec-chi_ref])

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
