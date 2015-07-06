try
  mpi_init;
end
more off

if ~exist('dat')
  dat_raw=read_npy('filtered.npy');
  %dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:));
  dat_q=squeeze(dat_raw(:,2,:));
  dat_u=squeeze(dat_raw(:,3,:));
  dat_v=squeeze(dat_raw(:,4,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=sum(isnan(dat))==0;
  dat=dat(:,good_chan);
  dat_q=dat_u(:,good_chan);
  dat_q=dat_u(:,good_chan);
  dat_v=dat_v(:,good_chan);
  freq_use=freq(good_chan);
  %crap=read_chains('chains/chain_tt_kiyo_test.txt',0.2,'');crap=crap(:,3:end);
  %crap=load('chain_circular_kiyo_test.txt_1.org');crap=crap(100:end,3:end);
  %crap=load('chains/chain_circular_kiyo_test2.txt_1.org');crap=crap(100:end,3:end);
  
  crap=load('chains/chain_tt_real_newdat.txt_1');crap=crap(100:end,3:end);


end


best_guess=mean(crap);best_guess(4)=1.5;
mycov=cov(crap);

myopts.nstep=300000;
%myopts.outroot='chains/chain_circular_kiyo_test2.txt';
myopts.outroot='chains/chain_circular_kiyo_powlaw_test.txt';
%myopts.func=@get_burst_chisq_linfit_cached_c;
[pp,ll,vft,nvec,imin]=mcmc_burst_intensity(dat_v,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);
return

if ~exist('vptr')
  model=floatcomplex2ptr(vft);
  vptr=floatcomplex2ptr(vft);
end

offset=(-2:0.1:-1);
slope=(1.5:0.1:2.5);
chimat=zeros(length(offset),length(slope));
gg=best_guess;gg(3)=0;gg(4)=0;chi_ref=get_burst_chisq_linfit_cached_c(gg,vptr,freq_use,length(dat),nvec,dt,imin,model);
for j=1:length(offset),
  for k=1:length(slope),
    gg=best_guess;gg(3)=slope(k);gg(4)=offset(j);
    chimat(j,k)=get_burst_chisq_linfit_cached_c(gg,vptr,freq_use,length(dat),nvec,dt,imin,model);
  end
end

[a,b]=min(chimat);[a,c]=min(a);b=b(c);disp([offset(b) slope(c) a-chi_ref chimat(b,c)-chi_ref])
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
