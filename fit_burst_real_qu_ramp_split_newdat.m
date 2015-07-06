try
  mpi_init;
end
more off
format short g

if ~exist('dat')
  dohigh=false;
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
  ishigh=freq>800;
  if (dohigh)
    good_chan(~ishigh)=false;
  else
    good_chan(ishigh)=false;
  end
  if dohigh
    tag='high';
  else
    tag='low';
  end


  dat=dat(:,good_chan);
  dat_q=-1*dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  freq_use=freq(good_chan);


  %crap=load('chains/chain_qu_ramp_conv_real_newdat.txt_1');crap=crap(round(0.2*end):end,3:end);
  crap=load(['chains/chain_qu_ramp_real_newdat_good_' tag 'chan.txt_1.trial']);crap=crap(round(0.2*end):end,3:end);
  best_guess=mean(crap);
  mycov=cov(crap);
  %mycov(end,end)=mycov(end,end);
   
end

myopts.nstep=500000;

myopts.outroot=['chains/chain_qu_ramp_real_newdat_good_' tag 'chan.txt'];
%myopts.func=@get_burst_real_chisq_qu_ramp_conv_cached_c;
myopts.func=@get_burst_real_chisq_qu_ramp_cached_c;
[pp,ll,nvec]=mcmc_burst_real_qu(dat_q,dat_u,myopts,freq_use,best_guess,0.5*mycov,dt);

return

DM=0;t0=0.2;dt=1e-3;n=round(1/dt);fwhm=5e-3;scat=5e-3;phi=0;omega=0;nu=764.2;sig=fwhm/sqrt(8*log(2));
vec=make_burst_model_real_qu_swing_c(dt,n,nu,fwhm,scat,0.0,t0,DM,omega);
omega=0.2;
vec2=make_burst_model_real_qu_swing_c(dt,n,nu,fwhm,scat,0.0,t0,DM,omega);
t=dt*(0:(n-1))';

%do the analytic comparison to make sure integral is behaving
v1=exp(-t/scat);v2=exp(-0.5*(t-t0).^2/sig^2);v3=exp(i*omega*(t-t0)/dt);vv=ifft(fft(v1).*fft(v2.*v3))/sum(v1)/sum(v2);

mm=[t real(vec) imag(vec) real(vec2) imag(vec2) real(vv) imag(vv)];
fid=fopen("burst_models_ramp.txt","w");fprintf(fid,'%8.4f %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n',mm');fclose(fid);







if (~exist('model'))
  model=floatcomplex2ptr(dat_q+i*dat_u);
  qptr=float2ptr(dat_q);
  uptr=float2ptr(dat_u);

  nvec=mean(dat_q.^2);nvec=1./nvec;
end
phivec=(0:0.1:2*pi)';chivec=phivec*0;

chi_ref=get_burst_real_chisq_qu_ramp_cached_c(best_guess,qptr,uptr,freq_use,length(dat_q),nvec,dt,model)
for j=1:length(phivec),
  gg=best_guess;gg(end)=phivec(j);chivec(j)=get_burst_real_chisq_qu_ramp_conv_cached_c(gg,qptr,uptr,freq_use,length(dat_q),nvec,dt,model);
end


return

phivec=(0:0.01:1)';chivec=0*phivec;
for j=1:length(phivec), gg=best_guess;gg(end)=phivec(j);chivec(j)=get_burst_real_chisq_qu_ramp_cached_c(gg,qptr,uptr,freq_use,length(dat_q),nvec,dt,model);end

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
