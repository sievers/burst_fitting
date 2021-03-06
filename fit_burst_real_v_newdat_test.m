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
  dat_v=squeeze(dat_raw(:,4,:));
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
  good_chan=sum(isnan(dat))==0;
  dat=dat(:,good_chan);
  dat_q=dat_q(:,good_chan);
  dat_u=dat_u(:,good_chan);
  dat_v=dat_v(:,good_chan);
  freq_use=freq(good_chan);

  %crap=read_chains('chains/chain_tt_kiyo_test.txt',0.2,'');crap=crap(:,3:end);
  crap2=load('chains/chain_tt_real_test3.txt_1');crap2=crap2(round(0.2*end):end,3:end);
  crap=load('chains/chain_circular_kiyo_linear.txt_1');crap=crap(round(0.2*end):end,3:end);
  best_guess=mean(crap);
  best_guess(5)=mean(crap2(:,5));
  %ref_freq= 764.2;
  %top_freq=max(freq_use);
  %tshift=4150*best_guess(6)*(1./ref_freq^2-1./top_freq^2);
  %best_guess(5)=best_guess(5)+tshift;

  mycov=cov(crap);

   
end


%dataptr=float2ptr(dat);
%model=float2ptr(dat);
%nvec=double(std(dat));nvec=1./nvec.^2;
%mm=make_burst_model_real_c(dt,length(dat),freq_use,best_guess(1),best_guess(2),best_guess(3),best_guess(5),best_guess(6));

%chi1=sum(   sum(dat.^2).*nvec);
%chi2=sum(   sum(  (dat-mm*best_guess(4)).^2).*nvec);
%dchi=chi2-chi1

%chisq=get_burst_real_chisq_cached_c(best_guess,dataptr,freq_use,length(dat),nvec,dt,model);



myopts.nstep=500000;
myopts.outroot='chains/chain_v_linfit_real_newdat_test.txt';
myopts.func=@get_burst_real_chisq_linfit_cached_c;
[pp,ll,nvec]=mcmc_burst_real(dat_v,myopts,freq_use,best_guess,0.5*mycov,dt,2.0);

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
