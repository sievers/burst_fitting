if ~exist('dat')
  dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:))+squeeze(dat_raw(:,4,:)); 
  clear dat_raw
  tvec=read_npy('time.npy');
  dt=median(diff(tvec));
end

sigs=std(dat);
x=(1:length(sigs));x=x-mean(x);x=x/max(x);
good_chan=(1:length(sigs));

x_org=x;
sigs_org=sigs;
ii=isfinite(sigs);
x=x(ii);
sigs=sigs(ii);
good_chan=good_chan(ii);

pp=polyfit(x,sigs,2);
mymod=polyval(pp,x);

delt=abs(sigs-mymod);
isbad=delt>4*median(delt);


xx=x(~isbad);
ss=sigs(~isbad);
good_chan=good_chan(~isbad);

pp=polyfit(xx,ss,3);
mymod2=polyval(pp,xx);

delt=abs(ss-mymod2);
isbad2=delt>4*median(delt);
good_chan=good_chan(~isbad2);


%clf;plot(x,sigs);hold on;plot(x,polyval(pp,x),'r'); plot(x(isbad),sigs(isbad),'k.');plot(x_org(good_chan),sigs_org(good_chan),'g');
freq_use=freq(good_chan);
dat_use=dat(:,good_chan);


%dd=dat_use;dm=625;dm0=4150;lags=dm*dm0*(1./freq_use.^2);lags=lags-min(lags);lags=round(lags/dt);
%for j=1:length(lags), dd(:,j)=circshift(dd(:,j),-lags(j));end;


scat=0.05;
alpha=-1;
dm=625;
fwhm=3e-3;
t0=23577*dt;
amp=1.5;


if (0)
  tic;flub=make_burst_model(dt,length(dat_use),freq_use([1 round(end/2) end]),fwhm,1/scat,alpha,t0,dm);toc;ff=flub;ff(end,:)=real(ff(end,:));ffr=real(ifft([ff;flipud(conj(ff(2:end-1,:)))]));
  [a,b]=max(ffr(:,2));myax=[b-200 b+200 -0.01 a*1.1];
  tic;flub=make_burst_model_c(dt,length(dat_use),freq_use,fwhm,scat,alpha,t0,dm);toc;
end



%guess=[fwhm 1./scat  alpha amp t0 dm];
guess=[fwhm scat  alpha amp t0 dm];
%sigs=[1e-3 0.01 0.05 0.1 2*dt 3]


sigs=[3e-3 3e-3 0.3 0.5 0.006 2];

%guess2=[0.0064756    0.010629     -3.6745      13.231       24.14 622.5];
%guess2=[ 0.002365    0.013541     -4.2628      13.742      24.142 622.11];
%guess2=[ 0.0065569   0.0051592      -4.681      12.057      24.139 624.58];

%guess2=[ 0.0026688    0.002345     -4.8478      9.7342      24.143      623.64];
%guess2=[3.0847000e-03   2.2877000e-03  -4.7066000e+00   1.0396000e+01   2.4143000e+01   6.2312000e+02 ];
%guess2=[ 3.0883000e-03   1.5139000e-03  -4.2911000e+00   9.4774000e+00   2.4143000e+01   6.2313000e+02];
guess2=[3.0989000e-03   9.0000000e-04  -4.4348000e+00   8.7676000e+00   2.4143000e+01   6.2332000e+02];
guess2=[2.4495946e-03   1.1714976e-03  -4.5047313e+00   8.5505905e+00   2.4143212e+01   6.2320896e+02 ];

%sigs=0.3*[4e-4 4e-4 0.3 0.5 0.006 2]; 
%sigs=0.3*[2e-4 2e-4 0.1 0.5 0.005 1]; 

guess2=[ 0.0024385974200246  0.00123751806476478  -7.10722299016074 10 24.1432294156555 623.208989011357 ];
sigs=0.3*[ 0.000265964168798882 0.000154557403917922  0.651994753250518 1 9.50335070713428e-05 0.0470434921823387];



best_guess=[2.0310062e-03   1.6932547e-03  -7.7203859e+00   8.0534911e+00   2.4143071e+01   6.2318296e+02];
sigmat=[  9.7173e-08  -5.4161e-08  5.2972e-05  1.5526e-05  1.6728e-08  4.4094e-07
  -5.4161e-08  7.8293e-08  -3.9634e-05  6.1869e-05  -2.1855e-08  -1.632e-06
  5.2972e-05  -3.9634e-05     0.53966     0.21587  2.7553e-05  -0.0051925
  1.5526e-05  6.1869e-05     0.21587      0.2051  -4.5197e-06  -0.0064241
  1.6728e-08  -2.1855e-08  2.7553e-05  -4.5197e-06  6.1204e-08  -2.4178e-05
  4.4094e-07  -1.632e-06  -0.0051925  -0.0064241  -2.4178e-05    0.012727
];



[dat_filt,nuvec,mylags]=stack_data(dat_use,freq_use,best_guess,dt,3);

dd1=mean(dat_filt(:,1:round(end/3)),2);
dd2=mean(dat_filt(:,round(end/3):round(2*end/3)),2);
dd3=mean(dat_filt(:,round(2*end/3):end),2);
myt=((1:length(dd1))-1000)*dt;
clf;hold on;h1=plot(myt,dd1/max(dd1),'b');h2=plot(myt,dd2/max(dd2),'g');h3=plot(myt,dd3/max(dd3),'r');axis([-0.03 0.08 -0.2 1.1])
set(h1,'linewidth',3);set(h2,'linewidth',3);set(h3,'linewidth',3);legend('high','middle','low');print -dpng burst_profile

clf;hold on;h1=plot(myt,dd1,'b');h2=plot(myt,dd2,'g');h3=plot(myt,dd3,'r');axis([-0.03 0.08 -1 5])
set(h1,'linewidth',3);set(h2,'linewidth',3);set(h3,'linewidth',3);legend('high','middle','low');print -dpng burst_profile_calib




%[dat,nuvec,mylags]=stack_data(dat_use,freq_use,best_guess,dt);
return
[pp,ll,dataft,nvec]=mcmc_burst_ampfit(dat_use,freq_use,guess2,sigs,dt,1.0);


ind=6;
guess3=guess2;vv=guess2(ind)*(1+(-0.2:0.01:0.2));
if ind==5, vv=guess2(ind)+(-0.1:0.01:0.1)*5*0.001;end;
if ind==6, vv=guess2(ind)+(-0.3:0.03:0.3);end;
chisq=zeros(size(vv));amps=zeros(size(vv));
for j=1:length(vv), guess3(ind)=vv(j);[chisq(j),amps(j)]=get_burst_chisq_ampfit_c(guess3,dataft,freq_use,length(dat_use),nvec,dt,57);end;
chisq=chisq-min(chisq);fitp=polyfit(vv,chisq,2);[-0.5*fitp(2)/fitp(1) sqrt(1/fitp(1))]
%[ 0.0024385974200246  0.000265964168798882]
%[ 0.00123751806476478  0.000154557403917922]
%[ -7.10722299016074     0.651994753250518]

%[24.1433040142732  0.000264918623491682]
%[24.1432294156555  9.50335070713428e-05]
%[623.208989011357    0.0470434921823387];