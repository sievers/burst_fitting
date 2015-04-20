function[dat,dat_q,dat_u,freq_use,nuvec,mylags]=get_trimmed_burst()

dat_raw=read_npy('calibrated.npy');
freq=read_npy('freq.npy');

dat=squeeze(dat_raw(:,1,:))+squeeze(dat_raw(:,4,:)); 
dat_q=(squeeze(dat_raw(:,1,:)-dat_raw(:,4,:)));
dat_u=2*squeeze(dat_raw(:,2,:));
clear dat_raw;
tvec=read_npy('time.npy');
dt=median(diff(tvec));


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
dat=dat(:,good_chan);
dat_q=dat_q(:,good_chan);
dat_u=dat_u(:,good_chan);

%dd=dat_use;dm=625;dm0=4150;lags=dm*dm0*(1./freq_use.^2);lags=lags-min(lags);lags=round(lags/dt);
%for j=1:length(lags), dd(:,j)=circshift(dd(:,j),-lags(j));end;

chains=load('chain_ampfit.txt.useful');
fitp=mean(chains(:,2:end));fwhm=fitp(1);scat=fitp(2);alpha=fitp(3);amp=fitp(4);t0=fitp(5);dm=fitp(6);


highpass=3;
[dat,nuvec,mylags]=stack_data(dat,freq_use,fitp,dt,highpass);
dat_q=stack_data(dat_q,freq_use,fitp,dt,highpass);
dat_u=stack_data(dat_u,freq_use,fitp,dt,highpass);

