try
  mpi_init;
end

if ~exist('dat')
  dat_raw=read_npy('calibrated.npy');
  freq=read_npy('freq.npy');
  dat=squeeze(dat_raw(:,1,:))+squeeze(dat_raw(:,4,:))/2; 
  dat_q=squeeze(dat_raw(:,1,:)-dat_raw(:,4,:))/2;
  dat_u=squeeze(dat_raw(:,2,:));
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
dat_q_use=dat_q(:,good_chan);
dat_u_use=dat_u(:,good_chan);



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


crap=load('chain_ampfit.txt.useful');crap=crap(:,2:end);
best_guess=mean(crap);best_guess(4)=best_guess(4)/4;
sigs=std(crap);sigs(4)=sigs(4)/4;


best_guess(end+1:end+2)=[-384 66*pi/180];
best_guess(4)=1.5;
sigs(end+1:end+2)=[2 2*pi/180];

crap=load('chain_polfit2.txt');best_guess=mean(crap(round(end/2):end,2:end));
sigs=std(crap(round(end/2):end,2:end));
best_guess(4)=1.837;sigs(4)=0.17;
best_guess(end-1)=-373.66;sigs(end-1)=5.2;
best_guess(end)=3.211;sigs(end)=0.093;

[pp,ll,qft,uft,nvec]=mcmc_burst_polfit(dat_q_use,dat_u_use,freq_use,best_guess,sigs*0.5,dt,1.0);

try 
  mpi_barrier;
end


%big_RM=(-450:-330)';big_amp=0*big_RM+best_guess(4);big_phi=0*big_RM*0+best_guess(end);

%big_RM=-500:500;big_amp=0*big_RM+best_guess(4);big_phi=0*big_RM*0+best_guess(end);
%long_chisq=get_burst_many_chisq_qu_c(best_guess,qft,uft,freq_use,length(dat_q_use),nvec,dt,57,big_amp,big_RM,big_phi);

return
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