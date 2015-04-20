scat=0.05;
alpha=-1;
dm=400;
freq=500:900;
fwhm=3e-3;
t0=1;
n=1e4;
amp=1;
tic;[flub,vec]=make_burst_model(1e-3,n,freq,fwhm,1/scat,alpha,t0,dm);toc
ff=flub;ff(end,:)=real(ff(end,:));ff=[ff;flipud(conj(flub(2:end-1,:)))];ff=amp*real(ifft(ff));
vv=[vec;flipud(conj(vec(2:end)))];vv=real(ifft(vv));


ffmax=max(max(ff));
ff_noisy=ff+0.2*ffmax*randn(size(ff));





guess=[fwhm 1./scat alpha amp t0 dm];
sigs=abs(0.25*guess);sigs(end-1)=fwhm;sigs(end)=2;



%to estimate sigmas:
%gg=guess;vec=399.5:0.01:400.5;chisq=0*vec;
%for j=1:length(vec), gg(6)=vec(j);chisq(j)=get_burst_chisq(gg,dataft,freq,length(ff),nvec);end;
%chisq=chisq-min(chisq);sum(chisq<1)
%sum(chisq<1)*mean(diff(vec))
%sigs(1)=0.00034
%sigs(2)=0.08
%sigs(3)=0.008
%sigs(4)=0.003
%sigs(5)=1e-4
%sigs(6)=0.02

%sigs=[0.00034 0.08 0.008 0.003 1e-4 0.02];
sigs=[0.00017211     0.12576    0.014941   0.0064253  6.5863e-05    0.023661];  %from chain


%from chain - 2d
%sigs=[
%  3.2138e-08  1.1649e-06  -8.1156e-08  -1.8006e-08  6.2044e-10  -3.2825e-07
%  1.1649e-06    0.013378  -4.8609e-05  -0.00032579  5.9499e-07  -1.0632e-05
%  -8.1156e-08  -4.8609e-05  0.00033141  0.00010683  4.0864e-08  -2.0282e-05
%  -1.8006e-08  -0.00032579  0.00010683  4.9596e-05  -5.384e-09  -4.6308e-06
%  6.2044e-10  5.9499e-07  4.0864e-08  -5.384e-09  3.9393e-09  -8.9239e-07
%  -3.2825e-07  -1.0632e-05  -2.0282e-05  -4.6308e-06  -8.9239e-07  0.00060747];

%from longer chain
sigs=[  3.3267e-08  1.4258e-06  9.8586e-08  5.8071e-08  -7.0193e-11  1.3763e-08
  1.4258e-06      0.0128  -0.00012439  -0.00034927  5.9783e-07  2.6698e-05
  9.8586e-08  -0.00012439  0.00030822  0.00010429  2.4538e-08  -2.6907e-05
  5.8071e-08  -0.00034927  0.00010429  5.1249e-05  -9.9741e-09  -7.1229e-06
  -7.0193e-11  5.9783e-07  2.4538e-08  -9.9741e-09  5.0873e-09  -1.2441e-06
  1.3763e-08  2.6698e-05  -2.6907e-05  -7.1229e-06  -1.2441e-06  0.00065939
];



%[pp,ll]=mcmc_burst(ff_noisy,freq,guess+[6e-3 0 0 0 0 0 ],sigs,[false true true true true true],guess);
%[pp,ll,dataft,nvec]=mcmc_burst(ff_noisy,freq,guess+2*randn(size(sigs)).*sigs,0.2*sigs,[false true true true true true],guess);
[pp,ll,dataft,nvec]=mcmc_burst(ff_noisy,freq,guess+0.2*randn(1,size(sigs,2)).*sqrt(diag(sigs)'),0.4*sigs);

return



dt=1e-3;
nu=700:900;
tmax=10;
tt=0:dt:tmax;
nn=floor(numel(tt)/2);
D0=4160;

datft=zeros(nn,numel(nu));
dnu=1/(tt(end)-tt(1));


%integral exp(-ax)exp(ikx)dx  x=[0,inf]
%         exp( (ik-a)x)dx     x=[0,inf]
%         1/(ik-a)exp(ik-a)xdx  
%         -1/(ik-a)
%dnu=0.1;k=[0:dnu:10 (-10):dnu:-dnu];

%integral sucks.  Try sum
%sum_0,n exp(2pi*i*k/n*x)*exp(-a*x)
%sum exp((2*pi*i*k/n-a)*x)
%k'=2*pi*k/n, have sum(exp(ik'-a *x))
%x'=x+dx, sum(exp( (ik'-a)(x+dx)))=exp(ik'-a *dx)sum

%after doing math, get sum is 1-exp(-adx)*exp(2*pi*i*k/n)

x=(0:1023)';
a=0.01;
y=exp(-a*x);
yft=fft(y);
%yft2=1-exp(-a*x).*exp(2*pi*i*x/numel(x));
aa=-a-2*pi*i*x/numel(x);yft2=(1-exp(numel(x)*aa))./(1-exp(aa));
%yft2=(1-exp(-./(1-exp(-a-2*pi*i*x/numel(x)));

sig=4;t0=17;scat=0.1;n=numel(x);

yft=ones(size(x));yft=yft.*exp(-0.5*(2*pi*sig*x).^2/n^2);
aa=-scat-2*pi*i*x/n;yft=yft.*(1-exp(n*aa))./(1-exp(aa));tmp=(1-exp(n*aa))./(1-exp(aa));
norm=(1-exp(-scat*n))/(1-exp(-scat));yft=yft/norm;
yft=yft.*exp(-2*pi*i*x*t0/n);


nn=ceil( (numel(x)+1)/2);
if iseven(n)
  yft=[yft(1:nn);conj(flipud(yft(2:nn-1)))];
  yft(nn)=real(yft(nn));
else
  yft=[yft(1:nn);conj(flipud(yft(2:nn)))];
end
%yft=[yft;conj(flipud(yft(2:end)))];

yy=real(ifft(yft));

if (1)
  %this is a check.  zz should be same as yy
  xx=x;xx(nn:end)=xx(nn:end)-n;
  z1=exp(-0.5*xx.^2/sig^2);z1=z1/sum(z1);
  z2=exp(-scat*x);z2=z2/sum(z2);
  zz=real(ifft(fft(z1).*fft(z2)));
  zz=circshift(zz,t0);
  sum(abs(zz-yy))/sum(abs(zz))
end



