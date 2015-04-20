function[mat,scatvec]=make_burst_model(dt,n,freq,fwhm,scat,alpha,t0,DM)
DM0=4150;
%freq_ref=1000; %reference frequency for scattering
freq_ref=800; %reference frequency for scattering & spectral index


nn=ceil((1+n)/2);
mat=zeros(nn,numel(freq));

scat=scat*dt;
sig=fwhm/sqrt(8*log(2))/dt;
t0=1+t0/dt;
dm_lags=DM0*DM./(freq.^2);dm_lags=dm_lags-min(dm_lags);
dm_lags=dm_lags/dt;
nfreq=numel(freq);
nu=(0:nn-1)';
sigvec=exp(-0.5*(2*pi*sig*nu).^2/n^2);
for j=1:nfreq,

  mylag=dm_lags(j)+t0;
  myscat=scat/(freq(j)/freq_ref)^-4;
  aa=-myscat-2*pi*i*nu/n;scatvec=(1-exp(n*aa))./(1-exp(aa));
  mynorm=(1-exp(-myscat*n))/(1-exp(-myscat))/ ((freq(j)/freq_ref)^alpha);

  %if j==nfreq, disp([j freq(j) myscat]);end;
  %scatvec=scatvec/scatvec(1);
  scatvec=scatvec/mynorm;
  lagvec=exp(-2*pi*i*nu*mylag/n);
  mat(:,j)=sigvec.*scatvec.*lagvec;
  if (j==2)
    disp([mylag myscat mynorm]);
    disp([aa(2) lagvec(2) scatvec(2)])
    %return
  end


end
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
aa=-scat-2*pi*i*x/n;yft=yft.*(1-exp(n*aa))./(1-exp(aa));
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



