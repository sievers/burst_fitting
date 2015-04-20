function[pp,ll,dataft,nvec]=mcmc_burst_ampfit(data,freqs,guess,sigs,dt,nu_min,fix_params,true_params)


if ~exist('dt')
  dt=1e-3;
end
if ~exist('nu_min')
  imin=1;
else
  dnu=1/(dt*length(data));
  imin=ceil(nu_min/dnu)
end

dataft=fft(data);
n=size(data,1);
nn=ceil((n+1)/2);
dataft=dataft(1:nn,:);


nstep=30000;
pp=zeros(nstep,numel(guess));
ll=zeros(nstep,1);
cur=guess;
nvec=0.5./mean(abs(dataft(imin:end,:)).^2);
nmat=repmat(nvec,[n 1]);


if exist('true_params')
  guess(fix_params)=true_params(fix_params);
end


fid=fopen('chain_ampfit.txt','w');

%fwhm=guess(1);
%scat=guess(2);
%alpha=guess(3);
%amp=guess(4);
%t0=guess(5);
%DM=guess(6);
%mod=amp*make_burst_model(dt,n,freqs,fwhm,scat,alpha,t0,DM);
%chisq=sum(abs(mod-dataft).^2);chisq=sum(chisq.*nvec)

tic;[chisq,amp]=get_burst_chisq_ampfit_c(guess,dataft,freqs,n,nvec,dt,imin);toc
guess(4)=amp;

pp(1,:)=guess;
ll(1)=chisq;


for j=2:nstep
  if rem(j,50)==0
    disp(j);
  end
  guess=cur+mc_step(sigs);
  if exist('true_params')
    guess(fix_params)=true_params(fix_params);
  end

  %fwhm=guess(1);
  %scat=guess(2);
  %alpha=guess(3);
  %amp=guess(4);
  %t0=guess(5);
  %DM=guess(6); 
  %mod=amp*make_burst_model(dt,n,freqs,fwhm,scat,alpha,t0,DM);
  %chisq=sum(abs(mod-dataft).^2);chisq=sum(chisq.*nvec);

  %if (guess(1)>0)
  %  [chisq,amp]=get_burst_chisq_ampfit_c(guess,dataft,freqs,n,nvec,dt,imin);
  %  guess(4)=amp;
  %else
  %  chisq=ll(j-1)+1e5;
  %end

  if ((guess(2)<0)|(guess(1)<0))
    chisq=ll(j-1)+1e5;
  else
    [chisq,amp]=get_burst_chisq_ampfit_c(guess,dataft,freqs,n,nvec,dt,imin);
    guess(4)=amp;    
  end
  
  
  if (chisq<ll(j-1))
    disp('should be accepting')
  end


  if (exp(-0.5*(chisq-ll(j-1)))>rand(1))
    ll(j)=chisq;
    pp(j,:)=guess;
    cur=guess;
  else
    ll(j)=ll(j-1);
    pp(j,:)=cur;
  end
  fprintf(fid,'%14.7f ',ll(j));
  fprintf(fid,'%14.8g ',pp(j,:));
  fprintf(fid,'\n');
  fflush(fid);

end
fclose(fid);
