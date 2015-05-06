function[pp,ll,qft,uft,nvec]=mcmc_burst_polfit_rmpow(dat_q,dat_u,freqs,guess,sigs,dt,nu_min,fix_params,true_params)


if ~exist('dt')
  dt=1e-3;
end
if ~exist('nu_min')
  imin=1;
else
  dnu=1/(dt*length(dat_q));
  imin=ceil(nu_min/dnu)
end

qft=fft(dat_q);
n=size(dat_q,1);
nn=ceil((n+1)/2);
qft=qft(1:nn,:);

uft=fft(dat_u);
uft=uft(1:nn,:);


nstep=300000;
pp=zeros(nstep,numel(guess));
ll=zeros(nstep,1);
cur=guess;
nvec=0.5./mean(abs(qft(imin:end,:)).^2);
nvec2=0.5./mean(abs(uft(imin:end,:)).^2);
disp(['error consistency is ' num2str([mean(abs(nvec-nvec2))/mean(abs(nvec))])]);
nmat=repmat(nvec,[n 1]);


if exist('true_params')
  guess(fix_params)=true_params(fix_params);
end


try
  myid=mpi_comm_rank;
  fid=fopen(['chain_polfit_rmpow_v2.txt_' num2str(myid)],'w');
catch
  fid=fopen('chain_polfit_rmpow_v2.txt','w');
end


%fwhm=guess(1);
%scat=guess(2);
%alpha=guess(3);
%amp=guess(4);
%t0=guess(5);
%DM=guess(6);
%mod=amp*make_burst_model(dt,n,freqs,fwhm,scat,alpha,t0,DM);
%chisq=sum(abs(mod-dataft).^2);chisq=sum(chisq.*nvec)

tic;chisq=get_burst_chisq_qu_rmpow_c(guess,qft,uft,freqs,n,nvec,dt,imin);toc


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

  if ((guess(2)<0)|(guess(1)<0))
    chisq=ll(j-1)+1e5;
  else
    chisq=get_burst_chisq_qu_rmpow_c(guess,qft,uft,freqs,n,nvec,dt,imin);
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
  if rem(j,10)==0
    fflush(fid);
  end

end
fclose(fid);
