function[pp,ll,qft,uft,nvec,imin]=mcmc_burst_pol2(dat_q,dat_u,myopts,freqs,guess,sigs,dt,nu_min,fix_params,true_params)
 
outroot=get_struct_mem(myopts,'outroot','chain_pol.txt');
nstep=get_struct_mem(myopts,'nstep',3);
chifun=get_struct_mem(myopts,'func',@get_burst_chisq_qu_cached_c);


if ~exist('dt')
  dt=1e-3;
end
if ~exist('nu_min')
  imin=1;
else
  dnu=1/(dt*length(dat_q));
  imin=ceil(nu_min/dnu);
end

qft=fft(dat_q);
n=size(dat_q,1);
nn=ceil((n+1)/2);
qft=qft(1:nn,:);
uft=fft(dat_u);
uft=uft(1:nn,:);



pp=zeros(nstep,numel(guess));
ll=zeros(nstep,1);
cur=guess;

noises=0.5*mean(abs(qft(imin:end,:)).^2);
nvec=1.0./noises;
noises2=0.5*mean(abs(uft(imin:end,:)).^2);
nvec2=1./noises2;
disp(['noise agreement is ' num2str(sqrt(mean((nvec2-nvec).^2)/mean(nvec.^2)))]);
unwind_protect

qptr=floatcomplex2ptr(qft);
uptr=floatcomplex2ptr(uft);
model=floatcomplex2ptr(qft);



if exist('true_params')
  guess(fix_params)=true_params(fix_params);
end


try
  myid=mpi_comm_rank+1;
  fid=fopen([outroot '_'  num2str(myid)],'w');
catch
  fid=fopen(outroot,'w');
end


noamp=guess;noamp(4)=0;
chi_ref=feval(chifun,noamp,qptr,uptr,freqs,n,nvec,dt,imin,model);
chi_ref=get_struct_mem(myopts,'chi_ref',chi_ref);disp(['chi_ref is ' sprintf('%14.6f',chi_ref)]);
tic;chisq=feval(chifun,guess,qptr,uptr,freqs,n,nvec,dt,imin,model)-chi_ref;toc



pp(1,:)=guess;
ll(1)=chisq;
nrep=1;

for j=2:nstep

  if rem(j,50)==0
    disp(j);
    fflush(fid);
  end
  guess=cur+mc_step(sigs);
  if exist('true_params')
    guess(fix_params)=true_params(fix_params);
  end

  chisq=feval(chifun,guess,qptr,uptr,freqs,n,nvec,dt,imin,model)-chi_ref;

  %if (chisq<ll(j-1))
  %  disp('should be accepting')
  %end


  if (exp(-0.5*(chisq-ll(j-1)))>rand(1))
    %disp('accepting');
    fprintf(fid,'%d ',nrep);
    fprintf(fid,'%14.7f ',ll(j-1));
    fprintf(fid,'%15.9g ',pp(j-1,:));
    fprintf(fid,'\n');
    %fflush(fid);


    ll(j)=chisq;
    pp(j,:)=guess;
    cur=guess;
    nrep=1;
  else
    %disp('not accepting')
    ll(j)=ll(j-1);
    pp(j,:)=cur;
    nrep=nrep+1;
  end

end
unwind_protect_cleanup
disp('cleaning up')
freeptr(model);
freeptr(qptr);
freeptr(uptr);
fclose(fid);
end_unwind_protect

