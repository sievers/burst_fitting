function[pp,ll,dataft,nvec,imin]=mcmc_burst_intensity(data,myopts,freqs,guess,sigs,dt,nu_min,fix_params,true_params)

outroot=get_struct_mem(myopts,'outroot','chain_tt.txt');
nstep=get_struct_mem(myopts,'nstep',3);
chifun=get_struct_mem(myopts,'func',@get_burst_chisq_cached_c);


if ~exist('dt')
  dt=1e-3;
end
if ~exist('nu_min')
  imin=1;
else
  dnu=1/(dt*length(data));
  imin=ceil(nu_min/dnu);
end

dataft=fft(data);
n=size(data,1);
nn=ceil((n+1)/2);
dataft=dataft(1:nn,:);

pp=zeros(nstep,numel(guess));
ll=zeros(nstep,1);
cur=guess;

noises=0.5*mean(abs(dataft(imin:end,:)).^2);
nvec=1.0./noises;

unwind_protect

dataptr=floatcomplex2ptr(dataft);
model=floatcomplex2ptr(dataft);



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
chi_ref=feval(chifun,noamp,dataptr,freqs,n,nvec,dt,imin,model);
chi_ref=get_struct_mem(myopts,'chi_ref',chi_ref);disp(['chi_ref is ' sprintf('%14.6f',chi_ref)]);
tic;chisq=feval(chifun,guess,dataptr,freqs,n,nvec,dt,imin,model)-chi_ref;toc




pp(1,:)=guess;
ll(1)=chisq;
nrep=1;

for j=2:nstep

  if rem(j,50)==0
    disp(j);
  end
  guess=cur+mc_step(sigs);
  if exist('true_params')
    guess(fix_params)=true_params(fix_params);
  end

  chisq=feval(chifun,guess,dataptr,freqs,n,nvec,dt,imin,model)-chi_ref;

  if (chisq<ll(j-1))
    disp('should be accepting')

  end


  if (exp(-0.5*(chisq-ll(j-1)))>rand(1))
    disp('accepting');
    fprintf(fid,'%d ',nrep);
    fprintf(fid,'%14.7f ',ll(j-1));
    fprintf(fid,'%15.9g ',pp(j-1,:));
    fprintf(fid,'\n');
    fflush(fid);


    ll(j)=chisq;
    pp(j,:)=guess;
    cur=guess;
    nrep=1;
  else
    disp('not accepting')
    ll(j)=ll(j-1);
    pp(j,:)=cur;
    nrep=nrep+1;
  end

end
unwind_protect_cleanup
disp('cleaning up')
freeptr(model);
freeptr(dataptr);
fclose(fid);
end_unwind_protect

