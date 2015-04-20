function[dat2,nuvec,mylags]=stack_data(dat,freqs,params,dt,highpass)
if exist('highpass')
  datft=fft(dat);
  dnu=1/(length(dat)*dt);
  nuvec=(1:length(dat))'-1;
  nuvec(floor(end/2+1):end)=nuvec(floor(end/2+1):end)-length(dat);
  nuvec=nuvec*dnu;
  tic;datft=fft(dat);toc;
  datft(abs(nuvec)<highpass,:)=0;
  dat=real(ifft(datft));
end

DM=params(6);DM0=4150;
mylags=DM*DM0./(freqs.^2);
mylags=(mylags-min(mylags));
mylags=round(mylags/dt);


targ=1+round(params(5)/dt);
pad=1000;
imin=targ-pad;
imax=targ+pad+max(mylags);
dat2=dat(imin:imax,:);


whos 


for j=1:length(mylags),
  dat2(:,j)=circshift(dat2(:,j),-round(mylags(j)));
end


