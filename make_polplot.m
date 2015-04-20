if ~exist('dat_q')
  aaa=now;[dat,dat_q,dat_u,freq,nuvec,mylags]=get_trimmed_burst();bbb=now;disp(86400*(bbb-aaa))
end


nu_ref=800;lref=2.9979e8/(1e6*nu_ref);
  
%rm=103:5:223;
rm=-[200:2:500 ];
qsum=zeros(size(dat,1),numel(rm));
usum=0*qsum;
nfreq=length(freq);
lambda=2.9979e8./(1e6*freq);
for jj=1:length(rm),
  qq=dat_q;
  uu=dat_u;
  for j=1:nfreq,
    qq(:,j)=dat_q(:,j)*cos(rm(jj)*(lambda(j)^2-lref^2))+dat_u(:,j)*sin(rm(jj)*(lambda(j)^2-lref^2));
    uu(:,j)=-dat_q(:,j)*sin(rm(jj)*(lambda(j)^2-lref^2))+dat_u(:,j)*cos(rm(jj)*(lambda(j)^2-lref^2));
  end
  qsum(:,jj)=sum(qq,2);
  usum(:,jj)=sum(uu,2);
end

bb=15;nb=size(dat_u,2)/bb;
u2=zeros(size(dat_u,1),nb);for j=1:nb, u2(:,j)=sum(uu(:,(j-1)*bb+1:j*bb),2);end;
q2=zeros(size(dat_q,1),nb);for j=1:nb, q2(:,j)=sum(qq(:,(j-1)*bb+1:j*bb),2);end;
i2=zeros(size(dat_q,1),nb);for j=1:nb, i2(:,j)=sum(dat(:,(j-1)*bb+1:j*bb),2);end;

