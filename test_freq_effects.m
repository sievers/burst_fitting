if ~exist('crap')
  crap=load('chains/chain_tt_real_test3.txt_1');crap=crap(round(0.2*end):end,3:end);
  gg=mean(crap);
end
dt=1e-3;n=50000;nu=900;dnu=200/4096;
mm=make_burst_model_real_c(dt,n,nu,gg(1),gg(2),0.0,gg(5),gg(6));
nuvec=nu+[-0.5 0 0.5]*dnu;

tvec=-2e-5:1e-6:2e-5;tvec=tvec';
tvec=0;
chivec=0*tvec;
for j=1:length(tvec),
  mm2=make_burst_model_real_c(dt,n,nuvec,gg(1),gg(2),0.0,gg(5)+tvec(j),gg(6));
  beta=0.23*(1-8/pi^2)+0.54/6;
  alpha=0.46*8/pi^2 + 0.54*2/3;
  mm3=mm2(:,2)*alpha+beta*(mm2(:,1)+mm2(:,3));
  chivec(j)=sum(abs(mm3-mm));
end
[tvec chivec]
