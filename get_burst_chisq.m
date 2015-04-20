function[chisq]=get_burst_chisq(guess,datft,freqs,n,nvec,dt,imin)

fwhm=guess(1);
scat=guess(2);
alpha=guess(3);
amp=guess(4);
t0=guess(5);
DM=guess(6);
if ~exist('dt')
  dt=1e-3;
end

mod=amp*make_burst_model_c(dt,n,freqs,fwhm,scat,alpha,t0,DM);

if exist('imin')
  v1=sum(real(mod(imin+1:end,:)).*conj(datft(imin+1:end,:)))
  v2=sum(real(mod(imin+1:end,:).*conj(mod(imin+1:end,:))));
  best_amp=sum(v1.*nvec)/sum(v2.*nvec)*amp;
  chisq=sum(abs(mod(imin+1:end,:)-datft(imin+1:end,:)).^2);disp(['chisq(11) is ' num2str(chisq(11))]);chisq=sum(chisq.*nvec);
else
  chisq=sum(abs(mod-datft).^2);chisq=sum(chisq.*nvec);
end

