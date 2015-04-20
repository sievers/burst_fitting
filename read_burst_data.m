function[dat]=read_burst_data(fname);
if ~exist('fname')
  fname='calibrated.npy';
end

fid=fopen(fname,'r');header=fread(fid,[1 96],'char=>char');
ii=find(header=='(');
jj=find(header==')');
dims=str2num([ '[' header(ii+1:jj-1) ']' ])
dat=fread(fid,inf,'float=>float');
assert(numel(dat)==prod(dims));
dat=reshape(dat,fliplr(dims));
fclose(fid);
