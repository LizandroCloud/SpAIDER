function wcoef = RWT2_BIO(x,M,ht,gt,sht,sgt)

% usage: wcoef = RWT2_BIO(x,M,ht,gt,sht,sgt)
% Non-Decimated 2D Wavelet Transform
% wcoef is a 4-dimensional tensor. wcoef(j,i,:,:) has the wavelet coefficients
% at scale j (j=1..M), orientation i (i=1,2,3 i.e. horizontal, vertical,
% diagonal)

[n m] = size(x);
w = zeros(M+1,3,n,m);
s = x;
for j=1:M
  whor = RWT_1step(s,j,gt,sgt);
  shor = RWT_1step(s,j,ht,sht);
  w(j,1,:,:) = RWT_1step(shor',j,gt,sgt)';
  w(j,2,:,:) = RWT_1step(whor',j,ht,sht)';
  w(j,3,:,:) = RWT_1step(whor',j,gt,sgt)';
  s = RWT_1step(shor',j,ht,sht)';
end
w(M+1,1,:,:) = s;
wcoef = w;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
