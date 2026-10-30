function wcoef = FWT2_BIO(x,L,ht,gt,sht,sgt)

% wcoef = FWT2_BIO(x,L,ht,gt,sht,sgt)
%
% 2D wavelet transform
% see FWT_BIO for more details
% dyadic square sized images only

w = x;
N = size(x,1);
M = size(x,2);
J = log2(N);
n = N;
m = M;
for j=J:-1:L+1
  w(1:n,1:m) = FWT2_1step(w(1:n,1:m),ht,gt,sht,sgt);
  n = n/2;
  m = m/2;
end
wcoef = w;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
