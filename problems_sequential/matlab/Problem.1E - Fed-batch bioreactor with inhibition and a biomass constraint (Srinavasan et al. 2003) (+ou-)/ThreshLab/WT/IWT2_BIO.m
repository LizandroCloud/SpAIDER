function f = IWT2_BIO(w,L,h,g,sh,sg)

% f = IWT2_BIO(w,L,h,g,sh,sg)
%
% inverse 2D wavelet transform
% see IWT_BIO for more details on parameters
% dyadic square sized images only

N = size(w,1);
M = size(w,2);
J = log2(N);
n = 2^L;
m = n*M/N;
for j=L+1:J
  n = 2*n;
  m = 2*m;
  w(1:n,1:m) = IWT2_1step(w(1:n,1:m),h,g,sh,sg);
end
f = w;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
