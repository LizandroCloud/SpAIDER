
function x = IFMF(w,L);

% usage: x = IFMF(w,L);
% Inverse Fast Multiscale Fisz transform for Poisson data
% L is coarsest resolution level
% 
% see also: IWT_BIO.m and FFMF.m

N = size(w,2);
if N == 1, N = size(w,1); end
J = log2(N);
x = w;
for j = L:J-1,
   d = x(:,2^j+1:2^(j+1));
   s = x(:,1:2^j);
   d = d.*sqrt(abs(s));
   x(:,2:2:2^(j+1)) = (s+d)/2;
   x(:,1:2:2^(j+1)) = (s-d)/2;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
