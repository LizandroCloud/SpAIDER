
function w = FFMF(x,L);

% usage: w = FFMF(x,L);
% x: input date
% L: coarsest scale
% Forward Fast Multiscale Fisz-transform for Poisson data
% L is coarsest resolution level
% cfr FWT_BIO

N = size(x,2);
if N == 1, N = size(x,1); end
J = log2(N);
w = x;
for j = J-1:-1:L,
   s = w(:,2:2:2^(j+1)) + w(:,1:2:2^(j+1));
   d = w(:,2:2:2^(j+1)) - w(:,1:2:2^(j+1));
   s1 = max(s,1);
   d = d./sqrt(s1);
   w(:,1:2^j) = s;
   w(:,2^j+1:2^(j+1)) = d;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
