
function w = FRMF(x,M);

% usage: w = FRMF(x,M);
% Forward Redundant Multiscale Fisz transform for Poisson data; 
% M is number of resolution levels
% cfr RWT_BIO

x = row(x);
N = length(x);
J = log2(N);
w = zeros(M+1,N);
w(M+1,:) = row(x);
for j = 0:M-1,
   x = w(M+1,:);
   x = [x x(1:2^j)];
   for k=1:2^j,
      s(k:2^j:N) = x(k+2^j:2^j:N+2^j) + x(k:2^j:N);
      d(k:2^j:N) = x(k+2^j:2^j:N+2^j) - x(k:2^j:N);
   end
   s1 = max(s,1);
   d = d./sqrt(s1);
   w(M+1,:) = s;
   w(j+1,:) = d;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
