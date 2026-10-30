
function x = IRMF(w,alfa);

% usage: x = IRMF(w,alfa);
% or   : x = IRMF(w);
% Inverse Redundant Multiscale Fisz transform for Poisson data; 
% alfa is weighting factor between different possible reconstructions
%
% see also IRT_BIO, FRMF, FRHFW

if nargin < 2,
   alfa = 0.5;
end
N = size(w,2);
M = size(w,1)-1;
J = log2(N);
x = w(M+1,:);
for j = M:-1:1,
   s = x;
   d = w(j,:);
   d = d.*sqrt(abs(s));
   for k=1:2^(j-1),
      sjk = s(k:2^j:N);
      djk = d(k:2^j:N);
      x1(k:2^j:N)          = (sjk - djk )/2;
      x1(k+2^(j-1):2^j:N)  = (sjk + djk )/2;
      sjk = s(k+2^(j-1):2^j:N);
      djk = d(k+2^(j-1):2^j:N);
      x2(k:2^j:N)          = (sjk - djk )/2;
      x2(k+2^(j-1):2^j:N)  = (sjk + djk )/2;
   end
   x2 = shift(x2,-2^(j-1));
   x = alfa*x1 + (1-alfa)*x2;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
