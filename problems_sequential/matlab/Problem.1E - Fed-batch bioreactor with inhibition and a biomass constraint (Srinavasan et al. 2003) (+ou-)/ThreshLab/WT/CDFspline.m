function [h, ht] = CDFspline(n,nt)

%   Usage
%     [h ht] = CDFspline(n);
%   Inputs
%     n    number of primal vanishing moments
%     nt   number of dual vanishing moments
%   Outputs
%     h    CDFspline filter coefficient
%     ht   dual CDFspline filter coefficient (bi-orthogonal transform)
%   Description
%     creates dual and primal scaling filter for bi-orthogonal spline wavelet
%     of Cohen-Daubechies-Feauveau (CDF)
%   Note
%     This file called by makeWTfilter
%
%     n+nt must be even
%
%     what is called "dual vanishing moments" and what "primal vanishing
%     moments" is matter of discussion.
%     In this routine nt refers to properties of ht (the dual two-scale
%     coefficients).
%     On the other hand, if the the "dual vanishing moments" are the vanishing
%     moments of the dual wavelet basis, then nt should correspond to the
%     design of gt (the dual wavelet equation), and hence to h. If that is what
%     you have in mind, then switch the roles of n and nt!!
%
%   See also
%     help makeWTfilter

nn = n + nt;

if rem(nn,2) 
  error('n + nt must be even');
end
nn = nn/2;

for i = 1:n+1
  h(i) = nchoosek(n,i-1);
end
h = h/2^n*sqrt(2);

gt2 = [-1 2 -1]/4;
ht2 = gt2;
ht = 1;
pfac = 1;
for p = 1:nn-1,
  pfac = pfac*(nn+p-1)/p;
  ht = [0,ht,0] + pfac*ht2;
  ht2 = conv(ht2,gt2);
end

for k = 1:nt
  ht = conv(ht,[0.5 0.5]);
end

ht = ht/sum(ht)*sqrt(2);

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
