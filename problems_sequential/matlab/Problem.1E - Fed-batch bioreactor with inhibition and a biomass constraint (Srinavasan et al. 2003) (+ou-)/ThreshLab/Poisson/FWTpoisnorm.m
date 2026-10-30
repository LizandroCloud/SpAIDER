
function [pn, first, last] = FWTpoisnorm(x,L,ht,gt,sht,sgt)

% usage: pn = FWTpoisnorm(x,L,ht,gt,sht,sgt);
% or: [pn, first, last] = FWTpoisnorm(x,L,ht,gt,sht,sgt);
%
% computes appropriate normalisation for wavelet coefficients of Poisson noise
% vector x. The parameters L,ht,gt,sht,sgt correspond to those in the function
% call wp = FWT_BIO(x,L,ht,gt,sht,sgt);
% wpn = wp./pn; is then the vector of Poisson normalised wavelet coefficients
%
% first is index into first sample in x participating in the normalisation
% factor; similar for last (see FWTdualsupport)
% 
% See also: FWTdualsupport

dotranspose = 0;
if size(x,2) == 1, x = x'; dotranspose = 1; end
[first, last] = FWTdualsupport(x(1,:),L,ht,gt,sht,sgt);
N = size(x,2);
pn = zeros(size(x));
for i = 1:N,
   if first(i) > last(i),
      pn(:,i) = sum(x(:,first(i):N),2) + sum(x(:,1:last(i)),2);
   else
      pn(:,i) = sum(x(:,first(i):last(i)),2);
   end
end
pn = sqrt(abs(pn));
if dotranspose, pn = pn'; end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
