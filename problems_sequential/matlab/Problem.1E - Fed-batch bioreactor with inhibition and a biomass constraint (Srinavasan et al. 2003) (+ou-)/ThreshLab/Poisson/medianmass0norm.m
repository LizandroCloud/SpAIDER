
function vhat = medianmass0norm(p,mu1,sig1);

% usage: vhat = medianmass0norm(p,mu1,sig1);
%
% Computes median of random variable with mixture density: probability mass p
% at origin + (1-p) normal density
%
% mu1: mean of random variable given that it is nonzero
% p: probability that this random variable is nonzero
% sig1: standard deviation of nonzero r.v. (heavy tail)
%
% see also CDFmass0norm

mu1sign = sign(mu1);
mu1 = abs(mu1);
v0 = 0;
v1 = mu1;
eps = 0.001*mu1;
F0 = CDFmass0norm(v0,p,mu1,sig1);
if (F0 > 0.5)
   vhat = 0;
else
   while (v1-v0)>eps
      v = (v0+v1)/2;
      F = CDFmass0norm(v,p,mu1,sig1);
      if (F > 0.5)
         v1 = v;
      else
         v0 = v;
      end
   end
vhat = mu1sign*(v0+v1)/2;
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
