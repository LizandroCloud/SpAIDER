
function Vsqrt = VsqrtPoisson(lambda,c);

% usage: Vsqrt = VsqrtPoisson(lambda,c);
% computes V(sqrt(X)), where X ~ Poisson(lambda)
% default for c is 0

if nargin == 1,
   c = 0;
end
Vsqrt = lambda + c - EsqrtPoisson(lambda,c).^2;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
