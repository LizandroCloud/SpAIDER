
function Vsqrt = VsqrtPoissonN0(lambda);

% usage: Vsqrt = VsqrtPoissonN0(lambda);
% computes V(sqrt(X) | X > 0), where X ~ Poisson(lambda)

Vsqrt = lambda./(1-exp(-lambda)) - EsqrtPoisson(lambda).^2./(1-exp(-lambda)).^2;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
