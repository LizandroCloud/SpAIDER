
function lnL = lnLmixednormal(p,x,sig0,sig1);

% usage: lnL = lnLmixednormal(p,x,sig0,sig1);
% log-likelikihood as a function of parameter p, for observations x from
% density (1-p)*phi0(x) + p*phi1(x)
% sig0 and sig1 are standard deviations of normal densities phi0 and phi1
% sig0 and/or sig1 are either scalars or vectors of length length(x)

if length(sig0) == 1,
   sig0 = sig0*ones(size(x));
end
if length(sig1) == 1,
   sig1 = sig1*ones(size(x));
end
for i = 1:length(p),
   lnL(i) = sum(log(p(i)*gauss(x./sig1)./sig1+(1-p(i))*gauss(x./sig0)./sig0));
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
