
function dlnLdp = dlnLmixednormaldp(p,x,sig0,sig1);

% usage: dlnLdp = dlnLmixednormaldp(p,x,sig0,sig1);
% derivative of the log-likelikihood as a function of parameter p, for 
% observations x from density (1-p)*phi0(x) + p*phi1(x)
% sig0 and sig1 are standard deviations of normal densities phi0 and phi1
% sig0 and/or sig1 are either scalars or vectors of length length(x)

if length(sig0) == 1,
   sig0 = sig0*ones(size(x));
end
if length(sig1) == 1,
   sig1 = sig1*ones(size(x));
end
maxx = 10*sig0;
xx = max(min(x,maxx),-maxx);
xinf = find(abs(xx)>maxx);
beta = gauss(xx./sig1) ./ gauss(xx./sig0) .* sig0./sig1 - 1;
for i = 1:length(p),
   dlnLdpterms = beta./(1+p(i)*beta);
   if p(i) < eps, 
      dlnLdpterms(xinf) = Inf;
   else
      dlnLdpterms(xinf) = 1/p(i);
   end
   dlnLdp(i) = sum(dlnLdpterms);
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
