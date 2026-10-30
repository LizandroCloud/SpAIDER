
function thr = minriskthr(v,stdev);

% usage: thr = minriskthr(v,stdev);
% finds minimum risk soft threshold for noise-free coefficients
% v is vector of uncorrupted coeff. Corresponding empirical coefficients have
% noise deviation equal to stdev.

if nargin == 1, stdev = 1; end
eps = 1.e-10;
maxiter = 5000;
thmax = max(abs(v));
N = length(v);
t1 = max(1,sqrt(abs(2*log(N))))*stdev;
t0 = 0;
iter = 0;
while (abs(t0 - t1) > eps)
 t0 = t1;
 iter = iter + 1;
 teller = gauss(t1 - v, stdev) + gauss(t1 + v,stdev);
 noemer = Rcumgauss(t1 - v,stdev) + Rcumgauss(t1 + v,stdev);
 if (iter < maxiter)
   t1 = stdev^2*sum(teller)/sum(noemer);
 end    
end
 iter
thr = t1;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
