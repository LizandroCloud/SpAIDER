
function thr = minriskhardthr(v,stdev);

% minriskhardthr -- finds hard threshold that minimizes expected L2-loss
%  Usage
%    thr = minriskhardthr(v,stdev);
%  Inputs
%    v      vector of noise-free data
%    stdev  standard deviation of normal noise (to be added to v)
%  Outputs
%    thr    optimal threshold w.r.t. expected L2-loss
%  Description
%    finds hard threshold that minimizes expected L2-loss in thresholding
%    w = v+noise where noise ~ N(0,stdev^2)
%  See also
%    help hardrisk


if nargin == 1, stdev = 1; end
eps = 1.e-10;
maxiter = 5000;
thmax = max(abs(v));
N = length(v);
t = sqrt(2*log(N))*stdev/2;
F = -1;
while F < -eps,
   t = 2*t;
   F = 2*sum(v.*(gauss(t - v,stdev) - gauss(t + v,stdev))) -...
       t*sum(gauss(t - v,stdev) + gauss(t + v,stdev));
end
t1 = t;
t0 = 0;
F0 = -1;
F1 = +1;
iter = 0;
while (abs(t0 - t1) > eps)
 t = (t0+t1)/2;
 iter = iter + 1;
 F = 2*sum(v.*(gauss(t - v,stdev) - gauss(t + v,stdev))) -...
     t*sum(gauss(t - v,stdev) + gauss(t + v,stdev));
 if F*F0 > eps,
    F0 = F;
    t0 = t;
 elseif F*F1 > eps,
    F1 = F;
    t1 = t;
 else
    t1 = t0;
 end
 if (iter == maxiter)
    t1 = t0;
 end    
end
 iter
thr = t1;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
