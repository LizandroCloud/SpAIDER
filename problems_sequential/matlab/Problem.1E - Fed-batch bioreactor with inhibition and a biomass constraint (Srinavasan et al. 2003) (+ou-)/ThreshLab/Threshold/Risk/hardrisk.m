
function Rh = hardrisk(t,v,stdev);

% hardrisk -- expected L2-loss in hard thresholding data with normal noise
%  Usage
%    Rh = hardrisk(t,v,stdev);
%  Inputs
%    t      threshold, or vector of thresholds
%    v      vector of noise-free data
%    stdev  standard deviation of normal noise (to be added to v) 
%  Outputs
%    Rh     hard threshold risk = expected L2-loss
%  Description
%    computes risk when thresholding w = v+noise where noise ~ N(0,stdev^2)
%  Note
%    This routine cannot be used for plotting the risk in a single coefficient
%    as a function of that coefficient's value: the function takes all values
%    of v as a data vector.
%  See also:
%    help softrisk

if nargin < 3, stdev = 1; end
if nargin < 2, v = 0; end
vcol = column(v);
trow = row(t);
N = length(v);
m = length(t);
t = repmat(trow,N,1);
v = repmat(vcol,1,m);
sigsig = stdev*stdev;
tv = t+v;
t_v = t-v;

Rh1 = sigsig*(tv).*gauss(tv,stdev)+sigsig*t_v.*gauss(t_v,stdev);
Rh2 = sigsig*(1-cumgauss(tv,stdev))+sigsig*(1-cumgauss(t_v,stdev));
Rh3 = v.*v.*(cumgauss(tv,stdev)+cumgauss(t_v,stdev)-1);
Rh = Rh1+Rh2+Rh3;

Rh = sum(Rh,1)/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
