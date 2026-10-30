
function Rs = softrisk(t,v,stdev);

% softrisk -- expected L2-loss in soft thresholding data with normal noise
%  Usage
%    Rs = softrisk(t,v,stdev);
%  Inputs
%    t      threshold, or vector of thresholds
%    v      vector of noise-free data
%    stdev  standard deviation of normal noise (to be added to v)
%  Outputs
%    Rs     soft threshold risk = expected L2-loss
%  Description
%    computes risk when thresholding w = v+noise where noise ~ N(0,stdev^2)
%  Note
%    This routine cannot be used for plotting the risk in a single coefficient
%    as a function of that coefficient's value: the function takes all values
%    of v as a data vector.
%  See also:
%    help hardrisk
%    help softbias
%    help minrisksoftthr


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

Rs1 = 2*(sigsig+t.*t)-v.^2;
Rs2 = (cumgauss(tv,stdev)+cumgauss(t_v,stdev)).*(v.^2-sigsig-t.*t);
Rs3 = -sigsig*(t_v.*gauss(tv,stdev)+tv.*gauss(t_v,stdev));
Rs = Rs1+Rs2+Rs3;

Rs = sum(Rs,1)/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
