
function Bs = softbias(l,v,stdev);

% softbias -- bias in soft thresholding data with normal noise
%  Usage
%    Bs = softbias(t,v,stdev);
%  Inputs
%    t      threshold, or vector of thresholds
%    v      vector of noise-free data
%    stdev  standard deviation of normal noise (to be added to v)
%  Outputs
%    Bs     soft threshold bias
%  Description
%    computes bias when thresholding w = v+noise where noise ~ N(0,stdev^2)
%    i.e., E(ST(w)) - v where ST = soft-threshold function
%  Note
%  See also:
%    help softrisk
%    help ST

vcol = column(v);
lrow = row(l);
N = length(v);
m = length(l);
l = repmat(lrow,N,1);
v = repmat(vcol,1,m);

Eom1 = stdev^2*(gauss(l-v,stdev) - gauss(l+v,stdev));
Eom2 = l.*(cumgauss(l-v,stdev) - cumgauss(l+v,stdev));
Eom3 = v.*(1 - cumgauss(l-v,stdev) - cumgauss(l+v,stdev));
Eom = Eom1 + Eom2 + Eom3;
Bs = sum(Eom.*Eom)/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
