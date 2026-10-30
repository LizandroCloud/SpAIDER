
function Bh = hardbias(l,v,stdev);

% hardbias -- bias in hard thresholding data with normal noise
%  Usage
%    Bh = hardbias(t,v,stdev);
%  Inputs
%    t      threshold, or vector of thresholds
%    v      vector of noise-free data
%    stdev  standard deviation of normal noise (to be added to v)
%  Outputs
%    Bh     hard threshold bias
%  Description
%    computes bias when thresholding w = v+noise where noise ~ N(0,stdev^2)
%    i.e., E(HT(w)) - v where HT = hard-threshold function
%  Note
%  See also:
%    help hardrisk
%    help HT


vcol = column(v);
lrow = row(l);
N = length(v);
m = length(l);
l = repmat(lrow,N,1);
v = repmat(vcol,1,m);

Eom1 = stdev^2*(gauss(l-v,stdev) - gauss(l+v,stdev));
Eom3 = v.*(1 - cumgauss(l-v,stdev) - cumgauss(l+v,stdev));
Eom = Eom1 + Eom3;
Bh = sum(Eom.*Eom)/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
