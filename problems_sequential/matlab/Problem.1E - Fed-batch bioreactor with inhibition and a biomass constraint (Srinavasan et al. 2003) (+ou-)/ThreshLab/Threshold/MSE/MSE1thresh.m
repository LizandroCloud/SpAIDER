
function mse = MSE1thresh(v,w,t,softORhard)

% MSE1thresh -- apply soft/hard threshold + compute mean squared error (MSE)
%  Usage
%    mse = MSE1thresh(v,w,t,softORhard)
%  Inputs
%    v      (1D) noise-free input coefficients
%    w      (1D) observed input coefficients
%    t      threshold (may be a vector)
%    softORhard boolean: 1 for soft; 0 for hard
%  Outputs
%    mse    the value of the mse
%  Description
%    Applies soft- or hard threshold t to all elements of w and computes MSE
%  Note
%    For minimizing MSEthresh, use minmsesoft or minmsehard

if nargin < 4, softORhard = 1; end
vcol = column(v);
wcol = column(w);
trow = row(t);
N = length(v);
m = length(t);
t = repmat(trow,N,1);
v = repmat(vcol,1,m);
w = repmat(wcol,1,m);
s = abs(w) - t; s = (s > 0);
if softORhard,
   wt = sign(w) .* (abs(w) - t) .* s;
else
   wt = w.*s;
end
mse = sum((wt-v).*(wt-v))/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
