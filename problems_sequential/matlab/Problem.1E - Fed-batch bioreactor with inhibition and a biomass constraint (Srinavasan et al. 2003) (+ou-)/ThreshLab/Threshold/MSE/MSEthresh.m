
function mse = MSEthresh(v,w,t,softORhard)

% MSEthresh -- apply soft/hard threshold + compute mean squared error (MSE)
%  Usage
%    mse = MSEthresh(v,w,t,softORhard)
%  Inputs
%    v      noise-free input coefficients
%    w      observed input coefficients
%    t      threshold (must be scalar)
%    softORhard boolean: 1 for soft; 0 for hard
%  Outputs
%    mse    the value of the mse
%  Description
%    Applies soft- or hard threshold t to all elements of w and computes MSE
%  Note
%    For making a plot of mse for a given vector of thresholds
%    use MSE1thresh (for one dimension)
%
%    For minimizing MSEthresh, use minmsethreshest

if nargin < 4, softORhard = 1; end
N = size(w,1)*size(w,2);
if softORhard,
   wt = ST(w,t);
else
   wt = HT(w,t);
end
mse = sum(sum((wt-v).*(wt-v)))/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
