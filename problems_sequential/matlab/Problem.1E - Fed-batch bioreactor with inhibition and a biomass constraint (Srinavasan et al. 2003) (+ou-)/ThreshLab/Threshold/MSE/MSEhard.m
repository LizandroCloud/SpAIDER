
function mse = MSEhard(v,w,t)

% MSEhard -- apply hard threshold + compute mean squared error (MSE)
%  Usage
%    mse = MSEhard(v,w,t)
%  Inputs
%    v      noise-free input coefficients
%    w      observed input coefficients
%    t      threshold (must be scalar)
%  Outputs
%    mse    the value of the mse
%  Description
%    Applies hard threshold t to all elements of w and computes MSE
%  Note
%    For making a plot of mse for a given vector of thresholds
%    use MSE1thresh (for one dimension)
%  See also
%    For minimizing MSEhard, use minmsehard
%    help MSEthresh
%    help MSE1thresh

mse = MSEthresh(v,w,t,0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
