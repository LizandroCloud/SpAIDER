
function mse = MSEsoft(v,w,t)

% MSEsoft -- apply soft threshold + compute mean squared error (MSE)
%  Usage
%    mse = MSEsoft(v,w,t)
%  Inputs
%    v      noise-free input coefficients
%    w      observed input coefficients
%    t      threshold (must be scalar)
%  Outputs
%    mse    the value of the mse
%  Description
%    Applies soft threshold t to all elements of w and computes MSE
%  Note
%    For making a plot of mse for a given vector of thresholds
%    use MSE1thresh (for one dimension)
%  See also
%    For minimizing MSEsoft, use minmsesoft 
%    help MSEthresh
%    help MSE1thresh

mse = MSEthresh(v,w,t,1);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
