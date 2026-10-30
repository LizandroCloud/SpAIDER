
function sure = SURE2soft(x,t,stdev)

% SURE2soft-- apply soft threshold + compute SURE in 2D
%  Usage
%    sure = SURE2soft(x,t,stdev)
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%    stdev  noise standard deviation
%  Outputs
%    sure    the value of the SURE, using soft thresholding
%  Description
%    Applies soft-threshold t to all elements of x and computes SURE
%    SURE = Stein's Unbiased Risk Estimator
%  See also
%    help SUREsoft

[m,n] = size(x);
N = m*n;
y = x;
s = abs(x) - t; s = (s > 0);
N1 = sum(s);
y = sign(x).*(abs(x) - t).*s;

sure = sum(sum((y-x).^2))/N - stdev^2 + 2*stdev^2*N1/N;



% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
