
function sure = SUREsoft(x,t,stdev)

% SUREsoft -- apply soft threshold + compute Stein's unbiased risk estimator 
%                                                 (SURE)
%  Usage
%    sure = SUREsoft(x,t,stdev);
%    sure = SUREsoft(x,t);
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%    stdev  noise standard deviation (default is 1)
%  Outputs
%    sure    the value of the SURE, using soft thresholding
%  Description
%    Applies soft-threshold t to all elements of x and computes SURE
%  Note
%    For making a plot of sure for a given vector of thresholds
%    use SURE1soft (for one dimension)
%
%    For minimizing SUREsoft, use minsuresoft (1D) 
%  See also:
%    help minsuresoft
%    help SUREhard
%    help GCVsoft
%    help SURE1soft

if nargin == 2,
   stdev = 1;
end

N = length(x);
s = abs(x) - t; s = (s > 0);
N1 = sum(s);
y = sign(x).*(abs(x) - t).*s;

sure = norm(y-x)^2/N - stdev^2 + 2*stdev^2*N1/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
