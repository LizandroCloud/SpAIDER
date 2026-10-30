
function sure = SUREhard(x,t,stdev)

% SUREhard -- apply hard threshold + compute Stein's unbiased risk estimator 
%                                                 (SURE)
%  Usage
%    sure = SUREhard(x,t,stdev);
%    sure = SUREhard(x,t);
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%    stdev  noise standard deviation (default is 1)
%  Outputs
%    sure    the value of the SURE, using hard thresholding
%  Description
%    Applies hard-threshold t to all elements of x and computes SURE
%  Note
%    For making a plot of sure for a given vector of thresholds
%    use SURE1hard (for one dimension)
%
%    For minimizing SUREhard, use minsurehard (1D)
%  See also:
%    help minsurehard
%    help SUREsoft
%    help GCVhard
%    help SURE1hard

if nargin == 2,
   stdev = 1;
end

N = length(x);
s = abs(x) - t; s = (s > 0);
N1 = sum(s);
y = sign(x).*(abs(x) - t).*s;

sure = norm(y-x)^2/N - stdev^2 + 2*stdev^2*N1/N ...
     + 4*stdev^2*t(1,:).*gauss(t(1,:),stdev);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
