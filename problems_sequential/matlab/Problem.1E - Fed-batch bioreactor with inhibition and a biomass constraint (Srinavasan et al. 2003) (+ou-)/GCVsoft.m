
function gcv = GCVsoft(x,t)

% GCVsoft -- apply soft threshold + compute Generalized Cross Validation (GCV)
%  Usage
%    gcv = GCVsoft(x,t)
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%  Outputs
%    gcv    the value of the gcv
%  Description
%    Applies soft-threshold t to all elements of x and computes GCV
%  Note
%    For making a plot of gcv for a given vector of thresholds
%    use GCV1soft (for one dimension)
%
%    For minimizing GCVsoft, use mingcvsoft (1D) / mingcv2soft (2D)
%  See also:
%    help mingcvsoft
%    help mingcv2soft
%    help GCV1soft
%    help GCV2soft


N = length(x);
s = abs(x) - t; s = (s > 0);
N0 = N - sum(s);
y = sign(x).*(abs(x) - t).*s;

if (N0 == 0)
   gcv = 0;
else
   gcv = N * (norm(y-x) / N0 )^2 ;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
