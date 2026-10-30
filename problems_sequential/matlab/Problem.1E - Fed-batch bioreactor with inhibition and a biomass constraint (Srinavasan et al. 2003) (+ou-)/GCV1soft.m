
function gcv = GCV1soft(x,t)

% GCV1soft -- apply soft threshold + compute Generalized Cross Validation (GCV)
%  Usage
%    gcv = GCV1soft(x,t)
%  Inputs
%    x      input coefficients
%    t      threshold (may be a vector)
%  Outputs
%    gcv    the value of the gcv
%  Description
%    Applies soft-threshold t to all elements of x and computes GCV
%  Note
%    A vector of thresholds may be interesting for plotting
%  See also
%    help GCVsoft

xcol = column(x);
trow = row(t);
N = length(x);
m = length(t);
t = repmat(trow,N,1);
x = repmat(xcol,1,m);

s = abs(x) - t; s = (s > 0);
N0 = N - sum(s);
y = sign(x).*(abs(x) - t).*s;

N00 = N0;
N0 = N0 + (N0 == 0);
gcv = N * sum((y-x).*(y-x)) ./ N0.^2 ;
gcv = gcv.*(N00>0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
