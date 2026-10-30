function gcv = GCV2soft(x,t)

% GCV2soft -- apply soft threshold + compute GCV in 2D
%  Usage
%    gcv = GCV2soft(x,t)
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%  Outputs
%    gcv    the value of the gcv
%  Description
%    Applies soft-threshold t to all elements of x and computes GCV
%    GCV = Generalized Cross Validation
%  See also
%    help GCVsoft
%    help mingcv2soft

[m,n] = size(x);
N = m*n;
y = x;
s = abs(x) - t; s = (s > 0);
N0 = N - sum(sum(s));
y = sign(x).*(abs(x) - t).*s;

gcv = n*m * sum(sum((y-x).^2)) / N0^2 ;


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
