function sure = SURE1soft(x,t,stdev)

% SURE1soft -- apply soft threshold + compute Stein's unbiased risk estimator
%                                         (SURE)
%  Usage
%    sure = SURE1soft(x,t,stdev)
%  Inputs
%    x      input coefficients
%    t      threshold (may be a vector)
%    stdev  noise standard deviation
%  Outputs
%    sure    the value of the sure
%  Description
%    Applies soft-threshold t to all elements of x and computes SURE
%  Note
%    A vector of thresholds may be interesting for plotting
%  See also
%    help SUREsoft

xcol = column(x);
trow = row(t);
N = length(x);
m = length(t);
t = repmat(trow,N,1);
x = repmat(xcol,1,m);

s = abs(x) - t; s = (s > 0);
N1 = sum(s);
y = sign(x).*(abs(x) - t).*s;

sure = sum((y-x).*(y-x))/N - stdev^2 + 2*stdev^2*N1/N;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
