function sure = SURE1hard(x,t,stdev)

% SURE1hard -- apply hard threshold + compute Stein's unbiased risk estimator
%                                         (SURE)
%  Usage
%    sure = SURE1hard(x,t,stdev)
%  Inputs
%    x      input coefficients
%    t      threshold (may be a vector)
%    stdev  noise standard deviation
%  Outputs
%    sure    the value of the sure
%  Description
%    Applies hard-threshold t to all elements of x and computes SURE
%  Note
%    A vector of thresholds may be interesting for plotting
%  See also
%    help SUREhard


xcol = column(x);
trow = row(t);
N = length(x);
m = length(t);
t = repmat(trow,N,1);
x = repmat(xcol,1,m);

s = abs(x) - t; s = (s > 0);
N1 = sum(s);

y = x.*s;

sure = sum((y-x).*(y-x))/N - stdev^2 + 2*stdev^2*N1/N ...
     + 4*stdev^2*t(1,:).*gauss(t(1,:),stdev);


% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
% www.cs.kuleuven.ac.be/~maarten/software/
