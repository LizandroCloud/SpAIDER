
function gcv = GCV1softhard(x,t)

% GCV1softhard -- apply hard threshold + evaluate expression of GCV for soft
%                 thresholding
%  Usage
%    gcv = GCV1softhard(x,t);
%  Inputs
%    x      input coefficients
%    t      threshold (may be a vector)
%  Outputs
%    gcv    value of gcv for soft thresholding
%  Description
%    This function applies hard thresholding but evaluates the quality of the
%    threshold by the expression of GCV for soft thresholding
%    This routine is for illustrative purpose only: the combination HT - GCV-ST
%    does not work properly.
%  Note
%    A vector of thresholds may be interesting for plotting
%  See also
%    help SURE1soft

xcol = column(x);
trow = row(t);
N = length(x);
m = length(t);
t = repmat(trow,N,1);
x = repmat(xcol,1,m);

s = abs(x) - t; s = (s > 0);
N0 = N - sum(s);
y = x.*s;

N00 = N0;
N0 = N0 + (N0 == 0);
gcv = N * sum((y-x).*(y-x)) ./ N0.^2 ;
gcv = gcv.*(N00>0);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
% www.cs.kuleuven.ac.be/~maarten/software/
