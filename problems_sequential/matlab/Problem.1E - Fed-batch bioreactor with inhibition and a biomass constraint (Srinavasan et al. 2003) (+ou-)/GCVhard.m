
function gcv = GCVhard(x,t,stdev)

% GCVhard -- apply hard threshold + compute Generalized Cross Validation (GCV)
%  Usage
%    gcv = GCVhard(x,t)
%  Inputs
%    x      input coefficients
%    t      threshold (must be scalar)
%    stdev  (optional) standard deviation of observational errors
%  Outputs
%    gcv    the value of the gcv
%  Description
%    Applies hard-threshold t to all elements of x and computes GCV
%    where GCV is defined as
%
%    GCV = (1/n)*[SS_E + thr^2*N1]/(1-N1/n)^2 + 2*dofprime*stdev^2 - thr^2*N1/n
%
%    where SS_E is the sum of squares of the residuals
%    and   dofprime = n*E[(1-epsilon^2/stdev^2)*I(epsilon^2>t^2)]
%    with  epsilon distributed like the observational errors
%    and   I(A) is the characteristic/identicator function for event A
%
%    in abscence of argument stdev, the standard deviation is estimated using
%    stdevest = sqrt(varestGCV(x))
%
%    This definition of GCV is stabilized version of GCV1hardDoFmirror
%  Note
%    For making a plot of gcv for a given vector of thresholds
%    use GCV1hard (for one dimension)
%
%    For minimizing GCVhard, use mingcvhard (1D) / mingcv2hard (2D)
%  See also:
%    help mingcvhard
%    help mingcv2hard
%    help GCV1hard
%    help GCV1hardDoFmirror
%    help GCVsoft


N = length(x);
s = abs(x) - t; s = (s > 0);
yh = x.*s;
N1 = sum(s);
SSe = norm(yh-x)^2;
% if nargin < 3, stdev = sqrt(SSe/(N-N1)); end
% if nargin < 3, stdev = MAD1(x)/invcumgauss(0.75); end
if nargin < 3, stdev = sqrt(varestGCV(x)); end
nuhatprime  = N*2*t*gauss(t,stdev);

if (N-N1 == 0)
   gcv = 0;
else
   gcv = SSe/N / (1-N1/N)^2 +2*stdev^2*nuhatprime/N - t^2*N1/N ;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
