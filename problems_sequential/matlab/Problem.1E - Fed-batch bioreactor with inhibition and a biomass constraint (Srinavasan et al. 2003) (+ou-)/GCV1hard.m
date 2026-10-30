function gcv = GCV1hard(x,t,stdev)

% GCV1hard -- apply hard threshold + compute 
%                      Generalized Cross Validation (GCV)
%  Usage
%    gcv = GCV1hard(x,t,stdev)
%  Inputs
%    x      input coefficients
%    t      threshold (may be a vector)
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
%    A vector of thresholds may be interesting for plotting
%  See also
%    help GCV1hardDoFmirror

xcol = column(x);
trow = row(t);
N = length(x);
m = length(t);
t = repmat(trow,N,1);
x = repmat(xcol,1,m);

s = abs(x) - t; s = (s > 0);
N1 = sum(s); N0 = N-N1;

N00 = N0;
N0 = N0 + (N0 == 0);

y = x.*s;

SSe = sum((y-x).*(y-x));
% if nargin < 3, stdev = sqrt(SSe./N0); end
% if nargin < 3, stdev = MAD1(x)/invcumgauss(0.75); end
if nargin < 3, stdev = sqrt(varestGCV(x)); end
nuhatprime  = 2*N*trow.*gauss(trow,stdev);
gcv = (SSe+trow.^2.*N1)/N ./ (1 - N1/N).^2  + 2*stdev^2*nuhatprime/N ...
                                            - trow.^2.*N1/N ;

% Copyright (c) Maarten Jansen
% 
% This software is part of PiefLab and is copyrighted material. 
