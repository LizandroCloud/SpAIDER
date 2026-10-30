function gcv = GCV1hardDoFmirror(x,t,stdev)

% GCV1hardDoFmirror -- apply hard threshold + compute 
%                      Generalized Cross Validation (GCV)
%  Usage
%    gcv = GCV1hardDoFmirror(x,t,stdev)
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
%    GCV = (1/n)*SS_E/(1-dof/n)^2
%    
%    where SS_E is the sum of squares of the residuals
%    and   dof = N1 + n*E[(1-epsilon^2/stdev^2)*I(epsilon^2>t^2)]
%    with  epsilon distributed like the observational errors
%    and   I(A) is the characteristic/identicator function for event A
%    
%    in abscence of argument stdev, the standard deviation is estimated using
%    stdevest = sqrt(varestGCV(x))
%
%    This definition of GCV is the `natural' extension of the soft-threshold
%    version, but it may be unstable; therefore, an alternative definition is
%    used in GCV1hard.
%  Note
%    A vector of thresholds may be interesting for plotting
%  See also
%    help GCV1hard

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
nuhat = N1+N*2*trow.*gauss(trow,stdev);
gcv = (SSe/N)./(1 - nuhat/N).^2;


% Copyright (c) Maarten Jansen
% 
% This software is part of PiefLab and is copyrighted material. 
