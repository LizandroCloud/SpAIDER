
function H = isingpot(X,tau)

% isingpot -- log probability of binary Ising model
%   Usage
%    H = isingpot(X,tau);
%   Inputs
%    X    binary (+/- 1) matrix
%    tau  rigidity/temperature constant
%   Outputs
%    H    potential/ log probability of X
%   Description
%   Note
%    This is a vector implementation, which means that update of one component
%    is not yet active in updating other components in this iteration step
%    (slower convergence rate per iteration, but faster iterations)
%   See also
%    help isingpotupdate
%    help condprob2

% H = isingpot(X,tau);
% X is binary image
% entries of X are +/- 1 (no zeros!)
% tau is rigidity factor (to be estimated by MPLEising.m)

[N M] = size(X);
positions = [-1 -1
             -1  0
             -1  1
              0 -1
              0  0
              0  1
              1 -1
              1  0
              1  1];
Uc = zeros(size(X(2:N-1,2:M-1)));
Ns = 0;
for s = 1:9
   sx = positions(s,1);
   sy = positions(s,2);
   if (s == 5)
       Uc = Uc + 0;
   else
       Uc = Uc - X(2+sx:N-1+sx,2+sy:M-1+sy).*X(2:N-1,2:M-1);
   end
end
H = tau*sum(sum(Uc));
% H = tau*Uc;

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
