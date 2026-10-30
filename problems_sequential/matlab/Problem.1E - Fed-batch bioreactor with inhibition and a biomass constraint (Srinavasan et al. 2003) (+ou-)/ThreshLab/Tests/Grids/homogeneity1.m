
function [gamma, i] = homogeneity1(x);

% homogeneity1 - ThreshLab/Tests/Grids/ --- homogeneity constant for 1D grids
%   Usage
%     [gamma, i] = homogeneity1(x);
%     gamma      = homogeneity1(x);
%   Inputs
%     x            data points
%   Outpouts
%     gamma        homogeneity constant (see definition below)
%     i            index into x of most "inhomogeneous" location
%   Description
%     This function computes the homogeneity of a grid x, defined as follows.
%     Let Delta = x(2:N) - x(1:N-1), then
%     gamma = max(Delta_i/Delta_{i-1},Delta_{i-1}/Delta_i), 
%     minimum value for gamma is 1, this is an equidistant grid
%     maximum value for gamma is Inf
%     homogeneity1 also returns the indices i of the points in x where the
%     inhomogeneity reaches its maximum.
%     In most cases, i will be a single number. In some degenerated cases,
%     however, i can be a vector. 
%     E.g.: x = [1 2 3 4 5 6 8 9 10 11]; will result in i being [6 7]. The gap 
%     between those 2 points is the only irregularity. Note that a 
%     `nearly-regular' grid can also have one point of maximum inhomogeneity, 
%     e.g.: x = [1 2 3 4 5 6 8 10 12 14]; Now i = 6. 
%   See also
%     help makextrmhomogeneousgrid
%     help dyadicbalance1
%   Reference
%     Ingrid Daubechies, Igor Guskov, and Wim Sweldens, Regularity of Irregular
%     Subdivision, Constructive Approximation, 15, pages 381-426, 1999.


x = sort(x);
N = length(x);
if length(x) < 3,
   gamma = NaN;
else
   delta = x(2:N) - x(1:N-1);
   if min(delta) < eps,
      gamma = Inf;
      if nargout > 1,
         i = find(abs(delta) < eps);
         i = unique([i i+1]);
      end
   else
      gamma1 = max(delta(2:N-1)./delta(1:N-2));
      gamma2 = max(delta(1:N-2)./delta(2:N-1));
      gamma = max([gamma1 gamma2]);
      if nargout > 1,
         i = find(delta(2:N-1)./delta(1:N-2) > gamma-2*eps)+1;
         i = [i find(delta(1:N-2)./delta(2:N-1) > gamma-2*eps)+1];
         i = unique(i);
      end
   end
end
