
function x = makenormalgrid(N,no_origin,stdev,realrandom)

% makecantorgrid - ThreshLab/Tests/Grids/ --- Grid with independent normal
% fluctuations on fixed design
%   Usage
%     x = makenormalgrid(N,no_origin,stdev,realrandom);
%   Inputs
%     N            number of gridpoints
%     no_origin    (optional; default = 1) first point is 1/N (default) or 0
%     stdev        (optional; default = 0) standard deviation of fluctuations
%     realrandom   (optional; default = 0) binary value; 0 for fixed seed
%   Outputs
%     x            data points on [0,1]
%   Description
%     typical use for stdev is stdev = 1/(6*N);
if nargin < 4, realrandom = 0; end
if nargin < 3, stdev = 0; end
if nargin < 2, no_origin = 1; end

if ~realrandom, randn('seed',931316785); end

if no_origin,
   x0 = (1:N)/N;
else
   x0 = (1:N-2)/(N-1);
end
epsilon = randn(size(x0))*stdev;
if no_origin,
   x = x0+epsilon;
else
   x = [0 x0+epsilon 1];
end
x = sort(x);
