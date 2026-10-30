
function x = makeuniformgrid(N,no_origin,xrandom,realrandom)

% makecantorgrid - ThreshLab/Tests/Grids/ --- Uniform grid (random or fixed
% design)
%   Usage
%     x = makeuniformgrid(N,no_origin,xrandom,realrandom);
%   Inputs
%     N            number of gridpoints
%     no_origin    (optional; default = 1) first point is 1/N (default) or 0
%     xrandom      (optional; default = 0) binary value for random design
%     realrandom   (optional; default = 0) binary value; 0 for fixed seed
%   Outputs
%     x            data points on [0,1]
%
if nargin < 4, realrandom = 0; end
if nargin < 3, xrandom = 0; end
if nargin < 2, no_origin = 1; end

if ~realrandom, rand('seed',931316785); end

if ~xrandom,
   if no_origin,
      x = (1:N)/N;
   else
      x = (0:N-1)/(N-1);
   end
else
   if no_origin,
      x = sort(rand(1,N));
   else
      x = [0 sort(rand(1,N-2)) 1];
   end
end
