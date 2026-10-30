
function x = makecantorgrid(N,gap,xrandom,realrandom)

% makecantorgrid - ThreshLab/Tests/Grids/ --- cantor's function (devil's
% staircase)
%   Usage
%     x = makecantorgrid(N,gap,xrandom,realrandom);
%   Inputs
%     N            number of gridpoints
%     gap          (optional; default = 1/3) paramater of function
%     xrandom      (optional; default = 0) binary value for random design
%     realrandom   (optional; default = 0) binary value; 0 for fixed seed
%   Outputs
%     x            data points on [0,1]
%
if nargin < 4, realrandom = 0; end
if nargin < 3, xrandom = 0; end
if nargin < 2, gap = 1/3; end
if nargin < 4 & (gap == 0 | gap == 1),
   realrandom = xrandom; xrandom = gap; gap = 1/3;
end

if ~realrandom, rand('seed',931316785); end

D = [0 (1-gap)/2 (1+gap)/2 1];
nn = 2;
x = (0:1);
while nn < N,
   for k = 1:2:nn-1,
      if xrandom, D = sort(rand(1,4)); end
      xx(2*k-1:2*k+2) = D*(x(k+1)-x(k))+x(k);
   end
   x = xx;
   nn = 2*nn;
end
x = x(1:N)/x(N);
