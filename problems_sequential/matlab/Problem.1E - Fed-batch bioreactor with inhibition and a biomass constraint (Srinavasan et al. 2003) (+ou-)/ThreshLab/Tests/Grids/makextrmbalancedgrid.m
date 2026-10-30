
function x = makextrmhomogeneousgrid(N,gamma)

% makextrmhomogeneousgrid - ThreshLab/Tests/Grids/ --- Grid that reaches most
% inhomogeneity for given homogeneity constant
% fluctuations on fixed design
%   Usage
%     x = makextrmhomogeneousgrid(N,gamma);
%   Inputs
%     N            number of gridpoints
%     gamma        (optional; default = 1) homogeneity constant
%   Outputs
%     x            data points on [0,1]
%   See also
%     help homogeneity1

if nargin < 2, gamma = 1; end

beta = 1/(1+gamma);
x = (0:2)/2;
% x = (0:4)/4;
% x = (0:8)/8;
nn = length(x);
while nn < N,
   xx(1:2:2*nn-1) = x;
   xx(2:4:2*nn-1) = beta*x(1:2:nn-2) + (1-beta)*x(2:2:nn-1);
   xx(4:4:2*nn-1) = (1-beta)*x(2:2:nn-1) + beta*x(3:2:nn);
   nn = length(xx);
   x = xx;
end
x = x(1:N)/x(N);

