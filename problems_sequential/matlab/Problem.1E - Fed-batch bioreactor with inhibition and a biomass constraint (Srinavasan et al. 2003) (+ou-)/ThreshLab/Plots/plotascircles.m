
function ans = plotascircles(x,y,r,varargin)

% plotascircles - ThreshLab/Plots -- plot z=f(x,y) as circles
%   Usage
%     plotascircles(x,y,z)
%     plotascircles(x,y,z,option)
%   Inputs
%     x,y    vectors of abscis values
%     z      vector of coordinate values z = f(x,y)
%     option plot options, as in matlab command plot
%   Outputs
%     none
%   Description
%     plot circles with radii proportional to values in z and centers in (x,y)
%   See also
%     help plotasspikes
%     help plotasblocks

if nargin < 3,
   error('Not enough input arguments.')
elseif nargin == 3,
   options = {'b'};
else
   options = varargin;
end

if length(r) == 1, r = r*ones(size(x)); end
h = ishold;
N = length(x);
m = 50;
t = (0:m)/m*2*pi;
for i = 1:N,
   u = x(i) + cos(t) * r(i);
   v = y(i) + sin(t) * r(i);
   plot(u,v,options{:});
   hold on
end
if ~h, hold off, end
ans = 1;

