function [xx yy] = plotasblocks(x,y,varargin);

% plotasblocks - ThreshLab/Plots -- plot vectors as blocks (piecewise const.)
%   Usage
%     plotasblocks(x,y)
%     plotasblocks(y)
%     plotasblocks(x,y,option)
%     [xx yy] = plotasblocks(x,y,option)
%   Inputs
%     x      vector of abscis values
%     y      vector of coordinate values
%     option plot options, as in matlab command plot
%   Outputs
%     void or [xx yy]
%     xx     vector corresponding to x, such that [xx,yy] is a piecewise
%            constant polyline
%     yy     vector corresponding to y, such that [xx,yy] is a piecewise
%            constant polyline
%            plotasblocks(x,y,option) = plot(xx,yy,option)
%   Description
%     plots vector as blocks, i.e., as piecewise constants, rather than as
%     piecewise linear (polylines), which is the matlab standard way
%   See also
%     help plotasspikes
%     help plotascircles
%     help pcwsconstant

if nargin == 1,
   options = {'b'};
   y = x;
   x = (1:length(y));
elseif nargin == 2,
   if isa(y,'numeric') | isa(y,'logical'),
      options = {'b'};
   else
      options = {y};
      y = x;
      x = (1:length(y));
   end
else
   if ~isa(y,'numeric'),
      options{1} = y;
      for a = 3:nargin, options{a-1} = varargin{a-2}; end
      y = x;
      x = (1:length(y));
   else
      options = varargin;
   end
end

n = length(x);
if n~=length(y)
   error('cannot plot vectors of unequal lengths')
end
xx = zeros(2*n-1,1);
  xx(1:2:2*n-1) = column(x);
  xx(2:2:2*n-1) = column(x(2:n));
yy = zeros(2*n-1,1);
  yy(1:2:2*n-1) = column(y);
  yy(2:2:2*n-1) = column(y(1:n-1));
plot(xx,yy,options{:})
if nargout == 0,
   clear xx yy;
end
