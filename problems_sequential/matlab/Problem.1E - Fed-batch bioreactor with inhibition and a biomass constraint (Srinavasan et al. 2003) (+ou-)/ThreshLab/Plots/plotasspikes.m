function [xx yy] = plotasspikes(x,y,varargin);

% plotasspikes - ThreshLab/Plots -- plot vectors as spikes from zero
%   Usage
%     plotasspikes(x,y)
%     plotasspikes(y)
%     plotasspikes(x,y,option)
%     [xx yy] = plotasspikes(x,y,option)
%   Inputs
%     x      vector of abscis values
%     y      vector of coordinate values
%     option plot options, as in matlab command plot
%   Outputs
%     void or [xx yy]
%     xx     vector corresponding to x, such that 
%                    plotasspikes(x,y,option) = plot(xx,yy,option)
%     yy     vector corresponding to y, such that
%                    plotasspikes(x,y,option) = plot(xx,yy,option)
%   Description
%     plots vector as spikes, i.e., with zeros in between values. This is
%     interesting for plotting sparse vectors
%   See also
%     help plotasblocks
%     help plotascircles

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
xx = zeros(3*n-2,1);
  xx(3:3:3*n-2) = column(x(2:n));
  xx(1:3:3*n-2) = column(x);
  xx(2:3:3*n-2) = column(x(1:n-1));
yy = zeros(3*n-2,1);
  yy(1:3:3*n-2) = column(y);
plot(xx,yy,options{:})
if nargout == 0,
   clear xx yy;
end
