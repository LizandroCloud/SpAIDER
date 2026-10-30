
function conv = convexhull(points,x,y);

% convexhull - ThreshLab2/Triangulations --- finds convex hull
%   Usage
%      conv = convexhull(points,x,y);
%   Inputs
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%      points     list of indices into x and y corresponding to vertices 
%   Outputs
%      conv       ordered list of vertices (= indices into x/y) that constitute
%                 the boundary of the convex hull of x(points),y(points)
%   Description
%      x([conv conv(1)]),y([conv conv(1)]) is the convex hull of the set
%      x(points),y(points)
%      Note that the closed polyline is represented by a list where the closure
%      is not made explicit.
%      convexhull relies on the delaunay triangulation. See convexhuull2 for a
%      different approach
%   See also
%      help convexhull2

if nargin == 2,
   y = x;
   x = points;
   points = (1:length(x));
end
if length(x) ~= length(y),
   warning('unequal lengths x and y coordinates')
   n = min(length(x),length(y));
   points = points(find(points<n+1));
end
if length(points) < 3,
   conv = points;
else
   xx = x(points);
   yy = y(points);
   tri = delaunay(xx,yy,{'QJ'});
   conv = triboundary(tri);
   conv = points(conv);
end
