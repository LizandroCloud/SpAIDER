
function conv = convexhull2(points,x,y);

% convexhull2 - ThreshLab2/Triangulations --- finds convex hull
%   Usage
%      conv = convexhull2(points,x,y);
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
%      This routine provides an alternative to convexhull. Unlike convexhull,
%      convexhull2 does not rely on delaunay triangulations
%      convexhull2 implements the gift wrapping algorithm (Jarvis march)
%   See also
%      help convexhull

if nargin == 2,
   y = x;
   x = points;
   N = length(x);
   points = (1:N);
else
   N = length(points);
end
x = row(x(points));
y = row(y(points));

[xmin p0] = min(x);
ymin = min(y); ymax = max(y);
x = [x xmin];
y = [y 2*ymax-ymin];
p1 = N+1;
conv = [p1 p0];

while p0 ~= conv(1),
   x0 = x(p0); y0 = y(p0);
   x1 = x(p1)-x0; y1 = y(p1)-y0;
   x2 = x - x0; y2 = y - y0;
   p1p2 = ((x1*x2)+(y1*y2))./(sqrt(x1^2+y1^2)*sqrt(x2.^2+y2.^2));
   [cosmaxangle p2] = min(p1p2);
   conv = [conv p2];
   if conv(1)~= conv(2), %i.e., in the first step
      conv(1) = conv(2); x = x(1:N); y = y(1:N);
   end
   p1 = p0;
   p0 = p2;
end
n = length(conv);
conv = conv(2:n-1);

conv = points(conv);
