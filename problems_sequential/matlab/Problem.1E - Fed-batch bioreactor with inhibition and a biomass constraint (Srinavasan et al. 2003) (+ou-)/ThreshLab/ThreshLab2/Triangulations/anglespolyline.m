
function angle = anglespolyline(points,x,y);

% anglespolyline - ThreshLab2/Triangulations --- find angles on polyline
%   Usage
%      angle = anglespolyline(polyline,x,y);
%   Inputs
%      polyline   ordered list of indices into x and y corresponding to
%                 vertices
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%   Outputs
%      angle      vector of size(polyline)-2 with angles
%   Description
%      This routine finds angles on polyline x(points),y(points).
%      Convention: the first angle is less than pi, next angles are measured on
%      the same side of the polyline as the first one.
%   See also
%      help anglesboundary (for inside angles on closed polylines)

n = length(points);
if n < 3,
   angle = zeros(0,0);
else
   points1 = points(1:n-2);
   points2 = points(2:n-1);
   points3 = points(3:n);
   angle = (x(points1) - x(points2)).*(x(points3) - x(points2)) ...
         + (y(points1) - y(points2)).*(y(points3) - y(points2));
   angle = angle ./ ...
           sqrt((x(points1) - x(points2)).^2 + (y(points1) - y(points2)).^2);
   angle = angle ./ ...
           sqrt((x(points3) - x(points2)).^2 + (y(points3) - y(points2)).^2);
   angle = acos(angle);
   for i = 3:n-1,
      ee = points([i-2 i]);
      ff = points([i-1 i+1]);
      if ~edgecrossing(ee,ff,x,y),
         angle(i-1:n-2) = 2*pi-angle(i-1:n-2);
      end
   end
end
