
function B = pointinside(x0,y0,points,x,y)

% pointinside - ThreshLab2/Triangulations --- checks if a point lies inside a
%               polygone
%   Usage
%      B = pointinside(x0,y0,points,x,y);
%   Inputs
%     x0,y0    coordinates of point
%     points   ordered vector of at least three indices into (x,y)
%     x        x-coordinates of 2d-scattered data
%     y        y-coordinates of 2d-scattered data
%   Outputs
%     B        Boolean: 1 if (x0,y0) inside area bordered by
%                       (x(points),y(points))
%   Description
%      The algorithm counts the number of intersections of the line x = x0 with
%      the boundary of the region. For a vector of points (x0,y0), use
%      pointsinside. pointsinside is slow when points is large.
%   Note
%      points on boundary are NOT inside according to this implementation
%   See also
%      help pointsinside
%      help pointonboundary

eps = 1.e-6;
n = length(points);
intersec = 0;
bound = 1;
xx = x(points);
yy = y(points);
for k = 1:n
   kk = mod(k,n)+1;
   x1 = xx(k);
   y1 = yy(k);
   x2 = xx(kk);
   y2 = yy(kk);
   if pointonline(x0,y0,x1,y1,x2,y2),
      bound = 0;
   end
   if (max(x1,x2) > x0)*(min(x1,x2) <= x0)
       r12 = (y2-y1)/(x2-x1);
       yint = y1+r12*(x0-x1);
       intersec = intersec + (yint < y0);
   end
end
disx = xx - x0;
disy = yy - y0;
dis = (abs(disx)+abs(disy))/2;
corner = all(dis>eps);
B = corner .* bound .* mod(intersec,2) ;
