
function inside = pointsinside(x0,y0,bound,x,y);

% pointsinside - ThreshLab2/Triangulations --- checks if points are inside a
%                                           triangle
%   Usage
%     B = pointsinside(x0,y0,points,x,y);
%   Inputs
%     x0,y0    two matrices of equal size containing coordinates of points
%     points   ordered vector of at least three indices into (x,y)
%     x        x-coordinates of 2d-scattered data
%     y        y-coordinates of 2d-scattered data
%   Outputs
%     B      boolean matrix with size equal to that of x0/y0
%   Description
%     B(i,j) = 1 <=> (x0(i,j),y0(i,j)) inside polygone defined by
%                                     (x(points),(points))
%     pointsinside triangulates (x(points),(points)) and then applies
%     pointsintri to all triangles. This approach is slow if points is a large
%     vector, i.e., if there are a lot of triangles to be checked. In that case
%     a for-loop of pointinside running for all points in (x0,y0) is preferable
%   Note
%     points on boundary are inside according to this implementation (unlike
%     pointinside). Use ~pointonboundary to mimick effect of pointinside
%   See also
%     help pointsintri
%     help pointinside
%     help pointonboundary

inside = zeros(size(x0));
tri = bndrycstrbndrydelaunay(bound,x,y);
for t = 1:size(tri,1),
   inside = inside | pointsintri(x0,y0,x(tri(t,1:3)),y(tri(t,1:3)));
end 
