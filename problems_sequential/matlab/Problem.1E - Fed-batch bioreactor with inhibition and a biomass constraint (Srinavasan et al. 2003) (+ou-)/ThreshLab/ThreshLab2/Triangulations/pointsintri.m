
function B = pointsintri(x0,y0,x3,y3);

% pointsintri - ThreshLab2/Triangulations --- checks if points are inside a
%                                           triangle
%   Usage
%     B = pointsintri(x0,y0,x3,y3);
%   Inputs
%     x0,y0  two matrices of equal size containing coordinates of points
%     x3,y3  two vectors of length 3, containing coordinates of triangle
%            vertices
%   Outputs
%     B      boolean matrix with size equal to that of x0/y0
%   Description
%     B(i,j) = 1 <=> (x0(i,j),y0(i,j)) inside triangle defined by x3 and y3
%   See also
%     help pointsinside
%     help pointinside

if size(x0)~=size(y0), error('inputs x0 and y0 must have equal size'), end

B = ones(size(x0));
c = [1 2 3 1 2];
for l = 1:3,
   c1 = c(l); c2 = c(l+1); c3 = c(l+2);
   line = (y3(c3)-y3(c1))*(x3(c2)-x3(c1)) - (y3(c2)-y3(c1))*(x3(c3)-x3(c1));
   line0 = (y0-y3(c1))*(x3(c2)-x3(c1)) - (y3(c2)-y3(c1))*(x0-x3(c1));
   B = B & (sign(line0) ~= -sign(line));
end
