
function y = shiftsegmentspolygone(x,r)

% shiftsegmentspolygone - ThreshLab2/Triangulations --- rescales closed polygone
%   Usage
%     y = shiftsegmentspolygone(x,r);
%   Inputs
%     x      matrix with two rows or two columns
%     r      rescaling factor (real number)
%   Outputs
%     y      matrix of same size as x
%   Description
%     This is a help function used in constructvoronoicell
%   Note
%     the use of r should be formalized; the current implementation is
%     heuristical
%   See also
%     help constructvoronoicell

transpose = 0;
sizex = size(x);
if sizex(1) == 2,
   n = sizex(2);
elseif sizex(2) == 2,
   n = sizex(1);
   x = x'; transpose = 1;
else
   error('incorrect size for input polygone')
end

cyc = 0;
if norm(x(1:2,1)-x(1:2,n)) > eps,
   xx = [x x(1:2,1)];
else
   xx = x;
   cyc = 1;
   n = n-1;
end

if n<3,
   error('insufficient number of vertices in polygone')
end

xx = [xx(1:2,n) xx];
X = x(1:2,1:n);

meanX = mean(X,2); 
d2 = 0; for i = 1:n, d2 = d2+norm(X(1:2,i)-meanX); end, d2 = d2/n;
r = (r-1)*d2;

angles = anglesboundary((1:n),xx(1,2:n+1),xx(2,2:n+1));
offsets = r./sin(angles/2);

CC = iscounterclockwise((1:n),xx(1,2:n+1),xx(2,2:n+1));

directions = xx(1:2,2:n+1)-xx(1:2,1:n);
normdirections = sqrt(directions(1,1:n).^2+directions(2,1:n).^2);
directions = directions./[normdirections;normdirections];
if CC,
   for i = 1:n,
      a = angles(i)/2;
      directions(1:2,i) = [cos(a) sin(a); -sin(a) cos(a)]*directions(1:2,i);
   end
else
   for i = 1:n,
      a = angles(i)/2;
      directions(1:2,i) = [cos(a) -sin(a); sin(a) cos(a)]*directions(1:2,i);
   end
end

% u = xx(1,1:n+2); v = xx(2,1:n+2);
% 
% tri = [(1:n); (2:n+1); (3:n+2)];
% x3 = u(tri); y3 = v(tri);
% C = triincenter(x3,y3);
% 
% sizeX = size(X)
% siceC = size(C)
% directions = X-C; 
% normdirections = sqrt(directions(1,1:n).^2+directions(2,1:n).^2);
% directions = directions./[normdirections;normdirections];

y = X + directions.*[offsets;offsets];
if cyc,
   y = [y y(1:2,1)];
end

if transpose,
   y = y';
end

