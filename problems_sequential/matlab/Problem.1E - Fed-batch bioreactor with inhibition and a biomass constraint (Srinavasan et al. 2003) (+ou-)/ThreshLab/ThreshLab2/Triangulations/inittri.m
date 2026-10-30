
function [tri, B] = inittri(x,y,anglethreshold,optionsfordelaunay)

% inittri - ThreshLab2/Triangulations --- initialise triangulation
%   Usage
%     [tri, B] = inittri(x,y,anglethreshold);
%     tri = inittri(x,y,anglethreshold);
%   Inputs
%     x     x-coordinates of 2d-scattered data
%     y     y-coordinates of 2d-scattered data
%     length(x) = length(y)
%     anglethreshold optional argument; default = pi/2 (see below)
%     optionsfordelaunay optional cell array used in function call
%                        tri = delaunay(x,y,optionsfordelaunay)
%   Outputs
%     tri   three-column matrix with indices into (x,y), the vertices of the
%           triangles
%     B     ordered list of indices into (x,y), the vertices on the
%           boundary of tri (see help boundary)
%   Description
%     delaunay triangulation + removal of all triangles on the boundary that
%     have inside angle larger than anglesthreshold. The default value of pi/2
%     removes all boundary triangles with obtuse inside angle.
%   See also
%     help initmultitri
%     help removenarrowtrinearboundary
%     help boundary

x = row(x);
y = row(y);

if nargin<3,
   anglethreshold = pi/2;
elseif nargin == 3,
   if iscell(anglethreshold),
      optionsfordelaunay = anglethreshold;
      anglethreshold = pi/2;
   end
end


% tri = delaunay(x,y,{'QJ'});
if exist('optionsfordelaunay') ~= 1,
   tri = delaunay(x,y);
else
   tri = delaunay(x,y,optionsfordelaunay);
end
[tri, B] = removenarrowtrinearboundary(tri,x,y,anglethreshold);
