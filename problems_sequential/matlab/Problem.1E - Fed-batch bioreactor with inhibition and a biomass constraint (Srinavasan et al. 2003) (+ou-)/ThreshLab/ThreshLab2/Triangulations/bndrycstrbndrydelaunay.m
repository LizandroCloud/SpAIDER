
function [tri, E] = bndrycstrbndrydelaunay(points,x,y);

% bndrycstrbndrydelaunay - ThreshLab2/Triangulations ---
%                       boundary constrained boundary delaunay triangulation
%   Usage
%      [tri E] = bndrycstrbndrydelaunay(points,x,y);
%   Inputs
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%      points     ordered list of indices into x and y corresponding to
%                 vertices
%   Outputs
%      tri        triangulation
%      E          edges of triangulation
%   Description
%      triangulation within boundary described by points. All vertices of
%      points are assumed to be on the boundary (if this is not the case, use
%      bndrycstrdelaunay).
%
%      This routine calls bndrycstrbndrydelaunay1 for technical details
%
%   See also
%      help bndrycstrbndrydelaunay1
%      help bndrycstrdelaunay

npoints = length(points);
if npoints < 3,
   tri = zeros(0,3);
   E = zeros(0,2);
else
   tri = bndrycstrbndrydelaunay1(points,x,y);
   E = listedges(tri);
end
