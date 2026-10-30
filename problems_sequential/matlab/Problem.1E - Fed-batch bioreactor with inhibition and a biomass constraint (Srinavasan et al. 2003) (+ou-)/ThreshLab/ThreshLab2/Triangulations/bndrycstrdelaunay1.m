
function tri = bndrycstrdelaunay1(x,y,boun);

% bndrycstrdelaunay1 - ThreshLab2/Triangulations ---
%                       boundary constrained delaunay triangulation
%   Usage
%      tri = bndrycstrdelaunay1(x,y,boun);
%   Inputs
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%      boun       ordered list of indices into x and y corresponding to
%                 vertices on the boundary (default = convex hull of (x,y))
%   Outputs
%      tri        triangulation
%   Description
%      triangulation within boundary described by boun. All vertices in x and y
%      take part in the triangulation, except those outside boun.
%
%      The routine proceeds in two basic steps:
%      First, the boundary constrained minimum spanning tree is established
%             The spanning tree branches act as a partition of the region
%             bounded by boun.
%      Second, the routine triangulates all subregions in the partition. The
%              partition of a subregion is easier since we know that all 
%              participating vertices lie on the boundary of the subregion. The
%              triangulation is performed by the routine
%              bndrycstrbndrydelaunay1
%      The routine is thus based on the property that a Delaunay triangulation
%      includes the minimum spanning tree.
%   See also
%      help bndrycstrbndrydelaunay1
%      help bndrycstrminspantree2d
%   Reported bug
%      the routine crashes when some elements in x,y are outside boun, because
%      this routine calls bndrycstrminspantree2d, which assumes that all
%      elements in (x,y) have to be linked to each other

if nargin<3,
   boun = convexhull(x,y);
end
if length(x) < 1,
   tri = zeros(0,3);
else
   E = bndrycstrminspantree2d(x,y,boun);
   n = length(boun);
   if ~iscounterclockwise(boun,x,y),
      boun = boun(n:-1:1);
   end
   if boun(1) == boun(n),
      n = n-1;
      boun = boun(1:n);
   end
   
   tri = zeros(0,3);
   p1 = boun(1);
   for b = 1:n,
      boun = shift(boun);
      p2 = boun(1);
      pend = boun(2);
      % we now try to connect p2 with pend using the MST, we then triangulate 
      % the domain enclosed by this MST section and the segment p2-pend. The 
      % point p1 only serves to select the segment of the MST starting in p2 
      % with the smallest angle wrt the boundary.
      pp = p2;
      while p2 ~= pend,
         points = findneighbours(p2,E);
         point = minanglepoint(p1,p2,points,x,y);
         pp = [pp point];
         p1 = p2;
         p2 = point;
      end
      trib = bndrycstrbndrydelaunay(pp,x,y);
      tri = [tri ; trib];
      p1 = boun(1);
   end
end
