
function [points, neighbourtri] = findneighbours(v0,tri);

% findneighbours - ThreshLab2/Triangulations --- find adjecent vertcies
%                  in triangulations, or sets of edges
%   Usage
%      [points, neighbourtri] = findneighbours(v0,tri);
%      points = findneighbours(v0,E);
%   Inputs
%      tri  triangulation: 3 times n matrix of indices
%      E    edges: 2 times n matrix of indices
%      v0   vertex: integer value, i.e., index
%   Outputs
%      points        adjacent vertices (vector of indices)
%      neighbourtri  triangles connecting v0 to its neighbours
%   Description
%   See also

if isempty(tri),
   tri = zeros(0,3);
end
if size(tri,2) == 2,
   % i.e., find neighbours by edges
   E = tri;
   TT = zeros(0,3);
   e = find(sum(E==v0,2));
   E = E(e,:);
   points = setdiff(row(unique(E)),v0);
else
   % i.e., find neighbours by triangles
   T1 = tri(find(tri(:,1)==v0),:);
   T2 = tri(find(tri(:,2)==v0),:);
   T3 = tri(find(tri(:,3)==v0),:);
   TT = [T1' T2' T3']';
   E = listedges(TT);
   e = find(1-sum(E==v0,2));
   E = E(e,:);
   points = row(unique(E));
   if length(points) < size(E,1)+2,
      % if neighbours can be linked in consistent order, do so
      points = linkedges(E);
   end
   neighbourtri = TT;
end
