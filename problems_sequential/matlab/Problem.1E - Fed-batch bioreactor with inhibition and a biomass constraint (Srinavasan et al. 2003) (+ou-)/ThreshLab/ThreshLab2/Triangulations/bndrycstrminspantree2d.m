
function E = bndrycstrminspantree2d(x,y,boundary);

% bndrycstrminspantree2d - ThreshLab2/Triangulations --- 
%                          boundary constrained minimum spanning tree
%   Usage
%      E = bndrycstrminspantree2d(x,y,boundary);
%   Inputs
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%      boundary   ordered list of indices into x and y corresponding to boundary
%   Outputs
%      E          edges of spanning tree
%   Description
%      minimum spanning tree, constrained by the fact that no branch/path must
%      connect vertices through zones outside the boundary. This is realized by
%      imposing an infinite weight (distance) on such paths
%   See also
%      help listedges
%      help linkedges

N = length(x);
E = nchoosek((1:N),2);
Elengths = sqrt((x(E(:,1)) - x(E(:,2))).^2 ...
              + (y(E(:,1)) - y(E(:,2))).^2);
[Elengths I] = sort(Elengths);
allE = E(I,:);
clear E;
ae = 0;
treecomponent = (1:N); % this variable indicates for every vertex to which
                       % tree component it belongs. A component is a set of
                       % connected vertices. At the end, all components are
                       % linked together into a single one.

for e = 1:N-1,
   ae = ae+1;
   ee = allE(ae,1:2);
   while treecomponent(ee(1)) == treecomponent(ee(2)) | ...
         ~(edginside(ee,boundary,x,y) | edgonboundary(ee,boundary)),
      ae = ae+1;
      ee = allE(ae,1:2);
   end
   E(e,1:2) = ee;
   tc2 = find(treecomponent == treecomponent(ee(2)));
   treecomponent(tc2) = treecomponent(ee(1));
end
