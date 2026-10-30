
function points = triboundary(tri);

% triboundary - ThreshLab2/Triangulations --- finds boundary of triangulation
%   Usage
%      points = triboundary(tri);
%   Inputs
%      tri three-column vector of indices (=vertices)
%   Outputs
%      points  ordered list of boundary vertices
%   See also
%      help listedges
%      help linkedges

[E etri] = listedges(tri);
e = find(~etri(:,2));
Ebound = E(e,:);
points = linkedges(Ebound);
