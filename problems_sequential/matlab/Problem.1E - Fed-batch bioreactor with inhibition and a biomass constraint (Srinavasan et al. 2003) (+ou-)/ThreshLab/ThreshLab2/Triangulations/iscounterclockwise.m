
function B = iscounterclockwise(boun,x,y);

% iscounterclockwise - ThreshLab2/Triangulations --- checks if a boundary is
%                                    ordered counterclockwise or clockwise
%   Usage
%      B = iscounterclockwise(boun,x,y);
%   Inputs
%      boun       ordered list of indices into x and y
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%   Outputs
%      B          boolean, 1 if boun has counterclockwise order
%                          0 if boun has clockwise order
%   Description
%      the polyline x(boun),y(boun) must not have self-intersections
%   See also
%      help anglesboundary

n = length(boun);
if boun(n) == boun(1),
   n = n-1;
   boun = boun(1:n);
end
[beta B] = anglesboundary(boun,x,y);
