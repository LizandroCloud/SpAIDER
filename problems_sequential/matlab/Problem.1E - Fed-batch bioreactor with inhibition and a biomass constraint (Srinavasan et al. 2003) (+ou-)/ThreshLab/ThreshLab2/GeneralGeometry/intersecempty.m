
function B = intersecempty(x1,y1,x2,y2,x3,y3,x4,y4,eps);

% intersecempty - ThreshLab2/GeneralGeometry --- checks if two line segments have
% intersection
%   Usage
%      B = intersecempty(x1,y1,x2,y2,x3,y3,x4,y4,eps);
%   Inputs
%      x1,y1,x2,y2,x3,y3,x4,y4 coordinates of (x1,y1), (x2,y2), etc.
%      eps    precision. If eps > 0, then B = 1 if one of the end points lies
%                        exactly on the other line segment
%   Outputs
%      B   boolean: 1 if the line segment between (x1,y1) and (x2,y2) has no
%      intersection with the line segment between (x3,y3) and (x4,y4)
%   Description
%      The routine checks if (x1,y1) and (x2,y2) are on the same side of the
%      line (x3,y3)--(x4,y4) and if (x3,y3) and (x4,y4) are on the same side of
%      the line (x1,y1)--(x2,y2)

if nargin < 9,
   eps = 1e-10;
end

z1 = (y1-y3).*(x4-x3) - (x1-x3).*(y4-y3);
z2 = (y2-y3).*(x4-x3) - (x2-x3).*(y4-y3);
z3 = (y3-y1).*(x2-x1) - (x3-x1).*(y2-y1);
z4 = (y4-y1).*(x2-x1) - (x4-x1).*(y2-y1);

z1 = sign(z1).*(abs(z1)>eps);
z2 = sign(z2).*(abs(z2)>eps);
z3 = sign(z3).*(abs(z3)>eps);
z4 = sign(z4).*(abs(z4)>eps);

B = (z1 == -z2) & (z3 == -z4);
if eps == 0,
   B = (z1 ~= z2) & (z3 ~= z4);
end
B = ~B;
