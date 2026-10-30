
function B = colinear(x1,y1,x2,y2,x3,y3,eps);

% colinear - ThreshLab2/GeneralGeometry --- checks if three points are colinear
%   Usage
%      B = colinear(x1,y1,x2,y2,x3,y3);
%   Inputs
%      x1,y1,x2,y2,x3,y3 coordinates of (x1,y1), (x2,y2), (x3,y3)
%                        all six vectors of the same size
%   Outputs
%      B   boolean: 1 if the points are colinear
%   Description

if nargin<7,
   eps = 1.e-10;
end

% A = [y1 x1 1
%      y2 x2 1
%      y3 x3 1];
% B = abs(det(A)) < eps;

determinant = x1.*y2+x2.*y3+x3.*y1-y1.*x2-y2.*x3-y3.*x1;
B = abs(determinant) < eps;
