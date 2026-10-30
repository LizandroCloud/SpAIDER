
function B = pointondisk(x,y,x0,y0,r);

% pointondisk - ThreshLab2/Triangulations --- checks if a point lies inside a
%                                           circle
% Usage
%    B = pointondisk(x,y,xC,yC,r);
% Inputs
%    x,y  coordinates of point (may be vectors of equal size)
%   xC,yC coordinates of circle center
%    r    circle radius
% Outputs
%    B    Boolean: 1 if (x,y) inside circel C(xC,yC;r)

eps = 1.e-8;
if r == Inf,
   B = true(size(x));
else
   d2 = (x-x0).^2 + (y-y0).^2;
   r2 = r^2;

   B = (d2 < r2*(1-eps));
end
