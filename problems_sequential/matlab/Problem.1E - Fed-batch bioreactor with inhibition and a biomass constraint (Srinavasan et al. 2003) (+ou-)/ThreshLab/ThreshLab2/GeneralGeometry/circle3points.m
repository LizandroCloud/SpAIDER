
function [C,r] = circle3points(P1,P2,P3);

% circle3points - ThreshLab2/GeneralGeometry --- finds circle from 3 points
%   Usage
%      [C,r] = circle3points(P1,P2,P3);
%   Inputs
%      P1, P2, P3: 2-dimensional vectors
%   Outputs
%      C   2-dimensional vector, center point
%      r   1-dimensional positive number: radius
%   See also
%      help triexcenter
%   Note
%      triexcenter works also for vectors of points

x1 = P1(1);
x2 = P2(1);
x3 = P3(1);
y1 = P1(2);
y2 = P2(2);
y3 = P3(2);

xx = mean([x1 x2 x3]); yy = mean([y1 y2 y3]);
x1 = x1-xx; x2 = x2-xx; x3 = x3-xx;
y1 = y1-yy; y2 = y2-yy; y3 = y3-yy;

A = [x1 y1 1
     x2 y2 1
     x3 y3 1];

D = -[x1^2+y1^2 y1 1
      x2^2+y2^2 y2 1
      x3^2+y3^2 y3 1];

E = [x1^2+y1^2 x1 1
     x2^2+y2^2 x2 1
     x3^2+y3^2 x3 1];

F = -[x1^2+y1^2 x1 y1
      x2^2+y2^2 x2 y2
      x3^2+y3^2 x3 y3];

a = det(A);
d = det(D);
e = det(E);
f = det(F);

if abs(a) > eps,
   C = -[d e]/(2*a);
   r = sqrt((d^2+e^2)/(4*a^2) - f/a);
else
   C = [NaN NaN];
   r = Inf;
end

C = C + [xx yy];
