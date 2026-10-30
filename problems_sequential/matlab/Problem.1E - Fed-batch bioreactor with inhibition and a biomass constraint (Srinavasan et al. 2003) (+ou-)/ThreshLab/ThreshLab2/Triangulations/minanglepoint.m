
function [point,alfa] = minanglepoint(p1,p2,points,x,y);

% minanglepoint - ThreshLab2/Triangulations --- finds vertex v in points with
%                 minimum angle p1-v-p2
%   Usage
%      [point,alfa] = minanglepoint(p1,p2,points,x,y);
%   Inputs
%      p1      index into x and y
%      p1      index into x and y
%      points  vector of indices into x and y
%      x       x-coordinates of data cloud (real vector)
%      y       y-coordinates of data cloud (real vector)
%   Outputs
%      point   element of points (minimizing angle p1-points-p2)
%      alfa    angle
%   Description
%   See also

i = sqrt(-1);
z21 = x(p1) - x(p2) + i * (y(p1) - y(p2));
alfa21 = angle(z21);
alfa21 = alfa21 + 2*pi *(alfa21<0);
z23 = x(points) - x(p2) + i * (y(points) - y(p2));
alfa23 = angle(z23);
alfa23 = alfa23 + 2*pi *(alfa23<0);

alfa = alfa23 - alfa21;
alfa = alfa + 2*pi *(alfa<=0);

[alfa, i] = min(alfa);
point = points(i);


