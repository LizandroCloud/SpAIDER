
function [vertexangle CC] = anglesboundary(B,x,y);

% anglesboundary - ThreshLab2/Triangulations --- find inside angles
%   Usage
%      [vertexangle CC] = anglesboundary(B,x,y);
%   Inputs
%      B                ordered list of indices into x and y corresponding to
%                       vertices
%      x                vector of x-coordinates of vertices
%      y                vector of y-coordinates of vertices
%   Outputs
%      vertexangle      vector of size(B) with angles 
%      CC               boolean; 1 if B is ordered counterclockwise
%                                0 if B is ordered clockwise
%   Description
%      this function finds the inside angles of all corners on boundary B

n = length(B);
cyc = 0;
if B(n) == B(1),
   B = B(1:n-1);
   n = n-1;
   cyc = 1;
end
B1 = shift(B,1);
xx = x(B1) - x(B); yy = y(B1) - y(B);
cosphases = xx./sqrt(xx.^2+yy.^2);
sinphases = yy./sqrt(xx.^2+yy.^2);
phases = acos(cosphases);
negsines = (sinphases<0);
phases(negsines) = 2*pi-phases(negsines);
% the above 5 lines can be replaced by matlab's function
%  phases = angle(xx+sqrt(-1)*yy);
vertexangle = pi - (phases-shift(phases,-1));
vertexangle = vertexangle - floor(vertexangle/(2*pi))*2*pi;

eps = 1.e-4;
o = (abs(vertexangle)<eps | abs(vertexangle-2*pi)<eps);
noto = ~o;
o = find(o);
   % vertices with angles near 0 or 2*pi will be eliminated from subsequent
   % steps, because the attributed angle value may be incorrect due to flipping
   % between 0 and 2*pi. The corrected value will be attributed afterwards
n0 = length(o);
n1 = length(noto);

sumangle = sum(vertexangle(noto));
if abs(sumangle-(n1+2)*pi) < eps,
   CC = 0;
   vertexangle = 2*pi-vertexangle;
elseif abs(sumangle-(n1-2)*pi) < eps,
   CC = 1;
else
   s = num2str(sumangle/pi);
   n = num2str(n1);
   error(['incorrect sum of angles: sum = ' s ' * pi; n = ' n ])
end

for i=1:n0,
   oi = o(i);
   vertexangle(oi) = 2*pi*pointinside(x(B(oi)),y(B(oi)),B(noto),x,y);
end
