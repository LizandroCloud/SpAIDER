
function [angle CC] = anglesboundary2(B,x,y);

% anglesboundary2 - ThreshLab2/Triangulations --- find inside angles
%   Usage
%      angle = anglesboundary2(B,x,y);
%   Inputs
%      B          matrix with on each row an ordered list of indices into x 
%                 and y corresponding to vertices
%      x          vector of x-coordinates of vertices
%      y          vector of y-coordinates of vertices
%   Outputs
%      angle      vector of size(B) with angles
%      CC         boolean vector; components 1 for rows that represent
%                 boundaries listed in counterclockwise order, components 0 for
%                 rows that represent boundaries listed in clockwise order
%   Description
%      this function finds the inside angles of all corners on boundary B
%      anglesboundary2 is the vector version of anglesboundary; computes angles
%      of several regions, e.g. in a triangulation.
%   See also
%      help anglesboundary

transpose = 0;
if size(B,2) == 1,
   B = B';
   transpose = 1;
end
n = size(B,2);

B1 = shift(B,1);
xx = x(B1) - x(B); yy = y(B1) - y(B);
cosphases = xx./sqrt(xx.^2+yy.^2);
sinphases = yy./sqrt(xx.^2+yy.^2);
phases = acos(cosphases);
negsines = (sinphases<0);
phases(negsines) = 2*pi-phases(negsines);
angle = pi - (phases-shift(phases,-1));
angle = angle - floor(angle/(2*pi))*2*pi;

sumangle = sum(angle,2);
eps = 1.e-4;
outer = find(abs(sumangle-(n+2)*pi) < eps);
inner = find(abs(sumangle-(n-2)*pi) < eps);
CC = zeros(size(B,1),1);
CC(inner) = 1;
angle(outer,1:n) = 2*pi-angle(outer,1:n);

if transpose, angle = angle'; end
