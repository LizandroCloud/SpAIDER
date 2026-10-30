
function yy = pcwsconstant(xx,x,y);

% pcwsconstant - ThreshLab2/Polylines --- evaluate piecewise constant defined by 
% knots
%   Usage
%     yy = pcwsconstant(xx,x,y);
%   Inputs
%     x     x-coordinates of knots
%     y     y-coordinates of knots
%     length(x) = length(y)
%     xx    x-coordinates of evaluation points
%   Outputs
%     yy    y-coordinates of evaluation points
%   Description
%     defines function values yy in xx of piecewise constant interpolating in 
%     (x,y)
%     y can be a matrix, the evaluation then takes place in a row-by-row or a
%     column-by-column fashion.
%   See also
%     help polyline

[r c] = size(y);
if c == size(x,2), y = y'; end
x = row(x);
[r c] = size(y);
if r == length(x), y = y'; end
[r c] = size(y);
if c ~= length(x),
   error('abscis-ordinate mismatch in polyline knots (x,y)')
end

[x o] = sort(x);
y = y(:,o);

p = length(x);
[rr cc] = size(xx);
xx = row(xx);
yy = zeros(size(xx));
for k = 1:p-1,
   ii = find(xx>=x(k) & xx<=x(k+1));
   for j = 1:r,
      yy(j,ii) = y(j,k);
   end
end
ii = find(xx>=x(p));
for j = 1:r,
   yy(j,ii) = y(j,p);
end
if cc == 1, yy = yy'; end
