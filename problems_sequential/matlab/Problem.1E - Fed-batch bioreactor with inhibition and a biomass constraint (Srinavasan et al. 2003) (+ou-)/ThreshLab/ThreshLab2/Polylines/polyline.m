
function yy = polyline(xx,x,y);

% polyline - ThreshLab2/Polylines --- evaluate polyline defined by knots
%   Usage
%     yy = polyline(xx,x,y);
%   Inputs
%     x     x-coordinates of knots
%     y     y-coordinates of knots
%     length(x) = length(y)
%     xx    x-coordinates of evaluation points
%   Outputs
%     yy    y-coordinates of evaluation points
%   Description
%     defines function values yy in xx of polyline interpolating in (x,y)
%     y can be a matrix, the evaluation then takes place in a row-by-row or a
%     column-by-column fashion.
%   See also
%     help polylineproduct

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
if cc == 1,
   xx = xx';
elseif rr>1 & r>1,
   error('cannot process multiple inputs and multiple outputs at once')
end
% xx = row(xx);
if r == 1,
   yy = zeros(size(xx));
   for k = 1:p-1,
      ii = find(xx>=x(k) & xx<=x(k+1));
      yy(ii) = (xx(ii)-x(k))/(x(k+1)-x(k))*y(k+1) + ...
                  (x(k+1)-xx(ii))/(x(k+1)-x(k))*y(k);
   end
else
   yy = zeros(r,length(xx));
   for k = 1:p-1,
      ii = find(xx>=x(k) & xx<=x(k+1));
      for j = 1:r,
         yy(j,ii) = (xx(ii)-x(k))/(x(k+1)-x(k))*y(j,k+1) + ...
                  (x(k+1)-xx(ii))/(x(k+1)-x(k))*y(j,k);
      end
   end
end
if cc == 1, yy = yy'; end
