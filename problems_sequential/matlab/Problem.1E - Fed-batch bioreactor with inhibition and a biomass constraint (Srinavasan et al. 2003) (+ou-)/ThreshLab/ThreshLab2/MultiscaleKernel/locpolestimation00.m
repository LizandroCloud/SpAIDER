
function yk = locpolestimation00(x,x0,y0,h,kerneltype,degree,side);

% locpolestimation00 - PiefLab2/Irregular/1D/ --- Local polynomial smoothing
%                                                 routine
%   Usage
%      for local polynomial smoothing
%         yk = kernelestimation00(x,x0,y0,h,kerneltype,degree,side);
%      for local polynomial density estimation
%         yk = kernelestimation00(x,X,h,kerneltype,degree,side);   
%   Inputs
%      x      x-coordinates in which kernel is evaluated
%      x0     x-coordinates of sample locations in smoothing problem
%      y0     sample observations in smoothing problem
%      X      observations/sample vector in desnity estimation problem
%      h      kernel bandwidth
%      kerneltype  string-type variable, see evalkernel
%      degree degree of polynomial (default 0)
%      side   -1,0,1: left, both sides, right
%   Outputs
%      yk     local polynomial smoothing evaluated in x
%   Notes
%      x, x0 and y0 are supposed to be arranged as rows.
%      if y0 is a square matrix of size equal to the length of x0, then the
%      smoothing operation is performed on all ROWS.
%   See also
%      help evalkernel
%      help kernelestimation01

x = row(x);
x0 = row(x0);

nobs = length(x0);
neval = length(x);
[nrows ny] = size(y0);
if nrows == 1 & ny == 1,
   if nargin >= 6, side = degree; else side = 0; end
   degree = kerneltype;
   kerneltype = h;
   h = y0;
   y0 = ones(1,nobs)/nobs;
   densityest = 1;
else
   if nargin < 7, side = 0; end
   densityest = 0;
   if ny ~= nobs, 
      if nrows ~= nobs, error('dimension mismatch x0 and y0')
      else nrows = ny; y0 = y0'; end
   end
end
if ischar(side),
   switch side(1)
   case {'R','r'},
      side = 1;
   case {'L','l'},
      side = -1;
   otherwise,
      side = 0;
   end
end

switch kerneltype
case {'Gaussian','gaussian','gauss','Gauss'}
   kernelsupport = Inf;
otherwise
   kernelsupport = h;
end

yk = zeros(nrows,neval);
sk = zeros(size(x));
for i = 1:neval,
   io = (abs(x0-x(i))<kernelsupport);
   if side == 1,
      io = io & (x0>=x(i));
   elseif side == -1,
      io = io & (x0<=x(i));
   end
   io = find(io);
   if ~isempty(io),
      xo = x0(io); 
      jo = find(max(abs(y0(1:nrows,io)),[],2)>eps);
      wo = evalkernel(x(i)-x0(io),h,kerneltype);
      for j = row(jo),
         yo = y0(j,io);
         yk(j,i) = weightedLSpol(x(i),xo,yo,wo,degree);
      end
   end
% plot(yk)
% pause
end
