
function yk = kernelestimation01(x,x0,y0,h,kerneltype,side);

% kernelestimation01 - PiefLab2/Irregular/1D/ --- Kernel smoothing routine
%   Usage
%      for kernel smoothing
%         yk = kernelestimation01(x,x0,y0,h,kerneltype,side);
%      for kernel density estimation
%         yk = kernelestimation01(x,X,h,kerneltype,side);   
%   Inputs
%      x      x-coordinates in which kernel is evaluated
%      x0     x-coordinates of sample locations in smoothing problem
%      y0     sample observations in smoothing problem
%      X      observations/sample vector in desnity estimation problem
%      h      kernel bandwidth
%      kerneltype  string-type variable, see evalkernel
%      side   -1,0,1: left, both sides, right
%   Outputs
%      yk     kernel smoothing evaluated in x
%   Notes
%      x, x0 and y0 are supposed to be arranged as rows.
%      if y0 is a square matrix of size equal to the length of x0, then the
%      smoothing operation is performed on all ROWS.
%      As a consequence
%      kernelestimation01(row(x),row(x),sparse(eye(N)),hj,kerneltype,side);
%      This version allows cut-off/one-sided kernels. (input variable side)
%      This version (01) is computationally more efficient than version 00, but
%      improvements are still possible.
%   See also
%      help evalkernel

x = row(x);
x0 = row(x0);

nobs = length(x0);
neval = length(x);
[nrows ny] = size(y0);
if nrows == 1 & ny == 1,
   if nargin >= 5, side = kerneltype; else side = 0; end
   kerneltype = h;
   h = y0;
   y0 = ones(1,nobs)/nobs;
   densityest = 1;
else
   if nargin < 6, side = 0; end
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
for i = 1:nobs,
   ii = (abs(x-x0(i))<kernelsupport);
   if side == 1,
      ii = ii & (x>=x0(i));
   elseif side == -1,
      ii = ii & (x<=x0(i));
   end
   ii = find(ii);
   if ~isempty(ii),
      jj = find(abs(y0(1:nrows,i))>eps);
      yk(jj,ii) = yk(jj,ii) + y0(jj,i)*evalkernel(x(ii)-x0(i),h,kerneltype);
      sk(ii) = sk(ii) + evalkernel(x(ii)-x0(i),h,kerneltype);
   end
% plot(yk)
% pause
end
if ~densityest,
   % lines below: points with no prediction point within the bandwidth take the
   % closest observation as prediction value
   oo = find(abs(sk)<eps);
   if ~isempty(oo) & ...
      ~ismember(kerneltype,{'none','dirac','Dirac','kronecker','Kronecker'})
      disp(['points without prediction point within bandwidth:' num2str(oo)])
   end
   for j = row(oo),
      [xi i] = min(abs(x(j)-x0(i)));
      yk(1:nrows,j) = y0(1:nrows,i);
      sk(j) = 1;
   end
   sk = repmat(sk,nrows,1);
   yk = yk./sk;
   % yk = yk/(nobs*h);
else
   yk = yk/h;
end
