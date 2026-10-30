
function y = weightedLSpol(x,x0,y0,w0,p);

% weightedLSpol - ThreshLab2/Irregular/MultiscaleKernel/ --- 
%                                          Weighted Least Squares Polynomial
%   Usage
%      y = weightedLSpol(x,x0,y0,w0,p);
%   Inputs
%      x      x-coordinates in which polynomial is evaluated
%      x0     x-coordinates of locations where polynonial is constructed
%      y0     y-coordinates of locations where polynonial is constructed
%      w0     weights corresponding to y0 (default = 1)
%      p      degree of polynomial
%   Outputs
%      y      polynomial evaluated in x
%   Notes
%   See also

if nargin < 5, p = 0; end

n = length(x);
n0 = length(x0);
if nargin < 5, w0 = ones(n0,1); end

if size(y0,1) ~= n0, y0 = y0'; end
if size(y0,1) ~= n0, error('dimension mismatch x0 and y0'); end
if size(w0,1) ~= n0, w0 = w0'; end
if size(w0,1) ~= n0, error('dimension mismatch x0 and w0'); end

p = min(p,n0);
x0 = column(x0);
xo = mean(x0);
X0 = zeros(n0,p+1); X0(1:n0,1) = ones(n0,1); 
for k = 1:p, X0(1:n0,k+1) = (x0-xo).^k; end
% working with orthogonal polynomials would be more stable
W0 = diag(w0);

coef = (X0'*W0*X0)\(X0'*W0*y0);

transposeresult = 0;
if size(x,2) == 1, transposeresult = 1; x = x'; end
X = zeros(n,p+1); X(1:n,1) = ones(n,1);
for k = 1:p, X(1:n,k+1) = (x-xo).^k; end

y = X*coef;
if transposeresult, y = y'; end
