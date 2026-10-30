
function f = laplaceplusnormalpdf(x,a,sigma);

% laplaceplusnormalpdf -- probablity density function of W = V+U
%                        with V ~ Laplace(a) and U ~ N(0,sigma^2)
%   Usage
%     f = laplaceplusnormalpdf(w,a,sigma);
%   Inputs
%     w     value in which pdf is evaluated
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     f     pdf of  W = V+U
%   Description
%     the definition of the Laplace distribution used here has pdf
%     f(v) = (a/2)*exp(-a*|v|)
%     i.e., V ~ Laplace(a), then E(V) = 0; var(V) = 2/a^2
%   See also
%     help laplaceplusnormalcdf
%     help JScauchyplusnormalpdf


if nargin < 3, sigma = 1; end
if nargin < 2, a = 0.5; end

if isinf(a), f = ones(size(x))/2; return, end

eps = 1.e-308;
% take eps such that 1/eps is still finite. It may happen that eps > 0, but
% 1/eps yields Inf.
s = sign(x);
oo = find(gauss(abs(x)+sigma^2*a,sigma)>eps);
xmax = max([min(abs(x)) abs(x(oo))]); 
     % if oo is empty, xmax will be min(abs(x))
xx = s.*min(abs(x),xmax);
f = (a/2)*gauss(xx,sigma).*...
           (cumgauss(xx-sigma^2*a,sigma)./gauss(xx-sigma^2*a,sigma)+...
           Rcumgauss(xx+sigma^2*a,sigma)./gauss(xx+sigma^2*a,sigma));
% alternative expression (possible numerically more stable)
% f = a*exp(sigma^2*a^2/2)/2* ...
%     (exp(-a*xx).*cumgauss(xx-a*sigma^2,sigma) + ...
%      exp(a*xx).*Rcumgauss(xx+a*sigma^2,sigma));

oo = find(gauss(abs(x)+sigma^2*a,sigma)<eps);
x = abs(x);
fapprox = (a/2)*exp(sigma^2*a^2/2)*exp(-a*x) + ...
          (a/2)*gauss(x,sigma)*sigma^2./(x+sigma^2*a);
f(oo) = fapprox(oo);
