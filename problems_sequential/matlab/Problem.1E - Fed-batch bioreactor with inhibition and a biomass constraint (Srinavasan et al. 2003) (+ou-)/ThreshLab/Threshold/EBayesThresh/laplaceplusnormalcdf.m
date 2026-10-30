
function F = laplaceplusnormalcdf(x,a,sigma);

% laplaceplusnormalcdf -- cumulative distribution function of W = V+U
%                         with V ~ Laplace(a) and U ~ N(0,sigma^2)
%   Usage
%     F = laplaceplusnormalcdf(w,a,sigma);
%   Inputs
%     w     value in which cdf is evaluated
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     F     cdf of  W = V+U
%   Description
%     the definition of the Laplace distribution used here has pdf
%     f(v) = (a/2)*exp(-a*|v|)
%     i.e., V ~ Laplace(a), then E(V) = 0; var(V) = 2/a^2
%   See also
%     help laplaceplusnormalpdf
%     help JScauchyplusnormalcdf


if nargin < 3, sigma = 1; end
if nargin < 2, a = 1/2; end

[n1 n2] = size(x);
if n1>1 & n2>1, 
   n = n1*n2;
   x = reshape(x,1,n);
end

eps = 1.e-308;
% take eps such that 1/eps is still finite. It may happen that eps > 0, but
% 1/eps yields Inf.
s = sign(x);
oo = find(gauss(abs(x)+sigma^2*a,sigma)>eps);
xmax = max([min(abs(x)) abs(x(oo))]);
     % if oo is empty, xmax will be min(abs(x))
xx = s.*min(abs(x),xmax);
F = cumgauss(xx,sigma) + ...
    (1/2)*gauss(xx,sigma).*...
           (-cumgauss(xx-sigma^2*a,sigma)./gauss(xx-sigma^2*a,sigma)+...
            Rcumgauss(xx+sigma^2*a,sigma)./gauss(xx+sigma^2*a,sigma));
oo = find(gauss(abs(x)+sigma^2*a,sigma)<eps);
x = abs(x);
Fapprox = cumgauss(xx,sigma) + ...
          (-1/2)*exp(sigma^2*a^2/2)*exp(-a*x) + ...
          (1/2)*gauss(x,sigma)*sigma^2./(x+sigma^2*a);
F(oo) = Fapprox(oo);
if n1>1 & n2>1,
   F = reshape(F,n1,n2);
end
