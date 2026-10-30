
function F = zeroinflatedlaplaceplusnormalcdf(x,p,a,sigma);

% laplaceplusnormalcdf -- cumulative distribution function of W = VX+U
%                         with X ~ Bernoulli(p)
%                              V ~ Laplace(a)
%                              U ~ N(0,sigma^2)
%   Usage
%     F = zeroinflatedlaplaceplusnormalcdf(w,p,a,sigma);
%   Inputs
%     w     value in which cdf is evaluated
%     p     parameter of Bernoulli distribution (default 0.05)
%     a     parameter of Laplace distribution (default 0.5)
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     F     cdf of  W = VX+U
%   Description
%     the definition of the Laplace distribution used here has pdf
%     f(v) = (a/2)*exp(-a*|v|)
%     i.e., V ~ Laplace(a), then E(V) = 0; var(V) = 2/a^2
%   See also
%     help laplaceplusnormalcdf
%     help JScauchyplusnormalcdf


if nargin < 4, sigma = 1; end
if nargin < 3, a = 1/2; end
if nargin < 2, p = 0.05; end

F = p*laplaceplusnormalcdf(x,a,sigma) + (1-p)*cumgauss(x,sigma);

