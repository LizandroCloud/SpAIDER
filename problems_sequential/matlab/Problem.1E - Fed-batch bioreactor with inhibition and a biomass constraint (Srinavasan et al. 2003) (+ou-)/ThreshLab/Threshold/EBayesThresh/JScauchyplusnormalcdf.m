
function F = JScauchyplusnormalcdf(x,sigma);

% JScauchyplusnormalcdf -- cumulative distribution function of W = V+U
%                        with V ~ JScauchy(sigma) and U ~ N(0,sigma^2)
%   Usage
%     F = JScauchyplusnormalcdf(w,sigma);
%   Inputs
%     w     value in which cdf is evaluated
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     F     cdf of  W = V+U
%   Reference
%     The Cauchy-like distribution JScauchy comes from:
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   See also
%     help laplaceplusnormalcdf
%     help JScauchycdf
%     help JScauchyplusnormalpdf

if (nargin == 1)
   sigma = 1;
end


F0 = cumgauss(x);
xx = x + (x == 0);
F1 = cumgauss(xx) + (exp(-xx.^2/(2*sigma^2))./xx - 1./xx)*sigma/sqrt(2*pi);
F = F0 .* (x == 0) + F1 .* (x ~= 0);

