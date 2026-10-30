
function f = JScauchy1cdf(u);

% JScauchy1cdf -- cumulative distribution function of Johnstone's and
%                 Silverman's Caucly-like prior
%   Usage
%     f = JScauchy1cdf(u);
%   Inputs
%     u     value in which cdf is evaluated
%   Outputs
%     f     cdf
%   Reference
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   See also
%     help JScauchycdf
%     help JScauchy1pdf
%     help JScauchyplusnormalcdf

uu = abs(u);
F0 = exp(uu.^2/2).*Rcumgauss(uu,1);
F = F0.*(u<0)+(1-F0).*(u>=0);
