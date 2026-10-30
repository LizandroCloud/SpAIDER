
function F = JScauchycdf(x,sigma);

% JScauchycdf -- cumulative distribution function of Johnstone's and Silverman's
%                Caucly-like prior
%   Usage
%     F = JScauchycdf(x,sigma);
%   Inputs
%     x     value in which cdf is evaluated
%     sigma parameter, default = 1
%   Outputs
%     F     cdf
%   Reference
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   Notes
%     This function calls JScauchy1cdf, which is the special case of sigma = 1;
%     for sigma =/= 1, we have F = JScauchy1cdf(x/sigma);
%   See also
%     help JScauchy1cdf
%     help JScauchypdf
%     help JScauchyplusnormalcdf

if (nargin == 1)
   sigma = 1;
end

F = JScauchy1cdf(x/sigma);


