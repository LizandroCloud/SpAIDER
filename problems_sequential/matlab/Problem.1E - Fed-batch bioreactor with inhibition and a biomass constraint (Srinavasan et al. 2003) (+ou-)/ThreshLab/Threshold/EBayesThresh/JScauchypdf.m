
function f = JScauchypdf(x,sigma);

% JScauchypdf -- probablity density function of Johnstone's and Silverman's
%                Caucly-like prior
%   Usage
%     f = JScauchypdf(x,sigma);
%   Inputs
%     x     value in which pdf is evaluated
%     sigma parameter, default = 1
%   Outputs
%     f     pdf
%   Reference
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   Notes
%     This function calls JScauchy1pdf, which is the special case of sigma = 1;
%     for sigma =/= 1, we have f = JScauchy1pdf(x/sigma)/sigma;
%   See also
%     help JScauchy1pdf
%     help JScauchycdf
%     help JScauchyplusnormalpdf

if (nargin == 1)
   sigma = 1;
end

f = JScauchy1pdf(x/sigma)/sigma;


