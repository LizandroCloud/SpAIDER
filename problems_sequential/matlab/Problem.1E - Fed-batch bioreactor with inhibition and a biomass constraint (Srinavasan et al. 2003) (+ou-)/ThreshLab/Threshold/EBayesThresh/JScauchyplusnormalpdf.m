
function f = JScauchyplusnormalpdf(x,sigma);

% JScauchyplusnormalpdf -- probablity density function of W = V+U
%                        with V ~ JScauchy(sigma) and U ~ N(0,sigma^2)
%   Usage
%     f = JScauchyplusnormalpdf(w,sigma);
%   Inputs
%     w     value in which pdf is evaluated
%     sigma st.dev. of normal distribution (default 1)
%   Outputs
%     f     pdf of  W = V+U
%   Reference
%     The Cauchy-like distribution JScauchy comes from:
%     I. Johnstone and B. Silverman, Needles and Straw in Haystacks, Annals of
%                                    Statistics, 2004.
%   See also
%     help laplaceplusnormalpdf
%     help JScauchypdf
%     help JScauchyplusnormalcdf

if (nargin == 1)
   sigma = 1;
end

xx = x+(x==0);
f1 = (sigma/sqrt(2*pi))*(1-exp(-xx.^2/(2*sigma)))./xx.^2;
f0 = 1/(2*sqrt(2*pi)*sigma);
f = f0.*(x==0)+f1.*(x~=0);

